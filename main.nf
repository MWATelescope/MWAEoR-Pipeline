#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Ensure the raw visibility files are present, or download via ASVO
process ensureRaw {
    input:
    val obsid
    output:
    tuple val(obsid), path("${obsid}.metafits"), path("${obsid}_2*.fits")

    // persist results in outdir, process will be skipped if files already present.
    storeDir "$params.outdir/$obsid/raw"
    // allow multiple retries
    maxRetries 5
    // exponential backoff: sleep for 2^attempt hours after each fail
    errorStrategy { sleep(Math.pow(2, task.attempt) * 60*60*1000 as long); 'retry' }

    script:
    """
    # echo commands, exit on any failures
    set -ex
    export MWA_ASVO_API_KEY="${params.asvo_api_key}"

    env | sort
    df \$TMPDIR

    # TODO: make this profile dependent
    export http_proxy="http://proxy.per.dug.com:3128"
    export https_proxy="http://proxy.per.dug.com:3128"
    export all_proxy="proxy.per.dug.com:3128"
    export ftp_proxy="http://proxy.per.dug.com:3128"

    function ensure_disk_space {
        local needed=\$1
        while read -r avail; do
            avail_bytes=\$((\$avail * 1000))
            if [[ \$avail_bytes -lt \$needed ]]; then
                echo "Not enough disk space available in \$PWD, need \$needed B, have \$avail_bytes B"
                exit 1
            fi
        done < <(df --output=avail . | tail -n 1)
    }

    # download any ASVO jobs for the obsid that are ready
    function get_first_ready_job {
        # extract id url and size from ready download vis jobs
        ${params.giant_squid} list -j --types download_visibilities --states ready -- $obsid \
            | tee /dev/stderr \
            | ${params.jq} -r '.[]|[.jobId,.files[0].fileUrl//"",.files[0].fileSize//"",.files[0].fileHash//""]|@tsv' \
            | tee ready.tsv
        read -r jobid url size hash < ready.tsv
        ensure_disk_space \$size || exit \$?
        # giant-squid download --keep-zip --hash -v \$jobid 2>&1 | tee download.log
        # if [ \${PIPESTATUS[0]} -ne 0 ]; then
        #     echo "Download failed, see download.log"
        #     exit 1
        # fi
        wget \$url -O \$jobid.tar --progress=dot:giga --wait=60 --random-wait
        if [ \$? -ne 0 ]; then
            echo "Download failed"
            exit 1
        fi
        sha1=\$(sha1sum \$jobid.tar | cut -d' ' -f1)
        if [ "\$sha1" != "\$hash" ]; then
            echo "Download failed, hash mismatch"
            exit 1
        fi
        [ -d "${task.storeDir}" ] || mkdir -p "${task.storeDir}"
        tar -xf \$jobid.tar -C "${task.storeDir}"
        return \$?
    }

    # download any ready jobs, exit if success, else try and submit a job,
    # if it's already there this will fail, and we can suppress the
    # warning, otherwise we can download and exit
    get_first_ready_job && exit 0 || ${params.giant_squid} submit-vis --delivery "acacia" $obsid -w \
        && get_first_ready_job && exit 0 || true

    # extract id and state from pending download vis jobs
    ${params.giant_squid} list -j --types download_visibilities --states queued processing -- $obsid \
        | ${params.jq} -r '.[]|[.jobId,.jobState]|@tsv' \
        | tee pending.tsv

    # for each pending job for the obsid, wait for completion and try downloading again
    while read jobid state; do
        echo "[obs:$obsid]: waiting until \$jobid with state \$state is Ready"
        ${params.giant_squid} wait -j -- \$jobid
        # a new job is ready, so try and get it
        get_first_ready_job
        exit \$?
    done <pending.tsv
    exit 1
    """
}

process metafitsStats {
    input:
    tuple val(obsid), path("${obsid}.metafits")
    output:
    tuple val(obsid), path("${obsid}.metafits.json")

    scratch false
    storeDir "$params.outdir/${obsid}/raw"

    // if this script returns failure, don't progress this vis
    errorStrategy 'ignore'

    // without this, nextflow checks for the output before fs has had time
    // to sync, causing it to think the job has failed.
    afterScript "sleep 20s; ls ${task.storeDir}"

    module 'python/3.9.7'

    script:
    """
    #!/usr/bin/env python

    from astropy.io import fits
    from json import dump as json_dump
    from collections import OrderedDict

    with \
        open('${obsid}.metafits.json', 'w') as out, \
        fits.open("${obsid}.metafits") as hdus \
    :
        data = OrderedDict()
        for key in hdus[0].header:
            if type(hdus[0].header[key]) not in [bool, int, str, float]:
                continue
            data[key] = hdus[0].header[key]
            if key in [ 'RECVRS', 'DELAYS', 'CHANNELS', 'CHANSEL' ]:
                data[key] = [*filter(None, map(lambda t: t.strip(), data[key].split(',')))]

        data['FLAGGED_INPUTS'] = [
            r['TileName'] + r['Pol']
            for r in hdus[1].data
            if r['Flag'] == 1
        ]

        json_dump(data, out, indent=4)
    """
}

// Ensure preprocessed observation is present, or preprocess raw with Birli
process ensurePrep {
    input:
    tuple val(obsid), path("${obsid}.metafits"), path("*") // <-preserves names of fits files
    output:
    tuple val(obsid), path("birli_${obsid}.uvfits"), path("${obsid}*.mwaf"), path("birli_prep.log")

    storeDir "$params.outdir/${obsid}/prep"

    // label jobs that need a bigger cpu allocation
    label "cpu"

    // this is necessary for singularity
    stageInMode "copy"
    module "singularity"

    script:
    """
    export flag_template="${obsid}_%%.mwaf"
    if [ ${obsid} -gt 1300000000 ]; then
        flag_template="${obsid}_ch%%%.mwaf"
    fi
    ${params.birli} \
        -u "birli_${obsid}.uvfits" \
        -f \$flag_template \
        -m "${obsid}.metafits" \
        --avg-time-res ${params.prep_time_res_s} \
        --avg-freq-res ${params.prep_freq_res_khz} \
        ${obsid}*.fits 2>&1 | tee birli_prep.log
    echo task.hash=$task.hash >> birli_prep.log
    """
}

// Simple check for vis, to catch "past end of file" errors
process prepStats {
    input:
    tuple val(obsid), path("birli_${obsid}.uvfits")
    output:
    tuple val(obsid), path("${obsid}_prep_stats.json")

    scratch false
    storeDir "$params.outdir/${obsid}/prep"

    // if this script returns failure, don't progress this vis
    errorStrategy 'ignore'

    // without this, nextflow checks for the output before fs has had time
    // to sync, causing it to think the job has failed.
    afterScript "sleep 20s; ls ${task.storeDir}"

    module 'python/3.9.7'

    script:
    """
    #!/usr/bin/env python
    # task.hash=$task.hash

    from astropy.io import fits
    from json import dump as json_dump
    from collections import OrderedDict

    with \
        open('${obsid}_prep_stats.json', 'w') as out, \
        fits.open("birli_${obsid}.uvfits") as hdus \
    :
        data = OrderedDict()
        first_shape=hdus[0].data[0].data.shape
        data['num_chans']=first_shape[2]
        all_times=hdus[0].data['DATE']
        data['num_times']=len(set(all_times))
        data['num_baselines']=len(all_times)//data['num_times']

        json_dump(data, out, indent=4)
    """
}

// QA tasks for flags.
process flagQA {
    input:
    tuple val(obsid), path("${obsid}.metafits"), path("*") // <-preserves names of mwaf files
    output:
    tuple val(obsid), path("${obsid}_occupancy.json")

    storeDir "$params.outdir/${obsid}/prep"

    // if this script returns failure, don't progress this obs to calibration, just ignore
    errorStrategy 'ignore'

    // without this, nextflow checks for the output before fs has had time
    // to sync, causing it to think the job has failed.
    afterScript "sleep 20s; ls ${task.storeDir}"

    module 'python/3.9.7'

    script:
    """
    #!/usr/bin/env python

    from astropy.io import fits
    from json import dump as json_dump
    from glob import glob
    from os.path import basename
    import re

    RE_MWAX_NAME = (
        r"\\S+_ch(?P<rec_chan>\\d{3}).mwaf"
    )
    RE_MWA_LGCY_NAME = (
        r"\\S+_(?P<gpubox_num>\\d{2}).mwaf"
    )

    def parse_filename(name, metafits_coarse_chans):
        result = {}
        if match:=re.match(RE_MWAX_NAME, name):
            result.update({
                'type': 'mwax',
                'rec_chan': int(match.group('rec_chan')),
            })
            result['gpubox_num'] = int(metafits_coarse_chans.index(result['rec_chan']))
        elif match:=re.match(RE_MWA_LGCY_NAME, name):
            result.update({
                'type': 'legacy',
                'gpubox_num': int(match.group('gpubox_num')),
            })
            result['rec_chan'] = int(metafits_coarse_chans[result['gpubox_num']-1])
        return result

    def split_strip_filter(str):
        return list(filter(None, map(lambda tok: tok.strip(), str.split(','))))

    paths = sorted(glob("*.mwaf"))
    total_occupancy = 0
    num_coarse_chans = len(paths)
    with \
        open('${obsid}_occupancy.json', 'w') as out, \
        fits.open("${obsid}.metafits") as meta \
    :
        data = {'obsid': ${obsid}, 'channels': {}}
        metafits_coarse_chans = [*map(int, split_strip_filter(meta[0].header['CHANNELS']))]
        for path in paths:
            filename_info = parse_filename(path, metafits_coarse_chans)
            if not data.get('type'):
                data['type'] = filename_info.get('type', '???')
            chan_id = filename_info.get('rec_chan', basename(path))
            data['channels'][chan_id] = {}
            with fits.open(path) as hdus:
                flag_data = hdus[1].data['FLAGS']
                occupancy = flag_data.sum() / flag_data.size
                data['channels'][chan_id]['occupancy'] = occupancy
                total_occupancy += occupancy
        total_occupancy /= num_coarse_chans
        data['total_occupancy'] = total_occupancy
        json_dump(data, out, indent=4)
    """
}

// ensure calibration solutions are present, or calibrate prep with hyperdrive
process ensureCalSol {
    input:
    tuple val(obsid), path("${obsid}.metafits"), path("${obsid}.uvfits"), path("cal_args.csv")
    output:
    tuple val(obsid), path("hyp_soln_${obsid}_{30l_src4k,50l_src4k}.fits"), path("hyp_di-cal_${obsid}_*.log")

    storeDir "$params.outdir/${obsid}/cal$params.cal_suffix"

    // TODO: figure out why this keeps failing
    errorStrategy 'ignore'

    // label jobs that need a bigger gpu allocation
    label "gpu"

    module "cuda/11.3.1:gcc-rt/9.2.0"
    stageInMode "copy"

    script:
    """
    export HYPERDRIVE_CUDA_COMPUTE=${params.cuda_compute}
    # use this control which GPU we're sending the job to
    export CUDA_VISIBLE_DEVICES=0

    set -ex
    ls -al
    export names=""
    export num_gpus="\$(nvidia-smi -L | wc -l)"
    if [ \$num_gpus -eq 0 ]; then
        echo "no gpus found"
        exit 1
    fi
    # for each calibration parameter set in cal_args.csv, run hyperdrive di-cal
    while IFS=, read -r name cal_args; do
        # on first iteration, this does nothing. on nth, wait for all background jobs to finish
        if [[ \$CUDA_VISIBLE_DEVICES -eq 0 ]]; then
            wait \$(jobs -rp)
        fi
        # hyperdrive di-cal backgrounded in a subshell
        (
            ${params.hyperdrive} di-calibrate \${cal_args} \
                --data "${obsid}.metafits" "${obsid}.uvfits" \
                --beam "${params.beam_path}" \
                --source-list "${params.sourcelist}" \
                --outputs "hyp_soln_${obsid}_\${name}.fits" \
                > hyp_di-cal_${obsid}_\${name}.log
        ) &
        names="\${names} \${name}"
        # increment the target device mod num_gpus
        CUDA_VISIBLE_DEVICES=\$(( CUDA_VISIBLE_DEVICES+1 % \${num_gpus} ))
    done < <(cut -d',' -f1,2 cal_args.csv | sort | uniq)

    # wait for all the background jobs to finish
    wait \$(jobs -rp)

    # print out important info from log
    grep -iE "err|warn|hyperdrive|chanblocks|reading|writing|flagged" *.log
    """
}

// ensure calibrated vis are present, or apply calsols to prep with hyperdrive
process ensureCalVis {
    input:
    tuple val(obsid), path("${obsid}.metafits"), path("${obsid}.uvfits"), path("*"), \
        val(cal_name), val(apply_name), val(apply_args)
    output:
    tuple val(obsid), path("hyp_${obsid}_${cal_name}_${apply_name}.uvfits"), path("hyp_apply_${cal_name}_${apply_name}.log")

    storeDir "$params.outdir/${obsid}/cal$params.cal_suffix"

    // TODO: figure out why this keeps failing
    errorStrategy 'ignore'

    // label jobs that need a bigger gpu allocation
    label "cpu"

    module "gcc-rt/9.2.0"
    stageInMode "copy"

    script:
    """

    # hyperdrive solutions apply
    ${params.hyperdrive} solutions-apply ${apply_args} \
        --data "${obsid}.metafits" "${obsid}.uvfits" \
        --solutions "hyp_soln_${obsid}_${cal_name}.fits" \
        --outputs "hyp_${obsid}_${cal_name}_${apply_name}.uvfits" \
        | tee hyp_apply_${cal_name}_${apply_name}.log

    # print out important info from log
    grep -iE "err|warn|hyperdrive|chanblocks|reading|writing|flagged" *.log
    """
    // ms:
    // ${params.hyperdrive} solutions-apply \${apply_args} \
    //   --outputs "hyp_${obsid}_\${name}.uvfits" "hyp_${obsid}_\${name}.ms" \
    //   zip -r "hyp_${obsid}_\${name}.ms.zip" "hyp_${obsid}_\${name}.ms"
}

// QA tasks that can be run on each visibility file.
process visQA {
    input:
    tuple val(obsid), val(name), path("*") // <- *.uvfits
    output:
    tuple val(obsid), path("${name}_vis_metrics.json")

    storeDir "$params.outdir/${obsid}/vis_qa"

    // TODO: figure out why this fails sometimes
    errorStrategy 'ignore'

    label 'cpu'

    // needed for singularity
    stageInMode "copy"
    module "singularity"

    // without this, nextflow checks for the output before fs has had time
    // to sync, causing it to think the job has failed.
    afterScript "sleep 20s; ls ${task.storeDir}"

    script:
    """
    set -ex
    singularity exec --bind \$PWD:/tmp --writable-tmpfs --pwd /tmp --no-home ${params.mwa_qa_sif} \
        run_visqa.py *.uvfits --out "${name}_vis_metrics.json"
    """
}

// QA tasks that can be run on each calibration parameter set
process calQA {
    input:
    tuple val(obsid), val(name), path("${obsid}.metafits"), path("*") // <- hyp_soln_*.fits
    output:
    tuple val(obsid), path("${name}_{X,Y}.json")

    storeDir "$params.outdir/${obsid}/cal_qa"

    // TODO: figure out why this fails sometimes
    errorStrategy 'ignore'
    // maxRetries 1

    stageInMode "copy"
    module "singularity"

    // without this, nextflow checks for the output before fs has had time
    // to sync, causing it to think the job has failed.
    afterScript "sleep 20s; ls ${task.storeDir}"

    script:
    """
    set -ex
    ls -al
    singularity exec --bind \$PWD:/tmp --writable-tmpfs --pwd /tmp --no-home ${params.mwa_qa_sif} \
        run_calqa.py *.fits *.metafits X --out "${name}_X.json"
    singularity exec --bind \$PWD:/tmp --writable-tmpfs --pwd /tmp --no-home ${params.mwa_qa_sif} \
        run_calqa.py *.fits *.metafits Y --out "${name}_Y.json"
    """
}

process plotSolutions {
    input:
    tuple val(obsid), val(name), path("${obsid}.metafits"), path("*") // <- hyp_soln_*.fits
    output:
    tuple val(obsid), path("${name}_{phases,amps}.png")

    storeDir "$params.outdir/${obsid}/cal_qa"
    errorStrategy 'ignore'

    beforeScript 'sleep 1s; ls -al'
    afterScript 'ls -al; sleep 1s'

    script:
    """
    hyperdrive solutions-plot -m "${obsid}.metafits" *.fits
    """
}

// create dirty iamges of xx,yy,v
process wscleanDirty {
    input:
    tuple val(obsid), val(name), path("${obsid}.metafits"), path("*") // <- hyp_*.uvfits
    output:
    tuple val(obsid), val(name), path("wsclean_${name}*-MFS-{XX,YY,V}-dirty.fits")

    storeDir "$params.outdir/${obsid}/img${params.img_suffix}"

    // TODO: figure out why this keeps failing
    errorStrategy 'ignore'

    label "cpu"
    label "img"
    maxRetries 1

    // this is necessary for singularity
    stageInMode "copy"

    script:
    """
    set -eux
    export uvfits="\$(ls -1 *.uvfits | head -1)"
    # convert uvfits to ms
    ${params.casa} -c "importuvfits('\${uvfits}', 'vis.ms')"
    # fix mwa ms
    singularity exec --bind \$PWD:/tmp --writable-tmpfs --pwd /tmp --no-home ${params.cotter_sif} \
        fixmwams vis.ms ${obsid}.metafits
    # imaging
    ${params.wsclean} \
        -weight briggs -1.0 \
        -name wsclean_${name} \
        -size ${params.img_size} ${params.img_size} \
        -scale ${params.img_scale} \
        -pol xx,yy,v \
        -abs-mem ${params.img_mem} \
        -channels-out ${params.img_channels_out} \
        vis.ms
    """
}

// create source list from deconvolved image of stokes I
// process wscleanSources {
//     maxRetries 1
//     // - shallow deconvolution, 1000 iter just enough for removing sidelobes
//     // - channels: need mfs
//     storeDir "$params.outdir/${obsid}/img${params.img_suffix}"
//     input:
//     tuple val(obsid), val(name), path("*") // *.ms.zip

//     output:
//     tuple val(obsid), val(name), path("wscleanSources_${name}*-I-image.fits"), path("wscleanSources_${name}*-sources.txt")

//     // label jobs that need a bigger cpu allocation
//     label "cpu"
//     label "img"

//     script:
//     """
//     unzip *.ms.zip
//     export name="\$(basename *.ms)"
//     name="\${name%.*}"
//     # todo: idg?
//     ${params.wsclean} \
//         -mgain 0.95 -weight briggs -1.0 -multiscale \
//         -name wscleanSources_\${name} \
//         -size ${params.img_size} ${params.img_size} \
//         -scale ${params.img_scale} \
//         -niter ${params.img_niter} \
//         -pol i \
//         -auto-threshold ${params.img_auto_threshold} \
//         -auto-mask ${params.img_auto_mask} \
//         -abs-mem ${params.img_mem} \
//         -channels-out ${params.img_channels_out} \
//         -save-source-list \
//         *.ms
//     """
// }

// power spectrum metrics via chips
process psMetrics {
    input:
    tuple val(obsid), val(name), path("*") // <- "*.uvfits"
    output:
    tuple val(obsid), path("output_metrics_${name}.dat"), path("${name}.log")

    storeDir "$params.outdir/${obsid}/ps_metrics"

    errorStrategy 'ignore'
    stageInMode "copy"

    // without this, nextflow checks for the output before fs has had time
    // to sync, causing it to think the job has failed.
    beforeScript 'sleep 5s; ls -al'
    afterScript 'ls -al; sleep 5s'

    module "cfitsio/3.470:gcc-rt/9.2.0"

    script:
    """
    set -eux
    export DATADIR="\$PWD"
    export OUTPUTDIR="\$PWD/"
    export BAND=0
    export NCHAN=384
    export UVFITS="\$(ls -1 *.uvfits | head -1)"
    export UVFITS="\$(basename \${UVFITS})"
    ${params.ps_metrics} "${name}" "\${BAND}" "\${NCHAN}" "\${UVFITS%.*}" 2>&1 | tee "${name}.log"
    ls -al "\${OUTPUTDIR}"
    cat "output_metrics_${name}.dat"
    """
}

// takes dirty V and clean XX, YY for a single vis
process imgQA {
    input:
    tuple val(obsid), val(name), path("*") // <- "*-dirty.fits"
    output:
    tuple val(obsid), path("${name}.json")

    storeDir "$params.outdir/${obsid}/img_qa"

    // TODO: figure out why this keeps failing
    errorStrategy 'ignore'

    stageInMode "copy"
    module "singularity"
    maxRetries 1

    // without this, nextflow checks for the output before fs has had time
    // to sync, causing it to think the job has failed.
    beforeScript 'sleep 1s; ls -al'
    afterScript 'ls -al; sleep 1s'

    script:
    """
    set -ex
    singularity exec --bind \$PWD:/tmp --writable-tmpfs --pwd /tmp --no-home ${params.mwa_qa_sif} \
        run_imgqa.py *.fits --out ${name}.json
    """
}

import groovy.json.JsonSlurper
def jslurp = new JsonSlurper()

// display a long list of ints
def displayRange = { s, e -> s == e ? "${s}," : s == e - 1 ? "${s},${e}," : "${s}-${e}," }
def displayInts = { l ->
    def sb, start, end
    (sb, start, end) = [''<<'', l[0], l[0]]
    for (i in l[1..-1]) {
        (sb, start, end) = i == end + 1 ? [sb, start, i] : [sb << displayRange(start, end), i, i]
    }
    (sb << displayRange(start, end))[0..-2].toString()
}

workflow {
    // get obsids from csv
    obsids = channel.fromPath(params.obsids_path).splitCsv().flatten()

    // ensure raw files are downloaded
    ensureRaw(obsids)

    // collect disk usage stats from ensureRaw stage
    ensureRaw.out
        // form row of tsv from obsid, metafits size, gpufits count, gpufits size
        .map { it -> [
            it[0], // obsid
            it[1].size(), // size of metafits [B]
            it[2].size(), // number of gpufits files
            it[2].collect(path -> path.size()).sum() // total gpufits size [B]
        ].join("\t") }
        .collectFile(
            name: "raw_stats.tsv", newLine: true, sort: true,
            seed: ["OBS","META SIZE","N GPUBOX","GPUBOX SIZE"].join("\t"),
            storeDir: "${projectDir}/results/"
        )
        | view

    // analyse metafits stats
    ensureRaw.out
        .map { it -> it[0..1] }
        | metafitsStats

    // collect metafits stats
    metafitsStats.out
        // for row of tsv from metafits json fields we care about
        .map { it ->
            def stats = jslurp.parse(it[1])
            def flagged_inps = stats.get("FLAGGED_INPUTS",[])
            [
                it[0],
                stats.get("DATE-OBS",""),
                stats.get("GRIDNUM",""),
                stats.get("CENTCHAN",""),
                stats.get("FINECHAN",""),
                stats.get("INTTIME",""),
                stats.get("NSCANS",""),
                stats.get("NINPUTS",""),
                flagged_inps.size(), // number of flagged inputs
                flagged_inps.join('\t') // flagged input names
            ].join("\t")
        }
        .collectFile(
            name: "metafits_stats.tsv", newLine: true, sort: true,
            seed: [
                "OBS","DATE","POINT","CENT CH","FREQ RES","TIME RES","N SCANS","N FLAG INPS","FLAG INPS"
            ].join("\t"),
            storeDir: "${projectDir}/results/",
        )
        | view

    // ensure preprocessed files have been generated from raw files
    ensurePrep(ensureRaw.out)

    // get shape of preprocessed visibilities
    ensurePrep.out
        .map { it -> it[0..1] }
        | prepStats

    // collect prep stats
    prepStats.out
        // join with uvfits, mwaf and birli log
        .join(ensurePrep.out)
        // form row of tsv from json fields we care about
        .map { it ->
            def stats = jslurp.parse(it[1])
            def missing_hdus = (
                it[4].getText() =~ /NoDataForTimeStepCoarseChannel \{ timestep_index: (\d+), coarse_chan_index: (\d+) \}/
            ).collect { m -> [m[1], m[2]] }
            // display a compressed version of missing hdus
            def missing_timesteps_by_gpubox = missing_hdus
                .groupBy { m -> m[1] }
                .collect { g ->
                    [g.key, displayInts(g.value.collect {p -> p[0] as int})].join(":")
                }
                .join("|");
            [
                it[0],
                // prep stats json
                stats.get('num_chans'),
                stats.get('num_times'),
                stats.get('num_baselines'),
                // disk usage
                it[2].size(), // size of uvfits [B]
                it[3].size(), // number of mwaf [B]
                it[3].collect(path -> path.size()).sum(), // total mwaf size [B]
                // parse birli log
                missing_hdus.size(),
                missing_timesteps_by_gpubox,
            ].join("\t")
        }
        .collectFile(
            name: "prep_stats.tsv", newLine: true, sort: true,
            seed: [
                "OBS","N CHAN","N TIME","N BL","UVFITS BYTES", "N MWAF","MWAF BYTES","MISSING HDUs"
            ].join("\t"),
            storeDir: "${projectDir}/results/"
        )
        | view

    // flag QA
    ensureRaw.out
        // get the obsid, metafits from ensureRaw
        .map { it -> it[0..1] }
        // join with the mwaf from ensurePrep
        .join(ensurePrep.out).map { it -> it[0..1] + [it[3]] }
        | flagQA

    // export flagQA results
    def sky_chans = (131..154).collect { ch -> "$ch".toString() }
    flagQA.out
        // form row of tsv from json fields we care about
        .map { it ->
            def stats = jslurp.parse(it[1])
            def chans = stats.get('channels',[:])
            def chan_occupancy = sky_chans
                .collect { ch -> chans.get(ch,[:]).get('occupancy', '') }
            ([
                it[0],
                stats.get('total_occupancy'),
            ] + chan_occupancy).join("\t")
        }
        .collectFile(
            name: "occupancy.tsv", newLine: true, sort: true,
            seed: (["OBS", "TOTAL OCCUPANCY"] + sky_chans).join("\t"),
            storeDir: "${projectDir}/results/"
        )
        | view

    // create a channel of a single csv for calibration args with columns:
    // - short name that appears in calibration solution filename
    // - args for hyperdrive soln apply args
    dical_args = channel
        .fromPath(params.dical_args_path)
        .collect()

    // create a channel of calibration parameter sets from a csv with columns:
    // - short name that appears in calibration solution filename
    // - short name that appears in the apply filename
    // - args for hyperdrive soln apply args
    apply_args = channel
        .fromPath(params.apply_args_path)
        .splitCsv(header: false)

    // calibration solutions for each obsid that passes QA
    ensureRaw.out
        // get obsid, metafits from ensureRaw
        .map { it -> it[0..1] }
        // join with uvfits from ensurePrep
        .join(ensurePrep.out).map { it -> it[0..2] }
        // filter out obsids which exceed flag occupancy threshold
        .join(
            flagQA.out
            .map { it -> [it[0], jslurp.parse(it[1]).get('total_occupancy')] }
            .filter(it -> it[1] < params.flag_occupancy_threshold)
            .map { it -> [it[0]] }
        )
        // combine with hyperdrive di-cal args file
        .combine(dical_args)
        | ensureCalSol

    // calibration QA
    ensureRaw.out
        // get obsid, metafits from ensureRaw
        .map { it -> it[0..1] }
        // join with hyp_soln from ensureCal
        .join(ensureCalSol.out).map( it -> it[0..2])
        // transpose to run once for each calibration solution
        .transpose()
        .map(it -> [
            it[0],
            // give each calibration a name from basename of uvfits
            it[2].getBaseName(),
            it[1],
            it[2]
        ])
        | (calQA & plotSolutions)

    // collect calQA results as .tsv
    calQA.out
        // separately process each cal solution
        .transpose()
        // form row of tsv from json fields we care about
        .map { it ->
            // give each cal a name from basename of metrics file
            def cal_name = it[1].getBaseName().split('_')[3..-1].join('_')
            // parse json
            // def stats = jslurp.parse(it[1])
            // TODO: fix nasty hack to deal with NaNs
            def stats = jslurp.parseText(it[1].getText().replace("NaN", '"NaN"'))
            [
                it[0],
                cal_name,
                stats.get("FLAGGED_BLS",""),
                stats.get("FLAGGED_CHS",""),
                stats.get("FLAGGED_ANTS",""),
                stats.get("NON_CONVERGED_CHS",""),
                (stats.get("CONVERGENCE_VAR",[])+[""])[0],
                stats.get("SKEWNESS_UCVUT", "")
            ].join("\t")
        }
        .collectFile(
            name: "cal_metrics.tsv", newLine: true, sort: true,
            seed: [
                "OBS", "CAL NAME","FLAG BLS", "FLAG CHS", "FLAG ANTS", "NON CONVG CHS", "CONVG VAR", "SKEW"
            ].join("\t"),
            storeDir: "${projectDir}/results/"
        ) \
        | view

    // apply calibration solutions
    ensureRaw.out
        // get obsid, metafits from ensureRaw
        .map { it -> it[0..1] }
        // join with uvfits from ensurePrep
        .join(ensurePrep.out).map { it -> it[0..2] }
        // join with hyp_solns from ensureCalSol
        .join(ensureCalSol.out).map { it -> it[0..3] }
        // combine with hyperdrive apply solution args file
        .combine(apply_args)
        | ensureCalVis

    // vis QA and ps_metrics
    ensureCalVis.out
        // get obsid, cal uvfits from ensureCalVis
        .map(it -> [
            it[0],
            // give each vis a name from basename of uvfits
            it[1].getBaseName(),
            it[1]
        ])
        | (visQA & psMetrics)

    // collect visQA results as .tsv
    visQA.out
        // form row of tsv from json fields we care about
        .map { it ->
            // give each vis a name from basename of metrics file
            def vis_name = it[1].getBaseName().split('_')[3..-3].join('_')
            // parse json
            // def stats = jslurp.parse(it[1])
            // TODO: fix nasty hack to deal with NaNs
            def stats = jslurp.parseText(it[1].getText().replace("NaN", '"NaN"'))
            def autos = stats.get("AUTOS",[:])
            [
                it[0],
                // give each image a name from basename of uvfits
                vis_name,
                // poor times can sometimes be [], this gives us a default
                (autos.get("XX",[:]).get("POOR_TIMES",[])+[''])[0],
                (autos.get("YY",[:]).get("POOR_TIMES",[])+[''])[0]
            ].join("\t")
        }
        .collectFile(
            name: "vis_metrics.tsv", newLine: true, sort: true,
            seed: ["OBS", "VIS NAME", "XX POOR TIMES", "YY POOR TIMES"].join("\t"),
            storeDir: "${projectDir}/results/"
        )
        | view

    // wsclean: make dirty images
    ensureRaw.out \
        // get obsid, metafits from ensureRaw
        .map { it -> it[0..1] }
        // with uvfits from ensureCalVis
        .join(ensureCalVis.out)
        .map(it -> [
            it[0],
            // give each image a name from basename of uvfits
            it[2].getBaseName(),
            it[1],
            it[2]
        ])
        | wscleanDirty

    // imgQA for all groups of images
    wscleanDirty.out
        // group images by the basename before "-${POL}-dirty.fits"
        .flatMap { it -> it[2]
            .collect { path -> [(path.getName() =~ /^(.*)-\w+-dirty.fits$/)[0][1], path] }
            .groupBy { tup -> tup[0] }
            .collect { entry -> [it[0], entry.key, entry.value.collect(item -> item[1])] }
        }
        | imgQA

    // collect imgQA results as .tsv
    imgQA.out
        // form row of tsv from json fields we care about
        .map { it ->
            // give each vis a name from basename of metrics file
            def img_name = it[1].getBaseName().split('_')[3..-1].join('_')
            // give each image a name from basename of metrics json
            // parse json
            // def stats = jslurp.parse(it[1])
            // TODO: fix nasty hack to deal with NaNs
            def stats = jslurp.parseText(it[1].getText().replace("NaN", '"NaN"'))
            def xx = stats.get("XX",[:])
            def yy = stats.get("YY",[:])
            def v = stats.get("V",[:])
            def v_xx = stats.get("V_XX",[:])
            [
                it[0],
                img_name,
                xx.get("RMS_ALL"),
                xx.get("RMS_BOX"),
                xx.get("PKS0023_026"),
                yy.get("RMS_ALL"),
                yy.get("RMS_BOX"),
                yy.get("PKS0023_026"),
                v.get("RMS_ALL"),
                v.get("RMS_BOX"),
                v.get("PKS0023_026"),
                v_xx.get("RMS_RATIO_ALL"),
                v_xx.get("RMS_RATIO_BOX")
            ].join("\t")
        }
        .collectFile(
            name: "img_metrics.tsv", newLine: true, sort: true,
            seed: [
                "OBS", "IMG NAME", "XX ALL", "XX BOX", "XX PKS", "YY ALL", "YY BOX", "YY PKS",
                "V ALL","V BOX","V PKS", "V:XX ALL", "V:XX BOX"
            ].join("\t"),
            storeDir: "${projectDir}/results/"
        )
        | view


    // collect psMetrics as a .dat
    psMetrics.out
        // read the content of each ps_metrics file including the trailing newline
        .map(it -> it[1].getText())
        .collectFile(
            name: "ps_metrics.dat",
            storeDir: "${projectDir}/results/"
        )
        | view

    // collect psMetrics as a .tsv
    psMetrics.out
        .map(it -> [it[0]] + it[1].getText().split('\n')[0].split(' ')[0..-1])
        // get calName from baseName
        .map(it -> [ it[0], it[1].split('_')[2..-1].join('_') ] + it[2..-1])
        // form each row of tsv
        .map(it -> it.join("\t"))
        .collectFile(
            name: "ps_metrics.tsv", newLine: true, sort: true,
            seed: [
                "OBS", "CAL NAME", "P_WEDGE", "NUM_CELLS", "P_WINDOW", "NUM_CELLS",
                "P_ALL", "D3"
            ].join("\t"),
            storeDir: "${projectDir}/results/"
        )
        | view
}