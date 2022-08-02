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
    errorStrategy {
        failure_reason = [
            5: "I/O error or hash mitch",
            28: "No space left on device",
        ][task.exitStatus]
        if (failure_reason) {
            println "task ${task.hash} failed with code ${task.exitStatus}: ${failure_reason}"
            return 'ignore'
        }
        retry_reason = [
            1: "general or permission",
            11: "Resource temporarily unavailable"
        ][task.exitStatus] ?: "unknown"
        println "retrying task ${task.hash} failed with code ${task.exitStatus}: ${retry_reason}"
        sleep(Math.pow(2, task.attempt) * 60*60*1000 as long)
        return 'retry'
    }

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
                exit 28  # No space left on device
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
        [ -z "\$size" ] && return 1
        ensure_disk_space \$size || exit \$?
        # giant-squid download --keep-zip --hash -v \$jobid 2>&1 | tee download.log
        # if [ \${PIPESTATUS[0]} -ne 0 ]; then
        #     echo "Download failed, see download.log"
        #     exit 1
        # fi
        wget \$url -O \$jobid.tar --progress=dot:giga --wait=60 --random-wait
        export wget_status=\$?
        if [ \$wget_status -ne 0 ]; then
            echo "Download failed. status=\$wget_status"
            exit \$wget_status
        fi
        sha1=\$(sha1sum \$jobid.tar | cut -d' ' -f1)
        if [ "\$sha1" != "\$hash" ]; then
            echo "Download failed, hash mismatch"
            exit 5  # Input/output error
        fi
        [ -d "${task.storeDir}" ] || mkdir -p "${task.storeDir}"
        tar -xf \$jobid.tar -C "${task.storeDir}"
        return \$?
    }

    # submit a job to ASVO, suppress failure if a job already exists.
    ${params.giant_squid} submit-vis --delivery "acacia" $obsid -w || true

    # download any ready jobs, exit if success, else suppress warning
    get_first_ready_job && exit 0 || true

    # extract id and state from pending download vis jobs
    ${params.giant_squid} list -j --types download_visibilities --states queued,processing -- $obsid \
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
process birliPrep {
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
process hypCalSol {
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

// ensure calibrated uvfits are present, or apply calsols to prep uvfits with hyperdrive
process hypApplyUV {
    input:
    tuple val(obsid), path("${obsid}.metafits"), path("${obsid}.uvfits"), path("*"), \
        val(cal_name), val(apply_name), val(apply_args)
    output:
    tuple val(obsid), val(vis_name), path("hyp_${obsid}_${vis_name}.uvfits"), \
        path("hyp_apply_${vis_name}.log")

    storeDir "$params.outdir/${obsid}/cal$params.cal_suffix"

    // TODO: figure out why this keeps failing
    errorStrategy 'ignore'

    // label jobs that need a bigger gpu allocation
    label "cpu"

    module "gcc-rt/9.2.0"
    stageInMode "copy"

    script:
    vis_name = "${cal_name}_${apply_name}"
    """
    # hyperdrive solutions apply uvfits
    ${params.hyperdrive} solutions-apply ${apply_args} \
        --data "${obsid}.metafits" "${obsid}.uvfits" \
        --solutions "hyp_soln_${obsid}_${cal_name}.fits" \
        --outputs "hyp_${obsid}_${vis_name}.uvfits" \
        | tee hyp_apply_${vis_name}.log

    # print out important info from log
    grep -iE "err|warn|hyperdrive|chanblocks|reading|writing|flagged" *.log
    """
    // ms:
    // ${params.hyperdrive} solutions-apply \${apply_args} \
    //   --outputs "hyp_${obsid}_\${name}.uvfits" "hyp_${obsid}_\${name}.ms" \
    //   zip -r "hyp_${obsid}_\${name}.ms.zip" "hyp_${obsid}_\${name}.ms"
}

// ensure calibrated ms are present, or apply calsols to prep uvfits with hyperdrive
process hypApplyMS {
    input:
    tuple val(obsid), path("${obsid}.metafits"), path("${obsid}.uvfits"), path("*"), \
        val(cal_name), val(apply_name), val(apply_args)
    output:
    tuple val(obsid), val(vis_name), path("hyp_${obsid}_${vis_name}.ms"), \
        path("hyp_apply_${vis_name}_ms.log")

    storeDir "$params.outdir/${obsid}/cal$params.cal_suffix"
    // storeDir "/data/curtin_mwaeor/FRB_hopper/"

    // TODO: figure out why this keeps failing
    errorStrategy 'ignore'

    // label jobs that need a bigger gpu allocation
    label "cpu"

    module "gcc-rt/9.2.0"
    stageInMode "copy"

    script:
    vis_name = "${cal_name}_${apply_name}"
    """
    # hyperdrive solutions apply ms
    ${params.hyperdrive} solutions-apply ${apply_args} \
        --data "${obsid}.metafits" "${obsid}.uvfits" \
        --solutions "hyp_soln_${obsid}_${cal_name}.fits" \
        --outputs "hyp_${obsid}_${vis_name}.ms" \
        | tee hyp_apply_${vis_name}_ms.log

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
    tuple val(obsid), val(name), path("hyp_${obsid}_${name}.uvfits")
    output:
    tuple val(obsid), val(name), path("hyp_${obsid}_${name}_vis_metrics.json")

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
        run_visqa.py *.uvfits --out "hyp_${obsid}_${name}_vis_metrics.json"
    """
}

// QA tasks that can be run on each calibration parameter set
process calQA {
    input:
    tuple val(obsid), val(name), path("${obsid}.metafits"), path("soln.fits")
    output:
    tuple val(obsid), val(name), path("hyp_soln_${obsid}_${name}_X.json")

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
    singularity exec --bind \$PWD:/tmp --writable-tmpfs --pwd /tmp --no-home ${params.mwa_qa_sif} \
        run_calqa.py soln.fits ${obsid}.metafits X --out "hyp_soln_${obsid}_${name}_X.json"
    """
}

process plotSolutions {
    input:
    tuple val(obsid), val(name), path("${obsid}.metafits"), path("hyp_soln_${obsid}_${name}.fits")
    output:
    tuple val(obsid), path("hyp_soln_${obsid}_${name}_{phases,amps}.png")

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
    tuple val(obsid), val(name), path("vis.ms")
    output:
    tuple val(obsid), val(name), path("wsclean_hyp_${obsid}_${name}-MFS-{XX,YY,V}-dirty.fits")

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
    # imaging
    ${params.wsclean} \
        -weight briggs -1.0 \
        -name wsclean_hyp_${obsid}_${name} \
        -size ${params.img_size} ${params.img_size} \
        -scale ${params.img_scale} \
        -pol xx,yy,v \
        -abs-mem ${params.img_mem} \
        -channels-out ${params.img_channels_out} \
        vis.ms
    """
}

// power spectrum metrics via chips
process psMetrics {
    input:
    tuple val(obsid), val(name), path("hyp_${obsid}_${name}.uvfits")
    output:
    tuple val(obsid), val(name), path("output_metrics_hyp_${obsid}_${name}.dat"), \
        path("hyp_${obsid}_${name}.log")

    storeDir "$params.outdir/${obsid}/ps_metrics"

    errorStrategy 'ignore'
    stageInMode "copy"

    // without this, nextflow checks for the output before fs has had time
    // to sync, causing it to think the job has failed.
    beforeScript 'sleep 5s; ls -al'
    afterScript 'ls -al; sleep 5s'

    module "cfitsio/3.470:gcc-rt/9.2.0"

    script:
    band = 0
    nchan = 384
    uv_base = "hyp_${obsid}_${name}"
    """
    set -eux
    export DATADIR="\$PWD"
    export OUTPUTDIR="\$PWD/"
    ${params.ps_metrics} "${uv_base}" "${band}" "${nchan}" "${uv_base}" 2>&1 \
        | tee "${uv_base}.log"
    """
}

// takes dirty V and clean XX, YY for a single vis
process imgQA {
    input:
    tuple val(obsid), val(name), path("*") // <- "*-dirty.fits"
    output:
    tuple val(obsid), val(name), path("wsclean_hyp_${obsid}_${name}-MFS.json")

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
        run_imgqa.py *.fits --out wsclean_hyp_${obsid}_${name}-MFS.json
    """
}

import groovy.json.JsonSlurper
def jslurp = new JsonSlurper()

// display a long list of ints, replace bursts of consecutive numbers with ranges
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
        // form row of tsv
        .map { def (obsid, metafits, gpufits) = it; [
            obsid,
            metafits.size(), // size of metafits [B]
            gpufits.size(), // number of gpufits files
            gpufits.collect(path -> path.size()).sum() // total gpufits size [B]
        ].join("\t") }
        .collectFile(
            name: "raw_stats.tsv", newLine: true, sort: true,
            // "seed" is the header of the tsv file
            seed: ["OBS","META SIZE","N GPUBOX","GPUBOX SIZE"].join("\t"),
            storeDir: "${projectDir}/results/"
        )
        // display output path and number of lines
        | view { [it, it.readLines().size()] }

    // channel of tuples (obsid, metafits)
    eachMetafits = ensureRaw.out.map { def (obsid, metafits) = it; [obsid, metafits] }

    // analyse metafits stats
    eachMetafits | metafitsStats

    // collect metafits stats
    metafitsStats.out
        // for row of tsv from metafits json fields we care about
        .map { def (obsid, json) = it;
            def stats = jslurp.parse(json)
            [
                obsid,
                stats."DATE-OBS",
                stats.GRIDNUM,
                stats.CENTCHAN,
                stats.FINECHAN,
                stats.INTTIME,
                stats.NSCANS,
                stats.NINPUTS,
                stats.FLAGGED_INPUTS?.size(), // number of flagged inputs
                stats.FLAGGED_INPUTS?.join('\t') // flagged input names
            ].join("\t")
        }
        .collectFile(
            name: "metafits_stats.tsv", newLine: true, sort: true,
            seed: [
                "OBS","DATE","POINT","CENT CH","FREQ RES","TIME RES","N SCANS", "N INPS",
                "N FLAG INPS","FLAG INPS"
            ].join("\t"),
            storeDir: "${projectDir}/results/",
        )
        // display output path and number of lines
        | view { [it, it.readLines().size()] }


    // ensure preprocessed files have been generated from raw files
    birliPrep(ensureRaw.out)

    // channel of tuples (obsid, prep uvfits)
    eachPrep = birliPrep.out.map { def (obsid, prep) = it; [obsid, prep] }

    // get shape of preprocessed visibilities
    eachPrep | prepStats

    // collect prep stats
    prepStats.out
        // join with uvfits, mwaf and birli log
        .join(birliPrep.out)
        // form row of tsv from json fields we care about
        .map { def (obsid, json, prepVis, prepMwafs, prepLog) = it;
            def stats = jslurp.parse(json)
            def missing_hdus = (
                prepLog.getText() =~ /NoDataForTimeStepCoarseChannel \{ timestep_index: (\d+), coarse_chan_index: (\d+) \}/
            ).collect { m -> [ts:m[1], cc:m[2]] }
            // display a compressed version of missing hdus
            def missing_timesteps_by_gpubox = missing_hdus
                .groupBy { m -> m.cc }
                .collect { cc, ts ->
                    [cc, displayInts(ts.collect {p -> p.ts as int})].join(":")
                }
                .join("|");
            [
                obsid,
                // prep stats json
                stats.num_chans,
                stats.num_times,
                stats.num_baselines,
                // disk usage
                prepVis.size(), // size of uvfits [B]
                prepMwafs.size(), // number of mwaf
                prepMwafs.collect(path -> path.size()).sum(), // total mwaf size [B]
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
        // display output path and number of lines
        | view { [it, it.readLines().size()] }

    // flag QA
    eachMetafits
        // join with the mwafs from birliPrep
        .join(birliPrep.out.map { def (obsid, _, prepMwafs) = it; [obsid, prepMwafs] })
        | flagQA

    // export flagQA results
    def sky_chans = (131..154).collect { ch -> "$ch".toString() }
    flagQA.out
        // form row of tsv from json fields we care about
        .map { def (obsid, json) = it;
            def stats = jslurp.parse(json)
            def chan_occupancy = sky_chans.collect { ch -> stats.channels?[ch]?.occupancy }
            ([obsid, stats.total_occupancy] + chan_occupancy).join("\t")
        }
        .collectFile(
            name: "occupancy.tsv", newLine: true, sort: true,
            seed: (["OBS", "TOTAL OCCUPANCY"] + sky_chans).join("\t"),
            storeDir: "${projectDir}/results/"
        )
        // display output path and number of lines
        | view { [it, it.readLines().size()] }

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

    // calibration solutions for each obsid that passes flagQA:
    // - do not calibrate obs where occupancy > params.flag_occupancy_threshold
    eachMetafits
        // join with uvfits from birliPrep
        .join(eachPrep)
        // filter out obsids which exceed flag occupancy threshold
        .join(flagQA.out
            .map { def (obsid, json) = it; [obsid, jslurp.parse(json).total_occupancy] }
            .filter { def (_, occ) = it; occ && occ < params.flag_occupancy_threshold }
            .map { def (obsid, _) = it; obsid }
        )
        // combine with hyperdrive di-cal args file
        .combine(dical_args)
        | hypCalSol

    // channel of tuples (obsid, calName, metafits, calSol)
    // hypCalSol gives multiple solutions, transpose gives 1 tuple per solution.
    eachCal = eachMetafits
        // join with hyp_soln from ensureCal
        .join(hypCalSol.out.map { def (obsid, solutions) = it; [obsid, solutions] })
        // transpose to run once for each calibration solution
        .transpose()
        .map { def (obsid, metafits, soln) = it
            // give each calibration a name from basename of solution fits
            def name = soln.getBaseName().split('_')[3..-1].join('_')
            [obsid, name, metafits, soln]
        }

    // calibration QA and plot solutions
    eachCal | (plotSolutions & calQA)

    // collect calQA results as .tsv
    calQA.out
        // form row of tsv from json fields we care about
        .map { def (obsid, name, json) = it;
            // parse json
            // def stats = jslurp.parse(json)
            // TODO: fix nasty hack to deal with NaNs
            def stats = jslurp.parseText(json.getText().replace("NaN", '"NaN"'))
            [
                obsid,
                name,
                stats.FLAGGED_BLS,
                stats.FLAGGED_CHS,
                stats.FLAGGED_ANTS,
                stats.NON_CONVERGED_CHS,
                stats.CONVERGENCE_VAR?[0],
                stats.SKEWNESS_UCVUT
            ].join("\t")
        }
        .collectFile(
            name: "cal_metrics.tsv", newLine: true, sort: true,
            seed: [
                "OBS", "CAL NAME", "FLAG BLS", "FLAG CHS", "FLAG ANTS",
                "NON CONVG CHS", "CONVG VAR", "SKEW"
            ].join("\t"),
            storeDir: "${projectDir}/results/"
        )
        // display output path and number of lines
        | view { [it, it.readLines().size()] }

    // apply calibration solutions
    eachMetafits
        // join with uvfits from birliPrep
        .join(eachPrep)
        // join with hyp_solns from hypCalSol
        .join(hypCalSol.out.map { def (obsid, solutions) = it; [obsid, solutions] })
        // combine with hyperdrive apply solution args file
        .combine(apply_args)
        | (hypApplyUV & hypApplyMS)

    // vis QA and ps_metrics
    hypApplyUV.out
        // get obsid, cal uvfits from hypApplyUV
        .map { def (obsid, name, vis) = it; [obsid, name, vis] }
        | (visQA & psMetrics)

    // collect visQA results as .tsv
    visQA.out
        // form row of tsv from json fields we care about
        .map { def (obsid, name, json) = it;
            // parse json
            // def stats = jslurp.parse(it[1])
            // TODO: fix nasty hack to deal with NaNs
            def stats = jslurp.parseText(json.getText().replace("NaN", '"NaN"'))
            view(stats)
            [
                obsid,
                name,
                stats.AUTOS?.XX?.POOR_TIMES?[0],
                stats.AUTOS?.YY?.POOR_TIMES?[0]
            ].join("\t")
        }
        .collectFile(
            name: "vis_metrics.tsv", newLine: true, sort: true,
            seed: ["OBS", "VIS NAME", "XX POOR TIMES", "YY POOR TIMES"].join("\t"),
            storeDir: "${projectDir}/results/"
        )
        // display output path and number of lines
        | view { [it, it.readLines().size()] }

    // wsclean: make dirty images
    // - get ms from hypApplyMS
    hypApplyMS.out
        .map { def (obsid, name, vis) = it; [obsid, name, vis] }
        | wscleanDirty

    // imgQA for all groups of images
    wscleanDirty.out
        // if doing multiple chans:
        // group dirty images by chan, where name is "wsclean_hyp_${obsid}_${name}-${chan}-${pol}-dirty.fits"
        // .flatMap { def (obsid, name, imgs) = it;
        //     imgs
        //         .groupBy { img -> (img.getBaseName().split(name)[1].split("-")[0] }
        //         .collect { chan, img -> [obsid, chan, img]}
        // }
        | imgQA

    // collect imgQA results as .tsv
    imgQA.out
        // form row of tsv from json fields we care about
        .map { def (obsid, name, json) = it;
            // parse json
            // def stats = jslurp.parse(json)
            // TODO: fix nasty hack to deal with NaNs
            def stats = jslurp.parseText(json.getText().replace("NaN", '"NaN"'))
            [
                obsid,
                name,
                stats.XX?.RMS_ALL, stats.XX?.RMS_BOX, stats.XX?.PKS0023_026,
                stats.YY?.RMS_ALL, stats.YY?.RMS_BOX, stats.YY?.PKS0023_026,
                stats.V?.RMS_ALL, stats.V?.RMS_BOX, stats.V?.PKS0023_026,
                stats.V_XX?.RMS_RATIO_ALL, stats.V_XX?.RMS_RATIO_BOX
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
        // display output path and number of lines
        | view { [it, it.readLines().size()] }


    // collect psMetrics as a .dat
    psMetrics.out
        // read the content of each ps_metrics file including the trailing newline
        .map { def (obsid, vis_name, dat) = it; dat.getText() }
        .collectFile(
            name: "ps_metrics.dat",
            storeDir: "${projectDir}/results/"
        )
        // display output path and number of lines
        | view { [it, it.readLines().size()] }

    // collect psMetrics as a .tsv
    psMetrics.out
        // form each row of tsv
        .map { def (obsid, vis_name, dat) = it;
            def dat_values = dat.getText().split('\n')[0].split(' ')[1..-1]
            ([obsid, vis_name] + dat_values).join("\t")
        }
        .collectFile(
            name: "ps_metrics.tsv", newLine: true, sort: true,
            seed: [
                "OBS", "CAL NAME", "P_WEDGE", "NUM_CELLS", "P_WINDOW", "NUM_CELLS",
                "P_ALL", "D3"
            ].join("\t"),
            storeDir: "${projectDir}/results/"
        )
        // display output path and number of lines
        | view { [it, it.readLines().size()] }
}