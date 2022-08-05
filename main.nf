#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Ensure the raw visibility files are present, or download via ASVO
process asvoRaw {
    input:
    val obsid
    output:
    tuple val(obsid), path("${obsid}.metafits"), path("${obsid}_2*.fits")

    // persist results in outdir, process will be skipped if files already present.
    storeDir "${params.outdir}/$obsid/raw"
    // TODO: can't fit in scratch because this writes both tar and contents to disk
    scratch false
    // tag to identify job in squeue and nf logs
    tag "${obsid}"

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
        wait_hours = Math.pow(2, task.attempt)
        println "sleeping for ${wait_hours} hours and retrying task ${task.hash}, which failed with code ${task.exitStatus}: ${retry_reason}"
        sleep(wait_hours * 60*60*1000 as long)
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
        if read -r avail < <(df --output=avail . | tail -n 1); then
            avail_bytes=\$((\$avail * 1000))
            if [[ \$avail_bytes -lt \$needed ]]; then
                echo "Not enough disk space available in \$PWD, need \$needed B, have \$avail_bytes B"
                exit 28  # No space left on device
            fi
            return 0
        fi
        echo "Could not determine disk space in \$PWD"
        exit 1
    }

    # download any ASVO jobs for the obsid that are ready
    function maybe_get_first_ready_job {
        # extract id url and size from ready download vis jobs
        ${params.giant_squid} list -j --types download_visibilities --states ready -- $obsid \
            | tee /dev/stderr \
            | ${params.jq} -r '.[]|[.jobId,.files[0].fileUrl//"",.files[0].fileSize//"",.files[0].fileHash//""]|@tsv' \
            | tee ready.tsv
        if read -r jobid url size hash < ready.tsv; then
            ensure_disk_space "\$size" || exit \$?
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
            tar -xf \$jobid.tar
            return \$?
        fi
        echo "no ready jobs"
        return 1
    }

    # submit a job to ASVO, suppress failure if a job already exists.
    ${params.giant_squid} submit-vis --delivery "acacia" $obsid -w || true

    # download any ready jobs, exit if success, else suppress warning
    maybe_get_first_ready_job && exit 0 || true

    # extract id and state from pending download vis jobs
    ${params.giant_squid} list -j --types download_visibilities --states queued,processing -- $obsid \
        | ${params.jq} -r '.[]|[.jobId,.jobState]|@tsv' \
        | tee pending.tsv

    # if pending jobs, wait for the first to complete and try downloading again
    if read -r jobid state < pending.tsv; then
        echo "[obs:$obsid]: waiting until \$jobid with state \$state is Ready"
        ${params.giant_squid} wait -j -- \$jobid
        # a new job is ready, so try and get it
        maybe_get_first_ready_job
        exit \$?
    fi

    # show errors if no pending jobs
    ${params.giant_squid} list -j --types download_visibilities --states error -- $obsid \
    | ${params.jq} -r '.[]|[.jobId,.jobState]|@tsv' \
    | tee errors.tsv
    """
}

process metaJson {
    input:
    tuple val(obsid), path("${obsid}.metafits")
    output:
    tuple val(obsid), path("${obsid}.metafits.json")

    scratch false
    storeDir "${params.outdir}/${obsid}/raw"

    tag "${obsid}"

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
        for key in hdus['PRIMARY'].header:
            if type(hdus['PRIMARY'].header[key]) not in [bool, int, str, float]:
                continue
            data[key] = hdus['PRIMARY'].header[key]
            if key in [ 'RECVRS', 'DELAYS', 'CHANNELS', 'CHANSEL' ]:
                data[key] = [*filter(None, map(lambda t: t.strip(), data[key].split(',')))]

        data['FLAGGED_INPUTS'] = [
            r['TileName'] + r['Pol']
            for r in hdus['TILEDATA'].data
            if r['Flag'] == 1
        ]

        data['DELAYS'] = hdus['TILEDATA'].data['Delays'].tolist()

        json_dump(data, out, indent=4)
    """
}

// Ensure preprocessed observation is present, or preprocess raw with Birli
process birliPrep {
    input:
    tuple val(obsid), val(spw), path("${obsid}.metafits"), path("*") // <-preserves names of fits files
    output:
    tuple val(obsid), val(spw), path("birli_${obsid}${spw}.uvfits"), path("${obsid}${spw}*.mwaf"), path("birli_prep.log")

    storeDir "${params.outdir}/${obsid}/prep"

    tag "${obsid}${spw}"
    // label jobs that need a bigger cpu allocation
    label "cpu"

    // this is necessary for singularity
    stageInMode "copy"
    module "singularity"

    script:
    flag_template = obsid as int > 1300000000 ? "${obsid}${spw}_ch%%%.mwaf" : "${obsid}${spw}_%%.mwaf"
    """
    ${params.birli} \
        -u "birli_${obsid}${spw}.uvfits" \
        -f "${flag_template}" \
        -m "${obsid}.metafits" \
        --avg-time-res ${params.prep_time_res_s} \
        --avg-freq-res ${params.prep_freq_res_khz} \
        ${obsid}*.fits 2>&1 | tee birli_prep.log
    """
}

// Simple check for vis, to catch "past end of file" errors
process prepStats {
    input:
    tuple val(obsid), path("birli_${obsid}.uvfits")
    output:
    tuple val(obsid), path("${obsid}_prep_stats.json")

    scratch false
    storeDir "${params.outdir}/${obsid}/prep"

    tag "${obsid}"

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

    storeDir "${params.outdir}/${obsid}/prep"

    tag "${obsid}"

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
    tuple val(obsid), val(spw), path("${obsid}.metafits"), path("${obsid}${spw}.uvfits"), path("dical_args.csv")
    output:
    tuple val(obsid), val(spw), path("hyp_soln_${obsid}${spw}_??l_src4k.fits"), path("hyp_di-cal_${obsid}${spw}_??l_src4k.log")

    storeDir "${params.outdir}/${obsid}/cal${params.cal_suffix}"

    tag "${obsid}${spw}"
    // label jobs that need a bigger gpu allocation
    label "gpu"

    errorStrategy 'ignore'

    module "cuda/11.3.1:gcc-rt/9.2.0"
    stageInMode "copy"

    script:
    """
    export HYPERDRIVE_CUDA_COMPUTE=${params.cuda_compute}
    # use this control which GPU we're sending the job to
    export CUDA_VISIBLE_DEVICES=0

    set -ex
    ls -al
    export num_gpus="\$(nvidia-smi -L | wc -l)"
    if [ \$num_gpus -eq 0 ]; then
        echo "no gpus found"
        exit 1
    fi
    # for each calibration parameter set in dical_args.csv, run hyperdrive di-cal
    while IFS=, read -r dical_name dical_args; do
        # on first iteration, this does nothing. on nth, wait for all background jobs to finish
        if [[ \$CUDA_VISIBLE_DEVICES -eq 0 ]]; then
            wait \$(jobs -rp)
        fi
        # hyperdrive di-cal backgrounded in a subshell
        (
            ${params.hyperdrive} di-calibrate \${dical_args} \
                --data "${obsid}.metafits" "${obsid}${spw}.uvfits" \
                --beam "${params.beam_path}" \
                --source-list "${params.sourcelist}" \
                --outputs "hyp_soln_${obsid}${spw}_\${dical_name}.fits" \
                > hyp_di-cal_${obsid}${spw}_\${dical_name}.log
        ) &
        # increment the target device mod num_gpus
        CUDA_VISIBLE_DEVICES=\$(( CUDA_VISIBLE_DEVICES+1 % \${num_gpus} ))
    done < <(cut -d',' -f1,2 dical_args.csv | sort | uniq)

    # wait for all the background jobs to finish
    wait \$(jobs -rp)

    # print out important info from log
    grep -iE "err|warn|hyperdrive|chanblocks|reading|writing|flagged" *.log
    """
}

// fit polynomial to calibration solution
process polyFit {
    input:
    tuple val(obsid), val(dical_name), path("hyp_soln_${obsid}_${dical_name}.fits")
    output:
    tuple val(obsid), val(name), path("hyp_soln_${obsid}_${name}.fits")

    storeDir "${params.outdir}/${obsid}/cal${params.cal_suffix}"
    tag "${obsid}.${dical_name}"
    errorStrategy 'ignore'

    // without this, nextflow checks for the output before fs has had time
    // to sync, causing it to think the job has failed.
    afterScript "sleep 20s; ls ${task.storeDir}"

    module 'python/3.9.7'

    script:
    name = "poly_${dical_name}"
    """
    #!/usr/bin/env python
    # this is nextflow-ified version of /astro/mwaeor/ctrott/poly_fit.py
    # TODO: this deserves its own version control

    from astropy.io import fits
    import numpy as np
    from numpy.polynomial import Polynomial


    def get_unflagged_indices(n_bands, n_chan, clip_width=0):
        # Returns indices of unflagged frequency channels
        edge_flags = 2 + clip_width
        centre_flag = n_chan // 2

        channels = np.arange(n_chan)
        channels = np.delete(channels, centre_flag)
        channels = channels[edge_flags:-edge_flags]

        all_channels = []
        for n in range(n_bands):
            for c in channels:
                all_channels.append(c+(n*n_chan))
        return all_channels


    def fit_hyperdrive_sols(sols, results, order=3, clip_width=0):
        # Fits polynomial to real and imaginary part of hyperdrive solutions
        # Remove flagged channels
        # Hyperdrive solutions have dimensions ((time),ant,freq,pols)

        n_bands = 24
        n_ants, n_freqs, n_pols = np.shape(sols)
        assert n_freqs % n_bands == 0
        n_chan = n_freqs // n_bands

        models_out = np.zeros_like(sols, dtype=complex)

        clipped_x = get_unflagged_indices(n_bands, n_chan, clip_width=clip_width)
        # filter any channels where result is NaN
        clipped_x = [ x for x in clipped_x if not np.isnan(results[x]) ]

        # Remove flagged tiles which are nan at first unflagged frequency and pol
        good_tiles = np.argwhere(~np.isnan(sols[:, clipped_x[0], 0]))

        freq_array = np.arange(n_freqs)

        for ant in good_tiles:
            for pol in range(n_pols):
                z_r = Polynomial.fit(clipped_x, np.real(
                    sols[ant, clipped_x, pol]), deg=order)
                z_i = Polynomial.fit(clipped_x, np.imag(
                    sols[ant, clipped_x, pol]), deg=order)
                models_out[ant, :, pol] = z_r(freq_array) + z_i(freq_array) * 1j

        return models_out

    infile = "hyp_soln_${obsid}_${dical_name}.fits"
    outfile = "hyp_soln_${obsid}_${name}.fits"

    with fits.open(infile) as hdus:
        # Reads hyperdrive solutions from fits file into numpy array
        data = hdus["SOLUTIONS"].data
        for i_timeblock in range(data.shape[0]):
            sols = data[i_timeblock, :, :, ::2] + data[i_timeblock, :, :, 1::2] * 1j
            results = hdus["RESULTS"].data[i_timeblock,:]

            fit_sols = fit_hyperdrive_sols(sols, results)

            # Write fitted solutions
            hdus["SOLUTIONS"].data[i_timeblock, :, :, ::2] = np.real(fit_sols)
            hdus["SOLUTIONS"].data[i_timeblock, :, :, 1::2] = np.imag(fit_sols)

        hdus.writeto(outfile, overwrite="True")
    """
}

// ensure calibrated uvfits are present, or apply calsols to prep uvfits with hyperdrive
process hypApplyUV {
    input:
    tuple val(obsid), path("${obsid}.metafits"), path("${obsid}.uvfits"), path("*"), \
        val(dical_name), val(apply_name), val(apply_args)
    output:
    tuple val(obsid), val(vis_name), path("hyp_${obsid}_${vis_name}.uvfits"), \
        path("hyp_apply_${vis_name}.log")

    storeDir "${params.outdir}/${obsid}/cal${params.cal_suffix}"

    tag "${obsid}.${dical_name}_${apply_name}"
    label "cpu"

    errorStrategy 'ignore'

    module "gcc-rt/9.2.0"
    stageInMode "copy"

    script:
    vis_name = "${dical_name}_${apply_name}"
    """
    # hyperdrive solutions apply uvfits
    ${params.hyperdrive} solutions-apply ${apply_args} \
        --data "${obsid}.metafits" "${obsid}.uvfits" \
        --solutions "hyp_soln_${obsid}_${dical_name}.fits" \
        --outputs "hyp_${obsid}_${vis_name}.uvfits" \
        | tee hyp_apply_${vis_name}.log

    # print out important info from log
    grep -iE "err|warn|hyperdrive|chanblocks|reading|writing|flagged" *.log
    """
}

// ensure calibrated ms are present, or apply calsols to prep uvfits with hyperdrive
process hypApplyMS {
    input:
    tuple val(obsid), path("${obsid}.metafits"), path("${obsid}.uvfits"), path("*"), \
        val(dical_name), val(apply_name), val(apply_args)
    output:
    tuple val(obsid), val(vis_name), path("hyp_${obsid}_${vis_name}.ms"), \
        path("hyp_apply_${vis_name}_ms.log")

    storeDir "${params.outdir}/${obsid}/cal${params.cal_suffix}"
    // storeDir "/data/curtin_mwaeor/FRB_hopper/"

    tag "${obsid}.${dical_name}_${apply_name}"
    label "cpu"

    errorStrategy 'ignore'

    module "gcc-rt/9.2.0"
    stageInMode "copy"

    script:
    vis_name = "${dical_name}_${apply_name}"
    """
    # hyperdrive solutions apply ms
    ${params.hyperdrive} solutions-apply ${apply_args} \
        --data "${obsid}.metafits" "${obsid}.uvfits" \
        --solutions "hyp_soln_${obsid}_${dical_name}.fits" \
        --outputs "hyp_${obsid}_${vis_name}.ms" \
        | tee hyp_apply_${vis_name}_ms.log

    # print out important info from log
    grep -iE "err|warn|hyperdrive|chanblocks|reading|writing|flagged" *.log
    """
}

// QA tasks that can be run on each visibility file.
process visQA {
    input:
    tuple val(obsid), val(name), path("hyp_${obsid}_${name}.uvfits")
    output:
    tuple val(obsid), val(name), path("hyp_${obsid}_${name}_vis_metrics.json")

    storeDir "${params.outdir}/${obsid}/vis_qa"

    tag "${obsid}.${name}"

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

    storeDir "${params.outdir}/${obsid}/cal_qa"

    tag "${obsid}.${name}"

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

process solJson {
    input:
    tuple val(obsid), val(dical_name), path("hyp_soln_${obsid}_${dical_name}.fits")
    output:
    tuple val(obsid), val(dical_name), path("hyp_soln_${obsid}_${dical_name}.fits.json")

    storeDir "${params.outdir}/${obsid}/cal_qa"
    tag "${obsid}.${dical_name}"
    errorStrategy 'ignore'

    // without this, nextflow checks for the output before fs has had time
    // to sync, causing it to think the job has failed.
    afterScript "sleep 20s; ls ${task.storeDir}"

    module 'python/3.9.7'

    script:
    """
    #!/usr/bin/env python
    # this is nextflow-ified version of /astro/mwaeor/ctrott/poly_fit.py
    # TODO: this deserves its own version control

    from astropy.io import fits
    from json import dump as json_dump
    from collections import OrderedDict
    import numpy as np

    with \
        open('hyp_soln_${obsid}_${dical_name}.fits.json', 'w') as out, \
        fits.open("hyp_soln_${obsid}_${dical_name}.fits") as hdus \
    :
        data = OrderedDict()
        # Reads hyperdrive results from fits file into json
        results = hdus["RESULTS"].data[:,:]
        data["RESULTS"] = hdus["RESULTS"].data[:,:].tolist()
        data["DIPOLE_GAINS"] = hdus['TILES'].data["DipoleGains"].tolist()
        data["TOTAL_CONV"] = np.nanprod(results)
        print(repr(data))
        json_dump(data, out, indent=4)
    """
}

process plotSols {
    input:
    tuple val(obsid), val(name), path("${obsid}.metafits"), path("hyp_soln_${obsid}_${name}.fits")
    output:
    tuple val(obsid), path("hyp_soln_${obsid}*_${name}_{phases,amps}.png")

    storeDir "${params.outdir}/${obsid}/cal_qa"

    tag "${obsid}.${name}"

    errorStrategy 'ignore'

    // without this, nextflow checks for the output before fs has had time
    // to sync, causing it to think the job has failed.
    afterScript "sleep 20s; ls ${task.storeDir}"

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

    storeDir "${params.outdir}/${obsid}/img${params.img_suffix}"

    tag "${obsid}.${name}"
    label "cpu"
    label "img"

    errorStrategy 'ignore'

    // this is necessary for singularity
    stageInMode "copy"

    // without this, nextflow checks for the output before fs has had time
    // to sync, causing it to think the job has failed.
    afterScript "sleep 20s; ls ${task.storeDir}"

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

    storeDir "${params.outdir}/${obsid}/ps_metrics"

    tag "${obsid}.${name}"

    errorStrategy 'ignore'
    stageInMode "copy"

    // without this, nextflow checks for the output before fs has had time
    // to sync, causing it to think the job has failed.
    afterScript "sleep 20s; ls ${task.storeDir}"

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

    storeDir "${params.outdir}/${obsid}/img_qa"

    tag "${obsid}.${name}"

    // TODO: figure out why this keeps failing
    errorStrategy 'ignore'

    stageInMode "copy"
    module "singularity"

    // without this, nextflow checks for the output before fs has had time
    // to sync, causing it to think the job has failed.
    afterScript "sleep 20s; ls ${task.storeDir}"

    script:
    """
    set -ex
    singularity exec --bind \$PWD:/tmp --writable-tmpfs --pwd /tmp --no-home ${params.mwa_qa_sif} \
        run_imgqa.py *.fits --out wsclean_hyp_${obsid}_${name}-MFS.json
    """
}

import groovy.json.JsonSlurper
import groovy.json.StringEscapeUtils
jslurp = new JsonSlurper()
def parseJson(path) {
    // TODO: fix nasty hack to deal with NaNs
    jslurp.parseText(path.getText().replaceAll(/(NaN|-?Infinity)/, '"NaN"'))
}

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

// collapse a list into its contiguous ranges
def contigRanges(l) {
    switch (l) {
        case { it.size() == 0 }: return [];
        case { it.size() == 1 }: return [(l[0]..l[0])]
        default:
            def sb, start, end
            (sb, start, end) = [[], l[0], l[0]]
            for (i in l.tail()) {
                (sb, start, end) = i == end + 1 ? [sb, start, i] : [sb << (start..end), i, i]
            }
            return (sb << (start..end))[0..-1]
    }
}

// ensure raw files are downloaded
workflow raw {
    take:
        // channel of obs ids
        obsids
    main:
        asvoRaw(obsids)

        // collect disk usage stats from asvoRaw stage
        asvoRaw.out
            // form row of tsv
            .map { def (obsid, metafits, gpuboxes) = it; [
                obsid,
                metafits.size(), // size of metafits [B]
                gpuboxes.size(), // number of gpubox files
                gpuboxes.collect(path -> path.size()).sum() // total gpubox size [B]
            ].join("\t") }
            .collectFile(
                name: "raw_stats.tsv", newLine: true, sort: true,
                // "seed" is the header of the tsv file
                seed: ["OBS","META SIZE","N GPUBOX","GPUBOX SIZE"].join("\t"),
                storeDir: "${params.outdir}/results${params.result_suffix}/"
            )
            // display output path and number of lines
            | view { [it, it.readLines().size()] }

        // analyse metafits stats
        asvoRaw.out
            .map { def (obsid, metafits, _) = it; [obsid, metafits] }
            | metaJson

        // collect metafits stats
        metaJson.out
            // for row of tsv from metafits json fields we care about
            .map { def (obsid, json) = it;
                def stats = jslurp.parse(json)
                // def row = [obsid]
                // row += [
                //     "DATE-OBS", "GRIDNUM", "CENTCHAN", "FINECHAN", "INTTIME",
                //     "NSCANS", "NINPUTS"
                // ].collect { stats.get(it) }
                def flagged_inputs = stats.FLAGGED_INPUTS?:[]
                def delays = (stats.DELAYS?:[]).flatten()
                [
                    obsid,
                    stats."DATE-OBS",
                    stats.GRIDNUM,
                    stats.CENTCHAN,
                    stats.FINECHAN,
                    stats.INTTIME,
                    stats.NSCANS,
                    stats.NINPUTS,
                    (delays.count { it == 32 } / delays.size()), // fraction of dead dipole
                    flagged_inputs.size(), // number of flagged inputs
                    flagged_inputs.join('\t') // flagged input names
                ].join("\t")
            }
            .collectFile(
                name: "metafits_stats.tsv", newLine: true, sort: true,
                seed: [
                    "OBS","DATE","POINT","CENT CH","FREQ RES","TIME RES","N SCANS", "N INPS",
                    "DEAD DIPOLE FRAC", "N FLAG INPS","FLAG INPS"
                ].join("\t"),
                storeDir: "${params.outdir}/results${params.result_suffix}/",
            )
            // display output path and number of lines
            | view { [it, it.readLines().size()] }

    emit:
        // channel of all raw files: tuple(obsid, metafits, gpuboxes)
        obsRawFiles = asvoRaw.out
}

// ensure preprocessed uvfits and flags have been generated from raw files
workflow prep {
    take:
        // channel of raw files: tuple(obsid, metafits, gpuboxes)
        obsRawFiles
    main:
        // preprocess raw files with Birli
        obsRawFiles
            // set spw to ""
            .map { def (obsid, metafits, gpuboxes) = it; [obsid, "", metafits, gpuboxes] }
            | birliPrep

        // get info about preprocessing stage
        birliPrep.out
            .map { def (obsid, _, prepUVFits, __, ___) = it; [obsid, prepUVFits] }
            | prepStats

        // collect prep stats
        prepStats.out
            // join with uvfits, mwaf and birli log
            .join( birliPrep.out.map { def (obsid, _, prepUVFits, prepMwafs, prepLog) = it;
                [obsid, prepUVFits, prepMwafs, prepLog]
            })
            // form row of tsv from json fields we care about
            .map { def (obsid, json, prepVis, prepMwafs, prepLog) = it;
                def stats = jslurp.parse(json)
                def missing_hdus = (
                    prepLog.getText() =~ /NoDataForTimeStepCoarseChannel \{ timestep_index: (\d+), coarse_chan_index: (\d+) \}/
                ).collect { m -> [tstep:m[1], cchan:m[2]] }
                // display a compressed version of missing hdus
                def missing_timesteps_by_gpubox = missing_hdus
                    .groupBy { m -> m.cchan }
                    .collect { cchan, pairs ->
                        [cchan, displayInts(pairs.collect {pair -> pair.tstep as int})].join(":")
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
                storeDir: "${params.outdir}/results${params.result_suffix}/"
            )
            // display output path and number of lines
            | view { [it, it.readLines().size()] }
    emit:
        // channel of preprocessed files: tuple(obsid, prepUVFits)
        obsPrepUVFits = birliPrep.out.map { def (obsid, _, prepUVFits, __) = it; [obsid, prepUVFits] }
        // channel of mwaf files for each obsid: tuple(obsid, mwafs)
        obsMwafs = birliPrep.out.map { def (obsid, _, __, prepMwafs) = it; [obsid, prepMwafs] }
}

// only emit obsids which pass flag occupancy threshold
workflow flagGate {
    take:
        // channel of required files: tuple(obsid, metafits, mwafs)
        flagFiles

    main:
        // analyse flag occupancy
        flagFiles | flagQA

        // collect flagQA results
        // TODO: collect sky_chans from flagQA
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
                storeDir: "${params.outdir}/results${params.result_suffix}/"
            )
            // display output path and number of lines
            | view { [it, it.readLines().size()] }

    emit:
        // channel of obsids which pass the flag gate
        pass = flagQA.out
            .map { def (obsid, json) = it; [obsid, jslurp.parse(json).total_occupancy] }
            .filter { def (_, occ) = it; occ && occ < params.flag_occupancy_threshold }
            .map { def (obsid, _) = it; obsid }
}

workflow cal {
    take:
        // channel of metafits and preprocessed uvfits: tuple(obsid, metafits, uvfits)
        obsMetaVis
    main:
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

        // hyperdrive di-calibrate on each obs
        obsMetaVis
            // set spw to ""
            .map { def (obsid, metafits, uvfits) = it; [obsid, "", metafits, uvfits] }
            .combine(dical_args)
            | hypCalSol

        // channel of solutions for each unflagged obsid: tuple(obsid, solutions)
        obsSolns = hypCalSol.out.map { def (obsid, _, solutions) = it; [obsid, solutions] }
        // channel of metafits for each obsid: tuple(obsid, metafits)
        obsMetafits = obsMetaVis.map { def (obsid, metafits, _) = it; [obsid, metafits] }

        // channel of individual dical solutions: tuple(obsid, name, soln)
        // - hypCalSol gives multiple solutions, transpose gives 1 tuple per solution.
        eachCal = obsSolns
            .transpose()
            .map { def (obsid, soln) = it
                // give each calibration a name from basename of solution fits
                def name = soln.getBaseName().split('_')[3..-1].join('_')
                [obsid, name, soln]
            }

        // do polyfit, generate json from solns
        eachCal | (polyFit & solJson)

        // collect solJson results as .tsv
        solJson.out
            // form row of tsv from json fields we care about
            .map { def (obsid, name, json) = it;
                def stats = parseJson(json)
                def results = stats.RESULTS?[0]?:[]
                def (nans, convergences) = results.split { it == "NaN" }
                [
                    obsid,
                    name,
                    nans.size() / results.size(),
                    results.join('\t')
                ].join("\t")
            }
            .collectFile(
                name: "cal_results.tsv", newLine: true, sort: true,
                seed: [
                    "OBS", "CAL NAME", "NAN FRAC", "RESULTS BY CH"
                ].join("\t"),
                storeDir: "${params.outdir}/results${params.result_suffix}/"
            )
            // display output path and number of lines
            | view { [it, it.readLines().size()] }

        // calibration QA and plot solutions
        obsMetafits
            // join with hyp_soln from ensureCal
            .cross(eachCal.mix(polyFit.out))
            .map { def (obsid, metafits) = it[0]; def (_, name, soln) = it[1]
                [obsid, name, metafits, soln]
            }
            | (plotSols & calQA)

        // collect calQA results as .tsv
        calQA.out
            // form row of tsv from json fields we care about
            .map { def (obsid, name, json) = it;
                def stats = parseJson(json)
                [
                    obsid,
                    name,
                    stats.FLAGGED_BLS,
                    stats.FLAGGED_CHS,
                    stats.FLAGGED_ANTS,
                    stats.NON_CONVERGED_CHS,
                    stats.NON_CONVERGED_CHS / 100,
                    stats.CONVERGENCE_VAR?[0],
                    stats.CONVERGENCE_VAR?[0] * 100000000000000,
                    stats.XX?.RMS_AMPVAR_FREQ?:'',
                    stats.XX?.RMS_AMPVAR_ANT?:'',
                    stats.XX?.RECEIVER_CHISQVAR?:'',
                    stats.XX?.DFFT_POWER?:'',
                    stats.YY?.RMS_AMPVAR_FREQ?:'',
                    stats.YY?.RMS_AMPVAR_ANT?:'',
                    stats.YY?.RECEIVER_CHISQVAR?:'',
                    stats.YY?.DFFT_POWER?:'',
                    (stats.SKEWNESS_UCVUT?:[]).flatten().join('\t')
                ].join("\t")
            }
            .collectFile(
                name: "cal_metrics.tsv", newLine: true, sort: true,
                seed: [
                    "OBS", "CAL NAME", "FLAG BLS", "FLAG CHS", "FLAG ANTS",
                    "NON CONVG CHS",
                    "NON CONVG CHS FRAC",
                    "CONVG VAR",
                    "CONVG VAR e14",
                    "XX RMS AMPVAR FREQ",
                    "XX RMS AMPVAR ANT",
                    "XX RECEIVER CHISQVAR",
                    "XX DFFT POWER",
                    "YY RMS AMPVAR FREQ",
                    "YY RMS AMPVAR ANT",
                    "YY RECEIVER CHISQVAR",
                    "YY DFFT POWER",
                    "SKEW UVCUT"
                ].join("\t"),
                storeDir: "${params.outdir}/results${params.result_suffix}/"
            )
            // display output path and number of lines
            | view { [it, it.readLines().size()] }

        // hyperdrive apply-solutions using apply solution args file
        // TODO: this is basically unreadable
        obsMetaVis
            .cross(
                apply_args
                    .cross(
                        eachCal
                            .mix(polyFit.out)
                            .map { def (obsid, name, soln) = it; [name.toString(), obsid, soln] }
                    )
                    .map {
                        def (dical_name, apply_name, apply_args) = it[0]
                        def (_, obsid, soln) = it[1]
                        [obsid, soln, dical_name, apply_name, apply_args]
                    }
            )
            .map {
                def (obsid, metafits, prepUVFits) = it[0]
                def (_, soln, dical_name, apply_name, apply_args) = it[1]
                [obsid, metafits, prepUVFits, soln, dical_name, apply_name, apply_args]
            }
            | (hypApplyUV & hypApplyMS)

    emit:
        // channel of calibrated uvfits: tuple(obsid, name, uvfits)
        obsNameUvfits = hypApplyUV.out.map { def (obsid, name, vis) = it; [obsid, name, vis] }
        // channel of calibrated measurement sets: tuple(obsid, name, ms)
        obsNameMS = hypApplyMS.out.map { def (obsid, name, vis) = it; [obsid, name, vis] }
}

// process uvfits visibilities
workflow uvfits {
    take:
        obsNameUvfits
    main:
        // vis QA and ps_metrics
        obsNameUvfits | (visQA & psMetrics)

        // collect visQA results as .tsv
        visQA.out
            // form row of tsv from json fields we care about
            .map { def (obsid, name, json) = it;
                def stats = parseJson(json)
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
                storeDir: "${params.outdir}/results${params.result_suffix}/"
            )
            // display output path and number of lines
            | view { [it, it.readLines().size()] }

        // collect psMetrics as a .dat
        psMetrics.out
            // read the content of each ps_metrics file including the trailing newline
            .map { def (obsid, vis_name, dat) = it; dat.getText() }
            .collectFile(
                name: "ps_metrics.dat",
                storeDir: "${params.outdir}/results${params.result_suffix}/"
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
                storeDir: "${params.outdir}/results${params.result_suffix}/"
            )
            // display output path and number of lines
            | view { [it, it.readLines().size()] }
}

// image measurement sets and QA images
workflow img {
    take:
        obsNameMS
    main:
        // wsclean: make dirty images
        obsNameMS | wscleanDirty

        // imgQA for all groups of images
        wscleanDirty.out | imgQA

        // collect imgQA results as .tsv
        imgQA.out
            // form row of tsv from json fields we care about
            .map { def (obsid, name, json) = it;
                // parse json
                // def stats = jslurp.parse(json)
                // TODO: fix nasty hack to deal with NaNs
                def stats = parseJson(json)
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
                storeDir: "${params.outdir}/results${params.result_suffix}/"
            )
            // display output path and number of lines
            | view { [it, it.readLines().size()] }
}

workflow {
    // get obsids from csv
    obsids = channel.fromPath(params.obsids_path).splitCsv().flatten()

    // download and preprocess raw obsids
    obsids | raw | prep

    // channel of metafits for each obsid: tuple(obsid, metafits)
    obsMetafits = raw.out.map { def (obsid, metafits, _) = it; [obsid, metafits] }

    // channel of obsids that pass the flag gate
    unflaggedObsids = obsMetafits.join(prep.out.obsMwafs) | flagGate

    // calibrate each obs that passes flag gate:
    unflaggedObsids.join(obsMetafits).join(prep.out.obsPrepUVFits) | cal

    // QA uvfits visibilities
    cal.out.obsNameUvfits | uvfits

    // QA images of measurementsets
    cal.out.obsNameMS | img
}