#!/usr/bin/env nextflow
nextflow.enable.dsl=2

def obsids_file = file(params.obsids_path)
if (params.obsids_suffix) {
    obsids_file = file("${obsids_file.getParent()}/${obsids_file.getSimpleName()}${params.obsids_suffix}.${obsids_file.getExtension()}")
}
def results_dir = "${params.outdir}/results${params.img_suffix}${params.cal_suffix}${params.obsids_suffix}${params.result_suffix}/"

// download observation metadata from webservices in json format
process wsMeta {
    input:
    val(obsid)
    output:
    tuple val(obsid), path("${obsid}_wsmeta.json"), path("${obsid}_files.json")

    maxForks 5

    // persist results in outdir, process will be skipped if files already present.
    storeDir "${params.outdir}/meta"
    // tag to identify job in squeue and nf logs
    tag "${obsid}"

    script:
    """
    #!/bin/bash -eux
    ${params.proxy_prelude} # ensure proxy is set if needed
    wget -O ${obsid}_wsmeta.json "http://ws.mwatelescope.org/metadata/obs?obs_id=${obsid}&extended=1&dict=1"
    wget -O ${obsid}_files.json "http://ws.mwatelescope.org/metadata/data_ready?obs_id=${obsid}"
    """
}

// download observation metadata from webservices in metafits format
process wsMetafits {
    input:
    val(obsid)
    output:
    tuple val(obsid), path("${obsid}.metafits")

    maxForks 5

    // persist results in outdir, process will be skipped if files already present.
    storeDir "${params.outdir}/${obsid}/raw"
    // tag to identify job in squeue and nf logs
    tag "${obsid}"

    script:
    """
    #!/bin/bash -eux
    ${params.proxy_prelude} # ensure proxy is set if needed
    wget -O ${obsid}.metafits "http://ws.mwatelescope.org/metadata/fits?obs_id=${obsid}&include_ppds=1"
    """
}

process wsSkyMap {
    input:
    val(obsid)
    output:
    tuple val(obsid), path("${obsid}_skymap.png")

    maxForks 1

    // persist results in outdir, process will be skipped if files already present.
    storeDir "${params.outdir}/${obsid}/meta"
    // tag to identify job in squeue and nf logs
    tag "${obsid}"

    script:
    """
    #!/bin/bash -eux
    ${params.proxy_prelude} # ensure proxy is set if needed
    wget -O ${obsid}_skymap.png "http://ws.mwatelescope.org/observation/skymap/?obs_id=${obsid}"
    """
}

// download preprocessed files from asvo
process asvoPrep {
    input:
    val obsid
    output:
    tuple val(obsid), path("birli_${obsid}_${params.prep_time_res_s}s_${params.prep_freq_res_khz}kHz.uvfits")

    storeDir "${params.outdir}/${obsid}/prep"
    // tag to identify job in squeue and nf logs
    tag "${obsid}"

    maxForks 10

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
            11: "Resource temporarily unavailable",
            75: "Temporary failure, try again"
        ][task.exitStatus] ?: "unknown"
        wait_hours = Math.pow(2, task.attempt)
        println "sleeping for ${wait_hours} hours and retrying task ${task.hash}, which failed with code ${task.exitStatus}: ${retry_reason}"
        sleep(wait_hours * 60*60*1000 as long)
        return 'retry'
    }

    script:
    """
    #!/bin/bash -eux

    export MWA_ASVO_API_KEY="${params.asvo_api_key}"

    ${params.proxy_prelude} # ensure proxy is set if needed

    # submit a job to ASVO, suppress failure if a job already exists.
    ${params.giant_squid} submit-conv -v \
        -p timeres=${params.prep_time_res_s},freqres=${params.prep_freq_res_khz},conversion=uvfits,preprocessor=birli \
        ${obsid} || true

    # list pending conversion jobs
    ${params.giant_squid} list -j --types conversion --states queued,processing,error -- ${obsid}

    # extract id url size hash from any ready download vis jobs for this obsid
    ${params.giant_squid} list -j --types conversion --states ready -- ${obsid} \
        | tee /dev/stderr \
        | ${params.jq} -r '.[]|[.jobId,.files[0].fileUrl//"",.files[0].fileSize//"",.files[0].fileHash//""]|@tsv' \
        | tee ready.tsv

    # download the first ready jobs
    if read -r jobid url size hash < ready.tsv; then
        if read -r avail < <(df --output=avail -B1 . | tail -n 1); then
            if [[ \$avail -lt \$size ]]; then
                echo "Not enough disk space available in \$PWD, need \$size B, have \$avail B"
                exit 28  # No space left on device
            fi
        else
            echo "Warning: Could not determine disk space in \$PWD"
        fi
        # wget into two pipes:
        # - first pipe untars the archive to disk
        # - second pipe validates the sha sum
        wget \$url -O- --progress=dot:giga --wait=60 --random-wait \
            | tee >(tar -x) \
            | sha1sum -c <(echo "\$hash -")
        ps=("\${PIPESTATUS[@]}")
        if [ \${ps[0]} -ne 0 ]; then
            echo "Download failed. status=\${ps[0]}"
            exit \${ps[0]}
        elif [ \${ps[1]} -ne 0 ]; then
            echo "Untar failed. status=\${ps[1]}"
            exit \${ps[1]}
        elif [ \${ps[2]} -ne 0 ]; then
            echo "Hash check failed. status=\${ps[2]}"
            exit \${ps[2]}
        fi
        mv ${obsid}.uvfits birli_${obsid}_${params.prep_time_res_s}s_${params.prep_freq_res_khz}kHz.uvfits
        exit 0 # success
    fi
    echo "no ready jobs"
    exit 75 # temporary
    """
}

// QA tasks for flags.
process flagQA {
    input:
    tuple val(obsid), path("${obsid}.metafits"), path("${obsid}.uvfits")
    output:
    tuple val(obsid), path("${obsid}_occupancy.json")

    storeDir "${params.outdir}/${obsid}/prep"

    tag "${obsid}"

    label 'python'

    script:
    """
    #!/usr/bin/env python

    from astropy.io import fits
    from json import dump as json_dump
    import numpy as np

    def split_strip_filter(str):
        return list(filter(None, map(lambda tok: tok.strip(), str.split(','))))

    total_occupancy = 0
    total_non_preflagged_bl_occupancy = 0
    num_cchans = 0
    num_unflagged_cchans = 0
    with \
        open('${obsid}_occupancy.json', 'w') as out, \
        fits.open("${obsid}.metafits") as meta, \
        fits.open("${obsid}.uvfits") as uv \
    :
        data = {'obsid': ${obsid}, 'channels': {}}
        td = meta['TILEDATA'].data
        preflagged_ants = np.unique(td[td['Flag'] > 0]['Antenna'])
        data['preflagged_ants'] = preflagged_ants.tolist()
        sky_chans = [*map(int, split_strip_filter(meta[0].header['CHANNELS']))]
        num_cchans = len(sky_chans)

        vis_hdu = uv['PRIMARY']
        baseline_array = np.int16(vis_hdu.data["BASELINE"])
        num_blts = len(baseline_array)
        num_bls = len(np.unique(baseline_array))
        data['num_baselines']=num_bls
        num_times = num_blts // num_bls
        data['num_times']=num_times
        vis_shape = vis_hdu.data.data.shape
        assert vis_shape[0] == num_blts, 'vis shape 0 != num_blts'

        num_chans = vis_hdu.header['NAXIS4']
        data['num_chans']=num_chans
        assert num_chans % num_cchans == 0, 'uvfits channels is not divisible by coarse chans'
        num_fchans = num_chans // num_cchans # number of fine chans per coarse
        assert num_fchans * num_cchans == vis_shape[3], 'vis shape 3 != NAXIS4'

        ant_2_array = baseline_array % 256 - 1
        ant_1_array = (baseline_array - ant_2_array) // 256 - 1
        antpair_array = np.stack((ant_1_array, ant_2_array), axis=1)
        preflagged_blts_mask = np.isin(antpair_array, preflagged_ants).any(axis=1)
        # unflagged_blt_idxs = np.where(np.logical_not(preflagged_blts_mask))[0]
        (unflagged_blt_idxs,) = np.where(np.logical_not(preflagged_blts_mask))
        num_unflagged_blts = len(unflagged_blt_idxs)
        num_unflagged_bls = num_unflagged_blts // num_times
        unflagged_fraction = num_unflagged_blts / num_blts

        # assumption: vis data axes 1,2 (ra/dec?) are not used
        # assumption: flags are the same for all pols, so only use pol=0
        # 4d weight array: [time, baseline, coarse channel, fine channel]
        weights = vis_hdu.data.data[unflagged_blt_idxs, 0, 0, :, 0, 2].reshape(
            (num_times, num_unflagged_bls, num_cchans, num_fchans))
        # flag count for unflagged bls by [time, coarse channel, fine channel]
        unflagged_bl_flag_count = np.sum(np.int64(weights <= 0), axis=1)
        # where unflagged_bl_flag_count is flagged for all baselines
        flagged_mask = unflagged_bl_flag_count == num_unflagged_bls
        data['flagged_timestep_idxs'] = np.where(flagged_mask.all(axis=(1,2)))[0].tolist()
        flagged_cchan_idxs = np.where(flagged_mask.all(axis=(0,2)))[0].tolist()
        data['flagged_cchan_idxs'] = flagged_cchan_idxs
        flagged_sky_chans = np.array(sky_chans)[flagged_cchan_idxs].tolist()
        data['flagged_sky_chans'] = flagged_sky_chans
        data['flagged_fchan_idxs'] = np.where(flagged_mask.all(axis=(0,1)))[0].tolist()

        # occupancy of non-preflagged baselines by coarse channel
        unflagged_bl_occupancy = np.sum(unflagged_bl_flag_count, axis=(0,2)) / (num_unflagged_blts * num_fchans)

        for chan_idx, cchan_unflagged_bl_occupancy in enumerate(unflagged_bl_occupancy):
            if chan_idx not in flagged_cchan_idxs:
                num_unflagged_cchans += 1
                total_non_preflagged_bl_occupancy += cchan_unflagged_bl_occupancy
            sky_chan = sky_chans[chan_idx]
            cchan_total_occupancy = (unflagged_fraction * cchan_unflagged_bl_occupancy) + (1 - unflagged_fraction)
            total_occupancy += cchan_total_occupancy
            data['channels'][sky_chan] = {
                'occupancy': cchan_total_occupancy,
                'non_preflagged_bl_occupancy': cchan_unflagged_bl_occupancy,
                'flagged': chan_idx in data['flagged_cchan_idxs']
            }
        if num_cchans:
            total_occupancy /= num_cchans
        data['total_occupancy'] = total_occupancy
        if num_unflagged_cchans:
            total_non_preflagged_bl_occupancy /= num_unflagged_cchans
        data['total_non_preflagged_bl_occupancy'] = total_non_preflagged_bl_occupancy
        json_dump(data, out, indent=4)
    """
}

def groovy2bashAssocArray(map, name) {
    // size = map.collect { _, v -> v.size() }.max()
    "declare -A ${name}=(" + map.collect { k, v -> "[${k}]=\"${v}\"".toString() }.join(" ") + ")".toString()
}

// ensure calibration solutions are present, or calibrate prep with hyperdrive
// do multiple calibration solutions for each obs, depending on dical_args
process hypCalSol {
    input:
    tuple val(obsid), val(dical_args), path("${obsid}.metafits"), path("${obsid}.uvfits"), val(tile_flags)
    output:
    tuple val(obsid), path("hyp_soln_${obsid}_${name_glob}.fits"), path("hyp_di-cal_${obsid}_${name_glob}.log")
    // todo: model subtract: path("hyp_model_${obsid}_${name_glob}.uvfits")

    storeDir "${params.outdir}/${obsid}/cal${params.cal_suffix}"

    tag "${obsid}"

    // label jobs that need a bigger gpu allocation
    label "hyperdrive"

    afterScript "ls ${task.storeDir}"

    script:
    dical_names = dical_args.keySet().collect()
    para = dical_names.size() > 1
    name_glob = para ? "{" + dical_names.join(',') + "}" : dical_names[0]
    flag_args = tile_flags.size() > 0 ? "--tile-flags ${tile_flags.join(' ')}" : ""
    """
    set -eux
    ${para ? "export CUDA_VISIBLE_DEVICES=0" : ""}
    export num_gpus="\$(nvidia-smi -L | wc -l)"
    if [ \$num_gpus -eq 0 ]; then
        echo "no gpus found"
        exit 1
    fi
    if [ \$num_gpus -ne ${params.num_gpus} ]; then
        echo "warning: expected \$num_gpus to be ${params.num_gpus}"
    fi
    ${groovy2bashAssocArray(dical_args, "dical_args")}
    if [ \${#dical_args[@]} -eq 0 ]; then
        echo "no dical args"
        exit 0
    fi
    if [ \${#dical_args[@]} -gt \$num_gpus ]; then
        echo "warning: more dical args than gpus";
    fi
    for name in \${!dical_args[@]}; do
        args=\${dical_args[\$name]}
        # hyperdrive di-cal backgrounded in a subshell
        (
            ${params.hyperdrive} di-calibrate \${args} \
                --data "${obsid}.metafits" "${obsid}.uvfits" \
                --beam "${params.beam_path}" \
                --source-list "${params.sourcelist}" \
                --outputs "hyp_soln_${obsid}_\${name}.fits" \
                ${flag_args} \
                > "hyp_di-cal_${obsid}_\${name}.log"
        ) &
        # TODO: model subtract: --model-filenames "hyp_model_${obsid}_\${name}.uvfits"
        # increment the target device mod num_gpus
        ${para ? "CUDA_VISIBLE_DEVICES=\$(( CUDA_VISIBLE_DEVICES+1 % \${num_gpus} ))" : ""}
    done

    # wait for all the background jobs to finish
    export result=0
    wait \$(jobs -rp) || result=\$?

    # print out important info from log
    if [ \$(ls *.log 2> /dev/null | wc -l) -gt 0 ]; then
        grep -iE "err|warn|hyperdrive|chanblocks|reading|writing|flagged" *.log
    fi

    exit \$result
    """
}

// fit polynomial to calibration solution
process polyFit {
    input:
    tuple val(obsid), val(dical_name), path("hyp_soln_${obsid}_${dical_name}.fits")
    output:
    tuple val(obsid), val(name), path("hyp_soln_${obsid}_${name}.fits"), path("polyfit_${obsid}_${name}.log")

    storeDir "${params.outdir}/${obsid}/cal${params.cal_suffix}"
    tag "${obsid}.${dical_name}"

    label 'python'

    script:
    name = "poly_${dical_name}"
    """
    #!/bin/bash -eux
    run_polyfit.py \
        hyp_soln_${obsid}_${dical_name}.fits \
        --outfile "hyp_soln_${obsid}_${name}.fits" \
        | tee polyfit_${obsid}_${name}.log
    """
}

// ensure calibrated uvfits are present, or apply calsols to prep uvfits with hyperdrive
process hypApplyUV {
    input:
    tuple val(obsid), path("${obsid}.metafits"), path("${obsid}.uvfits"), path("hyp_soln_${obsid}_${cal_name}.fits"), \
        val(cal_name), val(apply_name), val(apply_args)
    output:
    tuple val(obsid), val(vis_name), path("hyp_${obsid}_${vis_name}.uvfits"), path("hyp_apply_${vis_name}.log")

    storeDir "${params.outdir}/${obsid}/cal${params.cal_suffix}"

    tag "${obsid}.${cal_name}_${apply_name}"
    label "cpu"
    label "hyperdrive"

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
    if [ \$(ls *.log 2> /dev/null | wc -l) -gt 0 ]; then
        grep -iE "err|warn|hyperdrive|chanblocks|reading|writing|flagged" *.log
    fi
    """
}

// ensure calibrated ms are present, or apply calsols to prep uvfits with hyperdrive
process hypApplyMS {
    input:
    tuple val(obsid), path("${obsid}.metafits"), path("${obsid}.uvfits"), path("hyp_soln_${obsid}_${cal_name}.fits"), \
        val(cal_name), val(apply_name), val(apply_args)
    output:
    tuple val(obsid), val(vis_name), path("hyp_${obsid}_${vis_name}.ms"), path("hyp_apply_${vis_name}_ms.log")

    storeDir "${params.outdir}/${obsid}/cal${params.cal_suffix}"
    // storeDir "/data/curtin_mwaeor/FRB_hopper/"

    tag "${obsid}.${cal_name}_${apply_name}"
    label "cpu"
    label "hyperdrive"

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
    if [ \$(ls *.log 2> /dev/null | wc -l) -gt 0 ]; then
        grep -iE "err|warn|hyperdrive|chanblocks|reading|writing|flagged" *.log
    fi
    """
}

process hypSubUV {
    input:
    tuple val(obsid), val(name), path("${obsid}.metafits"), path("hyp_${obsid}_${name}.uvfits")
    output:
    tuple val(obsid), val(sub_name), path("hyp_${obsid}_${sub_name}.uvfits"), \
        path("hyp_vis-sub_${name}_uv.log")

    storeDir "${params.outdir}/${obsid}/cal${params.cal_suffix}"

    tag "${obsid}.${name}"
    label "hyperdrive"

    script:
    sub_name = "sub_${name}"
    """
    ${params.hyperdrive} vis-sub \
        --data "${obsid}.metafits" "hyp_${obsid}_${name}.uvfits" \
        --beam "${params.beam_path}" \
        --source-list "${params.sourcelist}" \
        --invert --num-sources 4000 \
        --outputs "hyp_${obsid}_${sub_name}.uvfits" \
        | tee hyp_vis-sub_${name}_uv.log
    # TODO: ^ num sources is hardcoded twice, would be better to re-use model from cal

    # print out important info from log
    if [ \$(ls *.log 2> /dev/null | wc -l) -gt 0 ]; then
        grep -iE "err|warn|hyperdrive|chanblocks|reading|writing|flagged" *.log
    fi
    """
}
process hypSubMS {
    input:
    tuple val(obsid), val(name), path("${obsid}.metafits"), path("hyp_${obsid}_${name}.ms")
    output:
    tuple val(obsid), val(sub_name), path("hyp_${obsid}_${sub_name}.ms"), \
        path("hyp_vis-sub_${name}_ms.log")

    storeDir "${params.outdir}/${obsid}/cal${params.cal_suffix}"

    tag "${obsid}.${name}"
    label "hyperdrive"

    script:
    sub_name = "sub_${name}"
    """
    ${params.hyperdrive} vis-sub \
        --data "${obsid}.metafits" "hyp_${obsid}_${name}.ms" \
        --beam "${params.beam_path}" \
        --source-list "${params.sourcelist}" \
        --invert --num-sources 4000 \
        --outputs "hyp_${obsid}_${sub_name}.ms" \
        | tee hyp_vis-sub_${name}_ms.log
    """
}

// QA tasks that can be run on preprocessed visibility files.
process prepVisQA {
    input:
    tuple val(obsid), path("birli_${obsid}.uvfits")
    output:
    tuple val(obsid), path("birli_${obsid}_prepvis_metrics.json")

    storeDir "${params.outdir}/${obsid}/vis_qa"

    tag "${obsid}"

    label 'cpu'
    label 'python'

    script:
    """
    #!/bin/bash -eux
    run_prepvisqa.py birli_${obsid}.uvfits --out "birli_${obsid}_prepvis_metrics.json"
    """
}

// QA tasks that can be run on calibrated visibility files.
process visQA {
    input:
    tuple val(obsid), val(name), path("hyp_${obsid}_${name}.uvfits")
    output:
    tuple val(obsid), val(name), path("hyp_${obsid}_${name}_vis_metrics.json")

    storeDir "${params.outdir}/${obsid}/vis_qa${params.cal_suffix}"

    tag "${obsid}.${name}"

    label 'cpu'
    label 'python'

    script:
    """
    #!/bin/bash -eux
    run_visqa.py *.uvfits --out "hyp_${obsid}_${name}_vis_metrics.json"
    """
}

// QA tasks that can be run on each calibration parameter set
process calQA {
    input:
    tuple val(obsid), val(name), path("${obsid}.metafits"), path("soln.fits")
    output:
    tuple val(obsid), val(name), path("hyp_soln_${obsid}_${name}_X.json")

    storeDir "${params.outdir}/${obsid}/cal_qa${params.cal_suffix}"

    tag "${obsid}.${name}"

    label 'python'

    script:
    """
    #!/bin/bash -eux
    run_calqa.py soln.fits ${obsid}.metafits --pol X --out "hyp_soln_${obsid}_${name}_X.json"
    """
}

// write info from solutions to json
process solJson {
    input:
    tuple val(obsid), val(dical_name), path("hyp_soln_${obsid}_${dical_name}.fits")
    output:
    tuple val(obsid), val(dical_name), path("hyp_soln_${obsid}_${dical_name}.fits.json")

    storeDir "${params.outdir}/${obsid}/cal_qa${params.cal_suffix}"
    tag "${obsid}.${dical_name}"

    label 'python'

    script:
    """
    #!/usr/bin/env python

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

process plotPrepVisQA {
    input:
    tuple val(obsid), path("birli_${obsid}_prepvis_metrics.json")
    output:
    tuple val(obsid), path("prepvis_metrics_${obsid}_rms.png")

    storeDir "${params.outdir}/${obsid}/prep"

    tag "${obsid}"

    label 'python'

    script:
    """
    #!/bin/bash -eux
    plot_prepvisqa.py "birli_${obsid}_prepvis_metrics.json" --out "prepvis_metrics_${obsid}_rms.png" --save
    """
}

process plotSols {
    input:
    tuple val(obsid), val(name), path("${obsid}.metafits"), path("hyp_soln_${obsid}_${name}.fits")
    output:
    tuple val(obsid), val(name), path("hyp_soln_${obsid}*_${name}_{phases,amps}.png")

    storeDir "${params.outdir}/${obsid}/cal_qa${params.cal_suffix}"

    tag "${obsid}.${name}"

    label "hyperdrive"

    script:
    """
    hyperdrive solutions-plot -m "${obsid}.metafits" *.fits
    """
}

process plotCalQA {
    input:
    tuple val(obsid), val(name), path("hyp_soln_${obsid}_${name}.json")
    output:
    tuple val(obsid), val(name), path("calmetrics_${obsid}_${name}_{fft,variance,dlyspectrum}.png")

    storeDir "${params.outdir}/${obsid}/cal_qa${params.cal_suffix}"

    tag "${obsid}.${name}"

    label 'python'

    script:
    """
    #!/bin/bash -eux
    plot_calqa.py "hyp_soln_${obsid}_${name}.json" --out "calmetrics_${obsid}_${name}.png" --save
    """
}

process plotVisQA {
    input:
    tuple val(obsid), val(name), path("hyp_${obsid}_${name}_vis_metrics.json")
    output:
    tuple val(obsid), val(name), path("hyp_${obsid}_${name}_vis_metrics_rms.png")

    storeDir "${params.outdir}/${obsid}/vis_qa${params.cal_suffix}"

    tag "${obsid}"

    label 'python'

    script:
    """
    #!/bin/bash -eux
    plot_visqa.py "hyp_${obsid}_${name}_vis_metrics.json" --out "hyp_${obsid}_${name}_vis_metrics_rms.png" --save
    """
}

// process plotImgQA {
//     input:
//     tuple val(obsid), val(name), path("wsclean_hyp_${obsid}_${name}-MFS.json")
//     output:
//     tuple val(obsid), val(name), \
//         path("wsclean_hyp_${obsid}_${name}-MFS_rms.png"), \
//         path("wsclean_hyp_${obsid}_${name}-MFS_pks.png")

//     storeDir "${params.outdir}/${obsid}/img_qa"

//     tag "${obsid}"

//     label 'python'
//     conda "${params.astro_conda}"

//     script:
//     """
//     #!/bin/bash
//     set -ex
//     plot_imgqa.py "wsclean_hyp_${obsid}_${name}-MFS.json" --out "wsclean_hyp_${obsid}_${name}-MFS" --save
//     """
// }

// create dirty iamges of xx,yy,v
process wscleanDirty {
    input:
    tuple val(obsid), val(name), path("vis.ms")
    output:
    tuple val(obsid), val(name), path("wsclean_hyp_${obsid}_${name}-MFS-{XX,YY,V}-dirty.fits"), path("wsclean_${name}.log")

    storeDir "${params.outdir}/${obsid}/img${params.img_suffix}${params.cal_suffix}"

    tag "${obsid}.${name}"
    label "wsclean"

    script:
    """
    set -eux
    # imaging
    ${params.wsclean} \
        ${params.wsclean_args} \
        -weight ${params.img_weight} \
        -name wsclean_hyp_${obsid}_${name} \
        -size ${params.img_size} ${params.img_size} \
        -scale ${params.img_scale} \
        -pol xx,yy,v \
        -channels-out ${params.img_channels_out} \
        -niter ${params.img_niter} \
        vis.ms | tee wsclean_${name}.log
    """
}

// power spectrum metrics via chips
process psMetrics {
    input:
    tuple val(obsid), val(name), path("hyp_${obsid}_${name}.uvfits")
    output:
    tuple val(obsid), val(name), path("output_metrics_hyp_${obsid}_${name}.dat"), \
        path("hyp_${obsid}_${name}.log")

    storeDir "${params.outdir}/${obsid}/ps_metrics${params.cal_suffix}"

    tag "${obsid}.${name}"

    label "chips"

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

    storeDir "${params.outdir}/${obsid}/img_qa${params.img_suffix}${params.cal_suffix}"

    tag "${obsid}.${name}"

    label 'python'

    script:
    """
    #!/bin/bash -eux
    run_imgqa.py *.fits --out wsclean_hyp_${obsid}_${name}-MFS.json
    """
}

// polarimetry composite raster
process polComp {
    input:
    tuple val(obsid), val(name), path("*") // <- "*-dirty.fits"
    output:
    tuple val(obsid), val(name), path("${obsid}_${name}_polcomp.png")
        // , path("${obsid}_${name}_{XX,YY,V}.png")

    storeDir "${params.outdir}/${obsid}/img_qa${params.img_suffix}${params.cal_suffix}"

    tag "${obsid}.${name}"

    label 'python'

    script:
    """
    #!/usr/bin/env python

    from astropy.io import fits
    from astropy.visualization import astropy_mpl_style, make_lupton_rgb, simple_norm, SqrtStretch
    from astropy.visualization import ImageNormalize, PowerDistStretch, SinhStretch, LogStretch, AsinhStretch
    from astropy.visualization.lupton_rgb import AsinhMapping, LinearMapping
    from astropy.visualization.wcsaxes import WCSAxesSubplot, WCSAxes
    from astropy.wcs import WCS
    from matplotlib.patches import Rectangle, Circle
    from matplotlib.colors import LinearSegmentedColormap
    import matplotlib.pyplot as plt
    from mpl_toolkits.axes_grid1.axes_divider import AxesDivider, HBoxDivider, make_axes_locatable, Size
    from mpl_toolkits.axes_grid1.axes_rgb import make_rgb_axes
    from mpl_toolkits.axes_grid1.inset_locator import inset_axes
    import numpy as np
    from os.path import basename
    from glob import glob

    plt.style.use(astropy_mpl_style)

    header = None
    data = {}
    for key in ["XX", "YY", "V"]:
        path = glob(f"*${obsid}*${name}*{key}*.fits")[0]
        with fits.open(path) as hdus:
            if not header:
                header = hdus[0].header
            data[key] = hdus[0].data[0,0,:,:]

    wcs = WCS(header)
    img_size = header['NAXIS1'], header['NAXIS2']
    # TODO: beam stuff
    # cell_size = header['CDELT1'], header['CDELT2']
    # beam_shape = header['BMAJ'], header['BMIN'], header['BPA']

    # normalize 0-percentile to 0-1, use same percentile for XX and YY
    i_percentile = np.percentile(np.stack((data["XX"], data["YY"])), 99.5)
    v_percentile = np.percentile(data["V"], 99.5)
    norm = {
        "XX": data["XX"] / i_percentile,
        "YY": data["YY"] / i_percentile,
        "V": data["V"] / v_percentile,
    }

    rgbMap = AsinhMapping(0., stretch=3., Q=5.)
    rgb = rgbMap.make_rgb_image(norm["XX"], norm["V"], norm["YY"])

    sources = np.array([
        [0, -27], # eor0
        [6.4549166666666675, -26.04], # PKS0023_026
        # [45, -26.04]
    ])

    # stretch=AsinhStretch(0.05)
    # stretch=PowerDistStretch(0.05)
    stretch=SinhStretch()

    dpi = 100
    with plt.style.context('dark_background'):
        # axis setup
        img_fig = plt.figure(dpi=dpi, figsize=(img_size[0]/dpi, 4*img_size[1]/3/dpi))

        axd = img_fig.subplot_mosaic(
            [
                ["Comp", "XX"],
                ["Comp", "YY"],
                ["Comp", "V"],
            ],
            subplot_kw={"projection": wcs, "slices": ('x', 'y', 0, 0)},
            gridspec_kw={ "width_ratios": [3, 1], },
        )

        for ax in axd.values():
            ax.set_ylabel("", visible=False)
            ax.set_xlabel("", visible=False)
            ax.set_label("")

        divider = AxesDivider(axd['Comp'])
        locator = divider.new_locator(nx=0, ny=0)
        axd['Comp'].set_axes_locator(locator)

        sub_xsize = (1/3) * Size.AxesX(axd['Comp'])
        cb_xsize = (1/10) * sub_xsize
        ysize = (1./3) * Size.AxesY(axd['Comp'])

        divider.set_horizontal([Size.AxesX(axd['Comp']), sub_xsize, cb_xsize])
        divider.set_vertical([ysize, ysize, ysize])

        axd['Comp'].set_axes_locator(divider.new_locator(0, 0, ny1=-1))
        axd["Comp"].tick_params(axis="y", direction="in", pad=-30, horizontalalighment="left")
        axd["Comp"].tick_params(axis="x", direction="in", pad=-20, verticalalignment="bottom")
        axd["Comp"].imshow(rgb)
        axd["Comp"].set_title("${obsid} - ${name}", y=0.95, fontdict={'verticalalignment':'top'})

        for key, ny, neg_col, pos_col, vmax in [
            ['XX', 2, 'cyan', 'red', i_percentile],
            ['YY', 1, 'yellow', 'blue', i_percentile],
            ['V', 0, 'magenta', 'green', v_percentile]
        ]:
            ax = axd[key]
            ax.set_title(key, y=0.95, fontdict={'verticalalignment':'top'})
            ax.set_axes_locator(divider.new_locator(nx=1, ny=ny))
            cmap = LinearSegmentedColormap.from_list(f'{neg_col}-bl-{pos_col}', [neg_col, 'black', pos_col], gamma=1)
            norm=ImageNormalize(data[key], vmin=-vmax, vmax=vmax, clip=False)
            im = ax.imshow(data[key], cmap=cmap, norm=norm)
            for coord in ax.coords:
                coord.set_ticklabel_visible(False)
            cax = inset_axes(
                ax,
                width="5%",
                height="100%",
                loc='lower left',
                bbox_to_anchor=(1.0, 0., 1, 1),
                bbox_transform=ax.transAxes,
                borderpad=0,
                axes_kwargs={ "axisbelow": False}
            )
            cbar = img_fig.colorbar(im, cax=cax)

        for ax in axd.values():
            ax.scatter(sources[:,0], sources[:,1], s=100, edgecolor='white', facecolor='none', transform=ax.get_transform('world'))
            ax.add_patch(Rectangle((0, 0), 100, 100, edgecolor='white', fill=False, zorder=999))

        #  hist_ax = inset_axes(
        #      axd["Comp"],
        #      width="100%",
        #      height="20%",
        #      loc='lower left',
        #      bbox_to_anchor=(0., -0.25, 1, 1),
        #      bbox_transform=axd['Comp'].transAxes,
        #      borderpad=0,
        #  )

        #  hist_ax.hist(data["XX"].flatten(), bins=1000, color='red', alpha=0.5, histtype='step')
        #  hist_ax.hist(data["YY"].flatten(), bins=1000, color='blue', alpha=0.5, histtype='step')
        #  hist_ax.hist(data["V"].flatten(), bins=1000, color='green', alpha=0.5, histtype='step')

        plt.savefig('${obsid}_${name}_polcomp.png', bbox_inches='tight', dpi=dpi)
    """
}

process plotCalJsons {
    input:
        val jsons
    output:
        path("cal_qa_{convergence,fft,rms,unused}.png")

    storeDir "${results_dir}"
    stageInMode "symlink"

    label 'python'

    script:
    """
    #!/bin/bash -eux
    plot_caljsons.py ${jsons.join(' ')} --save
    """
}

process ffmpeg {
    input:
        tuple val(name), path("??????????.png")

    output:
        path("${name}.mp4")

    storeDir "${results_dir}"
    stageInMode "symlink"

    tag "${name}"

    label 'ffmpeg'

    script:
    """
    #!/bin/bash -eux
    ffmpeg -y -framerate 5 \
        -pattern_type glob -i "??????????.png" \
        -vcodec libx264 \
        -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2" \
        -pix_fmt yuv420p \
        "${name}.mp4"
    """
}

process archive {
    input:
        tuple val(bucket_suffix), path(x)
    output:
        path("${x}.shadow")

    tag "$x"
    storeDir "${params.outdir}/.fakebuckets/${bucket}"

    label "rclone"
    stageInMode "symlink"

    script:
    bucket = "${params.bucket_prefix}.${bucket_suffix}"
    """
    #!/bin/bash -eux
    touch "${x}.shadow"
    ${params.proxy_prelude} # ensure proxy is set if needed
    rclone mkdir "${bucket}"
    rclone copyto --copy-links "$x" "${bucket}/$x"
    """
}

// collect ao quality metrics
process aoQuality {
    input:
    tuple val(obsid), val(name), path("vis.ms")

    output:
    tuple val(obsid), val(name), path("${obsid}_${name}_aoquality_{sum,rfi,b,t,f,q}.tsv")

    storeDir "${params.outdir}/${obsid}/vis_qa${params.cal_suffix}"
    stageInMode "symlink"
    script:
    """
    #!/bin/bash -eux
    ${params.aoquality} collect vis.ms
    ${params.aoquality} summarize vis.ms | tee ${obsid}_${name}_aoquality_sum.tsv
    ${params.aoquality} summarizerfi vis.ms | tee ${obsid}_${name}_aoquality_rfi.tsv
    for q in b t f q; do
        ${params.aoquality} liststats \
        | while read -r stat; do
            result=\$(${params.aoquality} query_\$q \$stat vis.ms) \
                && echo \$stat \$result \
                | tail -a ${obsid}_${name}_aoquality_\${q}.tsv;
        done
    done
    """

    // ${aoquality} collect vis.ms
    // ${aoquality} summarize vis.ms | tee ${obsid}_${name}_aoquality_sum.tsv
    // ${aoquality} summarizerfi vis.ms | tee ${obsid}_${name}_aoquality_rfi.tsv
    // for q in b t f q; do
    //     ${aoquality} liststats \
    //     | while read -r stat; do
    //         result=$(${aoquality} query_$q $stat vis.ms) \
    //             && echo $stat $result \
    //             | tee -a ${obsid}_${name}_aoquality_${q}.tsv;
    //     done
    // done
}

import groovy.json.JsonSlurper
import groovy.json.StringEscapeUtils
jslurp = new JsonSlurper()
def parseJson(path) {
    // TODO: fix nasty hack to deal with NaNs
    jslurp.parseText(path.getText().replaceAll(/(NaN|-?Infinity)/, '"$1"'))
}

import java.text.SimpleDateFormat
SimpleDateFormat logDateFmt = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss");

// display a long list of ints, replace bursts of consecutive numbers with ranges
def displayRange = { s, e -> s == e ? "${s}," : s == e - 1 ? "${s},${e}," : "${s}-${e}," }
max_ints_display = 50
def displayInts = { l ->
    switch (l) {
        case { it.size() == 0 }: return "";
        case { it.size() == 1 }: return displayRange(l[0],l[0]);
        default:
            def sb, start, end
            (sb, start, end) = [''<<'', l[0], l[0]]
            for (i in l[1..-1]) {
                (sb, start, end) = i == end + 1 ? [sb, start, i] : [sb << displayRange(start, end), i, i]
            }
            result = (sb << displayRange(start, end))[0..-2].toString()
            if (result.size() > max_ints_display) {
                return result.substring(0, (max_ints_display-1).intdiv(2))
                << "..."
                << result.substring(result.size() - (max_ints_display-1).intdiv(2))
            } else {
                return result
            }
    }
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

unit_conversions = [
    "s": 1,
    "ms": 1e-3,
    "µs": 1e-6,
    "ns": 1e-9,
]
def get_seconds(float time, String unit) {
    if (unit_conversions[unit] == null) {
        println "unknown duration unit ${unit} for time ${time}";
        time
    } else {
        time * unit_conversions[unit]
    }
}

// use mwa webservices to gate obsids by pointing and faults
workflow ws {
    take:
        // channel of obs ids
        obsids

    main:
        obsids | wsMeta

        quality_updates = file(params.quality_updates_path)
            .readLines()
            .collectEntries { line ->
                def (obsid, quality, comment) = line.split(',')
                [obsid, [dataquality: quality, dataqualitycomment: comment]]
            }

        wsSummary = wsMeta.out.map { obsid, json, filesJson ->
                // parse json
                def stats = parseJson(json)

                def ra_phase_center = stats.ra_phase_center;
                if (ra_phase_center == null) {
                    ra_phase_center = stats.metadata.ra_pointing
                }
                def dec_phase_center = stats.dec_phase_center;
                if (dec_phase_center == null) {
                    dec_phase_center = stats.metadata.dec_pointing
                }

                def pointing = stats.metadata.gridpoint_number
                def nscans = ((stats.stoptime?:0) - (stats.starttime?:0)) / (stats.int_time?:1)
                def delays = (stats.alldelays?:[:]).values().flatten()
                def quality = stats.quality?:[:]
                def tiles = stats.tdict?:[:]
                def bad_tiles = stats.bad_tiles?:[:]
                def n_bad_tiles = bad_tiles.size()
                def bad_tile_frac = n_bad_tiles / tiles.size()
                def n_dead_dipoles = delays.count { it == 32 }
                def dead_dipole_frac = n_dead_dipoles / delays.size()
                def quality_update = quality_updates[obsid]?:[:]
                def dataquality = Float.valueOf(quality_update.dataquality?:stats.dataquality?:0)
                def dataqualitycomment = quality_update.dataqualitycomment?:stats.dataqualitycomment?:''
                def faults = stats.faults?:[:]
                def badstates = (faults.badstates?:[:]).values().flatten()
                def badpointings = (faults.badpointings?:[:]).values().flatten()
                def badfreqs = (faults.badfreqs?:[:]).values().flatten()
                def badgains = (faults.badgains?:[:]).values().flatten()
                def badbeamshape = (faults.badbeamshape?:[:]).values().flatten()
                def fail_reasons = []

                def fileStats = parseJson(filesJson)
                if (!params.filter_pointings.contains(pointing)) {
                    fail_reasons += ["pointing=${pointing}"]
                }
                if (dataquality > params.filter_quality) {
                    fail_reasons += ["dataquality=${dataquality} (${dataqualitycomment})"]
                }
                if (bad_tile_frac > params.filter_bad_tile_frac) {
                    fail_reasons += ["bad_tiles(${bad_tiles.size()})=${displayInts(bad_tiles)}"]
                }
                if (dead_dipole_frac > params.filter_dead_dipole_frac) {
                    fail_reasons += ["dead_dipole_frac=${dead_dipole_frac}"]
                }
                def summary = [
                    fail_reasons: fail_reasons,
                    // obs metadata
                    nscans: nscans,
                    ra_pointing: stats.metadata.ra_pointing,
                    dec_pointing: stats.metadata.dec_pointing,
                    ra_phase_center: ra_phase_center,
                    dec_phase_center: dec_phase_center,
                    pointing: pointing,
                    lst: stats.metadata.local_sidereal_time_deg,
                    freq_res: stats.freq_res,
                    int_time: stats.int_time,
                    tiles: tiles,
                    bad_tiles: bad_tiles,
                    // iono quality
                    iono_magnitude: quality.iono_magnitude?:'',
                    iono_pca: quality.iono_pca?:'',
                    iono_qa: quality.iono_qa?:'',
                    // data quality
                    dataquality: dataquality,
                    dataqualitycomment: dataqualitycomment,
                    // fraction of flagged tiles
                    bad_tile_frac: bad_tile_frac,
                    n_bad_tiles: n_bad_tiles,
                    // fraction of dead dipoles
                    dead_dipole_frac: dead_dipole_frac,
                    n_dead_dipoles: n_dead_dipoles,
                    // faults
                    badstates: badstates.size(),
                    badpointings: badpointings.size(),
                    badfreqs: badfreqs.size(),
                    badgains: badgains.size(),
                    badbeamshape: badbeamshape.size(),
                    // significant_faults: significant_faults,
                    fault_str: faults.shortstring.replaceAll(/\n\s*/, '|'),
                    // files
                    // files: fileStats.files.toString(),
                    num_data_files: fileStats.num_data_files,
                    num_data_files_archived: fileStats.num_data_files_archived,
                ]
                [obsid, summary]
            }

        // display wsSummary
        wsSummary.map { obsid, summary ->
                [
                    obsid,
                    summary.fail_reasons.join('|'),
                    summary.ra_pointing,
                    summary.dec_pointing,
                    summary.ra_phase_center,
                    summary.dec_phase_center,
                    summary.pointing,
                    summary.lst,

                    summary.freq_res,
                    summary.int_time,
                    summary.nscans,

                    summary.num_data_files,
                    summary.num_data_files_archived,

                    summary.dataquality,
                    summary.dataqualitycomment,

                    summary.n_dead_dipoles,
                    summary.dead_dipole_frac,
                    summary.n_bad_tiles,
                    summary.bad_tile_frac,
                    displayInts(summary.bad_tiles),

                    summary.iono_magnitude,
                    summary.iono_pca,
                    summary.iono_qa,

                    summary.badstates,
                    summary.badpointings,
                    summary.badfreqs,
                    summary.badgains,
                    summary.badbeamshape,
                    summary.fault_str,
                ].join("\t")
            }
            .collectFile(
                name: "ws_stats.tsv", newLine: true, sort: true,
                seed: [
                    "OBS", "FAIL REASON", "RA POINT", "DEC POINT", "RA PHASE", "DEC PHASE", "POINT", "LST DEG",
                    "FREQ RES", "TIME RES","N SCANS",
                    "N FILES", "N ARCHIVED",
                    "QUALITY", "QUALITY COMMENT",
                    "N DEAD DIPOLES", "DEAD DIPOLE FRAC", "N FLAG TILES", "FLAG TILES FRAC", "FLAG TILES",
                    "IONO MAG", "IONO PCA", "IONO QA",
                    "STATE FAULTS", "POINTING FAULTS", "FREQ FAULTS", "GAIN FAULTS", "BEAM FAULTS", "FAULT STR"
                ].join("\t"),
                storeDir: "${results_dir}",
            )
            // display output path and number of lines
            | view { [it, it.readLines().size()] }

        pass = wsSummary
            .filter { _, summary -> summary.fail_reasons == [] }
            .map { obsid, _ -> obsid }

        pass | wsMetafits & wsSkyMap

    emit:
        // channel of good obsids with their metafits: tuple(obsid, metafits)
        obsMetafits = wsMetafits.out

        // channel of video name and frames to convert
        frame = wsSkyMap.out
            .map { _, png -> ["skymap", png] }
            .groupTuple()
}

// ensure preprocessed uvfits are downloaded
workflow prep {
    take:
        // channel of obsids with their metafits: tuple(obsid, metafits)
        obsMetafits
    main:
        // download preprocessed uvfits
        obsMetafits
            .map { obsid, _ -> obsid }
            | asvoPrep
            | prepVisQA
        // uncomment this to skip prep and flag QA
        // emit:
        //     obsMetaVis = obsMetafits.join(asvoPrep.out)
        //     obsFlags = asvoPrep.out.map { obsid, _ -> [obsid, []]}

        def tileUpdates = file(params.tile_updates_path)
            .readLines()
            .collect { line ->
                def (firstObsid, lastObsid, tileIdxs) = line.split(',') + ['', '']
                [firstObsid as int, lastObsid as int, (tileIdxs as String).split("\\|").collect {it as Integer} ]
            }

        // collect prepVisQA results as .tsv
        prepVisQA.out
            // form row of tsv from json fields we care about
            .map { obsid, json ->
                def stats = parseJson(json);
                [
                    obsid,
                    stats.NANTS?:'',
                    stats.NTIMES?:'',
                    stats.NFREQS?:'',
                    stats.NPOLS?:'',
                    stats.XX?.MXRMS_AMP_ANT?:'',
                    stats.XX?.MNRMS_AMP_ANT?:'',
                    stats.XX?.MXRMS_AMP_FREQ?:'',
                    stats.XX?.MNRMS_AMP_FREQ?:'',
                    stats.XX?.NPOOR_ANTENNAS?:'',
                    displayInts(stats.XX?.POOR_ANTENNAS?:[]),
                    stats.YY?.MXRMS_AMP_ANT?:'',
                    stats.YY?.MNRMS_AMP_ANT?:'',
                    stats.YY?.MXRMS_AMP_FREQ?:'',
                    stats.YY?.MNRMS_AMP_FREQ?:'',
                    stats.YY?.NPOOR_ANTENNAS?:'',
                    displayInts(stats.YY?.POOR_ANTENNAS?:[]),
                ].join("\t")
            }
            .collectFile(
                name: "prepvis_metrics.tsv", newLine: true, sort: true,
                seed: [
                    "OBS",
                    "NANTS",
                    "NTIMES",
                    "NFREQS",
                    "NPOLS",
                    "XX MXRMS_AMP_ANT",
                    "XX MNRMS_AMP_ANT",
                    "XX MXRMS_AMP_FREQ",
                    "XX MNRMS_AMP_FREQ",
                    "XX NPOOR_ANTENNAS",
                    "XX POOR_ANTENNAS",
                    "YY MXRMS_AMP_ANT",
                    "YY MNRMS_AMP_ANT",
                    "YY MXRMS_AMP_FREQ",
                    "YY MNRMS_AMP_FREQ",
                    "YY NPOOR_ANTENNAS",
                    "YY POOR_ANTENNAS",
                ].join("\t"),
                storeDir: "${results_dir}"
            )
            // display output path and number of lines
            | view { [it, it.readLines().size()] }

        // plot prepvisQA
        prepVisQA.out | plotPrepVisQA

        // analyse flag occupancy
        obsMetafits.join(asvoPrep.out) | flagQA

        // collect flagQA results
        // TODO: collect sky_chans from flagQA
        def sky_chans = (131..154).collect { ch -> "$ch".toString() }
        flagQA.out
            // form row of tsv from json fields we care about
            .map { obsid, json ->
                def stats = jslurp.parse(json)
                def flagged_sky_chans = stats.flagged_sky_chans?:[]
                def chan_occupancy = sky_chans.collect { ch -> stats.channels?[ch]?.non_preflagged_bl_occupancy?:'' }
                ([
                    obsid,
                    stats.num_chans?:'',
                    displayInts(flagged_sky_chans),
                    displayInts(stats.flagged_fchan_idxs?:[]),
                    flagged_sky_chans.size/chan_occupancy.size,
                    stats.num_times?:'',
                    stats.num_baselines?:'',
                    displayInts(stats.flagged_timestep_idxs?:[]),
                    displayInts(stats.preflagged_ants?:[]),
                    stats.total_occupancy,
                    stats.total_non_preflagged_bl_occupancy
                ] + chan_occupancy).join("\t")
            }
            .collectFile(
                name: "occupancy.tsv", newLine: true, sort: true,
                seed: ([
                    "OBS",
                    "N CHAN", "FLAGGED SKY CHANS", "FLAGGED FINE CHANS", "FLAGGED SKY CHAN FRAC",
                    "N TIME","N BL", "FLAGGED TIMESTEPS", "FLAGGED ANTS",
                    "TOTAL OCCUPANCY", "NON PREFLAGGED"
                ] + sky_chans).join("\t"),
                storeDir: "${results_dir}"
            )
            // display output path and number of lines
            | view { [it, it.readLines().size()] }

    emit:
        // channel of obsids which pass the flag gate: tuple(obsid, metafits, uvfits)
        obsMetaVis = flagQA.out
            .filter { obsid, json ->
                jslurp.parse(json).total_occupancy < params.flag_occupancy_threshold
            }
            .map { obsid, _ -> obsid }
            .join(obsMetafits)
            .join(asvoPrep.out)
        // channel of extra flags for each obsid
        obsFlags = prepVisQA.out
            .join(flagQA.out)
            .map { obsid, prepJson, flagJson ->
                def flagStats = parseJson(flagJson);
                def flagAntennas = (flagStats.preflagged_ants?:[]) as Set
                def prepStats = parseJson(prepJson);
                def prepAntennas = ((prepStats.XX?.POOR_ANTENNAS?:[]) + (prepStats.YY?.POOR_ANTENNAS?:[])) as Set
                def manualAntennas = ([]) as Set
                tileUpdates.each {
                    def (firstObsid, lastObsid, tileIdxs, comment) = it
                    if (obsid as int >= firstObsid && obsid as int <= lastObsid) {
                        manualAntennas.addAll(tileIdxs)
                    }
                }
                [
                    obsid,
                    flagAntennas as int[],
                    prepAntennas as int[],
                    manualAntennas as int[],
                    (prepAntennas + manualAntennas - flagAntennas) as int[]
                ]
            }
        // channel of video name and frames to convert
        frame = plotPrepVisQA.out
            .map { _, png -> ["prepvisqa_rms", png] }
            .groupTuple()
}

workflow cal {
    take:
        // channel of metafits and preprocessed uvfits: tuple(obsid, metafits, uvfits, antFlags)
        obsMetaVis
    main:
        // hyperdrive di-calibrate on each obs
        obsMetaVis
            .map { def (obsids, metafits, uvfits, antFlags) = it
                [obsids, params.dical_args, metafits, uvfits, antFlags]
            }
            | hypCalSol

        // hyperdrive dical log analysis
        hypCalSol.out
            .transpose()
            .map { obsid, soln, diCalLog ->
                def name = soln.getBaseName().split('_')[3..-1].join('_')
                def (convergedDurationSec, convergedNumerator, convergedDenominator) = ['', '', '']
                if (!diCalLog.isEmpty()) {
                    def startTime = logDateFmt.parse((diCalLog.getText() =~ /([\d: -]+) INFO  hyperdrive di-calibrate/)[0][1])
                    def convergedMatch = (diCalLog.getText() =~ /([\d: -]+) INFO  All timesteps: (\d+)\/(\d+)/)[0]
                    if (startTime && convergedMatch) {
                        def convergedTime = logDateFmt.parse(convergedMatch[1])
                        convergedNumerator = convergedMatch[2] as int
                        convergedDenominator = convergedMatch[3] as int
                        convergedDurationSec = (convergedTime.time - startTime.time) / 1000
                    }
                }

                [obsid, name, convergedDurationSec, convergedNumerator, convergedDenominator].join("\t")
            }
            .collectFile(
                name: "cal_timings${params.cal_suffix}.tsv", newLine: true, sort: true,
                seed: [
                    "OBS", "CAL NAME", "CAL DUR", "CHS CONVERGED", "CHS TOTAL"
                ].join("\t"),
                storeDir: "${results_dir}"
            )

        // channel of solutions for each obsid: tuple(obsid, solutions)
        obsSolns = hypCalSol.out.map { obsid, solutions, _ -> [obsid, solutions] }
        // channel of individual dical solutions: tuple(obsid, name, soln)
        // - hypCalSol gives multiple solutions, transpose gives 1 tuple per solution.
        eachCal = obsSolns
            .transpose()
            .map { obsid, soln ->
                // give each calibration a name from basename of solution fits
                def name = soln.getBaseName().split('_')[3..-1].join('_')
                [obsid, name, soln]
            }

        // generate json from solns
        eachCal | solJson

        // channel of all dical (and polyfit) solutions: tuple(obsid, name, soln)
        if (!params.nopoly) {
            // do polyFit if nopoly is not set
            eachCal | polyFit
            // channel of individual polyfit solutions: tuple(obsid, name, soln)
            eachPolyCal = polyFit.out.map { obsid, name, soln, _ -> [obsid, name, soln] }
            allCal = eachCal.mix(eachPolyCal)
        } else {
            allCal = eachCal
        }

        // collect solJson results as .tsv
        solJson.out
            // form row of tsv from json fields we care about
            .map { obsid, name, json ->
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
                name: "cal_results${params.cal_suffix}.tsv", newLine: true, sort: true,
                seed: [
                    "OBS", "CAL NAME", "NAN FRAC", "RESULTS BY CH"
                ].join("\t"),
                storeDir: "${results_dir}"
            )
            // display output path and number of lines
            | view { [it, it.readLines().size()] }

        // channel of metafits for each obsid: tuple(obsid, metafits)
        obsMetafits = obsMetaVis.map { obsid, metafits, _, __ -> [obsid, metafits] }
        // calibration QA and plot solutions
        obsMetafits
            // join with solutions from eachCal and eachPolyCal
            .cross(allCal)
            // .map { (obsid, metafits), (_, name, soln) ->
            .map { def (obsid, metafits) = it[0]; def (_, name, soln) = it[1]
                [obsid, name, metafits, soln]
            }
            | (plotSols & calQA)

        // plot each calQA result
        calQA.out | plotCalQA

        // plot aggregate calQA results
        calQA.out
            .map { _, __, json -> json }
            .collect(sort: true)
            | plotCalJsons

        // collect calQA results as .tsv
        calQA.out
            // form row of tsv from json fields we care about
            .map { obsid, name, json ->
                def stats = parseJson(json)
                def convg_var = stats.CONVERGENCE_VAR
                [
                    obsid,
                    name,
                    stats.STATUS?:'',
                    (stats.UNUSED_BLS?:0) / 100,
                    (stats.UNUSED_CHS?:0) / 100,
                    (stats.UNUSED_ANTS?:0) / 100,
                    (stats.NON_CONVERGED_CHS?:0) / 100,
                    convg_var,
                    (convg_var instanceof BigDecimal) ? convg_var * 1e14 : "NaN",
                    stats.XX?.SKEWNESS_UVCUT?:'',
                    stats.XX?.RMS_AMPVAR_ANT?:'',
                    stats.XX?.RMS_AMPVAR_FREQ?:'',
                    stats.XX?.DFFT_POWER?:'',
                    stats.XX?.RECEIVER_CHISQVAR?:'',
                    stats.YY?.SKEWNESS_UVCUT?:'',
                    stats.YY?.RMS_AMPVAR_ANT?:'',
                    stats.YY?.RMS_AMPVAR_FREQ?:'',
                    stats.YY?.DFFT_POWER?:'',
                    stats.YY?.RECEIVER_CHISQVAR?:'',
                    stats.FAILURE_REASON,
                ].join("\t")
            }
            .collectFile(
                name: "cal_metrics${params.cal_suffix}.tsv", newLine: true, sort: true,
                seed: [
                    "OBS", "CAL NAME", "STATUS",
                    "BL UNUSED FRAC", "CH UNUSED FRAC", "ANT UNUSED FRAC",
                    "NON CONVG CHS FRAC",
                    "CONVG VAR",
                    "CONVG VAR e14",
                    "XX SKEW UVCUT",
                    "XX RMS AMPVAR ANT",
                    "XX RMS AMPVAR FREQ",
                    "XX DFFT POWER",
                    "XX RECEIVER CHISQVAR",
                    "YY SKEW UVCUT",
                    "YY RMS AMPVAR ANT",
                    "YY RMS AMPVAR FREQ",
                    "YY DFFT POWER",
                    "YY RECEIVER CHISQVAR",
                    "FAILURE_REASON"
                ].join("\t"),
                storeDir: "${results_dir}"
            )
            // display output path and number of lines
            | view { [it, it.readLines().size()] }

    // channel of obsids and names that pass qa. tuple(obsid, name)
    // - take tuple(obsid, cal_name, json) from calQA.out
    // - filter on json.STATUS == "PASS"
    // - take obsid and name
    obsidNamePass = calQA.out
        .filter { _, __, json -> parseJson(json).STATUS == "PASS" }
        .map { obsid, name, _ -> [obsid, name] }

    emit:
        // channel of calibration solutions that pass qa. tuple(obsid, name, cal)
        // - take tuple(obsid, cal_name, soln) from allCal
        // - match with obsidNamePass on (obsid, cal_name)
        passCal = allCal
            .join(obsidNamePass, by: [0,1] )
        // channel of files to archive, and their buckets
        archive = calQA.out.map { _, __, json -> ["calqa", json]}
        // channel of video name and frames to convert
        frame = plotCalQA.out.mix(plotSols.out)
            .flatMap { _, name, pngs ->
                pngs.collect { png ->
                    def png_name = png.getSimpleName().split('_')[-1]
                    ["calqa_${name}_${png_name}", png]
                }
            }
            .groupTuple()
}

// process uvfits visibilities
workflow uvfits {
    take:
        obsNameUvfits
    main:
        // vis QA
        obsNameUvfits
            .filter { _, name, __ -> !name.toString().startsWith("sub_") }
            | visQA

        // collect visQA results as .tsv
        visQA.out
            // form row of tsv from json fields we care about
            .map { obsid, name, json ->
                def stats = parseJson(json);
                [
                    obsid,
                    name,
                    stats.AUTOS?.XX?.MXRMS_AMP_ANT?:'',
                    stats.AUTOS?.XX?.MNRMS_AMP_ANT?:'',
                    stats.AUTOS?.XX?.MXRMS_AMP_FREQ?:'',
                    stats.AUTOS?.XX?.MNRMS_AMP_FREQ?:'',
                    stats.AUTOS?.XX?.NPOOR_ANTENNAS?:'',
                    stats.AUTOS?.XX?.POOR_ANTENNAS?:'',
                    stats.AUTOS?.XX?.NPOOR_CHANNELS?:'',
                    stats.AUTOS?.XX?.POOR_CHANNELS?:'',
                    stats.AUTOS?.YY?.MXRMS_AMP_ANT?:'',
                    stats.AUTOS?.YY?.MNRMS_AMP_ANT?:'',
                    stats.AUTOS?.YY?.MXRMS_AMP_FREQ?:'',
                    stats.AUTOS?.YY?.MNRMS_AMP_FREQ?:'',
                    stats.AUTOS?.YY?.NPOOR_ANTENNAS?:'',
                    stats.AUTOS?.YY?.POOR_ANTENNAS?:'',
                    stats.AUTOS?.YY?.NPOOR_CHANNELS?:'',
                    stats.AUTOS?.YY?.POOR_CHANNELS?:'',
                    stats.REDUNDANT?.XX?.NPOOR_BLS?:'',
                    stats.REDUNDANT?.XX?.POOR_BLS?:'',
                    stats.REDUNDANT?.YY?.NPOOR_BLS?:'',
                    stats.REDUNDANT?.YY?.POOR_BLS?:'',
                ].join("\t")
            }
            .collectFile(
                name: "vis_metrics.tsv", newLine: true, sort: true,
                seed: [
                    "OBS", "VIS NAME",
                    "A:XX MXRMS_AMP_ANT",
                    "A:XX MNRMS_AMP_ANT",
                    "A:XX MXRMS_AMP_FREQ",
                    "A:XX MNRMS_AMP_FREQ",
                    "A:XX NPOOR_ANTENNAS",
                    "A:XX POOR_ANTENNAS",
                    "A:XX NPOOR_CHANNELS",
                    "A:XX POOR_CHANNELS",
                    "A:YY MXRMS_AMP_ANT",
                    "A:YY MNRMS_AMP_ANT",
                    "A:YY MXRMS_AMP_FREQ",
                    "A:YY MNRMS_AMP_FREQ",
                    "A:YY NPOOR_ANTENNAS",
                    "A:YY POOR_ANTENNAS",
                    "A:YY NPOOR_CHANNELS",
                    "A:YY POOR_CHANNELS",
                    "R:XX NPOOR_BLS",
                    "R:XX POOR_BLS",
                    "R:YY NPOOR_BLS",
                    "R:YY POOR_BLS",
                ].join("\t"),
                storeDir: "${results_dir}"
            )
            // display output path and number of lines
            | view { [it, it.readLines().size()] }

        // TODO: plotVisQA isn't working
        // channel.from([]) | plotVisQA
        visQA.out | plotVisQA

        // ps_metrics
        if (!params.nopsmetrics) {
            obsNameUvfits | psMetrics

            // collect psMetrics as a .dat
            psMetrics.out
                // read the content of each ps_metrics file including the trailing newline
                .map { obsid, vis_name, dat, _ -> dat.getText() }
                .collectFile(
                    name: "ps_metrics.dat",
                    storeDir: "${results_dir}"
                )
                // display output path and number of lines
                | view { [it, it.readLines().size()] }

            // collect psMetrics as a .tsv
            psMetrics.out
                // form each row of tsv
                .map { obsid, vis_name, dat, _ ->
                    def dat_values = dat.getText().split('\n')[0].split(' ')[1..-1]
                    ([obsid, vis_name] + dat_values).join("\t")
                }
                .collectFile(
                    name: "ps_metrics.tsv", newLine: true, sort: true,
                    seed: [
                        "OBS", "CAL NAME", "P_WEDGE", "NUM_CELLS", "P_WINDOW", "NUM_CELLS",
                        "P_ALL", "D3"
                    ].join("\t"),
                    storeDir: "${results_dir}"
                )
                // display output path and number of lines
                | view { [it, it.readLines().size()] }
        }

    emit:
        // channel of files to archive, and their buckets
        archive = visQA.out.map { _, __, json -> ["visqa", json]}
        // channel of video name and frames to convert
        frame = plotVisQA.out
            .map { _, name, png -> ["visqa_${name}_rms", png] }
            .groupTuple()
}

// image measurement sets and QA images
workflow img {
    take:
        obsNameMS
    main:
        // wsclean: make dirty images
        obsNameMS | wscleanDirty

        // imgQA for all groups of images
        wscleanDirty.out
            .map { obsid, name, imgs, _ -> [obsid, name, imgs] }
            | (imgQA & polComp)

        // collect imgQA results as .tsv
        imgQA.out
            // form row of tsv from json fields we care about
            .map { obsid, name, json ->
                def stats = parseJson(json)
                [
                    obsid,
                    name,
                    stats.XX?.RMS_ALL, stats.XX?.RMS_BOX, stats.XX?.PKS0023_026?.PEAK_FLUX, stats.XX?.PKS0023_026?.INT_FLUX,
                    stats.YY?.RMS_ALL, stats.YY?.RMS_BOX, stats.YY?.PKS0023_026?.PEAK_FLUX, stats.YY?.PKS0023_026?.INT_FLUX,
                    stats.V?.RMS_ALL, stats.V?.RMS_BOX, stats.V?.PKS0023_026?.PEAK_FLUX, stats.V?.PKS0023_026?.INT_FLUX,
                    stats.V_XX?.RMS_RATIO, stats.V_XX?.RMS_RATIO_BOX
                ].join("\t")
            }
            .collectFile(
                name: "img_metrics.tsv", newLine: true, sort: true,
                seed: [
                    "OBS", "IMG NAME",
                    "XX RMS ALL", "XX RMS BOX", "XX PKS0023_026 PEAK", "XX PKS0023_026 INT",
                    "YY RMS ALL", "YY RMS BOX", "YY PKS0023_026 PEAK", "YY PKS0023_026 INT",
                    "V RMS ALL","V RMS BOX", "V PKS0023_026 PEAK", "V PKS0023_026 INT" ,
                    "V:XX RMS RATIO", "V:XX RMS RATIO BOX"
                ].join("\t"),
                storeDir: "${results_dir}"
            )
            // display output path and number of lines
            | view { [it, it.readLines().size()] }

    emit:
        // channel of files to archive, and their buckets
        archive = wscleanDirty.out
            .transpose()
            .map { _, __, img, ___ -> ["img", img]}
            .mix( imgQA.out.map { _, __, json -> ["imgqa", json]} )
        // channel of video name and frames to convert
        frame = polComp.out
            .map { _, name, png -> ["imgqa_${name}_polcomp", png] }
            .groupTuple()
}

// separate entrypoint for moving unfiltered preprocessed uvfits from asvo accacia to mwaeor accacia
workflow archivePrep {
    obsids = channel.of(obsids_file)
        .splitCsv()
        .flatten()
        .filter { line -> !line.startsWith('#') }

    obsids | asvoPrep

    if (params.archive) {
        asvoPrep.out.map { _, vis -> ["prep", vis] } | archive
    }
}

workflow {
    // get obsids from csv
    obsids = channel.of(obsids_file)
        .splitCsv()
        .flatten()
        .filter { line -> !line.startsWith('#') }

    // analyse obsids with web services
    obsids | ws

    ws.out.obsMetafits
        .map { obsid, _ -> obsid }
        .collectFile(
            name: "filtered_obsids.csv", newLine: true, sort: true,
            storeDir: "${results_dir}"
        )
        | view { [it, it.readLines().size()] }

    // download preprocessed, unless noprep is set
    if (params.noprep) { println("params.noprep set, exiting before prep."); return }
    ws.out.obsMetafits | prep

    // channel of obsids that pass the flag gate
    prep.out.obsMetaVis
        .map { obsid, _, __ -> obsid }
        .collectFile(
            name: "unflagged_obsids.csv", newLine: true, sort: true,
            storeDir: "${results_dir}"
        )
        | view { [it, it.readLines().size()] }

    prep.out.obsFlags
        .filter { obsid, flagAnts, prepAnts, manualAnts, newAnts -> newAnts.size() > 0 }
        .map { obsid, flagAnts, prepAnts, manualAnts, newAnts ->
            ([
                obsid,
                displayInts(flagAnts),
                displayInts(prepAnts),
                displayInts(manualAnts),
                displayInts(newAnts),
            ]).join("\t")
        }
        .collectFile(
            name: "tile_flags.tsv", newLine: true, sort: true,
            seed: ([ "OBS", "ORIGINAL ANTS", "PREP ANTS", "MANUAL ANTS", "NEW ANTS" ]).join("\t"),
            storeDir: "${results_dir}"
        )
        // display output path and number of lines
        | view { [it, it.readLines().size()] }

    // calibrate each obs that passes flag gate unless nocal is set:
    if (params.nocal) {
        channel.from([]) | cal
    } else {
        prep.out.obsMetaVis.join(
            prep.out.obsFlags
                // .map {obsid, flagAnts, prepAnts, manualAnts, newAnts -> [obsid, []]}
                .map {obsid, flagAnts, prepAnts, manualAnts, newAnts -> [obsid, newAnts]}
        ) | cal
    }

    // channel of arguments for hypApply{UV,MS}
    // - take tuple(obsid, cal_name, soln) from obsNameSoln
    // - lookup [apply_name, apply_args] from params.apply_args on cal_name
    // - match with tuple(obsid, metafits, prepUVFits) by obsid
    if (params.noapply) {
        allApply = channel.from([])
    } else {
        allApply = prep.out.obsMetaVis
            .cross(
                cal.out.passCal
                    .filter { obsid, cal_name, _ -> params.apply_args[cal_name] != null }
                    .map { obsid, cal_name, soln ->
                        def (apply_name, apply_args) = params.apply_args[cal_name]
                        [obsid, soln, cal_name, apply_name, apply_args]
                    }
            )
            .map {
                def (obsid, metafits, prepUVFits) = it[0]
                def (_, soln, cal_name, apply_name, apply_args) = it[1]
                [obsid, metafits, prepUVFits, soln, cal_name, apply_name, apply_args]
            }
    }

    if (params.nouv) {
        channel.from([]) | hypApplyUV
    } else {
        allApply | hypApplyUV
    }

    // channel of calibrated (or subtracted) uvfits: tuple(obsid, name, uvfits)
    if (params.nosub) {
        obsNameUvfits = hypApplyUV.out
            .map { obsid, name, vis, _ -> [obsid, name, vis] }
    } else {
        // get subtracted uvfits vis
        prep.out.obsMetaVis
            .cross(hypApplyUV.out)
            .map {
                def (obsid, metafits, _) = it[0]
                def (__, name, vis) = it[1];
                [obsid, name, metafits, vis]
            }
            | hypSubUV

        obsNameUvfits = hypApplyUV.out
            .mix(hypSubUV.out)
            .map { obsid, name, vis, _ -> [obsid, name, vis] }
    }

    // QA uvfits visibilities
    obsNameUvfits | uvfits

    if (params.noms) {
        channel.from([]) | hypApplyMS
    } else {
        allApply | hypApplyMS
    }

    // channel of calibrated (or subtracted) measurement sets: tuple(obsid, name, ms)
    if (params.nosub) {
        obsNameMS = hypApplyMS.out
            .map { obsid, name, vis, _ -> [obsid, name, vis] }
    } else {
        // get subtracted ms vis
        prep.out.obsMetaVis
            .cross(hypApplyMS.out)
            .map {
                def (obsid, metafits, _) = it[0]
                def (__, name, vis) = it[1];
                [obsid, name, metafits, vis]
            }
            | hypSubMS

        obsNameMS = hypApplyMS.out
            .mix(hypSubMS.out)
            .map { obsid, name, vis, _ -> [obsid, name, vis] }
    }

    // image and qa measurementsets unless noimage is set
    if (params.noimg) {
        channel.from([]) | img
    } else {
        obsNameMS | img
    }

    // obsNameMS | aoQuality

    // aoQuality.out.view {it}
    //     .map { obsid, name,  -> [obsid, name, ms] } | img
    // }

    if (params.archive) {
        if (params.archive_prep) {
            prep_archive = prep.out.obsMetaVis.map { _, __, vis -> ["prep", vis] }
        } else {
            prep_archive = channel.from([])
        }
        if (params.archive_uvfits) {
            vis_archive = obsNameUvfits.map { _, __, vis -> ["uvfits", vis]}
        } else {
            vis_archive = channel.from([])
        }
        prep_archive.mix(vis_archive)
            .mix(cal.out.passCal.map { _, __, soln -> ["soln", soln] })
            .mix(cal.out.archive)
            .mix(uvfits.out.archive)
            .mix(img.out.archive)
            | archive
    }

    // make videos
    ws.out.frame
        .mix(prep.out.frame)
        .mix(cal.out.frame)
        .mix(uvfits.out.frame)
        .mix(img.out.frame)
        .map { name, frames -> [name, frames.sort()] }
        | ffmpeg
        | view { [it, it.size()] }
}