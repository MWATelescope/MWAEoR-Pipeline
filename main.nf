#!/usr/bin/env nextflow
nextflow.enable.dsl=2

def obsids_file = file(params.obsids_path)
if (params.obsids_suffix) {
    obsids_file = file("${obsids_file.parent}/${obsids_file.baseName}${params.obsids_suffix}.${obsids_file.extension}")
}
def results_dir = "${params.resultsdir}/results${params.obsids_suffix}${params.result_suffix}"

// whether imaging is configured for multiple channels
def img_channels_out = (params.img_channels_out instanceof String ? \
    params.img_channels_out.split(' ')[0] :\
    params.img_channels_out)
def multichannel = (img_channels_out as int > 1)
// whether imaging is configured for multiple intervals
def multiinterval = (params.img_intervals_out as int > 1) || params.img_split_intervals

// download observation metadata from webservices in json format
process wsMeta {
    input:
    val(obsid)
    output:
    tuple val(obsid), path(wsmeta), path(wsfiles)

    maxForks 1

    // persist results in outdir, process will be skipped if files already present.
    storeDir "${params.outdir}/meta"
    // tag to identify job in squeue and nf logs
    tag "${obsid}"

    script:
    wsmeta = "${obsid}_wsmeta.json"
    wsfiles = "${obsid}_files.json"
    """
    #!/bin/bash -eux
    ${params.proxy_prelude} # ensure proxy is set if needed
    wget -O "${wsmeta}" "http://ws.mwatelescope.org/metadata/obs?obs_id=${obsid}&extended=1&dict=1"
    wget -O "${wsfiles}" "http://ws.mwatelescope.org/metadata/data_ready?obs_id=${obsid}"
    """
}

// download observation metadata from webservices in metafits format
process wsMetafits {
    input:
    val(obsid)
    output:
    tuple val(obsid), path(metafits)

    maxForks 1

    // persist results in outdir, process will be skipped if files already present.
    storeDir "${params.outdir}/${obsid}/raw"
    // tag to identify job in squeue and nf logs
    tag "${obsid}"

    script:
    metafits = "${obsid}.metafits"
    """
    #!/bin/bash -eux
    ${params.proxy_prelude} # ensure proxy is set if needed
    wget -O "${metafits}" "http://ws.mwatelescope.org/metadata/fits?obs_id=${obsid}&include_ppds=${params.metafits_incl_ppds}"
    """
}

process wsSkyMap {
    input:
    val(obsid)
    output:
    tuple val(obsid), path(skymap)

    maxForks 1

    // persist results in outdir, process will be skipped if files already present.
    storeDir "${params.outdir}/${obsid}/meta"
    // tag to identify job in squeue and nf logs
    tag "${obsid}"

    script:
    skymap = "${obsid}_skymap.png"
    """
    #!/bin/bash -eux
    ${params.proxy_prelude} # ensure proxy is set if needed
    wget -O "${skymap}" "http://ws.mwatelescope.org/observation/skymap/?obs_id=${obsid}"
    """
}

process wsPPDs {
    input:
    val(obsid)
    output:
    tuple val(obsid), path(ppds)

    maxForks 1

    // persist results in outdir, process will be skipped if files already present.
    storeDir "${params.outdir}/${obsid}/meta"
    // tag to identify job in squeue and nf logs
    tag "${obsid}"

    script:
    ppds = "${obsid}_ppds.png"
    """
    #!/bin/bash -eux
    ${params.proxy_prelude} # ensure proxy is set if needed
    curl "http://ws.mwatelescope.org/observation/ppds/?replot=replot&obs_id=${obsid}&merge=on&corgains=on&adu=on&plotscale=1.0"
    wget -O "${ppds}" "http://ws.mwatelescope.org/observation/powerplot/?obs_id=${obsid}&group=False&plotscale=1.0&merge=1&corgains=1&adu=1&waterfall=False"
    """
}

// download preprocessed files from asvo
process asvoPrep {
    input:
    val obsid
    output:
    tuple val(obsid), path(uvfits)

    storeDir "${params.outdir}/${obsid}/prep"
    // tag to identify job in squeue and nf logs
    tag "${obsid}"

    if (params.pullPrep) {
        label "rclone"
    }

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
    uvfits = "birli_${obsid}_${params.prep_time_res_s}s_${params.prep_freq_res_khz}kHz.uvfits"
    """
    #!/bin/bash -eux

    export MWA_ASVO_API_KEY="${params.asvo_api_key}"

    ${params.proxy_prelude} # ensure proxy is set if needed

    """ + ( params.pullPrep ? """
    # download if available in accacia
    rclone copy "${params.bucket_prefix}.prep/${uvfits}" .
    if [ -f "${uvfits}" ]; then
        exit 0
    fi
    """ : "" ) + """

    # submit a job to ASVO, suppress failure if a job already exists.
    ${params.giant_squid} submit-conv -v \
        -p avg_time_res=${params.prep_time_res_s},avg_freq_res=${params.prep_freq_res_khz},output=uvfits \
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
        mv "${obsid}.uvfits" "${uvfits}"
        exit 0 # success
    fi
    echo "no ready jobs"
    exit 75 # temporary
    """
}

// QA tasks for flags.
process flagQA {
    input:
    tuple val(obsid), path(metafits), path(uvfits)
    output:
    tuple val(obsid), path(metrics)

    storeDir "${params.outdir}/${obsid}/prep"

    tag "${obsid}"

    label 'python'

    script:
    metrics = "${obsid}_occupancy.json"
    template "flagqa.py"
}

process ssins {
    input:
    tuple val(obsid), path(uvfits)
    output:
    tuple val(obsid),
        path("_SSINS_mask.h5"),
        path("{autos,cross,flagged}_SSINS*.png"),
        path("ssins_occ.json")
        // path(ssins_uvfits)
        // todo: path("ssins_VDH.png"),
        // todo: path("match_events.json"),

    storeDir "${params.outdir}/${obsid}/prep"

    tag "${obsid}"

    label "ssins"

    // errorStrategy "terminate"

    script:
    ssins_uvfits = "ssins_${obsid}_${params.prep_time_res_s}s_${params.prep_freq_res_khz}kHz.uvfits"
    guard_width = params.prep_freq_res_khz * 500
    title = "${obsid}"
    template "ssins.py"
}

process autoplot {
    input:
    tuple val(obsid), path(metafits), path(uvfits)
    output:
    tuple val(obsid), path(autoplot)

    storeDir "${params.outdir}/${obsid}/vis_qa"

    tag "${obsid}"

    label 'python'

    script:
    autoplot = "${obsid}_autoplot.png"
    title = "${obsid}"
    args = "${params.autoplot_args}"
    template "autoplot.py"
}

// make a reduced ao-format sourcelist with hyperdrive
process hypSrclistAO {
    input:
    tuple val(obsid), path(metafits)
    output:
    tuple val(obsid), path(reduced)

    storeDir "${params.outdir}/${obsid}/cal${params.cal_suffix}"

    tag "${obsid}"
    label "hyperdrive_cpu"
    time 15.minute

    script:
    reduced = "${obsid}_reduced_n${params.sub_nsrcs}.txt"
    """
    #!/bin/bash -eux

    # Reduce a sky-model source list to the top N brightest sources, given pointing information
    ${params.hyperdrive} srclist-by-beam \
        --metafits "${metafits}" \
        --number ${params.sub_nsrcs} \
        --beam-file "${params.beam_path}" \
        -o ao \
        "${params.sourcelist}" "${reduced}"
    """
}

// make a reduced yaml-format sourcelist with hyperdrive
process hypSrclistYaml {
    input:
    tuple val(obsid), path(metafits)
    output:
    tuple val(obsid), path(reduced)

    storeDir "${params.outdir}/${obsid}/cal${params.cal_suffix}"

    tag "${obsid}"
    label "hyperdrive_cpu"
    time 15.minute

    script:
    reduced = "${obsid}_reduced_n${params.sub_nsrcs}.yaml"
    """
    #!/bin/bash -eux

    # Reduce a sky-model source list to the top N brightest sources, given pointing information
    ${params.hyperdrive} srclist-by-beam \
        --metafits "${metafits}" \
        --number ${params.sub_nsrcs} \
        --beam-file "${params.beam_path}" \
        -o hyperdrive \
        "${params.sourcelist}" "${reduced}"
    """
}

// cluster a reduced sourcelist with mwa-reduce
process rexCluster {
    input:
    tuple val(obsid), path(reduced)
    output:
    tuple val(obsid), path(cluster)

    storeDir "${params.outdir}/${obsid}/cal${params.cal_suffix}"

    tag "${obsid}"
    label "mwa_reduce"

    script:
    cluster = "${obsid}_cluster_n${params.sub_nsrcs}_i${params.ionosub_nsrcs}.txt"
    nclusters = params.sub_nsrcs / params.ionosub_nsrcs
    """
    #!/bin/bash -eux

    ${params.mwa_reduce} cluster \
        "${reduced}" \
        "${cluster}" \
        ${nclusters}
    """
}

def groovy2bashAssocArray(map, name) {
    // size = map.collect { _, v -> v.size() }.max()
    "declare -A ${name}=(" + map.collect { k, v -> "[${k}]=\"${v}\"".toString() }.join(" ") + ")".toString()
}

// calibrate with mwa reduce
process rexCalSol {
    input:
    tuple val(obsid), val(dical_args), path(metafits), path(vis), val(tile_flags)
    output:
    tuple val(obsid), path("rex_soln_${obsid}_${name_glob}.bin"), path("rex_di-cal_${obsid}_${name_glob}.log")

    storeDir "${params.outdir}/${obsid}/cal${params.cal_suffix}"

    tag "${obsid}"

    label "mwa_reduce"

    script:
    dical_names = dical_args.keySet().collect()
    para = dical_names.size() > 1
    name_glob = para ? "{" + dical_names.join(',') + "}" : dical_names[0]
    flag_args = tile_flags.size() > 0 ? "--tile-flags ${tile_flags.join(' ')}" : ""
    vis_ms = "${uvfits.baseName}.ms"
    """
    #!/bin/bash -eux
    """ + groovy2bashAssocArray(dical_args, "dical_args") + """
    ${params.casa} -c "importuvfits('${vis}', '${vis_ms}')"
    singularity exec ${params.cotter_sif} fixmwams vis.ms ${metafits}

    for name in \${!dical_args[@]}; do
        export soln_name="rex_soln_${obsid}_\${name}.bin"
        export log_name="rex_di-cal_${obsid}_\${name}.log"
        export args=\${dical_args[\$name]}
        ${params.mwa_reduce} calibrate \${args} \
            -applybeam -mwa-path /astro/mwaeor/jline/software \
            -m "${params.sourcelist}" \
            -i 50 \
            -a 1e-4 1e-8 \
            ${vis_ms} \
            \${soln_name} | tee \${log_name}
    """
}

// ensure calibration solutions are present, or calibrate prep with hyperdrive
// do multiple calibration solutions for each obs, depending on dical_args
process hypCalSol {
    input:
    tuple val(obsid), val(dical_args), path(metafits), path(uvfits), val(tile_flags)
    output:
    tuple val(obsid), path("hyp_soln_${obsid}_${name_glob}.fits"), path("hyp_di-cal_${obsid}_${name_glob}.log")
    // todo: model subtract: path("hyp_model_${obsid}_${name_glob}.uvfits")

    storeDir "${params.outdir}/${obsid}/cal${params.cal_suffix}"

    tag "${obsid}"

    // label jobs that need a bigger gpu allocation
    label "hyperdrive"
    if (params.pullCalSol) {
        label "rclone"
    }

    script:
    dical_names = dical_args.keySet().collect()
    para = dical_names.size() > 1
    name_glob = para ? "{" + dical_names.join(',') + "}" : dical_names[0]
    flag_args = tile_flags.size() > 0 ? "--tile-flags ${tile_flags.join(' ')}" : ""

    """
    #!/bin/bash -eux
    """ + (para ? "export CUDA_VISIBLE_DEVICES=0" : "") + """
    export num_gpus="\$(nvidia-smi -L | wc -l)"
    if [ \$num_gpus -eq 0 ]; then
        echo "no gpus found"
        exit 1
    fi
    if [ \$num_gpus -ne ${params.num_gpus} ]; then
        echo "warning: expected \$num_gpus to be ${params.num_gpus}"
    fi
    """ + groovy2bashAssocArray(dical_args, "dical_args") + """
    if [ \${#dical_args[@]} -eq 0 ]; then
        echo "no dical args"
        exit 0
    fi
    if [ \${#dical_args[@]} -gt \$num_gpus ]; then
        echo "warning: more dical args than gpus";
    fi
    for name in \${!dical_args[@]}; do
        export soln_name="hyp_soln_${obsid}_\${name}.fits"
        export log_name="hyp_di-cal_${obsid}_\${name}.log"
    """ + (params.pullCalSol ? """
        # download if available in accacia
        rclone copy "${params.bucket_prefix}.soln/\${soln_name}" .
        if [ -f "\${soln_name}" ]; then
            touch \${log_name}
            continue
        fi
    """ : "") + """

        args=\${dical_args[\$name]}
        # hyperdrive di-cal backgrounded in a subshell
        (
            ${params.hyperdrive} di-calibrate \${args} \
                --data "${metafits}" "${uvfits}" \
                --beam "${params.beam_path}" \
                --source-list "${params.sourcelist}" \
                --outputs \$soln_name \
                ${flag_args} \
                | tee \$log_name
        ) &
        # TODO: model subtract: --model-filenames "hyp_model_${obsid}_\${name}.uvfits"
        # increment the target device mod num_gpus
        """ + (para ? "CUDA_VISIBLE_DEVICES=\$(( CUDA_VISIBLE_DEVICES+1 % \${num_gpus} ))" : "") + """
    done

    # wait for all the background jobs to finish
    export result=0
    wait \$(jobs -rp) || result=\$?

    # print out important info from log
    if [ \$(ls *.log 2> /dev/null | wc -l) -gt 0 ]; then
        grep -iE "err|warn|hyperdrive|chanblocks|reading|writing|flagged" *.log || echo no warnings
    fi

    exit \$result
    """
}

// fit polynomial to calibration solution
// meta is a map containing calibration solution metadata, fields:
// - name: name of the calibration (dical_name)
// - cal_prog: name of the calibration program (hyp or rex)
process polyFit {
    input:
    tuple val(obsid), val(meta), path(soln)
    output:
    tuple val(obsid), val(meta), path(poly_soln), path(logs)

    storeDir "${params.outdir}/${obsid}/cal${params.cal_suffix}"
    tag "${obsid}.${old_name}"

    label 'python'

    script:
    old_name = meta.name
    meta.name = "poly_${meta.name}"
    poly_soln = "${meta.cal_prog}_soln_${obsid}_${meta.name}.fits"
    logs = "polyfit_${obsid}_${meta.name}.log"
    """
    #!/bin/bash -eux
    run_polyfit.py "${soln}" --outfile "${poly_soln}" | tee "${logs}"
    """
}

// why this can't be a single process:
// - if nouv or noms changes, then this whole thing need to be re-run
// process hypApply {
//     input:
//     tuple val(obsid), path(metafits), path(uvfits), \
//         path("hyp_soln_${obsid}_${cal_prog}.fits"), \
//         val(cal_name), val(apply_name), val(apply_args)
//     output:
//     tuple val(obsid), val(name), \
//         path(params.nouv ? ".fakeuv" : "hyp_${obsid}_${name}.uvfits"), \
//         path(params.noms ? ".fakems" : "hyp_${obsid}_${name}.ms"), \
//         path("hyp_apply_${name}.log")
// }

// ensure calibrated uvfits are present, or apply calsols to prep uvfits with hyperdrive
// meta is a map containing calibration solution metadata, fields:
// - name: unique name of the calibration ((poly_)?$dical_name)
// - apply_name: name of the apply arguments used
// - time_res: time resolution to average when applying (seconds)
// - freq_res: frequency resolution to average when applying (kHz)
// - apply_args: any extra args for apply
process hypApplyUV {
    input:
    tuple val(obsid), val(meta), path(metafits), path(vis), path(soln)
    output:
    tuple val(obsid), val(meta), path(cal_vis), path(logs)

    storeDir "${params.outdir}/${obsid}/cal${params.cal_suffix}"

    tag "${obsid}.${meta.name}"
    label "cpu"
    label "hyperdrive"

    script:
    meta.name = "${meta.name}_${meta.apply_name}"
    cal_vis = "hyp_${obsid}_${meta.name}.uvfits"
    logs = "hyp_apply_${meta.name}.log"
    """
    ${params.hyperdrive} solutions-apply ${meta.apply_args} \
        --time-average=${meta.time_res}s \
        --freq-average=${meta.freq_res}kHz \
        --data "${metafits}" "${vis}" \
        --solutions "${soln}" \
        --outputs "${cal_vis}" \
        ${meta.nodut1 ? "--ignore-dut1" : ""} \
        | tee "${logs}"
    """
}

// ensure calibrated ms are present, or apply calsols to prep uvfits with hyperdrive
// meta is a map containing calibration solution metadata, fields:
// - name: unique name of the calibration ((poly_)?$dical_name)
// - apply_name: name of the apply arguments used
// - time_res: time resolution to average when applying (seconds)
// - freq_res: frequency resolution to average when applying (kHz)
// - apply_args: any extra args for apply
process hypApplyMS {
    input:
    tuple val(obsid), val(meta), path(metafits), path(vis), path(soln)
    output:
    tuple val(obsid), val(meta), path(cal_vis), path(logs)

    storeDir "${params.outdir}/${obsid}/cal${params.cal_suffix}"
    // storeDir "/data/curtin_mwaeor/FRB_hopper/"

    tag "${obsid}.${meta.name}"
    label "cpu"
    label "hyperdrive"

    script:
    meta.name = "${meta.name}_${meta.apply_name}"
    cal_vis = "hyp_${obsid}_${meta.name}.ms"
    logs = "hyp_apply_${meta.name}_ms.log"
    """
    # hyperdrive solutions apply ms
    ${params.hyperdrive} solutions-apply ${meta.apply_args} \
        --time-average=${meta.time_res}s \
        --freq-average=${meta.freq_res}kHz \
        --data "${metafits}" "${vis}" \
        --solutions "${soln}" \
        --outputs "${cal_vis}" \
        ${params.nodut1 ? "--ignore-dut1" : ""} \
        | tee "${logs}"
    """
}

// ensure subtracted uvfits are present, or subtract srclist from cal vis with hyperdrive
// meta is a map containing visibility metadata, fields:
// - name: unique name of the calibrated vis ((poly_)?${dical_name}_${apply_name})
// - sub_nsrcs: number of sources to subtract
// in output, meta will be updated:
// - name: add sub_ prefix
// - sub: set to "sub"
process hypSubUV {
    input:
    tuple val(obsid), val(meta), path(metafits), path(vis), path(srclist)
    output:
    tuple val(obsid), val(meta), path(sub_vis), path(logs)

    storeDir "${params.outdir}/${obsid}/cal${params.cal_suffix}"

    tag "${obsid}.${old_name}"
    label "hyperdrive"

    script:
    old_name = meta.name
    meta.sub = "sub"
    meta.name = "sub_${meta.name}"
    sub_vis = "hyp_${obsid}_${meta.name}.uvfits"
    logs = "hyp_vis-${meta.name}_uv.log"
    """
    ${params.hyperdrive} vis-sub \
        --data "${metafits}" "${vis}" \
        --beam "${params.beam_path}" \
        --source-list "${srclist}" \
        --invert --num-sources ${meta.sub_nsrcs} \
        --outputs "${sub_vis}" \
        | tee "${logs}"
    # TODO: ^ num sources is hardcoded twice, would be better to re-use model from cal
    """
}

// ensure subtracted ms are present, or subtract srclist from cal vis with hyperdrive
// meta is a map containing visibility metadata, fields:
// - name: unique name of the calibrated vis ((poly_)?${dical_name}_${apply_name})
// - sub_nsrcs: number of sources to subtract
// in output, meta will be updated:
// - name: add sub_ prefix
// - sub: set to "sub"
process hypSubMS {
    input:
    tuple val(obsid), val(meta), path(metafits), path(vis), path(srclist)
    output:
    tuple val(obsid), val(meta), path(sub_vis), path(logs)

    storeDir "${params.outdir}/${obsid}/cal${params.cal_suffix}"

    tag "${obsid}.${name}"
    label "hyperdrive"

    script:
    old_name = meta.name
    meta.sub = "sub"
    meta.name = "sub_${meta.name}"
    sub_vis = "hyp_${obsid}_${meta.name}.ms"
    logs = "hyp_vis-${meta.name}_ms.log"
    """
    ${params.hyperdrive} vis-sub \
        --data "${metafits}" "${vis}" \
        --beam "${params.beam_path}" \
        --source-list "${srclist}" \
        --invert --num-sources ${meta.sub_nsrcs} \
        --outputs "${sub_vis}" \
        | tee "${logs}"
    """
}

// ensure ionosubtracted uvfits are present, or ionosubtract srclist from cal vis with hyperdrive
// meta is a map containing visibility metadata, fields:
// - name: unique name of the calibrated vis ((poly_)?${dical_name}_${apply_name})
// - sub_nsrcs: total number of sources to subtract (including iono sources)
// - ionosub_nsrcs: number of sources to iono subtract
// in output, meta will be updated:
// - name: add ionosub_ prefix
// - sub: set to "ionosub"
process hypIonoSubUV {
    input:
    tuple val(obsid), val(meta), path(metafits), path(vis), path(srclist)
    output:
    tuple val(obsid), val(meta), path(sub_vis), path(json), path(logs)

    storeDir "${params.outdir}/${obsid}/cal${params.cal_suffix}"

    tag "${obsid}.${old_name}"
    label "hyperdrive"

    time 8.hour

    script:
    old_name = meta.name
    meta.sub = "ionosub"
    meta.name = "ionosub_${meta.name}"
    sub_vis = "hyp_${obsid}_${meta.name}.uvfits"
    logs = "hyp_vis-${meta.name}_uv.log"
    json = "hyp_peel_${obsid}_${meta.name}_uv.json"
    """
    ${params.hyperdrive} peel \
        --data "${metafits}" "${vis}" \
        --beam "${params.beam_path}" \
        --source-list "${srclist}" \
        --iono-sub ${meta.ionosub_nsrcs} \
        --sub ${meta.sub_nsrcs} \
        --outputs "${sub_vis}" "${json}" \
        | tee "${logs}"
    """
}

// ensure ionosubtracted ms are present, or ionosubtract srclist from cal vis with hyperdrive
// meta is a map containing visibility metadata, fields:
// - name: unique name of the calibrated vis ((poly_)?${dical_name}_${apply_name})
// - sub_nsrcs: total number of sources to subtract (including iono sources)
// - ionosub_nsrcs: number of sources to iono subtract
// in output, meta will be updated:
// - name: add ionosub_ prefix
// - sub: set to "ionosub"
process hypIonoSubMS {
    input:
    tuple val(obsid), val(meta), path(metafits), path(vis), path(srclist)
    output:
    tuple val(obsid), val(meta), path(sub_vis), path(json), path(logs)

    storeDir "${params.outdir}/${obsid}/cal${params.cal_suffix}"

    tag "${obsid}.${old_name}"
    label "hyperdrive"

    time 8.hour

    script:
    old_name = meta.name
    meta.sub = "ionosub"
    meta.name = "ionosub_${meta.name}"
    sub_vis = "hyp_${obsid}_${meta.name}.ms"
    logs = "hyp_vis-${meta.name}_ms.log"
    json = "hyp_peel_${obsid}_${meta.name}_ms.json"
    """
    echo nproc \$(nproc)
    ${params.hyperdrive} peel \
        --data "${metafits}" "${vis}" \
        --beam "${params.beam_path}" \
        --source-list "${srclist}" \
        --iono-sub ${meta.ionosub_nsrcs} \
        --sub ${meta.sub_nsrcs} \
        --outputs "${sub_vis}" "${json}" \
        -v \
        | tee "${logs}"
    """
}

process cthulhuPlot {
    input:
    tuple val(obsid), val(meta), path(srclist), path(offsets)
    output:
    tuple val(obsid), val(meta), path("cthuluplot_${title}*.png")

    storeDir "${params.outdir}/${obsid}/iono_qa${params.cal_suffix}"

    tag "${obsid}.${meta.name}"

    label 'python'

    script:
    title = "${obsid}_${meta.name}"
    cthulhuplot = "cthuluplot_${title}.png"
    extra = meta.time_res ? "--time_res=${meta.time_res}" : ""
    template "cthulhuplot.py"
}

// QA tasks that can be run on preprocessed visibility files.
process prepVisQA {
    input:
    tuple val(obsid), path(uvfits)
    output:
    tuple val(obsid), path(metrics)

    storeDir "${params.outdir}/${obsid}/vis_qa"

    tag "${obsid}"

    label 'python'

    script:
    metrics = "birli_${obsid}_prepvis_metrics.json"
    """
    #!/bin/bash -eux
    run_prepvisqa.py "${uvfits}" --out "${metrics}"
    """
}

// QA tasks that can be run on calibrated visibility files.
process visQA {
    input:
    tuple val(obsid), val(meta), path(uvfits)
    output:
    tuple val(obsid), val(meta), path(json)

    storeDir "${params.outdir}/${obsid}/vis_qa${params.cal_suffix}"

    tag "${obsid}.${meta.name}"

    label 'python'

    script:
    // TODO: json = "vis_metrics_${meta.cal_prog}_${obsid}_${meta.name}.json"
    json = "${meta.cal_prog}_${obsid}_${meta.name}_vis_metrics.json"
    """
    #!/bin/bash -eux
    run_visqa.py "${uvfits}" --out "${json}"
    """
}

process uvMeta {
    input:
    tuple val(obsid), val(meta), path(vis)
    output:
    tuple val(obsid), val(meta), path(uvmeta)

    storeDir "${params.outdir}/${obsid}/vis_qa${params.cal_suffix}"

    tag "${obsid}.${meta.name}"

    label 'python'
    time 15.minute

    script:
    uvmeta = "uvmeta_${meta.cal_prog}_${obsid}_${meta.name}.json"
    template "uvmeta.py"
}

// QA tasks that can be run on each calibration parameter set
process calQA {
    input:
    tuple val(obsid), val(meta), path(metafits), path(soln)
    output:
    tuple val(obsid), val(meta), path(metrics)

    storeDir "${params.outdir}/${obsid}/cal_qa${params.cal_suffix}"

    tag "${obsid}.${meta.name}"

    label 'python'

    script:
    metrics = "${meta.cal_prog}_soln_${obsid}_${meta.name}_X.json"
    """
    #!/bin/bash -eux
    run_calqa.py "${soln}" "${metafits}" --pol X --out "${metrics}"
    """
}

// write info from solutions to json
process solJson {
    input:
    tuple val(obsid), val(meta), path(soln)
    output:
    tuple val(obsid), val(meta), path(metrics)

    storeDir "${params.outdir}/${obsid}/cal_qa${params.cal_suffix}"
    tag "${obsid}.${meta.name}"

    label 'python'

    script:
    metrics = "${meta.cal_prog}_soln_${obsid}_${meta.name}.fits.json"
    template "soljson.py"
}

process plotPrepVisQA {
    input:
    tuple val(obsid), path(metrics)
    output:
    tuple val(obsid), path(img)

    storeDir "${params.outdir}/${obsid}/prep"

    tag "${obsid}"

    label 'python'

    script:
    img = "prepvis_metrics_${obsid}_rms.png"
    """
    #!/bin/bash -eux
    plot_prepvisqa.py "${metrics}" --out "${img}" --save
    """
}

process plotSols {
    input:
    tuple val(obsid), val(meta), path(metafits), path(soln)
    output:
    tuple val(obsid), val(meta), path(plots_glob)

    storeDir "${params.outdir}/${obsid}/cal_qa${params.cal_suffix}"

    tag "${obsid}.${meta.name}"

    label "hyperdrive"

    script:
    plots_glob = "${meta.cal_prog}_soln_${obsid}*_${meta.name}_{phases,amps}.png"
    """
    hyperdrive solutions-plot -m "${metafits}" ${soln}
    """
}

process plotCalQA {
    input:
    tuple val(obsid), val(meta), path(metrics)
    output:
    tuple val(obsid), val(meta), path(plots_glob)

    storeDir "${params.outdir}/${obsid}/cal_qa${params.cal_suffix}"

    tag "${obsid}.${meta.name}"

    label 'python'

    script:
    plots_glob = "calmetrics_${obsid}_${meta.name}_{fft,variance,dlyspectrum}.png"
    """
    #!/bin/bash -eux
    plot_calqa.py "${metrics}" --out "calmetrics_${obsid}_${meta.name}.png" --save
    """
}

process plotVisQA {
    input:
    tuple val(obsid), val(name), path(metrics)
    output:
    tuple val(obsid), val(name), path("${cal_prog}_${obsid}_${name}_vis_metrics_rms.png")

    storeDir "${params.outdir}/${obsid}/vis_qa${params.cal_suffix}"

    tag "${obsid}"

    label 'python'

    script:
    """
    #!/bin/bash -eux
    plot_visqa.py "${metrics}" --out "${cal_prog}_${obsid}_${name}_vis_metrics_rms.png" --save
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

process delaySpec {
    input:
    tuple val(obsid), val(name), path(uvfits)
    output:
    tuple val(obsid), val(name), path(dlyspec)

    storeDir "${params.outdir}/${obsid}/vis_qa${params.cal_suffix}"

    tag "${obsid}.${name}"

    label 'python'

    script:
    title = "${obsid}_${name}"
    dlyspec = "dlyspec_${title}.png"
    vmin = 3e13
    vmax = 1e15
    template "jline_delay_spec_from_uvfits.py"
}

// create dirty iamges of xx,yy,v
// TODO: merge img_params into meta?
process wscleanDirty {
    input:
    tuple val(obsid), val(meta), path(vis), val(img_params)
    output:
    tuple val(obsid), val(meta), path(img_glob)

    storeDir "${params.outdir}/${obsid}/img${params.img_suffix}${params.cal_suffix}"

    tag "${obsid}.${meta.name}"
    label "wsclean"

    time { 5.minute * (1 + (multiplier * pix_mult * chan_mult * inter_mult)) }

    beforeScript = 'pwd; hostname; df -h .'

    script:
    multiplier = vis.collect().size()
    mult_suffix = multiplier > 1 ? "_x${multiplier}" : ""
    // name_mult = "${name}${mult_suffix}"
    img_name = "${meta.cal_prog}_${obsid}${meta.subobs?:''}_${meta.name}${meta.inter_tok?:''}"

    // multipliers for determining compute resources
    pix_mult = 1 + (img_params.size / 1024) ** 2
    chan_mult = 1 + ("${img_params.channels_out}".split(' ')[0] as Double) / 25
    if (meta.interval) {
        inter_mult = meta.interval[1] - meta.interval[0]
    } else {
        inter_mult = 1 + ("${img_params.intervals_out}".split(' ')[0] as int) / 3
    }

    vis_ms = vis.collect {"${it.baseName}.ms"}
    vis = vis.collect()

    img_glob = "wsclean_${img_name}"
    if (multiinterval && !meta.inter_tok) {
        img_glob += "-t????"
    }
    if (multichannel) {
        img_glob += "-MFS"
    }
    img_glob += "-{XX,YY,XY,XYi,I,Q,U,V}-{dirty,uv-real,uv-imag}.fits"
    """
    #!/bin/bash -eux
    """ + (
            // convert any uvfits to ms
            [vis, vis_ms].transpose().collect { uv, ms ->
                (uv.extension == 'uvfits' ? \
                """${params.casa} -c "importuvfits('${uv}', '${ms}')" """ : "")
            }.join("\n")
        ) + """
        """ + (
        // run chgcentre if params.chgcentre_args specified.
        params.chgcentre_args ? \
            vis_ms.collect {"${params.chgcentre} ${params.chgcentre_args} ${it}"}.join("\n") : \
            ""
    ) + """
    ${params.wsclean} \
        ${img_params.args} \
        -weight ${img_params.weight} \
        -name wsclean_${img_name} \
        -size ${img_params.size} ${img_params.size} \
        -scale ${img_params.scale} \
        -pol ${img_params.pol} \
        -channels-out ${img_params.channels_out} \
        -save-uv \
        -niter 0 \
    """ + ((meta.interval) ? "-interval ${meta.interval[0]} ${meta.interval[1]}" : "") + """ \
    """ + vis_ms.join(' ')
}

// deconvolved images with wsclean
// note: I can't get `-reuse-dirty` to work with multi-interval data
process wscleanDConv {
    input:
    tuple val(obsid), val(meta), path(vis), val(img_params)  //, val(dirtyImgs)
    output:
    tuple val(obsid), val(meta), path(img_glob)
        // path("wsclean_${img_name}-sources.txt") <- only works for stokes I

    storeDir "${params.outdir}/${obsid}/img${img_params.suffix}${params.cal_suffix}"

    tag "${obsid}.${meta.name}"
    label "wsclean"

    time { 20.min * (1 + (multiplier * pix_mult * chan_mult * iter_mult * inter_mult)) }

    beforeScript = 'pwd; hostname; df -h .'

    script:
    multiplier = vis.collect().size()
    img_name = "${meta.cal_prog}_${obsid}${meta.subobs?:''}_${meta.name}${meta.inter_tok?:''}"

    // multipliers for determining compute resources
    pix_mult = 1 + (img_params.size / 1024) ** 2
    chan_mult = 1 + ("${img_params.channels_out}".split(' ')[0] as int) / 25
    if (meta.interval) {
        inter_mult = meta.interval[1] - meta.interval[0]
    } else {
        inter_mult = 1 + ("${img_params.intervals_out}".split(' ')[0] as int) / 3
    }
    iter_mult = 1 + Math.sqrt(img_params.niter as Double) / 100
    vis_ms = vis.collect {"${it.baseName}.ms"}
    vis = vis.collect()
    img_glob = "wsclean_${img_name}"
    if (multiinterval && !meta.inter_tok) {
        img_glob += "-t????"
    }
    if (multichannel) {
        img_glob += "-MFS"
    }
    img_glob += "-{XX,YY,XY,XYi,I,Q,U,V}-image.fits"
    """
    #!/bin/bash -eux
    """ + (
            // convert any uvfits to ms
            [vis, vis_ms].transpose().collect { uv, ms ->
                (uv.extension == 'uvfits' ? \
                """${params.casa} -c "importuvfits('${uv}', '${ms}')" """ : "")
            }.join("\n")
    ) + """
    """ + (
            // run chgcentre if params.chgcentre_args specified.
            params.chgcentre_args ? \
                vis_ms.collect {"${params.chgcentre} ${params.chgcentre_args} ${it}"}.join("\n") : \
                ""
    ) + """
    ${params.wsclean} \
        ${img_params.args} \
        -name wsclean_${img_name} \
        """ + // (dirtyImgs?"-reuse-dirty wsclean_${img_name}":"") +
        """ \
        -weight ${img_params.weight} \
        -size ${img_params.size} ${img_params.size} \
        -scale ${img_params.scale} \
        -pol ${img_params.pol} \
        -channels-out ${img_params.channels_out} \
        -intervals-out ${img_params.intervals_out} \
        -niter ${img_params.niter} \
        -mgain ${img_params.major_clean_gain} -gain ${img_params.minor_clean_gain} \
        -auto-threshold ${img_params.auto_threshold} -auto-mask ${img_params.auto_mask} \
        -mwa-path ${img_params.mwa_path} \
        -circular-beam \
    """ + ((meta.interval) ? "-interval ${meta.interval[0]} ${meta.interval[1]}" : "") + """ \
    """ + vis_ms.join(' ')
}

process imgQuantiles {
    input:
    tuple val(obsid), val(meta), path(fits)
    output:
    tuple val(obsid), val(meta), path(csv)

    storeDir "${params.outdir}/${obsid}/img_qa${params.img_suffix}${params.cal_suffix}"

    tag "${obsid}.${meta.name}.${meta.suffix}"

    label 'python'

    time {5.min * task.attempt}

    script:
    csv = "quantile_${meta.name}_${meta.suffix}.csv"
    """
    #!/usr/bin/env python

    from astropy.io import fits
    import numpy as np

    with fits.open("${fits}") as hdus:
        data = hdus[0].data[0,0,:,:]

    quantiles = [0.90, 0.95, 0.99, 1.0]
    values = np.quantile(data, quantiles)
    np.savetxt(
        "${csv}", np.array([quantiles, values]).T,
        delimiter=',', fmt='%.6e', header='quantile,value'
    )

    # TODO: histogram?
    # import matplotlib.pyplot as plt
    # plt.style.use([astropy_mpl_style, 'dark_background'])
    """
}

// power spectrum metrics via chips
process psMetrics {
    input:
    tuple val(obsid), val(meta), path("vis.uvfits")
    output:
    tuple val(obsid), val(meta), path(out_metrics)

    storeDir "${params.outdir}/${obsid}/ps_metrics${params.cal_suffix}"

    tag "${obsid}.${meta.name}"

    label "chips"

    script:
    band = 0
    uv_base = "vis"
    out_metrics = "output_metrics_${meta.cal_prog}_${obsid}_${meta.name}.dat"
    """
    #!/bin/bash -eux
    export DATADIR="\$PWD"
    export OUTPUTDIR="\$PWD/"
    ${params.ps_metrics} "${uv_base}" "${band}" "${meta.nchan?:384}" "${uv_base}" 2>&1 \
        | tee "${uv_base}.log"
    mv "output_metrics_vis.dat" "${out_metrics}"
    """
}

// analyse images of V,XX,YY
process imgQA {
    input:
    tuple val(obsid), val(name), path(fits)
    output:
    tuple val(obsid), val(name), path(json)

    storeDir "${params.outdir}/${obsid}/img_qa${params.img_suffix}${params.cal_suffix}"

    tag "${obsid}.${name}"

    label 'python'

    script:
    json = "wsclean_hyp_${obsid}_${name}-MFS.json"
    """
    #!/bin/bash -eux
    run_imgqa.py ${fits.join(' ')} --out ${json}
    """
}

process uvPlot {
    input:
    tuple val(obsid), val(name), path(uvfits)
    output:
    tuple val(obsid), val(name), path("uvplot_${title}_{XX,YY,XY,YX}.png")

    storeDir "${params.outdir}/${obsid}/vis_qa${params.cal_suffix}"

    tag "${obsid}.${name}"

    label 'python'

    script:
    title = "${obsid}_${name}"
    uvplot = "uvplot_${title}"
    template "uvplot_3d.py"
}

// make a thumbnail png from a fits image
process thumbnail {
    input:
    tuple val(obsid), val(meta), path(img)
    output:
    tuple val(obsid), val(meta), path(thumb)

    storeDir "${params.outdir}/${obsid}/img_qa${params.img_suffix}${params.cal_suffix}"

    tag "${obsid}.${meta.name}.${meta.suffix}"

    label 'python'
    time {10.min * task.attempt}

    script:
    thumb = "${obsid}_${meta.name}_${meta.suffix}.png"
    title = "${obsid} ${meta.name} ${meta.suffix}"
    args = "${params.thumbnail_args}"
    template "thumbnail.py"
}

// polarimetry composite raster
// meta contains info about how to the composite is to be rastered.
// - name: name of the composite image (for grouping of frames)
// - subobs: (optional) suffix after obsid
// - prod: type of image (dirty, image, etc.), should all be the same
// - order: polarization order to use (e.g. ['XX', 'V', 'YY'])
// - limits: brightness limits for each polarization image.
process polComp {
    input:
    tuple val(obsid), val(meta), path(fits)
    output:
    tuple val(obsid), val(meta), path(polcomp)

    storeDir "${params.outdir}/${obsid}/img_qa${params.img_suffix}${params.cal_suffix}"

    tag "${obsid}${meta.subobs?:''}.${meta.name}.${meta.prod}.${meta.orderName}"

    label 'python'

    script:
    meta.orderName = meta.order.join('')
    polcomp = "${obsid}${meta.subobs?:''}_${meta.name}_polcomp_${meta.prod}_${meta.orderName}.png"
    title = "${obsid}${meta.subobs?:''} ${meta.name} ${meta.prod} ${meta.orderName}"
    args = ""
    if (meta.order) {
        args += " --pol_order ${meta.order.join(' ')}"
    }
    if (meta.limits) {
        args += " --limits ${meta.limits.join(' ')}"
    }
    template "polcomp.py"
}

// montage multiple thumbnails into a single image
// meta contains info about how to montage the images
// - name: name of montage (for grouping of frames)
// - subobs: (optional) suffix after obsid
process polMontage {
    input:
    tuple val(obsid), val(meta), path(thumbs)
    output:
    tuple val(obsid), val(meta), path(montage)

    storeDir "${params.outdir}/${obsid}/img_qa${params.img_suffix}${params.cal_suffix}"

    tag "${obsid}${meta.subobs?:''}.${meta.name}"

    label "imagemagick"

    script:
    montage = "${obsid}${meta.subobs?:''}_${meta.name}_montage.png"
    """
    montage \
        -font /usr/share/fonts/truetype/droid/DroidSansFallbackFull.ttf \
        ${thumbs.join(' ')} \
        -geometry +0+0 \
        -background none \
        ${montage}
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

    storeDir "${results_dir}${params.img_suffix}${params.cal_suffix}"
    stageInMode "symlink"

    tag "${name}"

    label 'ffmpeg'

    script:
    """
    #!/bin/bash -eux
    ffmpeg -y -framerate 5 \
        -pattern_type glob -i "??????????.png" \
        -vcodec libx264 \
        -vf "scale='min(3840,iw)':'min(2160,ih)':force_original_aspect_ratio=decrease,pad=ceil(iw/2)*2:ceil(ih/2)*2" \
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

def parseCsv(path, header = true, skipBytes = 0) {
    allLines = path.readLines()
    if (!header) {
        return allLines.collect { it.split(',') }
    } else {
        head = allLines[0][skipBytes..-1].split(',')
        body = allLines[1..-1]
        return body.collect { row ->
            [head, row.split(',')].transpose().collectEntries()
        }
    }
}

import java.text.SimpleDateFormat
import java.util.LinkedHashMap
import groovy.transform.Synchronized

@Synchronized
def deepcopy(orig) {
     bos = new ByteArrayOutputStream()
     oos = new ObjectOutputStream(bos)
     oos.writeObject(orig); oos.flush()
     bin = new ByteArrayInputStream(bos.toByteArray())
     ois = new ObjectInputStream(bin)
     return ois.readObject()
}

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

// decompose an image filename into:
// - interval if multiple intervals or -1 if combined
// - channel of multiple channels or -1 if combined
// - polarization
// - image product name (dirty, image, model, uv-{real,imag})
def decomposeImg = { img ->
    def tokens = img.baseName.split('-').collect()
    // defaults:
    def meta = [chan: -1, inter: -1, chan_tok: "MFS"]
    // this handles the case where product is "uv-{real,imag}"
    if (tokens[-2] == "uv") {
        tokens = tokens[0..-3] + ["uv-${tokens[-1]}"]
    }
    meta.prod = tokens.removeLast()
    meta.pol = tokens.removeLast()
    // channel is only present in multi-frequency imaging
    if (multichannel && (chan_tok = tokens.removeLast()) != "MFS") {
        meta.chan_tok = chan_tok
        meta.chan = (chan_tok as int)
    }
    // suffix without interval
    meta.inter_suffix = [meta.chan_tok, meta.pol, meta.prod].join('-')
    meta.suffix = meta.inter_suffix
    if (multiinterval && (inter_tok = tokens.removeLast()) =~ /t\d{4}/) {
        meta.inter_tok = inter_tok
        meta.inter = inter_tok[1..-1] as int
        meta.suffix = "${inter_tok}-${meta.inter_suffix}"
    }
    meta
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

                def obs_name = stats.obsname;
                def groupid = stats.groupid;

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
                if (params.filter_pointings && !params.filter_pointings.contains(pointing)) {
                    fail_reasons += ["pointing=${pointing}"]
                }
                if (params.filter_quality && dataquality > params.filter_quality) {
                    fail_reasons += ["dataquality=${dataquality} (${dataqualitycomment})"]
                }
                if (params.filter_bad_tile_frac && bad_tile_frac > params.filter_bad_tile_frac) {
                    fail_reasons += ["bad_tiles(${bad_tiles.size()})=${displayInts(bad_tiles)}"]
                }
                if (params.filter_dead_dipole_frac && dead_dipole_frac > params.filter_dead_dipole_frac) {
                    fail_reasons += ["dead_dipole_frac=${dead_dipole_frac}"]
                }
                // if (ra_phase_center != 0.0) {
                //     fail_reasons += ["ra_phase_center=${ra_phase_center}"]
                // }
                // if (dec_phase_center != -27.0) {
                //     fail_reasons += ["dec_phase_center=${dec_phase_center}"]
                // }
                def summary = [
                    fail_reasons: fail_reasons,
                    // obs metadata
                    obs_name: obs_name,
                    groupid: groupid,
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

                    summary.groupid,
                    summary.ra_pointing,
                    summary.dec_pointing,
                    summary.ra_phase_center,
                    summary.dec_phase_center,
                    summary.pointing,
                    summary.lst,
                    summary.obs_name,

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
                    "OBS", "FAIL REASON", "GROUP ID", "RA POINT", "DEC POINT", "RA PHASE", "DEC PHASE", "POINT", "LST DEG", "OBS NAME",
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

        pass | wsMetafits & wsSkyMap & wsPPDs

    emit:
        // channel of good obsids with their metafits: tuple(obsid, metafits)
        obsMetafits = wsMetafits.out

        // channel of video name and frames to convert
        frame = wsSkyMap.out
            .map { _, png -> ["skymap", png] }
            .mix( wsPPDs.out.map { _, png -> ["ppd", png] } )
            .groupTuple()

        // channel of (obsid, groupid, pointing)
        obsGroupPoint = wsSummary.map { obsid, summary ->
            [obsid, summary.groupid, summary.pointing]
        }
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

        if (params.noprepqa) {
            channel.from([]) | prepVisQA
        } else {
            obsMetafits.join(asvoPrep.out) | prepVisQA
        }

        if (params.nossins) {
            channel.from([]) | ssins
        } else {
            asvoPrep.out | ssins
        }

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
        if (params.noprepqa) {
            channel.from([]) | flagQA & autoplot
        } else {
            obsMetafits.join(asvoPrep.out) | flagQA & autoplot
        }

        // collect flagQA results
        // TODO: collect sky_chans from flagQA
        def sky_chans = (params.sky_chans).collect { ch -> "$ch".toString() }
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
        obsMetaVis = (params.noprepqa ? \
            // no filter unless prepqa
            obsMetafits.join(asvoPrep.out) : \
            // filter
            flagQA.out
                .filter { obsid, json ->
                    params.flag_occupancy_threshold && jslurp.parse(json).total_occupancy < params.flag_occupancy_threshold
                }
                .map { obsid, _ -> obsid }
                .join(obsMetafits)
                .join(asvoPrep.out))
        // channel of extra flags for each obsid
        obsFlags = (params.noprepqa ? \
            // flags empty unless prepqa
            asvoPrep.out.map { obsid, _ -> [obsid, [], [], [], []] } : \
            // calculate extra flags if prepqa
            prepVisQA.out
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
                    def newAntennas = ([]) as Set
                    if (!params.noPrepFlags) {
                        newAntennas.addAll(prepAntennas)
                    }
                    if (!params.noManualFlags) {
                        newAntennas.addAll(manualAntennas)
                    }
                    [
                        obsid,
                        flagAntennas as int[],
                        prepAntennas as int[],
                        manualAntennas as int[],
                        (newAntennas - flagAntennas) as int[]
                    ]
                })
        // channel of video name and frames to convert
        frame = plotPrepVisQA.out.map { _, png -> ["prepvisqa_rms", png] }
            .mix(ssins.out.flatMap { _, __, imgs, ___ ->
                imgs.collect { img ->
                    tokens = img.baseName.split('_')
                    prefix = tokens[0]
                    suffix = tokens[-1]
                    if (suffix == "SSINS") {
                        ["ssins_${prefix}", img]
                    } else {
                        ["ssins_${prefix}_${suffix}", img]
                    }
                }
            })
            .mix(autoplot.out.map {_, img -> ["prepvisqa_autoplot", img]})
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
                    def startMatch = diCalLog.getText() =~ /([\d: -]+) INFO  hyperdrive di-calibrate/
                    def startTime = null
                    if (startMatch) {
                        startTime = logDateFmt.parse(startMatch[0][1])
                    }
                    def convergedMatch = (diCalLog.getText() =~ /([\d: -]+) INFO  All timesteps: (\d+)\/(\d+)/)
                    if (startTime && convergedMatch) {
                        def convergedTime = logDateFmt.parse(convergedMatch[0][1])
                        convergedNumerator = convergedMatch[0][2] as int
                        convergedDenominator = convergedMatch[0][3] as int
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
                storeDir: "${results_dir}${params.cal_suffix}"
            )

        // channel of individual dical solutions: tuple(obsid, meta, soln)
        // - hypCalSol gives multiple solutions, transpose gives 1 tuple per solution.
        eachCal = hypCalSol.out
            .map { obsid, solutions, _ -> [obsid, solutions] } // ignore logs
            .transpose()
            .map { obsid, soln ->
                // give each calibration a name from basename of solution fits.
                // this is everything after the obsid
                def meta = [
                    dical_name: soln.baseName.split('_')[3..-1].join('_'),
                    cal_prog: "hyp"
                ]
                meta.name = meta.dical_name
                [obsid, meta, soln]
            }

        // generate json from solns
        eachCal | solJson

        // do polyFit unless --nopoly is set
        if (params.nopoly) {
            channel.from([]) | polyFit
        } else {
            eachCal | polyFit
        }
        // channel of all dical (and polyfit) solutions: tuple(obsid, meta, soln)
        allCal = eachCal.mix(
            polyFit.out.map { obsid, meta, soln, _ -> [obsid, meta, soln] }
        )

        // collect solJson results as .tsv
        solJson.out
            // form row of tsv from json fields we care about
            .map { obsid, meta, json ->
                def name = meta.name
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
                storeDir: "${results_dir}${params.cal_suffix}"
            )
            // display output path and number of lines
            | view { [it, it.readLines().size()] }

        // channel of metafits for each obsid: tuple(obsid, metafits)
        obsMetafits = obsMetaVis.map { obsid, metafits, _, __ -> [obsid, metafits] }
        // calibration QA and plot solutions
        obsMetafits
            // join with solutions from eachCal and eachPolyCal
            .cross(allCal)
            .map { def (obsid, metafits) = it[0]; def (_, meta, soln) = it[1]
                [obsid, meta, metafits, soln]
            }
            | (plotSols & calQA)

        // plot each calQA result
        if (params.noplotcalqa) {
            channel.from([]) | plotCalQA
        } else {
            calQA.out | plotCalQA
        }

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
                storeDir: "${results_dir}${params.cal_suffix}"
            )
            // display output path and number of lines
            | view { [it, it.readLines().size()] }

    // channel of obsids and names that pass qa. tuple(obsid, name)
    // - take tuple(obsid, cal_name, json) from calQA.out
    // - filter on json.STATUS == "PASS"
    // - take obsid and name
    obsidNamePass = calQA.out
        .filter { _, __, json -> parseJson(json).STATUS == "PASS" }
        .map { obsid, meta, _ -> [obsid, meta] }

    emit:
        // channel of calibration solutions that pass qa. tuple(obsid, name, cal)
        // - take tuple(obsid, meta, soln) from allCal
        // - match with obsidNamePass on (obsid, cal_name)
        passCal = allCal
            .cross(obsidNamePass) {it -> [it[0], it[1].name]}
            .map{ allCal, _ -> allCal }
        // channel of files to archive, and their buckets
        archive = calQA.out.map { _, __, json -> ["calqa", json]}
        // channel of video name and frames to convert
        frame = plotCalQA.out.mix(plotSols.out)
            .flatMap { _, meta, pngs ->
                pngs.collect { png ->
                    def png_name = png.baseName.split('_')[-1]
                    ["calqa_${meta.name}_${png_name}", png]
                }
            }
            .groupTuple()
}

// process uvfits visibilities
workflow uvfits {
    take:
        obsMetaVis
    main:
        // vis QA
        obsMetaVis
            .filter { _, meta, __ -> !meta.sub }
            | visQA
        channel.from([]) | uvPlot

        // collect visQA results as .tsv
        visQA.out
            // form row of tsv from json fields we care about
            .map { obsid, meta, json ->
                def stats = parseJson(json);
                [
                    obsid,
                    meta.name,
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
                storeDir: "${results_dir}${params.cal_suffix}"
            )
            // display output path and number of lines
            | view { [it, it.readLines().size()] }

        // TODO: plotVisQA isn't working
        channel.from([]) | plotVisQA
        // visQA.out | plotVisQA

        // ps_metrics
        if (params.nopsmetrics) {
            passVis = obsMetaVis
                .map { obsid, meta, __ -> [obsid, meta.name] }
                .unique()
        } else {
            obsMetaVis | psMetrics

            // collect psMetrics as a .dat
            psMetrics.out
                // read the content of each ps_metrics file including the trailing newline
                .map { obsid, vis_name, dat -> dat.getText() }
                .collectFile(
                    name: "ps_metrics.dat",
                    storeDir: "${results_dir}${params.cal_suffix}"
                )
                // display output path and number of lines
                | view { [it, it.readLines().size()] }

            // collect psMetrics as a .tsv
            psMetrics.out
                // form each row of tsv
                .map { obsid, meta, dat ->
                    def dat_values = dat.getText().split('\n')[0].split(' ')[1..-1]
                    def vis_name = meta.name
                    ([obsid, vis_name] + dat_values).join("\t")
                }
                .collectFile(
                    name: "ps_metrics.tsv", newLine: true, sort: true,
                    seed: [
                        "OBS", "CAL NAME", "P_WEDGE", "NUM_CELLS", "P_WINDOW", "NUM_CELLS",
                        "P_ALL", "D3"
                    ].join("\t"),
                    storeDir: "${results_dir}${params.cal_suffix}"
                )
                // display output path and number of lines
                | view { [it, it.readLines().size()] }

            // filter obs with big ps_metrics
            // TODO: this is a dirty hack until we get proper vis metrics filtering
            passVis = psMetrics.out
                // form each row of tsv
                .map { obsid, vis_name, dat ->
                    def (p_wedge, num_cells, p_window) = dat.getText().split('\n')[0].split(' ')[1..-1]
                    [obsid, vis_name, p_wedge, p_window]
                }
                // .filter { obsid, vis_name, p_wedge, p_window ->
                //     if (vis_name.contains("sub")) {
                //         return p_wedge.toFloat() < 100 && p_window.toFloat() < 50
                //     } else {
                //         return p_wedge.toFloat() < 500 && p_window.toFloat() < 500
                //     }
                // }
                .map { obsid, vis_name, p_wedge, p_window -> [obsid, vis_name] }
                .unique()
        }

        // delayspectrum
        // obsMetaVis | delaySpec
        channel.from([]) | delaySpec

    emit:
        // channel of files to archive, and their buckets
        archive = visQA.out.map { _, __, json -> ["visqa", json]}
        // channel of video name and frames to convert
        frame = plotVisQA.out
            .map { _, name, png -> ["visqa_${name}_rms", png] }
            .mix(uvPlot.out.flatMap {_, name, pngs -> pngs.collect { png ->
                def pol = png.baseName.split('_')[-2..-1].join('_')
                ["visqa_${name}_${pol}", png]
            }})
            .mix(delaySpec.out.map { _, name, png -> ["visqa_dlyspec_${name}", png] })
            .groupTuple()

        passVis = passVis
}

// image measurement sets and QA images
workflow img {
    take:
        obsMetaVis
    main:

        // wsclean: make deconvolved images
        wscleanParams = [
            suffix: params.img_suffix,
            weight: params.img_weight,
            size: params.img_size,
            scale: params.img_scale,
            channels_out: params.img_channels_out,
            intervals_out: params.img_intervals_out,
            split_intervals: params.img_split_intervals,
            pol: params.img_pol,
            args: params.wsclean_args,
        ]
        wscleanDConvParams = deepcopy(wscleanParams) + [
            args: params.wsclean_dconv_args,
            niter: params.img_niter,
            minor_clean_gain: params.img_minor_clean_gain,
            major_clean_gain: params.img_major_clean_gain,
            auto_threshold: params.img_auto_threshold,
            auto_mask: params.img_auto_mask,
            mwa_path: params.img_mwa_path,
        ]
        if (params.img_split_intervals) {
            // add -t???? suffix to name, to split wsclean over multiple jobs
            splitObsMetaVis = obsMetaVis.flatMap { obsid, meta, vis ->
                (0..meta.ntimes).collect { i ->
                    def newMeta = deepcopy(meta)
                    newMeta.interval = [i, i+1]
                    // interval token
                    newMeta.inter_tok = sprintf("-t%04d", i)
                    [obsid, newMeta, vis, deepcopy(wscleanParams)]
                }
            }
        } else {
            splitObsMetaVis = obsMetaVis.map { obsid, meta, vis ->
                [obsid, meta, vis, deepcopy(wscleanParams)]
            }
        }
        if (params.nodeconv) {
            splitObsMetaVis | wscleanDirty
            channel.from([]) | wscleanDConv
        } else {
            // splitObsMetaVis.cross(wscleanDirty.out) {
            //         def (obsid, meta) = it;
            //         [obsid, meta.subobs, meta.name, meta.inter_tok]
            //     }
            //     .map { splitObsMetaVis_, wscleanDirty_ ->
            //         def (obsid, _, vis, imgParams) = splitObsMetaVis_
            //         def (__, meta, imgs) = wscleanDirty_
            //         // def dirtyImgs = imgs.findAll {img -> img.baseName.split('-')[-1] == 'dirty'}
            //         def newImgParams = deepcopy(imgParams + wscleanDConvParams)
            //         [obsid, meta, vis, newImgParams, dirtyImgs]
            //     }
            channel.from([]) | wscleanDirty
            splitObsMetaVis.map { obsid, meta, vis, imgParams ->
                    [obsid, meta, vis, deepcopy(wscleanDConvParams) + imgParams]
                }
            | wscleanDConv
        }

        obsMetaImg = wscleanDConv.out.mix(wscleanDirty.out)
            .flatMap { obsid, meta, imgs ->
                imgs.collect { img ->
                    imgMeta = deepcopy(meta)
                    imgMeta += decomposeImg(img)
                    [obsid, imgMeta, img] }}
        // obsMetaImg but only MFS images unless thumbnail_all_chans
        obsMetaImgMfs = obsMetaImg.filter { _, meta, img ->
                meta.prod !=~ /uv-.*/ \
                && (params.thumbnail_all_chans || meta.chan?:-1 == -1)
            }
        // obsMetaImg but only gridded uv images
        if (params.thumbnail_uvs) {
            obsMetaGrid = obsMetaImg.filter { _, meta, img -> meta.prod ==~ /uv-.*/ }
        } else {
            obsMetaGrid = channel.from([])
        }

        // calculate quantiles (what values are at nth percentile)
        obsMetaImgMfs \
            // .mix(obsMetaGrid) \
            | imgQuantiles

        // limits are used to set the color scale of each type of image.
        imgLimits = imgQuantiles.out
            .map { obsid, meta, hist ->
                limit = parseCsv(hist, true, 2)
                    // get 90th percentile
                    .find { row -> Float.compare(row.quantile as float, 0.90) == 0 }
                [obsid, meta.name, meta.inter_suffix, limit.value as float]
            }
            .groupTuple(by: 0..2)
            .map { obsid, name, suff, limits ->
                [obsid, name, suff, limits.max()]
            }
            .groupTuple(by: 0..1)
            .map { obsid, name, suffs, limits ->
                suffLimits = [suffs, limits].transpose().collectEntries()
                [obsid, name, suffLimits]
            }

        allSuffs = imgLimits.flatMap { obsid, name, limits -> limits.keySet()}
            .unique()
            .collect(sort: true)
            .toList()

        // write img limits to file
        allSuffs.subscribe { suffs ->
            imgLimits.map { obsid, name, limits ->
                    ([obsid, name] + suffs.collect { limits[it]?:'' }).join("\t")}
            .collectFile(
                name: "img_limits.tsv", newLine: true, sort: true,
                seed: (["OBSID", "IMG NAME"] + suffs).join("\t"),
                storeDir: "${results_dir}${params.img_suffix}${params.cal_suffix}"
            )
        }

        // value channel containing a map from img suffix to max limit
        suffLimits = imgLimits.map{ _, __, limits -> limits }
            .reduce([:]) { acc, limits ->
                limits.collectEntries { k, v -> [k, [v, acc[k]?:v].max()] }
            }
            .collect()
            .view { "suffLimits: ${it}" }

        // all valid pol orders for polComp
        polOrders = [
            // ["I", "V", "Q"],
            ["XX", "V", "YY"],
            // ["XX", "XY", "YY"]
        ]

        // do a polcomp for each prod in ['dirty', 'image']
        obsMetaImgMfs.filter {_, meta, __ -> meta.prod in ['dirty', 'image']}
            .map { obsid, meta, img -> [obsid, meta.name, meta.prod, meta, img] }
            .groupTuple(by: 0..2)
            .combine(suffLimits)
            .flatMap { obsid, name, prod, metas, imgs, limits ->
                polLimits = metas.collect { meta -> [meta.pol, limits[meta.inter_suffix]] }
                    .findAll { pol, limit -> limit != null }
                    .collectEntries()
                // this line is so damn salient
                polOrders.findAll { order -> polLimits.keySet().containsAll(order) }
                    .collect { order ->
                        limits = order.collect { pol -> polLimits[pol] }
                        meta_imgs = [metas, imgs].transpose()
                        imgs_ordered = order.collect {
                            pol -> meta_imgs.find { meta, img -> meta.pol == pol }[1]
                        }
                        [obsid, [name:name, prod:prod, order:order, limits:limits], imgs_ordered]
                    }
            }
            | polComp

        // make thumbnails
        obsMetaImgMfs \
            .mix(obsMetaGrid) \
            | thumbnail

        // imgQA and polComp for all groups of images
        // obsNameMfsImgs | (imgQA & polComp)
        channel.from([]) \
            | imgQA

        // montage of polarizations
        thumbnail.out.flatMap { obs, meta, png ->
                // visibility name, without sub*
                meta.vis_name = meta.name
                meta.sub = meta.sub?:'nosub'
                if (m = meta.vis_name =~ /(sub|ionosub)_(.*)/) {
                    meta.sub = m.group(1).toString()
                    meta.vis_name = m.group(2).toString()
                }

                // suffix without interval (-t????)
                // meta.subobs = meta.subobs?:''
                // meta.inter_suffix = meta.suffix
                // if (m = meta.inter_suffix =~ /(t????)-(.*)/) {
                //     meta.subobs = m.group(1).toString() + meta.subobs
                //     meta.inter_suffix = m.group(2).toString()
                // }
                // def montage_name = [sub: suffix, pol: sub_name].get(params.montage_by)
                subobs = "${meta.subobs?:''}${meta.inter_tok?:''}".toString()
                montages = []
                if (params.montage_by_sub) {
                    montages << [obsid, subobs, "${meta.vis_name}-${meta.chan_tok}-${meta.pol}-${meta.prod}".toString(), png]
                }
                if (params.montage_by_pol) {
                    montages << [obsid, subobs, "${meta.name}-${meta.chan_tok}-${meta.prod}".toString(), png]
                }
                montages
            }
            .groupTuple(by: 0..2)
            .map { obs, subobs, name, pngs ->
                [obs, [subobs:subobs, name:name], pngs.sort()]
            }
            .view {"polmontage ${it}"}
            | polMontage
        // thumbnail.out.map { obs, meta, thumb -> [obs, meta.name, meta, thumb]}
        //     .groupTuple(by: [0, 1])
        //     .map { obs, name, _, thumbs ->
        //         [obs, name, thumbs.sort()]
        //     }

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
                storeDir: "${results_dir}${params.img_suffix}${params.cal_suffix}"
            )
            // display output path and number of lines
            | view { [it, it.readLines().size()] }

    emit:
        // channel of files to archive, and their buckets
        archive = obsMetaImgMfs.filter { obsid, _, __ -> (!obsid.startsWith("e")) }
            .transpose()
            .map { _, __, img -> ["img", img]}
            .mix( imgQA.out.map { _, __, json -> ["imgqa", json]} )
        // channel of video name and frames to convert
        frame = polComp.out.map { _, meta, png ->
                ["imgqa_${meta.name}_polcomp_${meta.prod}_${meta.orderName}", png]
            }
            .mix(polMontage.out.map { _, meta, png ->
                ["imgqa_${meta.name}_polmontage", png]
            })
            .mix(thumbnail.out.flatMap { obsid, meta, png ->
                def frames_ = []
                if (meta.chan?:-1 == -1) {
                    frames_.push(["imgqa_${meta.name}_${meta.inter_suffix}", png])
                } else {
                    // scan through channels for each obsid, pol
                    if (params.frame_chan_scan) {
                        frames_.push(["imgqa_${obsid}${meta.subobs}_${meta.name}_${meta.pol}-${meta.prod}", png])
                    }
                    // scan through obsids for each channel, pol, prod
                    if (params.frame_obs_scan) {
                        frames_.push(["imgqa_${meta.name}_${meta.inter_suffix}", png])
                    }
                }
                frames_
            })
            .groupTuple()
}

// entrypoint: move unfiltered preprocessed uvfits from asvo accacia to mwaeor accacia
workflow archivePrep {
    obsids = channel.of(obsids_file)
        .splitCsv()
        .flatten()
        .filter { line -> !line.startsWith('#') }
        .unique()

    obsids | asvoPrep

    if (params.archive) {
        asvoPrep.out.map { _, vis -> ["prep", vis] } | archive
    }
}

// default entrypoint: get preprocessed vis from asvo and run qa
workflow {
    // get obsids from csv
    obsids = channel.of(obsids_file)
        .splitCsv()
        .flatten()
        .filter { line -> !line.startsWith('#') }
        .unique()

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
    if (params.noprep) {
        channel.from([]) | prep
    } else {
        ws.out.obsMetafits | prep
    }

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

    obsFlags = prep.out.obsFlags
        .map {obsid, flagAnts, prepAnts, manualAnts, newAnts -> [obsid, newAnts]}
    frame = ws.out.frame
            .mix(prep.out.frame)
    qaPrep(prep.out.obsMetaVis, obsFlags, frame)
}

// given a channel of tuple (obs, metafits, vis), calibrate and analyse
workflow qaPrep {
    take:
        obsMetafitsVis
        obsFlags
        frame
    main:
        // get sourcelists for each obs (only currently used in subtraction, not calibration)
        obsMetafitsVis.map { obsid, metafits, _ -> [obsid, metafits ] }
            | (hypSrclistAO & hypSrclistYaml)

        // obsMetaSrclist: channel of tuple(obsid, metafits, srclist)
        // cluster unless --nocluster
        if (params.nocluster) {
            channel.from([]) | rexCluster
            obsMetaSrclist = obsMetafitsVis.join(hypSrclistAO.out)
                .map { obsid, metafits, _, srclist ->
                    [obsid, metafits, srclist] }
        } else {
            hypSrclistAO.out | rexCluster
            obsMetaSrclist = obsMetafitsVis.join(rexCluster.out)
                .map { obsid, metafits, _, cluster ->
                    [obsid, metafits, cluster] }
        }

        // calibrate each obs that passes flag gate unless --nocal:
        if (params.nocal) {
            // empty channel disables a process
            channel.from([]) | cal
        } else {
            obsMetafitsVis.join(obsFlags)
                .filter { _, metafits, vis, __ -> metafits != null && vis != null}
                .map { obs, metafits, vis, flags -> [obs, metafits, vis, flags ? flags : [] ]}
                | cal
        }

        // channel of arguments for hypApply{UV,MS}
        // - take tuple(obsid, cal_name, soln) from cal.out.passCal
        // - match with tuple(obsid, metafits, prepUVFits) by obsid
        if (params.noapply) {
            allApply = channel.from([])
        } else {
            allApply = obsMetafitsVis
                .join(cal.out.passCal)
                .map { obsid, metafits, prepUVFits, meta, soln ->
                    def newMeta = deepcopy(meta)
                    newMeta += [
                        time_res: params.apply_time_res,
                        freq_res: params.apply_freq_res,
                        nodut1: params.nodut1,
                        apply_args: params.apply_args,
                    ]
                    newMeta.apply_name = "${newMeta.time_res}s_${newMeta.freq_res}kHz"
                    [obsid, newMeta, metafits, prepUVFits, soln]
                }
        }

        // apply calibration solutions to uvfits and ms unless --nouv or --noms
        if (params.nouv) {
            channel.from([]) | hypApplyUV
        } else {
            allApply | hypApplyUV
        }
        if (params.noms) {
            channel.from([]) | hypApplyMS
        } else {
            allApply | hypApplyMS
        }
        // get uvfits subtraction arguments from apply output.
        subArgsUV = obsMetaSrclist.cross(hypApplyUV.out)
            .map { obsMetaSrclist_, hypApplyUV_ ->
                def (obsid, metafits, srclist) = obsMetaSrclist_
                def (__, meta, vis) = hypApplyUV_;
                def newMeta = deepcopy(meta)
                newMeta.sub_nsrcs = params.sub_nsrcs
                [obsid, newMeta, metafits, vis, srclist]
            }
        subArgsMS = obsMetaSrclist.cross(hypApplyMS.out)
            .map { obsMetaSrclist_, hypApplyUV_ ->
                def (obsid, metafits, srclist) = obsMetaSrclist_
                def (__, meta, vis) = hypApplyUV_;
                def newMeta = deepcopy(meta)
                newMeta.sub_nsrcs = params.sub_nsrcs
                [obsid, newMeta, metafits, vis, srclist]
            }
        if (params.nosub) {
            channel.from([]) | (hypSubUV & hypSubMS)
        } else {
            subArgsUV | hypSubUV
            subArgsMS | hypSubMS
        }
        if (params.noionosub) {
            channel.from([]) | (hypIonoSubUV & hypIonoSubMS)
        } else {
            subArgsUV.map {obsid, meta, metafits, vis, srclist ->
                def newMeta = deepcopy(meta)
                newMeta.ionosub_nsrcs = params.ionosub_nsrcs
                [obsid, newMeta, metafits, vis, srclist]}
                | hypIonoSubUV
            subArgsMS.map {obsid, meta, metafits, vis, srclist ->
                def newMeta = deepcopy(meta);
                newMeta.ionosub_nsrcs = params.ionosub_nsrcs
                [obsid, newMeta, metafits, vis, srclist]}
                | hypIonoSubMS
        }
        // make cthulhu plots
        hypIonoSubUV.out.cross(hypSrclistYaml.out) {it[0]}
            .map { hypIonoSubUV_, hypSrclistYaml_ ->
                def (obsid, meta, _, offsets, __) = hypIonoSubUV_;
                def (___, srclist) = hypSrclistYaml_;
                [obsid, meta, srclist, offsets]
            }
            | cthulhuPlot

        // channel of calibrated, subtracted and ionosubtracted uvfits: tuple(obsid, name, uvfits)
        obsMetaUV = hypApplyUV.out.map { obsid, meta, vis, _ -> [obsid, meta, vis] }
            .mix(hypIonoSubUV.out.map { obsid, meta, vis, _, __ -> [obsid, meta, vis] })
            .mix(hypSubUV.out.map { obsid, meta, vis, _ -> [obsid, meta, vis] })

        // improve uvfits meta
        obsMetaUV | uvMeta
        obsMetaUVSmart = obsMetaUV.cross(uvMeta.out) {[it[0], it[1].name]}
            .map { obsMetaUV_, uvMeta_ ->
                def (obsid, meta, vis) = obsMetaUV_; def (_, __, json) = uvMeta_;
                def jsonMeta = parseJson(json)
                def meta2 = [
                    nchans:jsonMeta.freqs.size(),
                    ntimes:jsonMeta.times.size()
                ];
                [obsid, meta + meta2, vis] }

        // QA uvfits visibilities
        obsMetaUVSmart | uvfits

        // channel of calibrated, subtracted and ionosubtracted ms: tuple(obsid, name, vis)
        obsMetaMS = hypApplyMS.out.map { obsid, meta, vis, _ -> [obsid, meta, vis] }
            .mix(hypIonoSubMS.out.map { obsid, meta, vis, _, __ -> [obsid, meta, vis] })
            .mix(hypSubMS.out.map { obsid, meta, vis, _ -> [obsid, meta, vis] })

        // image and qa measurementsets or uvfits unless --noimage
        if (params.noimg) {
            channel.from([]) | img
        } else {
            // visibilities for imaging
            // prefer measurement sets if availabe, otherwise importuvfits
            // filter by visibilities that pass uvfits QA if --nouv is not set
            if (params.nouv) {
                obsMetaVisPass = obsMetaMS
            } else {
                obsMetaVisPass = obsMetaUVSmart
            }
            // if (params.noms) {
            //     obsMetaVisPass = uvfits.out.passVis
            //         .cross(obsMetaUVSmart.map {
            //             def (obs, meta, vis) = it
            //             [obs, meta.name, meta, vis]
            //         }) { it[0..1] }
            //         .map { def (obs, _, meta, vis) = it[1]; [obs, meta, vis] }
            //         .view { it -> "obsMetaVisPass ${it}"}
            // } else {
            //     obsMetaVisPass = uvfits.out.passVis
            //         .cross(obsMetaMS.map {
            //             def (obs, meta, vis) = it
            //             [obs, meta.name, meta, vis]
            //         }) { it[0..1] }
            //         .map { def (obs, _, meta, vis) = it[1]; [obs, meta, vis] }
            // }
            // group obsids by groupid and pointing for imaging
            imgEpochs = obsMetaVisPass
                .map { obs, meta, vis ->
                    def epoch = obs[0..5]
                    ["e${epoch}X", "e_${meta.name}", meta, vis]
                }
                .groupTuple(by: [0, 1])
                .filter { _, __, ___, viss -> viss.size() > 1 }
                .map { group, name, metas, viss ->
                    def ntimes = metas.collect { it.ntimes?:0 }.sum()
                    def newMeta = deepcopy(metas[0])
                    newMeta.name = name
                    newMeta.ntimes = ntimes
                    [group, newMeta, viss.unique()] }
            obsMetaVisPass
                .map { obs, meta, vis -> [obs, meta, [vis]]}
                // .mix( imgEpochs )
                | img
            imgEpochs
                .map { group, meta, viss -> ([group, meta.name, viss.size()]).join("\t") }
                .collectFile(
                    name: "img_epoch.tsv", newLine: true, sort: true,
                    seed: ([ "EPOCH", "NAME", "VISS" ]).join("\t"),
                    storeDir: "${results_dir}${params.img_suffix}${params.cal_suffix}"
                )
        }

        // archive data to object store
        if (params.archive) {
            if (params.archive_prep) {
                prep_archive = obsMetafitsVis.map { _, __, vis -> ["prep", vis] }
            } else {
                prep_archive = channel.from([])
            }
            if (params.archive_uvfits) {
                vis_archive = obsMetaUV.map { _, __, vis -> ["uvfits", vis]}
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
        frame.mix(cal.out.frame)
            .mix(uvfits.out.frame)
            .mix(img.out.frame)
            .mix(cthulhuPlot.out.flatMap {_, meta, pngs -> pngs.collect{ png -> ["cthulhuplot_${meta.name}", png]}}.groupTuple())
            .map { name, frames -> [name, frames.sort()] }
            | ffmpeg
            | view { [it, it.size()] }
}

// entrypoint: get externally preprocessed uvfits files from csv file and run qa
workflow extPrep {
    obsVis = channel.of(file("external${params.external_suffix}.csv"))
        .splitCsv(header: ["obsid", "vis"])
        .filter { !it.obsid.startsWith('#') }
        .map { [it.obsid, file(it.vis)] }

    obsids = obsVis.map { obsid, _ -> obsid }
    obsids | ws

    qaPrep(ws.out.obsMetafits.join(obsVis), channel.from([]), ws.out.frame)
}

