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

def coerceList(x) {
    if (x instanceof List) {
        x
    } else {
        [x]
    }
}

// download observation metadata from webservices in json format
process wsMeta {
    input:
    val(obsid)
    output:
    tuple val(obsid), path(wsmeta), path(wsfiles)

    label "rate_limit"

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

process tapMeta {
    input:
    val(obsid)
    output:
    tuple val(obsid), path(tapmeta)

    label "tap"

    // persist results in outdir, process will be skipped if files already present.
    storeDir "${params.outdir}/meta"
    // tag to identify job in squeue and nf logs
    tag "${obsid}"

    script:
    tapmeta = "${obsid}_tapmeta.json"
    obs_fields = [
        'starttime_mjd',
        'starttime_utc',
        'mwa_array_configuration'
    ]
    template "tapmeta.py"
}

// download observation metadata from webservices in metafits format
process wsMetafits {
    input:
    val(obsid)
    output:
    tuple val(obsid), path(metafits)

    label "rate_limit"

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

    label "rate_limit"

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

    label "rate_limit"

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

process metaJson {
    input:
        tuple val(obsid), path(metafits)
    output:
        tuple val(obsid), path(json), path(tsv)

    storeDir "${params.outdir}/${obsid}/meta"
    tag "${obsid}"

    errorStrategy 'terminate'

    label "python"

    script:
    // metrics = "${obsid}_occupancy.json"
    json = "${obsid}_meta.json"
    tsv = "${obsid}_inputs.tsv"
    txt = "${obsid}_inputs.txt"
    template "metajson.py"
}

// Download preprocessed files from asvo
//
// If the preprocessed files are not present, the `asvoPrep` process schedules
// a conversion job on ASVO using [Giant Squid](github.com/mwaTelescope/giant-squid).
// An exponential backoff is used to wait until the conversion job is ready to
// download from [Acacia](https://pawsey.org.au/systems/acacia/), then the
// archive is downloaded with `wget`, and hash-validated as it is inflated with
// `tar`.
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

    label "rate_limit"
    label "nvme"
    label "mem_half"

    time { 1.5.hour * task.attempt }

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
    prefix = "birli_"
    suffix = "_${params.prep_time_res_s}s_${params.prep_freq_res_khz}kHz"
    args = [
        output: "uvfits"
    ]
    if (params.prep_time_res_s != null) {
        args.avg_time_res = "${params.prep_time_res_s}"
    }
    if (params.prep_freq_res_khz != null) {
        args.avg_freq_res = "${params.prep_freq_res_khz}"
    }
    if (params.prep_rfi != null && !params.prep_rfi) {
        args.no_rfi = "true"
        suffix += '_norfi'
    }
    argstr = args.collect { k, v -> "${k}=${v}" }.join(',');
    uvfits = "${prefix}${obsid}*${suffix}.uvfits"

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
    ${params.giant_squid} submit-conv -v -p ${argstr} ${obsid} || true

    # list pending conversion jobs
    ${params.giant_squid} list -j --types conversion --states queued,processing,error -- ${obsid}

    # extract id url size hash from any ready download vis jobs for this obsid
    ${params.giant_squid} list -j --types conversion --states ready -- ${obsid} \
        | tee /dev/stderr \
        | ${params.jq} -r '.[]|[.jobId,.files[0].fileUrl//"",.files[0].fileSize//"",.files[0].fileHash//""]|@tsv' \
        | sort -r \
        | tee ready.tsv

    # download the most recent ready job
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
        for f in ${obsid}*.uvfits; do
            mv "\$f" "${prefix}\${f%.uvfits}${suffix}.uvfits"
        done
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

    storeDir "${params.outdir}/${obsid}/prep_qa"

    tag "${obsid}"

    label "python"
    label "nvme"
    label "mem_half"

    script:
    // metrics = "${obsid}_occupancy.json"
    names = coerceList(uvfits).collect { f -> f.name }
    metrics = names.size > 1 ? "occupancy_{${names.join(',')}}.json" : "occupancy_${names[0]}.json"
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

    storeDir "${params.outdir}/${obsid}/prep_qa"

    tag "${obsid}"

    label "ssins"
    label "nvme"
    label "mem_half"

    // errorStrategy "terminate"

    script:
    ssins_uvfits = "ssins_${obsid}_${params.prep_time_res_s}s_${params.prep_freq_res_khz}kHz.uvfits"
    guard_width = params.prep_freq_res_khz * 500
    title = "${obsid}"
    template "ssins.py"
}

process autoplot {
    input:
    tuple val(obsid), val(meta), path(metafits), path(uvfits)
    output:
    tuple val(obsid), val(meta), path(autoplot)

    storeDir "${params.outdir}/${obsid}/prep_qa"

    tag "${obsid}${suffix}"

    label "python"
    label "nvme"
    label "mem_half"

    script:
    suffix = meta.suffix?:""
    autoplot = "autoplot_${obsid}${suffix}.png"
    args = "${meta.autoplot_args}"
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
    label "hyperdrive"
    label "cpu_full"
    time 15.minute

    script:
    reduced = "${obsid}_reduced_n${params.sub_nsrcs}.txt"
    """
    #!/bin/bash -eux

    # Reduce a sky-model source list to the top N brightest sources, given pointing information
    ${params.hyperdrive} srclist-by-beam ${params.hyp_srclist_args} \
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
    label "hyperdrive"
    label "cpu_full"
    time 15.minute

    script:
    reduced = "${obsid}_reduced_n${params.sub_nsrcs}.yaml"
    """
    #!/bin/bash -eux

    # Reduce a sky-model source list to the top N brightest sources, given pointing information
    ${params.hyperdrive} srclist-by-beam ${params.hyp_srclist_args} \
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
    label "nvme"
    label "mem_full"
    label "cpu_full"

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
        ps=("\${PIPESTATUS[@]}")
        if [ \${ps[0]} -ne 0 ]; then
            echo "mwa_reduce failed. status=\${ps[0]}"
            exit \${ps[0]}
        fi
    fi
    """
}

// ensure calibration solutions are present, or calibrate prep with hyperdrive
// do multiple calibration solutions for each obs, depending on dical_args
// meta is a hashmap containing extra arguments, think `**kwargs` in Python
process hypCalSol {
    input:
    tuple val(obsid), val(meta), path(metafits), path(uvfits), val(dical_args)
    output:
    tuple val(obsid), val(meta), path("hyp_soln_${obsid}_${name_glob}.fits"), path("hyp_di-cal_${obsid}_${name_glob}.log")
    // todo: model subtract: path("hyp_model_${obsid}_${name_glob}.uvfits")

    storeDir "${params.outdir}/${obsid}/cal${params.cal_suffix}"

    tag "${obsid}"

    // label jobs that need a bigger gpu allocation
    label "hyperdrive"
    label "mem_full"
    label "cpu_half"
    label "gpu_nvme"
    if (params.pullCalSol) {
        label "rclone"
    }

    time { 2.hour }

    script:
    dical_names = dical_args.keySet().collect()
    para = dical_names.size() > 1
    name_glob = para ? "{" + dical_names.join(',') + "}" : dical_names[0]
    flag_args = ""
    prepFlags = meta.prepflags?:[]
    fineChanFlags = meta.fineChanFlags?:[]
    if (prepFlags.size() > 0) {
        flag_args += " --tile-flags ${prepFlags.join(' ')}"
    }
    if (fineChanFlags.size() > 0) {
        flag_args += " --fine-chan-flags-per-coarse-chan ${fineChanFlags.join(' ')}"
    }

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
                --data "${metafits}" ${uvfits} \
                --beam-file "${params.beam_path}" \
                --source-list "${params.sourcelist}" \
                --outputs \$soln_name \
                ${flag_args} \
                | tee \$log_name
            ps=("\${PIPESTATUS[@]}")
            if [ \${ps[0]} -ne 0 ]; then
                echo "hyperdrive failed. status=\${ps[0]}"
                exit \${ps[0]}
            fi
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
    tag "${obsid}.${meta.dical_name}"

    label "python"

    script:
    meta = deepcopy(meta) + [name: "poly_${meta.name}"]
    poly_soln = "${meta.cal_prog}_soln_${obsid}_${meta.name}.fits"
    logs = "polyfit_${obsid}_${meta.name}.log"
    """
    #!/bin/bash -eux
    run_polyfit.py "${soln}" --outfile "${poly_soln}" | tee "${logs}"
    ps=("\${PIPESTATUS[@]}")
    if [ \${ps[0]} -ne 0 ]; then
        echo "run_polyfit.py failed. status=\${ps[0]}"
        exit \${ps[0]}
    fi
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
    tuple val(obsid), val(newMeta), path(cal_vis), path(logs)

    storeDir "${params.outdir}/${obsid}/cal${params.cal_suffix}"

    tag "${obsid}.${meta.name}"
    label "hyperdrive"
    label "cpu_half"
    label "nvme"
    label "mem_half"

    script:
    newMeta = deepcopy(meta) + [name: "${meta.name}_${meta.apply_name}"]
    cal_vis = "hyp_${obsid}_${newMeta.name}.uvfits"
    logs = "hyp_apply_${newMeta.name}.log"
    """
    #!/bin/bash -eux
    ${params.hyperdrive} solutions-apply ${newMeta.apply_args} \
        --time-average=${newMeta.time_res}s \
        --freq-average=${newMeta.freq_res}kHz \
        --data "${metafits}" "${vis}" \
        --solutions "${soln}" \
        --outputs "${cal_vis}" \
        ${newMeta.nodut1 ? "--ignore-dut1" : ""} \
        | tee "${logs}"
    ps=("\${PIPESTATUS[@]}")
    if [ \${ps[0]} -ne 0 ]; then
        echo "Hyperdrive failed. status=\${ps[0]}"
        exit \${ps[0]}
    fi
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
    tuple val(obsid), val(newMeta), path(cal_vis), path(logs)

    storeDir "${params.outdir}/${obsid}/cal${params.cal_suffix}"
    // storeDir "/data/curtin_mwaeor/FRB_hopper/"

    tag "${obsid}.${meta.name}"
    label "hyperdrive"
    label "cpu_half"
    label "gpu"

    script:
    newMeta = deepcopy(meta) + [name: "${meta.name}_${meta.apply_name}"]
    cal_vis = "hyp_${obsid}_${newMeta.name}.ms"
    logs = "hyp_apply_${newMeta.name}_ms.log"
    """
    #!/bin/bash -eux
    # hyperdrive solutions apply ms
    ${params.hyperdrive} solutions-apply ${newMeta.apply_args} \
        --time-average=${newMeta.time_res}s \
        --freq-average=${newMeta.freq_res}kHz \
        --data "${metafits}" "${vis}" \
        --solutions "${soln}" \
        --outputs "${cal_vis}" \
        ${params.nodut1 ? "--ignore-dut1" : ""} \
        | tee "${logs}"
    ps=("\${PIPESTATUS[@]}")
    if [ \${ps[0]} -ne 0 ]; then
        echo "Hyperdrive failed. status=\${ps[0]}"
        exit \${ps[0]}
    fi
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
    tuple val(obsid), val(newMeta), path(sub_vis), path(logs)

    storeDir "${params.outdir}/${obsid}/cal${params.cal_suffix}"

    tag "${obsid}.${old_name}"
    label "hyperdrive"
    label "cpu_half"
    label "gpu"

    script:
    old_name = meta.name
    newMeta = deepcopy(meta) + [sub: "sub", name: "sub_${meta.name}"]
    sub_vis = "hyp_${obsid}_${newMeta.name}.uvfits"
    logs = "hyp_vis-${newMeta.name}_uv.log"
    """
    #!/bin/bash -eux
    ${params.hyperdrive} vis-sub \
        --data "${metafits}" "${vis}" \
        --beam-file "${params.beam_path}" \
        --source-list "${srclist}" \
        --invert --num-sources ${newMeta.sub_nsrcs} \
        --outputs "${sub_vis}" \
        | tee "${logs}"
    ps=("\${PIPESTATUS[@]}")
    if [ \${ps[0]} -ne 0 ]; then
        echo "Hyperdrive failed. status=\${ps[0]}"
        exit \${ps[0]}
    fi
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
    tuple val(obsid), val(newMeta), path(sub_vis), path(logs)

    storeDir "${params.outdir}/${obsid}/cal${params.cal_suffix}"

    tag "${obsid}.${old_name}"
    label "hyperdrive"
    label "cpu_half"
    label "gpu"

    script:
    old_name = meta.name
    newMeta = deepcopy(meta) + [sub: "sub", name: "sub_${meta.name}"]
    sub_vis = "hyp_${obsid}_${newMeta.name}.ms"
    logs = "hyp_vis-${newMeta.name}_ms.log"
    """
    #!/bin/bash -eux
    ${params.hyperdrive} vis-sub \
        --data "${metafits}" "${vis}" \
        --beam-file "${params.beam_path}" \
        --source-list "${srclist}" \
        --invert --num-sources ${newMeta.sub_nsrcs} \
        --outputs "${sub_vis}" \
        | tee "${logs}"
    ps=("\${PIPESTATUS[@]}")
    if [ \${ps[0]} -ne 0 ]; then
        echo "Hyperdrive failed. status=\${ps[0]}"
        exit \${ps[0]}
    fi
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
    tuple val(obsid), val(newMeta), path(sub_vis), path(json), path(logs)

    storeDir "${params.outdir}/${obsid}/cal${params.cal_suffix}"

    tag "${obsid}.${old_name}"
    label "hyperdrive"
    label "cpu_half"
    label "gpu"

    time 45.minute

    script:
    old_name = meta.name
    newMeta = deepcopy(meta) + [sub: "ionosub", name: "ionosub_${meta.name}"]
    sub_vis = "hyp_${obsid}_${newMeta.name}.uvfits"
    logs = "hyp_vis-${newMeta.name}_uv.log"
    json = "hyp_peel_${obsid}_${newMeta.name}_uv.json"
    """
    #!/bin/bash -eux
    ${params.hyperdrive} peel \
        --data "${metafits}" "${vis}" \
        --beam-file "${params.beam_path}" \
        --source-list "${srclist}" \
        --iono-sub ${newMeta.ionosub_nsrcs} \
        --sub ${newMeta.sub_nsrcs} \
        --outputs "${sub_vis}" "${json}" \
        | tee "${logs}"
    ps=("\${PIPESTATUS[@]}")
    if [ \${ps[0]} -ne 0 ]; then
        echo "Hyperdrive failed. status=\${ps[0]}"
        exit \${ps[0]}
    fi
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
    tuple val(obsid), val(newMeta), path(sub_vis), path(json), path(logs)

    storeDir "${params.outdir}/${obsid}/cal${params.cal_suffix}"

    tag "${obsid}.${old_name}"
    label "hyperdrive"
    label "cpu_half"
    label "gpu"

    script:
    old_name = meta.name
    newMeta = deepcopy(meta) + [sub: "ionosub", name: "ionosub_${meta.name}"]
    sub_vis = "hyp_${obsid}_${newMeta.name}.ms"
    logs = "hyp_vis-${newMeta.name}_ms.log"
    json = "hyp_peel_${obsid}_${newMeta.name}_ms.json"
    """
    #!/bin/bash -eux
    ${params.hyperdrive} peel \
        --data "${metafits}" "${vis}" \
        --beam-file "${params.beam_path}" \
        --source-list "${srclist}" \
        --iono-sub ${newMeta.ionosub_nsrcs} \
        --sub ${newMeta.sub_nsrcs} \
        --outputs "${sub_vis}" "${json}" \
        -v \
        | tee "${logs}"
    ps=("\${PIPESTATUS[@]}")
    if [ \${ps[0]} -ne 0 ]; then
        echo "Hyperdrive failed. status=\${ps[0]}"
        exit \${ps[0]}
    fi
    """
}

process cthulhuPlot {
    input:
    tuple val(obsid), val(meta), path(srclist), path(offsets)
    output:
    tuple val(obsid), val(meta), path("cthuluplot_${title}*.png"), path(csv), path(json)

    storeDir "${params.outdir}/${obsid}/iono_qa${params.cal_suffix}"

    tag "${obsid}.${meta.name}"

    label "python"

    time 10.minute

    script:
    title = "${obsid}_${meta.name}"
    plot = "cthuluplot_${title}.png"
    csv = "ionoqa_${title}.csv"
    json = "ionoqa_${title}.json"
    extra = meta.time_res ? "--time_res=${meta.time_res}" : ""
    template "cthulhuplot.py"
}

// QA tasks that can be run on preprocessed visibility files.
process prepVisQA {
    input:
    tuple val(obsid), path(metafits), path(uvfits)
    output:
    tuple val(obsid), path(metrics)

    storeDir "${params.outdir}/${obsid}/prep_qa"

    tag "${obsid}"

    label "python"
    label "nvme"
    label "mem_half"

    script:
    // metrics = "birli_${obsid}_prepvis_metrics.json"
    basenames = coerceList(uvfits).collect { f -> f.baseName }
    metrics = basenames.size > 1 ? "{${basenames.join(','')}}_prepvis_metrics.json" : "${basenames[0]}_prepvis_metrics.json"
    """
    #!/bin/bash -eux
    for f in ${uvfits}; do
        run_prepvisqa.py \$f "${metafits}" --out "\${f%.uvfits}_prepvis_metrics.json"
    done
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

    label "python"

    script:
    // TODO: json = "vis_metrics_${meta.cal_prog}_${obsid}_${meta.name}.json"
    json = "${meta.cal_prog}_${obsid}_${meta.name}_vis_metrics.json"
    """
    #!/bin/bash -eux
    run_visqa.py ${uvfits} --out "${json}"
    """
}

process uvMeta {
    input:
    tuple val(obsid), val(meta), path(vis)
    output:
    tuple val(obsid), val(meta), path(uvmeta)

    storeDir "${params.outdir}/${obsid}/vis_qa${params.cal_suffix}"

    tag "${obsid}.${meta.name}"

    label "python"
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

    label "python"

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

    label "python"

    script:
    metrics = "${meta.cal_prog}_soln_${obsid}_${meta.name}.fits.json"
    template "soljson.py"
}

process plotPrepVisQA {
    input:
    tuple val(obsid), path(metrics)
    output:
    tuple val(obsid), path("${plot_base}_{rms,modz}.png")

    storeDir "${params.outdir}/${obsid}/prep_qa"

    tag "${obsid}"

    label "python"

    script:
    plot_base = "prepvis_metrics_${obsid}"
    """
    #!/bin/bash -eux
    plot_prepvisqa.py "${metrics}" --out "${plot_base}.png" --save
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
    hyperdrive solutions-plot ${params.hyp_sols_plot_args} -m "${metafits}" ${soln}
    """
}

process plotCalQA {
    input:
    tuple val(obsid), val(meta), path(metrics)
    output:
    tuple val(obsid), val(meta), path("${plot_base}_{rms,fft}.png")

    storeDir "${params.outdir}/${obsid}/cal_qa${params.cal_suffix}"

    tag "${obsid}.${meta.name}"

    label "python"

    script:
    plot_base = "calmetrics_${obsid}_${meta.name}"
    """
    #!/bin/bash -eux
    plot_calqa.py "${metrics}" --out "${plot_base}" --save
    """
}

process plotVisQA {
    input:
    tuple val(obsid), val(meta), path(metrics)
    output:
    tuple val(obsid), val(meta), path("${meta.cal_prog}_${obsid}_${meta.name}_vis_metrics_rms.png")

    storeDir "${params.outdir}/${obsid}/vis_qa${params.cal_suffix}"

    tag "${obsid}"

    label "python"

    script:
    """
    #!/bin/bash -eux
    plot_visqa.py "${metrics}" --out "${meta.cal_prog}_${obsid}_${meta.name}_vis_metrics_rms.png" --save
    """
}

process plotImgQA {
    input:
    tuple val(name), path(jsons)
    output:
    tuple val(name), path("${base}_*.png")

    storeDir "${results_dir}"
    stageInMode "symlink"

    tag "${name}"

    label "python"

    script:
    base = "wsclean_hyp_${name}-MFS"
    """
    #!/bin/bash
    set -ex
    plot_imgqa.py --out "${base}" --save ${jsons.join(' ')}
    """
}

process delaySpec {
    input:
    tuple val(obsid), val(meta), path(uvfits)
    output:
    tuple val(obsid), val(meta), path(dlyspec)

    storeDir "${params.outdir}/${obsid}/vis_qa${params.cal_suffix}"

    tag "${obsid}.${meta.name}"

    label "python"

    script:
    title = "${obsid}_${meta.name}"
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

    tag "${obsid}${meta.inter_tok?:''}.${meta.name}"
    label "wsclean"
    label "cpu_half"
    label "mem_half"
    label "nvme"

    time { 1.8.minute * (1 + (multiplier * pix_mult * chan_mult * inter_mult)) }

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
        img_glob += "-{MFS,????}"
        // img_glob += "-*"
    }
    img_glob += "-{XX,YY,XY,XYi,I,Q,U,V}-{dirty,uv-real,uv-imag}.fits"
    img_args = img_params.args
    if (img_params.weight) {
        img_args += " -weight ${img_params.weight}"
    }
    if (img_params.pol) {
        img_args += " -pol ${img_params.pol}"
    }
    """
    #!/bin/bash -eux
    """ + (
            // convert any uvfits to ms
            [vis, vis_ms].transpose().collect { uv, ms ->
                (uv.extension == 'uvfits' ? \
                """${params.casa} -c "importuvfits('${uv}', '${ms}')" ;
                rm -rf ${uv} ; """ : "")
            }.join("\n")
        ) + """
        """ + (
        // run chgcentre if params.chgcentre_args specified.
        params.chgcentre_args ? \
            vis_ms.collect {"${params.chgcentre} ${params.chgcentre_args} ${it}"}.join("\n") : \
            ""
    ) + """
    ${params.wsclean} \
        ${img_args} \
        -name wsclean_${img_name} \
        -size ${img_params.size} ${img_params.size} \
        -scale ${img_params.scale} \
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
    tuple val(obsid), val(meta), path(vis), val(img_params), path(dirtyImgs)
    output:
    tuple val(obsid), val(meta), path(img_glob)
        // path("wsclean_${img_name}-sources.txt") <- only works for stokes I

    storeDir "${params.outdir}/${obsid}/img${img_params.suffix}${params.cal_suffix}"

    tag "${obsid}${meta.inter_tok?:''}.${meta.name}"
    label "wsclean"
    label "cpu_full"
    label "mem_full"
    label "nvme"

    time { 15.min * (1 + (multiplier * pix_mult * chan_mult * iter_mult * inter_mult)) }

    script:
    multiplier = Math.sqrt(vis.collect().size())
    img_name = img_params.name?:''
    if (img_name == '') {
        img_name = "${meta.cal_prog}_${obsid}${meta.subobs?:''}_${meta.name}${meta.inter_tok?:''}"
    }

    // multipliers for determining compute resources
    pix_mult = 1 + (img_params.size / 1024)
    chan_mult = 1 + ("${img_params.channels_out}".split(' ')[0] as int) / 25
    if (meta.interval) {
        inter_mult = meta.interval[1] - meta.interval[0]
    } else {
        inter_mult = 1 + ("${img_params.intervals_out}".split(' ')[0] as int) / 3
    }
    iter_mult = 1 + Math.sqrt(img_params.niter as Double) / 1000
    vis_ms = vis.collect {"${it.baseName}.ms"}
    vis = vis.collect()
    img_glob = img_params.glob?:''
    if (img_glob == '') {
        img_glob = "wsclean_${img_name}"
        if (multiinterval && !meta.inter_tok) {
            img_glob += "-t????"
        }
        if (multichannel) {
            img_glob += "-MFS"
        }
        if (img_params.pol) {
            img_glob += "-{XX,YY,XY,XYi,I,Q,U,V}"
        }
        img_glob += "-image.fits"
    }
    img_args = img_params.args
    if (img_params.weight) {
        img_args += " -weight ${img_params.weight}"
    }
    if (img_params.pol) {
        img_args += " -pol ${img_params.pol}"
    }
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
        ${img_args} \
        -name wsclean_${img_name} \
        """ + (dirtyImgs?"-reuse-dirty wsclean_${img_name}":"") +
        """ \
        -size ${img_params.size} ${img_params.size} \
        -scale ${img_params.scale} \
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

    label "python"

    time {5.min * task.attempt}

    script:
    csv = "quantile_${meta.name}_${meta.suffix}.csv"
    template "img_meta.py"
}

// power spectrum metrics via chips
process psMetrics {
    input:
    tuple val(obsid), val(meta), path("${meta.cal_prog}_${obsid}_${meta.name}.uvfits")
    output:
    tuple val(obsid), val(meta), path(out_metrics)

    storeDir "${params.outdir}/${obsid}/ps_metrics${params.cal_suffix}"

    tag "${obsid}.${meta.name}"

    label "chips"
    label "cpu_half"

    time 20.minute

    script:
    uvbase = "${meta.cal_prog}_${obsid}_${meta.name}"
    nchans = meta.nchans
    eorband = meta.eorband
    out_metrics = "output_metrics_${uvbase}.dat"
    """
    #!/bin/bash -eux
    export DATADIR="\$PWD"
    export OUTPUTDIR="\$PWD/"
    export OMP_NUM_THREADS=${task.cpus}
    ${params.ps_metrics} "${uvbase}" "${eorband}" "${nchans}" "${uvbase}" 2>&1 \
        | tee "${uvbase}.log"
    ps=("\${PIPESTATUS[@]}")
    if [ \${ps[0]} -ne 0 ]; then
        echo "ps_metrics failed. status=\${ps[0]}"
        exit \${ps[0]}
    fi
    """
}

process storeManifest {
    input:
    tuple val(group), val(obsids)
    output:
    path(manifest)

    tag "${group}"

    time 5.minutes

    storeDir "${params.outdir}/${group}"

    script:
    manifest = "manifest_${group}.csv"
    """
    cat <<EOF > ${manifest}
${obsids.join('\n')}
EOF"""
}

// power spectrum with chips
process chipsGrid {
    input:
    tuple val(group), val(meta), val(obsids), path(viss)
    output:
    tuple val(group), val(newMeta), \
        path("{vis_tot,vis_diff,noise_tot,noise_diff,weights}_{xx,yy}.${ext}.dat")

    storeDir "${params.outdir}/${group}/ps_metrics${params.cal_suffix}/${meta.name}"

    tag "${group}.${meta.name}"

    label "chips"
    label "cpu_full"
    label "mem_full"

    time { 30.minute * obsids.size() }

    script:
    ext = meta.ext
    nchans = meta.nchans
    eorband = meta.eorband
    eorfield = meta.eorfield
    lowfreq = meta.lowfreq
    period = "8.0"
    freq_res = "${(meta.freq_res?:80) * 1000}"
    // pol = meta.pol?:"xx"
    freq_idx_start = 0

    newMeta = deepcopy(meta) + [lowfreq: lowfreq]

    """
    #!/bin/bash -eux

    echo "meta=${meta}"
    """ + (eorfield == null || eorband == null || !nchans || !ext || lowfreq == null || !freq_res || freq_idx_start == null || period == null ? "exit 68" : "" ) + """

    export DATADIR="\$PWD"
    export INPUTDIR="\$PWD/"
    export OUTPUTDIR="\$PWD/"
    export OBSDIR="\$PWD/"
    export OMP_NUM_THREADS=${task.cpus}

    # copy data files
    for f in present_vals.dat missing_vals.dat krig_weights.dat; do
        [ -f \$f ] || cp "/astro/mwaeor/ctrott/output/\$f" .
    done

    """ + (
        [coerceList(obsids), coerceList(viss)].transpose().collect { obsid, vis ->
            """gridvisdiff "${vis}" "${obsid}" "${ext}" "${eorband}" -f "${eorfield}" \
                2>&1 > syslog_gridvisdiff_${obsid}.txt """
        }.join("\n")
    ) + """

    for pol in xx yy; do
        prepare_diff "${ext}" "${nchans}" "${freq_idx_start}" "\${pol}" "${ext}" "${eorband}" -p "${period}" -c "${freq_res}" -n "${lowfreq}" \
            2>&1 > syslog_prepare_diff_\${pol}.txt
    done
    """
}

process chipsCombine {
    input:
    tuple val(group), val(meta), val(exts), path(grids)

    output:
    tuple val(group), val(meta), path("combine_{xx,yy}.${combined_ext}.txt"),
        path("{vis_tot,vis_diff,noise_tot,noise_diff,weights}_{xx,yy}.${combined_ext}.dat")

    storeDir "${params.outdir}/${group}/ps_metrics${params.cal_suffix}/${meta.name}"

    tag "${group}.${meta.name}"

    label "chips"
    label "cpu_full"
    label "mem_full"
    label "nvme"

    time { 1.hour }

    script:
    // ext = "${group}_${meta.name}"
    combined_ext = meta.ext
    nchans = meta.nchans
    """
    #!/bin/bash -eux

    echo "meta=${meta}"
    """ + (!nchans || !combined_ext || exts.size() == 0 ? "exit 68" : "" ) + """

    export DATADIR="\$PWD"
    export INPUTDIR="\$PWD/"
    export OUTPUTDIR="\$PWD/"
    export OBSDIR="\$PWD/"
    export OMP_NUM_THREADS=${task.cpus}

    # copy data files
    for f in present_vals.dat missing_vals.dat krig_weights.dat; do
        [ -f \$f ] || cp "/astro/mwaeor/ctrott/output/\$f" .
    done

    for pol in xx yy; do
        echo -n "" > combine_\${pol}.${combined_ext}.txt
        """ + (
            exts.collect { ext ->
                """echo "\${pol}.${ext}" >> combine_\${pol}.${combined_ext}.txt"""
            }.join("\n")
        ) + """

        combine_data "combine_\${pol}.${combined_ext}.txt" ${nchans} "\${pol}.${combined_ext}" 1

        # todo: check first bytes not null
    done
    """
}

process chipsLssa {
    input:
    tuple val(group), val(meta), path(grid)

    output:
    tuple val(group), val(newMeta),
        path("{crosspower,residpower,residpowerimag,totpower,flagpower,fg_num,outputweights}_{xx,yy}_${bias_mode}.iter.${ext}.dat")

    storeDir "${params.outdir}/${group}/ps_metrics${params.cal_suffix}/${meta.name}"

    tag "${group}.${meta.name}"

    label "chips"
    label "cpu_full"
    label "mem_full"
    label "nvme"

    time 1.hour

    script:
    ext = meta.ext
    nchans = meta.nchans
    eorband = meta.eorband
    // pol = meta.pol?:"xx"
    maxu = 300
    nbins = 80
    freq_idx_start = 0
    bias_mode = 0

    newMeta = deepcopy(meta) + [nbins: nbins, maxu: maxu]

    """
    #!/bin/bash -eux
    echo "meta=${meta}"
    """ + (eorband == null || !nchans || nbins == null || !ext || maxu == null || bias_mode == null ? "exit 68" : "" ) + """

    export DATADIR="\$PWD"
    export INPUTDIR="\$PWD/"
    export OUTPUTDIR="\$PWD/"
    export OBSDIR="\$PWD/"
    export OMP_NUM_THREADS=${task.cpus}

    for pol in xx yy; do
        lssa_fg_simple "${ext}" "${nchans}" "${nbins}" "\${pol}" "${maxu}" "${ext}" "${bias_mode}" "${eorband}" \
            2>&1 > syslog_lssa_simple_\${pol}.txt

        export cross_size="\$(stat -c%s crosspower_\${pol}_${bias_mode}.iter.${ext}.dat)"
        if (( cross_size < 4097 )); then
            echo "crosspower_\${pol}_${bias_mode}.iter.${ext}.dat is too small (\$cross_size), exiting"
            exit 69
        fi
    done
    """
}

// TODO: process chipsCombine

process chipsPlot {
    input:
    tuple val(group), val(meta), path(grid)
    output:
    tuple val(group), val(meta), path("chips${dims}D_${pol}_${suffix}.png")

    storeDir "${params.outdir}/${group}/ps_metrics${params.cal_suffix}/${meta.name}"

    tag "${group}.${meta.name}.${ptype}"

    label "python"

    time 15.minute

    script:
    // nchans = meta.nchans
    // eorband = meta.eorband
    // eorfield = meta.eorfield
    // ext = meta.name
    // diff_ext = ""
    // lowfreq = "${meta.lowfreq?:167075000}"
    // freq_res = "${(meta.freq_res?:80) * 1000}"
    pol = meta.pol?:"both"
    if (pol == "both") {
        pol = "xx+yy"
    }
    ptype = meta.ptype?:""
    dims = ptype[0]?:2
    suffix = ""
    if (ptype =~ /.*comp/) {
        suffix = "_comparison"
    } else if (ptype =~ /.*diff/) {
        suffix = "_diff"
    } else if (ptype =~ /.*ratio/) {
        suffix = "_ratio"
    } else if (ptype[0] == "2") {
        suffix = "_crosspower"
    }
    if ((meta.tags?:[])[1]) {
        suffix = "_${meta.tags[1]}${suffix}"
    }
    if ((meta.tags?:[])[0]) {
        suffix = "${meta.tags[0]}${suffix}"
    } else {
        suffix = "${meta.ext}${suffix}"
    }
    basedir = "./"
    // title = "${group}\\n${ext}"
    // ptypes = "1D 2D"
    args = [
        title: meta.title,
        // file group
        basedir: "./",
        chips_tag: meta.ext,
        chips_tag_one: (meta.tags?:[])[0],
        chips_tag_two: (meta.tags?:[])[1],
        chips_tag_three: (meta.tags?:[])[2],
        // plot group
        plot_type: ptype,
        polarisation: meta.pol,
        min_power: meta.min_power,
        max_power: meta.max_power,
        // plot group 2D
        colourscale: meta.colourscale,
        max_neg_power: meta.max_neg_power,
        min_neg_power: meta.min_neg_power,
        // plot group 1D
        wedge_factor: meta.wedge_factor,
        low_k_edge: meta.low_k_edge,
        high_k_edge: meta.high_k_edge,
        num_k_edges: meta.nbins,
        kperp_max: meta.kperp_max,
        kperp_min: meta.kperp_min,
        kparra_min: meta.kparra_min,
        plot_wedge_cut_2D: meta.plot_wedge_cut_2D,
        chips_tag_one_label: (meta.labels?:[])[0],
        chips_tag_two_label: (meta.labels?:[])[1],
        chips_tag_three_label: (meta.labels?:[])[2],
        // chips group
        // N_kperp: meta.N_kperp,
        // N_chan: meta.N_chan,
        lowerfreq: meta.lowfreq,
        chan_width: (meta.freq_res * 1e3),
        umax: meta.maxu,
        // density_correction: meta.density_correction,
        // omega_matter: meta.omega_matter,
        // omega_baryon: meta.omega_baryon,
        // omega_lambda: meta.omega_lambda,
        // hubble: meta.hubble,

    ].findAll { _, v -> v != null }
        .collect { k, v -> """--${k} "${v}" """ }
        .join(" ")
    if (meta.plot_delta) {
        args += " --plot_delta"
    }
    if (meta.plot_wedge_cut) {
        args += " --plot_wedge_cut_2D"
    }
    if (meta.nchans) {
        args += " --N_chan ${meta.nchans}"
    }
    if (meta.kperp) {
        args += " --K_perp ${meta.kperp}"
    }

    template "jline_plotchips.py"
}

// analyse images of V,XX,YY
process imgQA {
    input:
    tuple val(obsid), val(meta), path(fits)
    output:
    tuple val(obsid), val(meta), path(json)

    storeDir "${params.outdir}/${obsid}/img_qa${params.img_suffix}${params.cal_suffix}"

    tag "${obsid}.${meta.name}"

    label "python"

    script:
    json = "wsclean_hyp_${obsid}_${meta.name}-MFS.json"
    """
    #!/bin/bash -eux
    run_imgqa.py ${fits.join(' ')} --out ${json}
    """
}

process uvPlot {
    input:
    tuple val(obsid), val(meta), path(uvfits)
    output:
    tuple val(obsid), val(meta), path("uvplot_${title}_{XX,YY,XY,YX}.png")

    storeDir "${params.outdir}/${obsid}/vis_qa${params.cal_suffix}"

    tag "${obsid}.${meta.name}"

    label "python"

    script:
    title = "${obsid}_${meta.name}"
    uvplot = "uvplot_${title}"
    template "uvplot_2d.py"
}

// make a thumbnail png from a fits image
process thumbnail {
    input:
    tuple val(obsid), val(meta), path(img)
    output:
    tuple val(obsid), val(meta), path(thumb)

    storeDir "${params.outdir}/${obsid}/img_qa${params.img_suffix}${params.cal_suffix}"

    tag "${obsid}.${meta.name}.${meta.suffix}"

    label "python"
    time {10.min * task.attempt}

    script:
    thumb = "${obsid}_${meta.name}_${meta.suffix}.png"
    args = [
        fits: img,
        title: meta.title?:"${obsid} ${meta.name} ${meta.suffix}",
        thumb: thumb,
    ] + meta
    // argstr = "${params.thumbnail_args}"
    argtokens = args.findAll { k, v ->
            v != null && [
                'fits', 'thumb', 'title', 'dpi', 'limit', 'vmin', 'vmax', 'vmin_quantile',
                'vmax_quantile', 'cmap', 'norm_args'
            ].contains(k)
        }
        .collect { k, v -> ["--${k}"] + coerceList(v) }
        .flatten()
    argtokens += args.findAll { k, v ->
            v && ['transparent', 'symmetric'].contains(k)
        }
        .collect { k, _ -> "--${k}" }

    argstr = argtokens
        .collect { token -> "\"${token}\"" }
        .join(' ')

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

    label "python"

    script:
    meta = deepcopy(meta) + ['orderName': meta.order.join('')]
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
        path(jsons)
    output:
        path("cal_qa.png")

    storeDir "${results_dir}"
    stageInMode "symlink"

    label "python"

    script:
    """
    #!/bin/bash -eux
    plot_calqa.py --out cal_qa --save ${jsons.join(' ')}
    """
}

process stackImgs {
    input:
        tuple val(chunk), val(meta), path(imgs)
    output:
        tuple val(chunk), val(meta), path(stack)
    storeDir "${params.outdir}/${chunk}/img${params.img_suffix}${params.cal_suffix}"

    tag "${chunk}${meta.subobs?:''}.${meta.name}"

    label 'imagemagick'

    script:
    stack = "stack_${meta.name}_${meta.pol}.fits"
    """
    #!/bin/bash -eux
    convert ${imgs.join(' ')} -average -auto-gamma -auto-level ${stack}
    """
}

process stackThumbnail {
    input:
    tuple val(chunk), val(meta), path(img)
    output:
    tuple val(chunk), val(meta), path("${chunk}_${meta.name}_${meta.suffix}.png")

    storeDir "${params.outdir}/${chunk}/img_qa${params.img_suffix}${params.cal_suffix}"

    tag "${chunk}.${meta.name}.${meta.suffix}"

    label "python"
    time {10.min * task.attempt}

    script:
    title = "${chunk} ${meta.name} ${meta.suffix}"
    thumb = "stack_${chunk}_${meta.name}_${meta.suffix}.png"
    args = [
        fits: img,
        title: title,
        thumb: thumb,
        vmin_quantile: 0.3,
        vmax_quantile: 0.99
    ]
    argstr = (args + meta).findAll { k, v ->
            v != null && [
                'fits', 'thumb', 'title', 'dpi', 'limit', 'vmin_quantile',
                'vmax_quantile', 'transparent', 'symmetric', 'cmap'
            ].contains(k)
        }
        .collect { k, v -> ["--${k}"] + coerceList(v) }
        .flatten()
        .collect { token -> "\"${token}\"" }
        .join(' ')

    template "thumbnail.py"
}

process tarchive {
    input:
        tuple val(name), path(files), val(cachebust)
    output:
        tuple path(out), path("${dot_cachebust}")

    storeDir "${results_dir}${params.img_suffix}${params.cal_suffix}"
    stageInMode "copy"

    tag "${name}"

    script:
    dot_cachebust = ".${name}.${cachebust}.cachebust"
    out = "${name}.tar.gz"
    """
    #!/bin/bash -eux
    touch ${dot_cachebust}
    tar cvzf ${name}.tar.gz ${files.join(' ')}
    """
}

process ffmpeg {
    input:
        tuple val(name), path("??????????.png"), val(cachebust)

    output:
        tuple path("${name}.mp4"), path("${dot_cachebust}")

    storeDir "${results_dir}${params.img_suffix}${params.cal_suffix}"
    stageInMode "symlink"

    tag "${name}"

    label "ffmpeg"

    time 1.hour

    script:
    dot_cachebust = ".${name}.${cachebust}.cachebust"
    """
    #!/bin/bash -eux
    touch ${dot_cachebust}
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

    // label "rclone"
    label "datamover"

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

def parseFloatOrNaN(s) {
    try {
        return Float.parseFloat(s) // .toFloat()?
    } catch (NumberFormatException e) {
        return Float.NaN
    }
}

SimpleDateFormat logDateFmt = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss");

def isNaN(n) {
    if (n==null) {
        return true
    }
    def cls = null
    try {
        cls = n.getClass();
        return n.isNaN()
    } catch (groovy.lang.MissingMethodException e) {
        print("WARN isNaN(${n})<${cls}>: ${e}")
        org.codehaus.groovy.runtime.StackTraceUtils.sanitize(e).printStackTrace()
        return true
    }
}

// display a long list of ints, replace bursts of consecutive numbers with ranges
def displayRange(Integer s, Integer e) {
    return s == e ? "${s}," : s == e - 1 ? "${s},${e}," : "${s}-${e},"
}
max_ints_display = 50
mid_ints_display = (max_ints_display-1).intdiv(2)

def displayInts(l_) {
    def l = (l_ as ArrayList).sort(false).unique()
    switch (l) {
        case { l.size == 0 }: return "";
        case { l.size == 1 }: return "${l[0]}";
        default:
            def sb, start, end
            (sb, start, end) = [''<<'', l[0], l[0]]
            for (i in l[1..-1]) {
                (sb, start, end) = i == end + 1 ? [sb, start, i] : [sb << displayRange(start, end), i, i]
            }
            def result = (sb << displayRange(start, end))[0..-2].toString()
            def size = result.size()
            if (size > max_ints_display) {
                try {
                    return result.substring(0, mid_ints_display)
                    << "..."
                    << result.substring(size - mid_ints_display)
                } catch (StringIndexOutOfBoundsException e) {
                    println("error: ${e}\n${result} ${size} ${mid_ints_display}")
                    throw e
                }
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


// POWER	P_win_sub, P_win	< 20	Small window power overall
def cmt_ps_metrics_pass(meta) {
    if (meta.p_window != null) {
        if (meta.p_window > params.filter_max_ps_window) {
            return "${meta.sub?:""}_p_win(${meta.p_window}) > max_ps_window(${params.filter_max_ps_window})"
        }
    }
    return null
}

// POWER	Normal P_win/P_wg	< 0.1	Window power a small fraction of wedge power
// POWER	P_wg_sub/P_wg	< 0.3	More than 70% wedge power subtracted
// POWER	P_win_sub/P_win	< 1.0, >0.1	Window power not crazy after subtraction
def cmt_ps_metrics_pass_sub(nosubMeta, subMeta) {
    def subReason = cmt_ps_metrics_pass(subMeta)
    if (subReason != null) {
        return subReason
    }
    if (nosubMeta.p_window != null && nosubMeta.p_wedge != null) {
        def p_win_p_wg = nosubMeta.p_window / nosubMeta.p_wedge
        if (p_win_p_wg > params.filter_max_ps_ratio) {
            return "p_win:p_wg(${p_win_p_wg}) > max_ps_ratio(${params.filter_max_ps_ratio})"
        }
    }
    if (subMeta.p_wedge != null && nosubMeta.p_wedge != null) {
        def sub_p_wg = subMeta.p_wedge / nosubMeta.p_wedge
        if (sub_p_wg > params.filter_max_ps_wedge_sub) {
            return "sub_p_wg:p_wg(${sub_p_wg}) > max_ps_wedge_sub{${params.filter_max_ps_wedge_sub}}"
        }
    }
    if (subMeta.p_window != null && nosubMeta.p_window != null) {
        def sub_p_win = subMeta.p_window / nosubMeta.p_window
        if (sub_p_win < 0.1) {
            return "sub_p_win:p_win(${sub_p_win}) < min_ps_win_sub(${params.filter_min_ps_win_sub})"
        }
        if (sub_p_win < 0.1) {
            return "sub_p_win:p_win(${sub_p_win}) > max_ps_win_sub(${params.filter_max_ps_win_sub})"
        }
    }
    return null
}

def cmt_imgqa_pass(meta) {
    if (params.filter_max_vrms_box!= null && meta.v_rms_box != null) {
        if (meta.v_rms_box > params.filter_max_vrms_box) {
            return "${meta.sub?:""}_v_rms_box(${meta.v_rms_box}) > max_vrms_box(${params.filter_max_vrms_box})"
        }
    }
    // if (meta.xx_pks_int != null && meta.yy_pks_int != null && meta.v_pks_int != null) {
    //     def pks_int_v_ratio = meta.v_pks_int / (meta.xx_pks_int + meta.yy_pks_int)
    //     if (pks_int_v_ratio > params.filter_max_pks_int_v_ratio) {
    //         return "${meta.sub?:""}_pks_int v/(xx+yy) (${pks_int_v_ratio}) > max_pks_int_v_ratio(${params.filter_max_pks_int_v_ratio})"
    //     }
    // }
    if (params.filter_max_pks_int_diff != null && meta.xx_pks_int != null && meta.yy_pks_int != null) {
        def pks_int_diff = (meta.xx_pks_int - meta.yy_pks_int).abs()
        if (pks_int_diff > params.filter_max_pks_int_diff) {
            return "${meta.sub?:""}_pks_int |xx+yy| (${pks_int_diff}) > max_pks_int_diff(${params.filter_max_pks_int_diff})"
        }
    }
    return null
}

// IMG	Vrms box	< 0.05	RMS V should be small
// IMG	V/(XX+YY) PKS int	< 0.001	V should be small compared with XX and YY
// IMG	PKS XX and YY	|XX-YY| < 10.0	XX and YY integrated should be similar
// IMG	XX_sub/XX integ	< 0.2	Most flux subtracted
// IMG	YY_sub/YY integ	< 0.2
// IMG	XX_sub integ	< 0.5	Integrated remaining flux after subtraction is small
// IMG	YY_sub integ	< 0.5	Integrated remaining flux after subtraction is small
def cmt_imgqa_pass_sub(nosubMeta, subMeta) {
    if (params.filter_max_pks_int_v_ratio != null && nosubMeta.xx_pks_int != null && nosubMeta.yy_pks_int != null && nosubMeta.v_pks_int != null) {
        def pks_int_v_ratio = nosubMeta.v_pks_int / (nosubMeta.xx_pks_int + nosubMeta.yy_pks_int)
        if (pks_int_v_ratio > params.filter_max_pks_int_v_ratio) {
            return "pks_int v/(xx+yy) (${pks_int_v_ratio}) > max_pks_int_v_ratio(${params.filter_max_pks_int_v_ratio})"
        }
    }
    if (params.filter_max_pks_int_sub_ratio != null && nosubMeta.xx_pks_int != null && subMeta.xx_pks_int != null) {
        def pks_int_sub_ratio = subMeta.xx_pks_int / nosubMeta.xx_pks_int
        if (pks_int_sub_ratio > params.filter_max_pks_int_sub_ratio) {
            return "${subMeta.sub?:""}_xx_pks_int:xx_pks_int (${pks_int_sub_ratio}) > max_pks_int_sub_ratio(${params.filter_max_pks_int_sub_ratio})"
        }
    }
    if (params.filter_max_pks_int_sub_ratio != null && nosubMeta.yy_pks_int != null && subMeta.yy_pks_int != null) {
        def pks_int_sub_ratio = subMeta.yy_pks_int / nosubMeta.yy_pks_int
        if (pks_int_sub_ratio > params.filter_max_pks_int_sub_ratio) {
            return "${subMeta.sub?:""}_yy_pks_int:yy_pks_int (${pks_int_sub_ratio}) > max_pks_int_sub_ratio(${params.filter_max_pks_int_sub_ratio})"
        }
    }
    if (params.filter_max_pks_int_sub != null && subMeta.xx_pks_int != null) {
        if (subMeta.xx_pks_int > params.filter_max_pks_int_sub) {
            return "${subMeta.sub?:""}_xx_pks_int(${subMeta.v_rms_box}) > max_pks_int_sub(${params.filter_max_pks_int_sub})"
        }
    }
    if (params.filter_max_pks_int_sub != null && subMeta.yy_pks_int != null) {
        if (subMeta.yy_pks_int > params.filter_max_pks_int_sub) {
            return "${subMeta.sub?:""}_yy_pks_int(${subMeta.v_rms_box}) > max_pks_int_sub(${params.filter_max_pks_int_sub})"
        }
    }
    return null
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
    if (tokens.size > 1 && tokens[-2] == "uv") {
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
        obsids | wsMeta & tapMeta

        quality_updates = file(params.quality_updates_path)
            .readLines()
            .findAll { !it.startsWith('#') && it.length() > 13 }
            .collectEntries { line ->
                def (obsid, quality, comment) = line.split(',')
                [obsid, [dataquality: quality, dataqualitycomment: comment]]
            }

        wsSummary = wsMeta.out.join(tapMeta.out).map { obsid, json, filesJson, tapJson ->
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

                def eorfield = null;
                if (ra_phase_center.round() == 0) {
                    def nearest_dec = dec_phase_center.round()
                    if (dec_phase_center.round() == -27) {
                        eorfield = 0
                    } else if (dec_phase_center.round() == 30) {
                        eorfield = 1
                    }
                }

                def coarse_chans = ((stats.rfstreams?:[:])["0"]?:[:]).frequencies?:[]
                def center_chan = coarse_chans[coarse_chans.size()/2]
                def eorband = null;
                if (center_chan == 143) {
                    eorband = 1
                } else if (center_chan == 108) {
                    eorband = 0
                }

                def pointing = stats.metadata.gridpoint_number
                def nscans = ((stats.stoptime?:0) - (stats.starttime?:0)) / (stats.int_time?:1)
                def delays = (stats.alldelays?:[:]).values().flatten()
                def quality = stats.quality?:[:]
                def tiles = stats.tdict?:[:]
                def tile_nums = tiles.collect { k, _ -> k as int }
                def tile_names = tiles.collect { _, v -> v[0] }
                def tile_rxs = tiles.collect { _, v -> v[1] }
                def n_tiles = tile_names.size()
                def n_lb = tile_names.count { it =~ /(?i)lb/ }
                def n_hex = tile_names.count { it =~ /(?i)hex/ }
                def config = "ph1"
                if (n_tiles > 128) {
                    config = "ph3"
                } else if (n_lb > 50) {
                    config = "ph2b"
                } else if (n_hex > 50) {
                    config = "ph2a"
                }
                config += "-${n_tiles}T"

                def bad_tiles = stats.bad_tiles?:[:]
                def n_bad_tiles = bad_tiles.size()
                def bad_tile_frac = n_bad_tiles / n_tiles
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
                // TODO: filter_sweets
                // if (params.filter_sweets && !params.filter_sweets.contains()) {
                //     fail_reasons += ["pointing=${pointing}"]
                // }
                if (params.filter_quality != null && dataquality > params.filter_quality) {
                    fail_reasons += ["dataquality=${dataquality} (${dataqualitycomment})"]
                }
                if (params.filter_bad_tile_frac && bad_tile_frac > params.filter_bad_tile_frac) {
                    fail_reasons += ["bad_tiles(${bad_tiles.size()})=${displayInts(bad_tiles)}"]
                }
                if (params.filter_dead_dipole_frac && dead_dipole_frac > params.filter_dead_dipole_frac) {
                    fail_reasons += ["dead_dipole_frac=${dead_dipole_frac}"]
                }
                if (params.filter_eorfield != null && eorfield != params.filter_eorfield) {
                    fail_reasons += ["phase_radec=(${ra_phase_center}, ${dec_phase_center})"]
                }
                if (params.filter_eorband != null && eorband != params.filter_eorband) {
                    fail_reasons += ["center_chan=${center_chan}"]
                }
                if (params.filter_ionoqa && quality.iono_qa != null && quality.iono_qa > params.filter_ionoqa) {
                    fail_reasons += ["iono_qa>${params.filter_ionoqa}"]
                }
                // if (ra_phase_center != 0.0) {
                //     fail_reasons += ["ra_phase_center=${ra_phase_center}"]
                // }
                // if (dec_phase_center != -27.0) {
                //     fail_reasons += ["dec_phase_center=${dec_phase_center}"]
                // }

                def tapStats = parseJson(tapJson)

                def summary = [
                    fail_reasons: fail_reasons,
                    // obs metadata
                    obs_name: obs_name,
                    groupid: groupid,
                    corrmode: stats.mode,
                    delaymode: stats.delaymode_name,
                    starttime_mjd: tapStats.starttime_mjd,
                    starttime_utc: tapStats.starttime_utc,

                    // pointing
                    ra_pointing: stats.metadata.ra_pointing,
                    dec_pointing: stats.metadata.dec_pointing,
                    az_pointing: stats.metadata.azimuth_pointing,
                    el_pointing: stats.metadata.elevation_pointing,
                    ra_phase_center: ra_phase_center,
                    dec_phase_center: dec_phase_center,
                    pointing: pointing,
                    lst: stats.metadata.local_sidereal_time_deg,
                    eorfield: eorfield,

                    // channels
                    freq_res: stats.freq_res,
                    coarse_chans: coarse_chans,
                    eorband: eorband,

                    // times
                    int_time: stats.int_time,
                    nscans: nscans,

                    // tiles
                    config: config,
                    tile_nums: tile_nums,
                    tile_rxs: tile_rxs,
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
                    summary.starttime_mjd?:'',
                    summary.starttime_utc?:'',
                    summary.groupid,
                    summary.corrmode,
                    summary.delaymode?:'',

                    // pointing
                    summary.ra_pointing,
                    summary.dec_pointing,
                    summary.az_pointing,
                    summary.el_pointing,
                    summary.ra_phase_center,
                    summary.dec_phase_center,
                    summary.pointing,
                    summary.lst,
                    summary.obs_name,
                    summary.eorfield?:'',

                    // channels
                    summary.freq_res,
                    displayInts(summary.coarse_chans),
                    summary.eorband?:'',

                    // times
                    summary.int_time,
                    summary.nscans,

                    // config
                    summary.config,
                    summary.n_dead_dipoles,
                    summary.dead_dipole_frac,
                    summary.n_bad_tiles,
                    summary.bad_tile_frac,
                    displayInts(summary.tile_nums),
                    displayInts(summary.bad_tiles),
                    displayInts(summary.tile_rxs),

                    // iono quality
                    summary.iono_magnitude,
                    summary.iono_pca,
                    summary.iono_qa,

                    // archive
                    summary.num_data_files,
                    summary.num_data_files_archived,

                    // errors
                    summary.dataquality,
                    summary.dataqualitycomment,
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
                    "OBS", "FAIL REASON", "START MJD", "START UTC",
                    "GROUP ID", "CORR MODE", "DELAY MODE",
                    "RA POINT", "DEC POINT", "AZ POINT", "EL POINT", "RA PHASE", "DEC PHASE", "POINT", "LST DEG", "OBS NAME", "EOR FIELD",
                    "FREQ RES", "COARSE CHANS", "EOR BAND",
                    "TIME RES", "N SCANS",
                    "CONFIG", "N DEAD DIPOLES", "DEAD DIPOLE FRAC", "N FLAG TILES", "FLAG TILES FRAC", "TILE NUMS", "FLAG TILES", "TILE RXS",
                    "IONO MAG", "IONO PCA", "IONO QA",
                    "N FILES", "N ARCHIVED",
                    "QUALITY", "QUALITY COMMENT",
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
        if (params.noppds) {
            channel.empty() | wsPPDs
        } else {
            pass | wsPPDs
        }

        wsMetafits.out | metaJson

    emit:
        // channel of good obsids with their metafits: tuple(obsid, metafits)
        obsMetafits = wsMetafits.out

        // channel of video name and frames to convert
        frame = wsSkyMap.out
            .map { _, png -> ["skymap", png] }
            .mix( wsPPDs.out.map { _, png -> ["ppd", png] } )
            .groupTuple()

        // channel of (obsid, metadata hashmap)
        obsMeta = wsSummary.map { obsid, summary ->
            def meta = [:]
            [
                "groupid", "pointing", "obs_name", "starttime_utc", "starttime_mjd",
                "bad_tiles", "eorband", "eorfield", "lst",
            ].each { key ->
                if (summary[key] != null) {
                    meta[key] = summary[key]
                }
            }
            [obsid, meta]
        }
}

// ensure preprocessed uvfits are downloaded
workflow prep {
    take:
        // channel of obsids with their metafits: tuple(obsid, metafits)
        obsMetafits
        obsMeta
    main:
        // download preprocessed uvfits
        obsMetafits
            .map { obsid, _ -> obsid }
            | asvoPrep

        if (params.noprepqa) {
            channel.empty() | prepVisQA
        } else {
            obsMetafits.join(asvoPrep.out) | prepVisQA
        }

        if (params.nossins) {
            channel.empty() | ssins
        } else {
            asvoPrep.out | ssins
        }

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
                bad_ants = stats.BAD_ANTS?:[]
                [
                    obsid,
                    stats.STATUS?:'',
                    stats.NANTS?:'',
                    stats.NTIMES?:'',
                    stats.NCHAN?:'',
                    stats.NPOLS?:'',
                    bad_ants.size(),
                    displayInts(bad_ants),
                ].join("\t")
            }
            .collectFile(
                name: "prepvis_metrics.tsv", newLine: true, sort: true,
                seed: [
                    "OBS",
                    "STATUS",
                    "NANTS",
                    "NTIMES",
                    "NCHAN",
                    "NPOLS",
                    "N_BAD_ANTS",
                    "BAD_ANTS",
                ].join("\t"),
                storeDir: "${results_dir}"
            )
            // display output path and number of lines
            | view { [it, it.readLines().size()] }

        [
            ["XX_RMS", { stats -> (stats.XX?:[:]).RMS?:[] }],
            ["YY_RMS", { stats -> (stats.YY?:[:]).RMS?:[] }],
            ["XX_MODZ0", { stats -> ((stats.XX?:[:]).MODZ_SCORE?:[:])["0"]?:[] }],
            ["XX_MODZ1", { stats -> ((stats.XX?:[:]).MODZ_SCORE?:[:])["1"]?:[] }],
            ["YY_MODZ0", { stats -> ((stats.YY?:[:]).MODZ_SCORE?:[:])["0"]?:[] }],
            ["YY_MODZ1", { stats -> ((stats.YY?:[:]).MODZ_SCORE?:[:])["1"]?:[] }],
        ].each { metric, getMetric ->
            prepVisQA.out
                // form row of tsv from json fields we care about
                .map { obsid, json ->
                    def stats = parseJson(json);
                    ([ obsid ] + getMetric(stats)).join("\t")
                }
                .collectFile(
                    name: "prepvis_${metric}.tsv", newLine: true, sort: true,
                    seed: [
                        "OBS",
                        metric,
                    ].join("\t"),
                    storeDir: "${results_dir}"
                )
                // display output path and number of lines
                | view { [it, it.readLines().size()] }
        }

        // plot prepvisQA
        prepVisQA.out | plotPrepVisQA

        // analyse aoflagger occupancy
        if (params.noprepqa) {
            channel.empty() | flagQA
        } else {
            obsMetafits.join(asvoPrep.out) | flagQA
        }
        if (params.noautoplot) {
            channel.empty() | autoplot
        } else {
            obsMetafits.join(asvoPrep.out).join(flagQA.out).flatMap { obsid, metafits, uvfits, flagJson ->
                    def flagStats = parseJson(flagJson)
                    def rxTiles = (flagStats.INPUTS?:[])
                        .findAll { it['Pol'] == "X" }
                        .collect { [it['Rx'], it['Antenna']] }
                        .groupBy { rx, tile -> rx }
                        .collect { rx, rx_tiles -> [rx, rx_tiles.collect { it[1] }] };
                    rxTiles.collect { rx, tiles ->
                        def suffix = sprintf("_rx%02d", rx);
                        def plotMeta = [
                            suffix: suffix,
                            "autoplot_args": "${params.autoplot_args} --sel_ants ${tiles.join(' ')} --plot_title '${obsid}${suffix}'"
                        ]
                        [obsid, plotMeta, metafits, uvfits]
                    }
                }
                | autoplot
        }

        // collect flagQA results
        // TODO: collect sky_chans from flagQA
        def sky_chans = (params.sky_chans).collect { ch -> "$ch".toString() }
        flagQA.out
            // form row of tsv from json fields we care about
            .map { obsid, json ->
                def stats = parseJson(json)
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

        // analyse ssins occupancy
        def ssinsOcc = ssins.out.map { obsid, _, __, occ_json ->
                [obsid, parseJson(occ_json)]
            }
        ssinsOcc.map { obsid, occ ->
                [
                    obsid,
                    occ['streak']?:'',
                    occ['total']?:'',
                ].join("\t")
            }
            .collectFile(
                name: "ssins_occupancy.tsv", newLine: true, sort: true,
                seed: ["OBS", "STREAK", "TOTAL"].join("\t"),
                storeDir: "${results_dir}"
            )
            | view { [it, it.readLines().size()] }

    emit:
        // channel of obsids which pass the flag gate: tuple(obsid, metafits, uvfits)
        obsMetaVis = (params.noprepqa ? \
            // no filter unless prepqa
            obsMetafits.join(asvoPrep.out) : \
            // filter
            flagQA.out
                .map { obsid, json -> [obsid, parseJson(json)] }
                .filter { obsid, flagStats ->
                    params.flag_occupancy_threshold == null || flagStats.total_occupancy < params.flag_occupancy_threshold
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
                    def prepStats = parseJson(prepJson);
                    def flagStats = parseJson(flagJson);
                    [obsid, prepStats, flagStats]
                }
                .filter { obsid, prepStats, _ ->
                    params.noprepqafilter || prepStats.STATUS == "GOOD"
                }
                .map { obsid, prepStats, flagStats ->
                    def flagAntennas = (flagStats.preflagged_ants?:[]) as Set
                    def prepAntennas = (prepStats.BAD_ANTS?:[]) as Set
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
                    def flagged_fchan_idxs = flagStats.flagged_fchan_idxs?:[]
                    [
                        obsid,
                        (flagAntennas as ArrayList).sort(false),
                        (prepAntennas as ArrayList).sort(false),
                        (manualAntennas as ArrayList).sort(false),
                        ((newAntennas - flagAntennas) as ArrayList).sort(false),
                        flagged_fchan_idxs,
                    ]
                })
        // channel of video name and frames to convert
        frame = plotPrepVisQA.out.flatMap { _, imgs ->
                imgs.collect { img ->
                    suffix = img.baseName.split('_')[-1]
                    ["prepvisqa_${suffix}", img]
                }
            }
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
            .mix(autoplot.out.map {_, meta, img -> ["prepvisqa_autoplot${meta.suffix?:''}", img]})
            .groupTuple()

        zip = prepVisQA.out.map { _, json -> ["prepvisqa", json]}
            .mix(flagQA.out.map { _, json -> ["flagqa", json]})
            .groupTuple()
}

workflow cal {
    take:
        // channel of metafits and preprocessed uvfits: tuple(obsid, meta, metafits, uvfits)
        obsMetaVis
    main:
        // hyperdrive di-calibrate on each obs
        obsMetaVis
            .map { def (obsids, meta, metafits, uvfits) = it
                [obsids, meta, metafits, uvfits, params.dical_args]
            }
            | hypCalSol

        // hyperdrive dical log analysis
        // hypCalSol.out
        //     .flatMap { obsid, meta, solns, logs ->
        //         [
        //             (solns instanceof List ? solns : [solns]),
        //             (logs instanceof List ? logs : [logs])
        //         ].transpose().collect { soln, log ->
        //             def name = soln.getBaseName().split('_')[3..-1].join('_')
        //             def (convergedDurationSec, convergedNumerator, convergedDenominator) = ['', '', '']
        //             if (!diCalLog.isEmpty()) {
        //                 def startMatch = diCalLog.getText() =~ /([\d: -]+) INFO  hyperdrive di-calibrate/
        //                 def startTime = null
        //                 if (startMatch) {
        //                     startTime = logDateFmt.parse(startMatch[0][1])
        //                 }
        //                 def convergedMatch = (diCalLog.getText() =~ /([\d: -]+) INFO  All timesteps: (\d+)\/(\d+)/)
        //                 if (startTime && convergedMatch) {
        //                     def convergedTime = logDateFmt.parse(convergedMatch[0][1])
        //                     convergedNumerator = convergedMatch[0][2] as int
        //                     convergedDenominator = convergedMatch[0][3] as int
        //                     convergedDurationSec = (convergedTime.time - startTime.time) / 1000
        //                 }
        //             }
        //             [obsid, name, convergedDurationSec, convergedNumerator, convergedDenominator].join("\t")
        //         }
        //     }
        //     .collectFile(
        //         name: "cal_timings${params.cal_suffix}.tsv", newLine: true, sort: true,
        //         seed: [
        //             "OBS", "CAL NAME", "CAL DUR", "CHS CONVERGED", "CHS TOTAL"
        //         ].join("\t"),
        //         storeDir: "${results_dir}${params.cal_suffix}"
        //     )
        //     | view { [it, it.readLines().size()] }

        // channel of individual dical solutions: tuple(obsid, meta, soln)
        // - hypCalSol gives multiple solutions, transpose gives 1 tuple per solution.
        eachCal = hypCalSol.out
            .flatMap { obsid, meta, solns, _ ->
                coerceList(solns).collect { soln ->
                    // give each calibration a name from basename of solution fits.
                    // this is everything after the obsid
                    def dical_name = soln.baseName.split('_')[3..-1].join('_');
                    def newMeta = [
                        dical_name: dical_name,
                        name: dical_name,
                        cal_prog: "hyp"
                    ]
                    [obsid, deepcopy(meta) + newMeta, soln]
                }
            }

        // generate json from solns
        eachCal | solJson

        // do polyFit unless --nopoly is set
        if (params.nopoly) {
            channel.empty() | polyFit
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
        obsMetafits = obsMetaVis.map { obsid, _, metafits, __ -> [obsid, metafits] }
        // calibration QA and plot solutions
        obsMetafits.cross(allCal)
            // cross with solutions from eachCal and eachPolyCal
            // TODO: eliminate the need for metafits in plotSols and calQA
            .map { obsMetafits_, allCal_ ->
                def (obsid, metafits) = obsMetafits_
                def (_, meta, soln) = allCal_
                [obsid, meta, metafits, soln]
            }
            | (plotSols & calQA)

        // plot each calQA result
        if (params.noplotcalqa) {
            channel.empty() | plotCalQA
        } else {
            calQA.out | plotCalQA
        }

        // collect calQA results as .tsv
        calQA.out
            // form row of tsv from json fields we care about
            .map { obsid, meta, json ->
                def stats = parseJson(json)
                def bad_ants = stats.BAD_ANTS?:[]
                // def convg_var = stats.CONVERGENCE_VAR
                [
                    obsid,
                    meta.name,
                    stats.STATUS?:'',
                    bad_ants.size,
                    displayInts(bad_ants),
                    (stats.PERCENT_UNUSED_BLS?:0) / 100,
                    (stats.PERCENT_NONCONVERGED_CHS?:0) / 100,
                    stats.RMS_CONVERGENCE?:'',
                    stats.SKEWNESS?:'',
                    stats.RECEIVER_VAR?:'',
                    stats.DFFT_POWER?:'',
                    stats.FAILURE_REASON,
                ].join("\t")
            }
            .collectFile(
                name: "cal_metrics${params.cal_suffix}.tsv", newLine: true, sort: true,
                seed: [
                    "OBS", "CAL NAME", "STATUS",
                    "N BAD ANTS", "BAD ANTS",
                    "BL UNUSED FRAC", "NON CONVG CHS FRAC",
                    "RMS CONVG",
                    "SKEW",
                    "RECV VAR",
                    "DFFT POW",
                    // "CONVG VAR",
                    // "CONVG VAR e14",
                    // "XX SKEW",
                    // "XX DFFT POWER",
                    // "YY SKEW",
                    // "YY DFFT POWER",
                    "FAILURE_REASON"
                ].join("\t"),
                storeDir: "${results_dir}${params.cal_suffix}"
            )
            // display output path and number of lines
            | view { [it, it.readLines().size()] }

        [
            ["XX_RMS", { stats -> (stats.XX?:[:]).RMS?:[] }],
            ["YY_RMS", { stats -> (stats.YY?:[:]).RMS?:[] }],
            ["XX_MODZ", { stats -> (stats.XX?:[:]).RMS_MODZ?:[] }],
            ["YY_MODZ", { stats -> (stats.YY?:[:]).RMS_MODZ?:[] }],
        ].each { metric, getMetric ->
            calQA.out
                // form row of tsv from json fields we care about
                .map { obsid, meta, json ->
                    def stats = parseJson(json);
                    ([ obsid, meta.name ] + getMetric(stats)).join("\t")
                }
                .collectFile(
                    name: "calqa_${metric}.tsv", newLine: true, sort: true,
                    seed: [
                        "OBS",
                        metric,
                    ].join("\t"),
                    storeDir: "${results_dir}"
                )
                // display output path and number of lines
                | view { [it, it.readLines().size()] }
        }

        // channel of obsids and names that pass qa. tuple(obsid, name)
        // - take tuple(obsid, cal_name, json) from calQA.out
        // - filter on json.STATUS == "PASS"
        // - take obsid and name
        obsMetaPass = calQA.out
            .map { obsid, meta, json ->
                def stats = parseJson(json);
                def newMeta = [:];
                if (stats.BAD_ANTS) {
                    def prepflags = (meta.prepflags?:[]) as Set
                    def calflags = ([]) as Set
                    if (!params.noCalFlags) {
                        calflags.addAll(stats.BAD_ANTS?:[])
                    }
                    def newflags = (calflags - prepflags) as ArrayList
                    if (newflags) {
                        newMeta.calflags = deepcopy(newflags.sort(false))
                    }
                }
                [obsid, deepcopy(meta) + newMeta, stats]
            }
            // TODO: reintroduce status filter
            .filter { _, __, stats -> stats.STATUS == null || stats.STATUS == "PASS" }
            .map { obsid, meta, stats -> [obsid, meta] }

    emit:
        // channel of calibration solutions that pass qa. tuple(obsid, name, cal)
        // - take tuple(obsid, meta, soln) from allCal
        // - match with obsMetaPass on (obsid, cal_name)
        obsMetaCalPass = allCal
            .cross(obsMetaPass) {it -> def (obsid, meta) = it; [obsid, meta.name]}
            .map{ allCal_, obsMetaPass_ ->
                def (obsid, _, soln) = allCal_
                def (__, meta) = obsMetaPass_
                [obsid, meta, soln]
            }
        // channel of files to archive, and their buckets
        archive = calQA.out.map { _, __, json -> ["calqa", json]}
        // channel of video name and frames to convert
        frame = plotCalQA.out.mix(plotSols.out)
            .flatMap { _, meta, pngs ->
                coerceList(pngs).collect { png ->
                    def suffix = png.baseName.split('_')[-1]
                    ["calqa_${meta.name}_${suffix}", png]
                }
            }
            .groupTuple()
        // channel of files to zip
        zip = calQA.out.map { _, meta, json -> ["calqa_${meta.name}", json] }
            .mix(solJson.out.map { obsid, meta, json -> ["soljson_${meta.name}", json] })
            .mix(allCal.map { obsid, meta, soln -> ["${meta.cal_prog}_soln_${meta.name}", soln] })
            .groupTuple()
}

// process uvfits visibilities
workflow uvfits {
    take:
        obsMetaUV
    main:
        // vis QA on each phase2 obsid, no subtractions
        obsMetaUV
            .filter { obsid, meta, __ -> meta.sub == null && meta.config =~ /phase2a.*/ }
            | visQA

        if (params.nouvplot) {
            channel.empty() | uvPlot
        } else {
            obsMetaUV | uvPlot
        }

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

        // write bands and fields to a file
        obsMetaUV.map { obsid, meta, vis ->
            [
                obsid, meta.name, meta.config?:'', meta.eorband?:'', meta.eorfield?:'',
                meta.nchans?:'', meta.lowfreq?:'', meta.freq_res?:'', vis
            ].join('\t') }
            .collectFile(
                name: "eor_params.tsv", newLine: true, sort: true,
                seed: [
                    "OBSID", "NAME", "CONFIG", "BAND", "FIELD", "NCHAN", "LOWFREQ", "FREQRES", "VIS"
                ].join('\t'),
                storeDir: "${results_dir}${params.cal_suffix}"
            )
            // display output path and number of lines
            | view { [it, it.readLines().size()] }

        if (params.noplotvisqa) {
            channel.empty() | plotVisQA
        } else {
            visQA.out | plotVisQA
        }

        if (params.noeor) {
            obsMetaUVEoR = obsMetaUV
        } else {
            obsMetaUVEoR = obsMetaUV
                .filter { _, meta, __ -> (meta.eorband != null && meta.eorfield != null) }
        }

        // ps_metrics
        if (params.nopsmetrics) {
            obsMetaUVPass_ = obsMetaUVEoR
        } else {
            obsMetaUVEoR | psMetrics

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

            psMeta = psMetrics.out
                // add ps_metrics values to vis meta
                .map { obsid, meta, dat ->
                    def (p_wedge, num_cells, p_window) = dat.getText().split('\n')[0].split(' ')[1..-1]
                    def newMeta = [
                        p_window: parseFloatOrNaN(p_window),
                        p_wedge: parseFloatOrNaN(p_wedge),
                        num_cells: num_cells
                    ]
                    [obsid, deepcopy(meta) + newMeta]
                }

            passReasons = psMeta
                // get cmt reduced metrics failures
                .groupTuple(by: 0)
                .flatMap { obsid, metas ->
                    def nosubMeta = metas.find { it.sub == null} ?: [:]
                    def nosubReason = cmt_ps_metrics_pass(nosubMeta)
                    def subMetas = metas.findAll { it.sub != null } ?: []
                    return [
                        [obsid, nosubMeta, nosubReason]
                    ] + subMetas.collect { subMeta ->
                        [obsid, subMeta, cmt_ps_metrics_pass_sub(nosubMeta, subMeta)?:nosubReason]
                    }
                }

            passReasons.map { obsid, meta, reason ->
                    [obsid, meta.name, meta.p_window, meta.p_wedge, meta.num_cells, reason].join("\t")
                }
                .collectFile(
                    name: "ps_reasons.tsv", newLine: true, sort: true,
                    seed: [
                        "OBSID", "NAME", "P_WINDOW", "P_WEDGE", "NUM_CELLS", "REASON"
                    ].join("\t"),
                    storeDir: "${results_dir}${params.cal_suffix}"
                )
                // display output path and number of lines
                | view { [it, it.readLines().size()] }

            obsMetaPass = passReasons
                .filter { obsid, meta, reason -> reason == null }
                .map { obsid, meta, _ -> [obsid, meta] }

            obsMetaUVPass_ = obsMetaUV.cross(obsMetaPass) { def (obsid, meta) = it; [obsid, meta.name] }
                .map { obsMetaUV_ , obsMetaPass_ ->
                    def (obsid, _, vis) = obsMetaUV_
                    def (__, meta) = obsMetaPass_
                    [obsid, meta, vis]
                }
        }

        // delay spectrum
        if (params.nodelayspec) {
            channel.empty() | delaySpec
        } else {
            obsMetaUVEoR | delaySpec
        }

    emit:
        // channel of files to archive, and their buckets
        archive = visQA.out.map { _, __, json -> ["visqa", json]}
        // channel of video name and frames to convert
        frame = plotVisQA.out
            .map { _, meta, png -> ["visqa_${meta.name}_rms", png] }
            .mix(uvPlot.out.flatMap {_, name, pngs -> pngs.collect { png ->
                def pol = png.baseName.split('_')[-2..-1].join('_')
                ["visqa_${name}_${pol}", png]
            }})
            .mix(delaySpec.out.map { _, meta, png -> ["visqa_dlyspec_${meta.name}", png] })
            .groupTuple()
        // channel of files to zip
        zip = visQA.out.map { _, meta, json -> ["visqa_${meta.name}", json] }
            .groupTuple()
        obsMetaUVPass = obsMetaUVPass_
}

def wscleanParams = [
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
def wscleanDConvParams = deepcopy(wscleanParams) + [
    args: "${params.wsclean_args} ${params.wsclean_dconv_args}",
    // args: params.wsclean_dconv_args
    niter: params.img_niter,
    minor_clean_gain: params.img_minor_clean_gain,
    major_clean_gain: params.img_major_clean_gain,
    auto_threshold: params.img_auto_threshold,
    auto_mask: params.img_auto_mask,
    mwa_path: params.img_mwa_path,
]

// image measurement sets and QA images
workflow img {
    take:
        obsMetaVis
    main:

        // wsclean: make deconvolved images
        if (params.img_split_intervals) {
            // add -t???? suffix to name, to split wsclean over multiple jobs
            splitObsMetaVis = obsMetaVis.flatMap { obsid, meta, vis ->
                (0..meta.ntimes).collect { i ->
                    def newMeta = deepcopy(meta)
                    newMeta.interval = [i, i+1]
                    // interval token
                    newMeta.inter_tok = sprintf("-t%04d", i)
                    [obsid, newMeta, vis]
                }
            }
        } else {
            splitObsMetaVis = obsMetaVis
        }

        if (params.nodeconv) {
            splitObsMetaVis.map {obsid, meta, vis ->
                    def imgParams = deepcopy(wscleanParams)
                    imgParams.img_name = "${meta.cal_prog}_${obsid}${meta.subobs?:''}_${meta.name}${meta.inter_tok?:''}"
                    [obsid, meta, vis, imgParams]
                }
                | wscleanDirty
            channel.empty() | wscleanDConv
        } else {
            channel.empty() | wscleanDirty
            // splitObsMetaVis.cross(wscleanDirty.out) {
            //         def (obsid, meta) = it;
            //         [obsid, meta.subobs, meta.name, meta.inter_tok]
            //     }
            //     .map { splitObsMetaVis_, wscleanDirty_ ->
            //         def (obsid, _, vis, imgParams) = splitObsMetaVis_
            //         def (__, meta, imgs) = wscleanDirty_
            //         def dirtyImgs = imgs.findAll {img -> img.baseName.split('-')[-1] == 'dirty'}
            //         def newImgParams = deepcopy(imgParams + wscleanDConvParams)
            //         [obsid, meta, vis, newImgParams, dirtyImgs]
            //     }

            splitObsMetaVis.map { obsid, meta, vis ->
                    def imgParams = deepcopy(wscleanDConvParams)
                    imgParams.img_name = "${meta.cal_prog}_${obsid}${meta.subobs?:''}_${meta.name}${meta.inter_tok?:''}"
                    [obsid, meta, vis, imgParams, []]
                }
                .flatMap {obsid, meta, vis, imgParams_, dirtyImgs ->
                    [
                        ["-{XX,YY}", "xx,yy -join-polarizations"],
                        ["-{I,V}", "i,v -join-polarizations"],
                        // ["-{XX,YY,XY,XYi}", "xx,xy,yx,yy -link-polarizations xx,yy"],
                        // ["-{Q,U}", "q,u -join-polarizations -squared-channel-joining"],
                    ].collect { polGlob, polArg ->
                        def imgParams = deepcopy(imgParams_)
                        imgParams.pol = polArg
                        imgParams.glob = "wsclean_${imgParams.img_name}"
                        if (multiinterval && !meta.inter_tok) {
                            imgParams.glob += "-t????"
                        }
                        if (multichannel) {
                            imgParams.glob += "-MFS"
                        }
                        imgParams.glob += polGlob
                        imgParams.glob += "-image.fits"
                        [obsid, meta, vis, imgParams, dirtyImgs]
                    }
                }
                | wscleanDConv
        }

        obsMetaImg = wscleanDirty.out.mix(wscleanDConv.out)
            .flatMap { obsid, meta, imgs ->
                coerceList(imgs).collect { img ->
                    [obsid, deepcopy(meta) + decomposeImg(img), img] }}
            .branch { obsid, meta, img ->
                // target product is image unless deconv is disabled, separate MFS
                imgMfs: meta.prod == (params.nodeconv ? "dirty" : "image") && meta.chan == -1
                imgNoMfs: meta.prod == (params.nodeconv ? "dirty" : "image")
                // target product is uv grids
                gridMfs: meta.prod ==~ /uv-.*/ && meta.chan == -1
                gridNoMfs: meta.prod ==~ /uv-.*/
                // other potential products: psf, residual, model
            }

        if (!params.nopolcomp || params.thumbnail_limits) {
            // calculate quantiles (what values are at nth percentile)
            obsMetaImg.imgMfs \
                // .mix(obsMetaGrid) \
                | imgQuantiles

            // limits are used to set the color scale of each type of image.
            imgLimits = imgQuantiles.out
                .map { obsid, meta, hist ->
                    high = parseCsv(hist, true, 2)
                        .find { row -> Float.compare(parseFloatOrNaN(row.quantile), params.thumbnail_quantile as Float) == 0 }
                    high = parseFloatOrNaN(high.value)
                    low = parseCsv(hist, true, 2)
                        .find { row -> Float.compare(parseFloatOrNaN(row.quantile), (1-params.thumbnail_quantile) as Float) == 0 }
                    low = parseFloatOrNaN(low.value)
                    // subobs = meta.subobs?:''
                    [obsid, meta.inter_tok, meta.name, meta.inter_suffix, high, low]
                }
                .filter { _obsid, _interval, _name, _suff, high, low -> !(isNaN(high) || isNaN(low)) }
                // get max limit for each obs, interval, name, suffix (deals with multiple channels)
                .groupTuple(by: 0..3)
                .map { obsid, interval, name, suff, highs, lows ->
                    high = coerceList(highs).max()
                    low = coerceList(lows).max()
                    [obsid, interval, name, suff, high, low]
                }
                // for each obs, interval, name, produce a mapping from suff to limit
                .groupTuple(by: 0..2)
                .map { obsid, interval, name, suffs, highs, lows ->
                    _suffLimits = [
                            coerceList(suffs),
                            coerceList(highs),
                            coerceList(lows),
                        ]
                        .transpose()
                        .collect { suff, high, low -> [suff, [high, low]]}
                        .collectEntries()
                    [obsid, interval, name, _suffLimits]
                }
                // .view { "imgLimits ${it}" }

            obsMetaImgMfsPass = obsMetaImg.imgMfs.map { obsid, meta, img ->
                    [obsid, meta.inter_tok, meta.name, meta.inter_suffix, meta, img]
                }
                .cross(imgLimits.flatMap { obsid, interval, name, _suffLimits ->
                    _suffLimits.collect { k, _ -> [ obsid, interval, name, k ]}}
                ) { it[0..3] }
                .map { obsMetaImgMfs_, imgLimits_ ->
                    def (obsid, _, __, ___, meta, img) = obsMetaImgMfs_
                    [obsid, meta, img]
                }
                // .view { "obsMetaImgMfs ${it}" }

            // channel of all suffixes in
            suffs_ = obsMetaImg.imgMfs
                .map { obsid, meta, img -> meta.inter_suffix }
                .unique()
                .toSortedList()
                .map{it -> [it]}

            // write img limits to file
            imgLimits.combine(suffs_)
                .map { obsid, interval, name, limits, suffs ->
                        ([obsid, interval, name] + suffs.collect { suff ->
                                hilo = ['', '']
                                if (limits[suff]) {
                                    hilo = limits[suff]
                                }
                                [suff] + hilo
                            }
                            .flatten()
                        ).join("\t")
                }
                .collectFile(
                    name: "img_limits.tsv", newLine: true, sort: true,
                    seed: (["OBSID", "interval", "IMG NAME", "suff", "hi", "lo"]).join("\t"),
                    storeDir: "${results_dir}${params.img_suffix}${params.cal_suffix}"
                )
                // display output path and number of lines
                | view { [it, it.readLines().size()] }

            // value channel containing a map from img suffix to max limit excluding outliers
            // suffLimits = channel.from([:])
            suffLimits = imgLimits.combine(suffs_)
                .flatMap{ obsid, interval, name, limits, suffs ->
                    suffs.collect { suff ->
                        hilo = [Float.NaN, Float.NaN]
                        if (limits[suff]) {
                            hilo = limits[suff]
                        }
                        [new Tuple(name, suff)] + hilo
                    }
                }
                .groupTuple(by: 0)
                .map { group, highs_, lows_ ->
                    highs = coerceList(highs_)
                        .findAll { !isNaN(it) }
                        .sort(false)
                    high_index = ((highs.size() - 1) * params.thumbnail_quantile).round() as Integer

                    lows = coerceList(lows_)
                        .findAll { !isNaN(it) }
                        .sort(false)
                    low_index = (lows.size() * (1-params.thumbnail_quantile)).round() as Integer
                    [group, [highs[high_index], lows[low_index]]]
                }
                .toList()
                .map { it.collectEntries() }
        } else {
            suffLimits = channel.from([:])
            obsMetaImgMfsPass = obsMetaImg.imgMfs
        }

        suffLimits.view { "suffLimits \n${it.collect().join('\n')}"}

        // do a polcomp for each allowed prod
        if (params.nopolcomp) {
            channel.empty() | polComp
        } else {

            // all valid pol orders for polComp
            polOrders = [
                // ["I", "V", "Q"],
                ["XX", "V", "YY"],
                // ["XX", "XY", "YY"]
            ]

            obsMetaImgMfsPass
                // look up the sufflimit for each suffix
                .combine(suffLimits)
                .map { obsid, meta, img, _suffLimits ->
                    def name_suff = new Tuple(meta.name, meta.inter_suffix)
                    def hilo = [Float.NaN, Float.NaN]
                    if (_suffLimits && _suffLimits[name_suff]) {
                        hilo = _suffLimits[name_suff]
                    }
                    [obsid, meta, img] + hilo
                }
                // filter out NaN limits
                .filter { obsid, meta, img, high, low ->
                    !(isNaN(high) || isNaN(low))
                }
                // group (pol, img, limit) by (obs, interval, name, prod).
                .map { obsid, meta, img, high, low ->
                    [obsid, meta.inter_tok, meta.name, meta.prod, meta.pol, img, high, low]
                }
                .groupTuple(by: 0..3)
                // make a hashmap of pol -> (img, limit) for each (obs, interval, name, prod)
                .map { obsid, interval, name, prod, pols, imgs, highs, lows ->
                    polImgLimits = [
                        coerceList(pols),
                        coerceList(imgs),
                        coerceList(highs),
                        coerceList(lows),
                    ].transpose()
                        .collect { pol, img, high, low -> [pol, [img, high, low]] }
                        .collectEntries()
                    [obsid, interval, name, prod, polImgLimits]
                }
                // for every pol order where every pol is present, do a polComp
                .flatMap { obsid, interval, name, prod, polImgLimits ->
                    polOrders.findAll { order -> polImgLimits.keySet().containsAll(order) }
                        .collect { order ->
                            def (imgs, highs, _)  = order.collect { pol -> polImgLimits[pol] }.transpose()
                            polcompMeta = [
                                name:name,
                                prod:prod,
                                order:order,
                                limits:highs,
                                subobs:interval,
                            ]
                            [obsid, polcompMeta, imgs]
                        }
                }
                | polComp
        }

        // make thumbnails
        if (params.nothumbnail) {
            channel.empty() | thumbnail
        } else {
            obsMetaImgMfsPass
                .mix(params.thumbnail_uvs ? obsMetaGrid : channel.empty())
                // look up the sufflimit for each suffix
                .combine(params.thumbnail_limits ? suffLimits : channel.from([:]))
                // .view { "thumbnail limit ${it}"}
                .map { obsid, meta, img, _suffLimits ->
                    def name_suff = new Tuple(meta.name, meta.inter_suffix)
                    def hilo = [Float.NaN, Float.NaN]
                    if (_suffLimits && _suffLimits[name_suff]) {
                        hilo = _suffLimits[name_suff]
                    }
                    [obsid, meta, img] + hilo
                }
                // filter out NaN limits
                // .filter { obsid, meta, img, high, low -> !isNaN(limit) }
                .map { obsid, meta, img, high, low  ->
                    def newMeta = [:]
                    if (!isNaN(high)) {
                        newMeta.vmax = high
                    } else {
                        newMeta.vmax_quantile = params.thumbnail_quantile
                    }
                    if (!isNaN(low)) {
                        newMeta.vmin = low
                    } else {
                        newMeta.vmin_quantile = (1-params.thumbnail_quantile)
                    }
                    if (['Q', 'V', 'U', 'XY', 'YX'].contains(meta.pol)) {
                        newMeta.symmetric = true
                        newMeta.norm_args = '{\\\\"stretch\\\\":\\\\"asinh\\\\",\\\\"asinh_a\\\\":0.8}'
                    } else {
                        if (newMeta.vmin as Float < 0 && newMeta.vmax as Float > 0) {
                            ratio = (-newMeta.vmin) / (newMeta.vmax - newMeta.vmin)
                        } else {
                            ratio = 0.1
                        }
                        // newMeta.norm_args = '{\\\\"stretch\\\\":\\\\"asinh\\\\",\\\\"asinh_a\\\\":'+ "${ratio}" + '}'
                        newMeta.norm_args = '{\\\\"stretch\\\\":\\\\"power\\\\",\\\\"power\\\\":2,\\\\"clip\\\\":true}'
                    }
                    [obsid, deepcopy(meta) + newMeta, img]
                }
                .view { "thumbnail ${it}" }
                | thumbnail
        }

        // montage of polarizations
        if (!params.nomontage) {
            thumbnail.out.flatMap { obs, meta, png ->
                // visibility name, without sub*
                def newMeta = deepcopy(meta) + [vis_name: meta.name, sub: meta.sub?:'nosub']
                if (m = newMeta.vis_name =~ /(sub|ionosub)_(.*)/) {
                    newMeta.sub = m.group(1).toString()
                    newMeta.vis_name = m.group(2).toString()
                }

                // suffix without interval (-t????)
                // newMeta.subobs = newMeta.subobs?:''
                // newMeta.inter_suffix = newMeta.suffix
                // if (m = newMeta.inter_suffix =~ /(t????)-(.*)/) {
                //     newMeta.subobs = m.group(1).toString() + newMeta.subobs
                //     newMeta.inter_suffix = m.group(2).toString()
                // }
                // def montage_name = [sub: suffix, pol: sub_name].get(params.montage_by)
                subobs = "${newMeta.subobs?:''}${newMeta.inter_tok?:''}".toString()
                montages = []
                if (params.montage_by_sub) {
                    montages << [obsid, subobs, "${newMeta.vis_name}-${newMeta.chan_tok}-${newMeta.pol}-${newMeta.prod}".toString(), png]
                }
                if (params.montage_by_pol) {
                    montages << [obsid, subobs, "${newMeta.name}-${newMeta.chan_tok}-${newMeta.prod}".toString(), png]
                }
                montages
            }
            .groupTuple(by: 0..2)
            .map { obs, subobs, name, pngs ->
                [obs, [subobs:subobs, name:name], pngs.sort(false)]
            }
            // .view {"polmontage ${it}"}
            | polMontage
        } else {
            channel.empty() | polMontage
        }


        // each vis name can have multiple images in obsMetaImgMfs
        // group by obsid, vis name using original meta from obsMetaVis
        obsMetaImgGroup = obsMetaVis.cross(obsMetaImgMfsPass) { def (obsid, meta) = it; [obsid, meta.name] }
            .map { obsMetaVis_, obsMetaImgMfs_ ->
                def (obsid, meta) = obsMetaVis_
                def (_, imgMeta, img) = obsMetaImgMfs_
                [obsid, meta.name, meta, imgMeta, img]
            }
            .filter { obsid, name, meta, imgMeta, img -> ['XX', 'YY', 'V'].contains(imgMeta.pol)}
            .groupTuple(by: 0..1)
            .map { obsid, _, metas, imgMetas, imgs -> [obsid, metas[0], imgMetas, imgs] }
            // | view { def (obsid, meta) = it; "obsMetaImgGroup ${obsid}, ${meta.name}"}

        // imgQA for MFS images and groups of images
        //  - imgs need to be deconvolved for QA
        //  - can't handle multiple intervals
        if (params.noimgqa || params.nodeconv || params.img_split_intervals) {
            channel.empty() | imgQA
            obsMetaImgPass_ = obsMetaImgGroup
        } else {
            obsMetaImgGroup
                .map { obsid, meta, imgMetas, imgs -> [obsid, meta, imgs] }
                .filter { _, meta, __ ->
                    // imgQA looks at pks flux, which needs to be deconvolved, only works with eor0
                    // meta.prod == 'image' &&
                    meta.eorfield == 0
                }
                // .map { obsid, meta, img -> [obsid, meta.name, meta, img] }
                // .groupTuple(by: 0..1)
                // .map {obsid, name, metas, imgs -> [obsid, metas[0], imgs]}
                | imgQA

            imgQA.out.map { _, meta, json -> [meta.name, json] }
                .groupTuple(by: 0)
                | plotImgQA

            // collect imgQA results as .tsv
            imgQA.out
                // form row of tsv from fields we care about
                .map { obsid, meta, json ->
                    def stats = parseJson(json)
                    def xx = stats.XX?:[:]
                    def xx_pks = xx.PKS0023_026?:[:]
                    def yy = stats.YY?:[:]
                    def yy_pks = yy.PKS0023_026?:[:]
                    def v = stats.V?:[:]
                    def v_pks = v.PKS0023_026?:[:]
                    [
                        obsid,
                        meta.name,
                        xx.RMS_ALL, xx.RMS_BOX, xx_pks.PEAK_FLUX, xx_pks.INT_FLUX,
                        yy.RMS_ALL, yy.RMS_BOX, yy_pks.PEAK_FLUX, yy_pks.INT_FLUX,
                        v.RMS_ALL, v.RMS_BOX, v_pks.PEAK_FLUX, v_pks.INT_FLUX,
                        // stats.V_XX?.RMS_RATIO, stats.V_XX?.RMS_RATIO_BOX
                    ].join("\t")
                }
                .collectFile(
                    name: "img_metrics.tsv", newLine: true, sort: true,
                    seed: [
                        "OBS", "IMG NAME",
                        "XX RMS ALL", "XX RMS BOX", "XX PKS0023_026 PEAK", "XX PKS0023_026 INT",
                        "YY RMS ALL", "YY RMS BOX", "YY PKS0023_026 PEAK", "YY PKS0023_026 INT",
                        "V RMS ALL","V RMS BOX", "V PKS0023_026 PEAK", "V PKS0023_026 INT" ,
                        // "V:XX RMS RATIO", "V:XX RMS RATIO BOX"
                    ].join("\t"),
                    storeDir: "${results_dir}${params.img_suffix}${params.cal_suffix}"
                )
                // display output path and number of lines
                | view { [it, it.readLines().size()] }

            obsMetaReasons = imgQA.out
                .map { obsid, meta, json ->
                    def stats = parseJson(json)
                    def xx = stats.XX?:[:]
                    def xx_pks = xx.PKS0023_026?:[:]
                    def yy = stats.YY?:[:]
                    def yy_pks = yy.PKS0023_026?:[:]
                    def v = stats.V?:[:]
                    def v_pks = v.PKS0023_026?:[:]
                    def newMeta = [
                        xx_pks_int: xx_pks.INT_FLUX,
                        yy_pks_int: yy_pks.INT_FLUX,
                        v_pks_int: v_pks.INT_FLUX,
                        v_rms_box: v.RMS_BOX,
                    ]
                    [obsid, deepcopy(meta) + newMeta]
                }
                .groupTuple(by: 0)
                .map { obsid, metas -> [obsid, coerceList(metas)] }
                .flatMap { obsid, metas ->
                    def nosubMeta = metas.find { meta -> meta.sub == null }
                    def nosubReason = cmt_imgqa_pass(nosubMeta)
                    def subMetas = metas.findAll { it.sub != null } ?: []
                    def subReasons = subMetas.collect { subMeta -> cmt_imgqa_pass_sub(nosubMeta, subMeta) }
                    def asubReason = subReasons.find { it != null }
                    if (nosubReason == null && asubReason) {
                        nosubReason = asubReason
                    }
                    return [
                        [obsid, nosubMeta, nosubReason]
                    ] + [subMetas, subReasons].transpose().collect { subMeta, subReason ->
                        [obsid, subMeta, subReason?:nosubReason]
                    }
                }

            obsMetaReasons.map { obsid, meta, reason ->
                    [obsid, meta.name, meta.xx_pks_int, meta.yy_pks_int, meta.v_pks_int, meta.v_rms_box, reason].join("\t")
                }
                .collectFile(
                    name: "imgqa_reasons.tsv", newLine: true, sort: true,
                    seed: [
                        "OBSID", "NAME", "XX PKS INT", "YY PKS INT", "V PKS INT", "V RMS BOX", "REASON"
                    ].join("\t"),
                    storeDir: "${results_dir}${params.img_suffix}${params.cal_suffix}"
                )
                // display output path and number of lines
                | view { [it, it.readLines().size()] }

            obsMetaPass = obsMetaReasons
                .filter { obsid, meta, reason -> reason == null }
                .map { obsid, meta, _ -> [obsid, meta] }

            obsMetaImgPass_ = obsMetaImgGroup.cross(obsMetaPass) { def (obsid, meta) = it; [obsid, meta.name] }
                .map { obsMetaImgGroup_, obsMetaPass_ ->
                    def (obsid, _, imgMetas, imgs) = obsMetaImgGroup_
                    def (__, meta) = obsMetaPass_
                    [obsid, meta, imgMetas, imgs]
                }
        }

    emit:
        // channel of files to archive, and their buckets
        archive = obsMetaImgMfsPass.filter { obsid, _, __ -> !obsid.startsWith("e") }
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
        // channel of files to zip
        zip = imgQA.out.map { _, meta, json -> ["imgqa_${meta.name}", json] }
            .groupTuple()
        // channel of tuple (obsid, imgMeta, img) that pass qa
        obsMetaImgPass = obsMetaImgPass_
}

workflow imgCombine {
    take:
        // tuple of (obsid, meta, uvfits)
        obsMetaUV
        // tuple of (chunk, chunkMeta, obsids) for grouping
        chunkMetaObs

    main:

        // combined images
        groupMetaVisImgParams = obsMetaUV.cross(
                chunkMetaObs.flatMap { chunk, chunkMeta, obsids ->
                    obsids.collect { obsid -> [obsid, deepcopy(chunkMeta), chunk] }
                }
            ) { def (obsid, meta) = it; [obsid, meta.name] }
            .map { obsMetaUV_, chunkMetaObs_ ->
                def (obsid, obsMeta, vis) = obsMetaUV_
                def (__, chunkMeta, chunk) = chunkMetaObs_
                [chunk, chunkMeta.name, chunkMeta, vis]
            }
            .groupTuple(by: 0..1)
            // .view { it -> "\n -> imgCombine obsMetaUV x chunkMetaObs ${it}\n" }
            .map { chunk, _, chunkMetas, viss ->
                [chunk, deepcopy(chunkMetas[0]), viss.flatten(), deepcopy(wscleanParams)]
            }
            // .view { it -> "\n -> imgCombine before filter ${it}\n" }
            .filter { chunk, chunkMeta, viss, imgParams -> chunkMeta.name }

        if (params.nodeconv) {
            groupMetaVisImgParams | wscleanDirty
            channel.empty() | wscleanDConv
        } else {
            channel.empty() | wscleanDirty
            groupMetaVisImgParams.map { chunk, meta, vis, imgParams ->
                    [chunk, meta, vis, deepcopy(wscleanDConvParams) + imgParams, []]
                }
                | wscleanDConv
        }

        // make thumbnails
        // - only MFS images unless thumbnail_all_chans
        // - exclude dirty if deconv is enabled
        wscleanDConv.out.mix(wscleanDirty.out)
            .flatMap { obsid, meta, imgs ->
                imgs.collect { img ->
                    [obsid, deepcopy(meta) + decomposeImg(img), img] }}
            .filter { _, meta, img ->
                meta.prod !=~ /uv-.*/ \
                && (meta.chan?:-1 == -1 || params.thumbnail_all_chans) \
                && (meta.prod == (params.nodeconv ? "dirty" : "image"))
            }
            | thumbnail
}

// make combined grids and images using chips and wsclean
workflow chips {
    take:
        // tuple of (obsid, meta, uvfits)
        obsMetaUV
        // tuple of (chunk, chunkMeta, obsids) for grouping
        chunkMetaObs

    main:
        // power spectrum
        if (params.nopowerspec) {
            channel.empty() | chipsGrid
        } else {
            obsMetaUV
                .map { obsid, meta, viss ->
                    [obsid, deepcopy(meta) + [nobs: 1, ext: "${obsid}_${meta.name}"], [obsid], viss]
                }
                | chipsGrid
        }

        chipsGrid.out.cross(
                chunkMetaObs.flatMap { chunk, chunkMeta, obsids ->
                    obsids.collect { obsid -> [obsid, deepcopy(chunkMeta), chunk] }
                }
            ) { def (obsid, meta) = it; [obsid, meta.name] }
            .map { chipsGrid_, chunkMetaObs_ ->
                def (obsid, obsMeta, grid) = chipsGrid_
                def (__, chunkMeta, chunk) = chunkMetaObs_
                [chunk, chunkMeta.name, chunkMeta, obsMeta.ext, grid]
            }
            .groupTuple(by: 0..1)
            // .view { it -> "\n -> chipsGrid x chunkMetaObs ${it}\n" }
            .map { chunk, _, chunkMetas, exts, grids ->
                def chunkMeta = chunkMetas[0]
                def newMeta = [ ext: "${chunk}_${chunkMeta.name}"]
                // print("\n -> chunkMetas 0: ${chunkMeta}\n")
                [chunk, deepcopy(chunkMeta) + newMeta, exts, grids.flatten()]
            }
            // .view { it -> "\n -> before chipsCombine ${it}\n" }
            .filter { chunk, chunkMeta, exts, grid -> chunkMeta.name }
            // .take(10)
            | chipsCombine

        // chipsGrid.out.mix(chipsCombine.map { group, meta, _, grid -> group, meta,grid }) | chipsLssa
        chipsCombine.out.map { group, meta, _, grid -> [group, meta, grid] } | chipsLssa

        def singles1D = chipsLssa.out.map { chunk, meta, grid ->
                def newMeta = [
                    ptype: '1D',
                    pol: 'both',
                    title: "crosspower\\n${chunk}\\n${meta.name}",
                    plot_name: 'chips1d',
                    max_power: 1e15,
                    min_power: 1e3,
                    tags: [],
                ]
                [chunk, deepcopy(meta) + newMeta, grid]
            }

        def singles2D = chipsLssa.out.map { chunk, meta, grid ->
                def newMeta = [
                    ptype: '2D',
                    pol: 'both',
                    title: "crosspower\\n${chunk}\\n${meta.name}",
                    plot_name: 'chips2d',
                    max_power: 1e15,
                    min_power: 1e3,
                    tags: [],
                ]
                [chunk, deepcopy(meta) + newMeta, grid]
            }

        def comps1D = chipsLssa.out
            // group grids from vis name together
            .groupTuple(by: 0)
            // .view { chunk, metas, _ -> "\n -> comps1D metas: \n${metas.join('\n')}"}
            .filter { _, metas, __ -> metas.size() >= 3 }
            .map { chunk, metas, grids ->
                def nosubMeta = metas.find { it.sub == null} ?: [:]
                def subMeta = metas.find { it.sub == 'sub'} ?: [:]
                def ionosubMeta = metas.find { it.sub == 'ionosub'} ?: [:]
                def newMeta = [
                    ptype: '1D_comp',
                    pol: 'both',
                    title: "crosspower\\n${chunk}",
                    plot_name: 'chips1d_comp',
                    name: nosubMeta.name,
                    ext: nosubMeta.ext,
                    tags: [
                        nosubMeta.ext,
                        subMeta.ext,
                        ionosubMeta.ext,
                    ],
                    labels: [
                        nosubMeta.name,
                        subMeta.name,
                        ionosubMeta.name,
                    ],
                    max_power: 1e15,
                    min_power: 1e3,
                ]
                [chunk, deepcopy(nosubMeta) + newMeta, grids.flatten()]
            }
        def diffs2D_sub = chipsLssa.out
            // group grids from vis name together
            .groupTuple(by: 0)
            .filter { _, metas, __ -> metas.size() >= 2 }
            .map { chunk, metas, grids ->
                def nosubMeta = metas.find { it.sub == null} ?: [:]
                def subMeta = metas.find { it.sub == 'sub'} ?: [:]
                def newMeta = [
                    ptype: '2D_diff',
                    pol: 'both',
                    title: "crosspower diff (nosub-sub)\\n${chunk}",
                    plot_name: 'chips2d_diff_nosub_sub',
                    name: subMeta.name,
                    ext: subMeta.ext,
                    tags: [
                        nosubMeta.ext,
                        subMeta.ext,
                    ],
                    labels: [
                        nosubMeta.name,
                        subMeta.name,
                    ],
                ]
                [chunk, deepcopy(nosubMeta) + newMeta, grids.flatten()]
            }
        def diffs2D_iono = chipsLssa.out
            // group grids from vis name together
            .groupTuple(by: 0)
            .filter { _, metas, __ -> metas.size() >= 2 }
            .map { chunk, metas, grids ->
                def nosubMeta = metas.find { it.sub == null} ?: [:]
                def ionosubMeta = metas.find { it.sub == 'ionosub'} ?: [:]
                def newMeta = [
                    ptype: '2D_diff',
                    pol: 'both',
                    title: "crosspower diff (nosub-ionosub)\\n${chunk}",
                    plot_name: 'chips2d_diff_nosub_ionosub',
                    name: ionosubMeta.name,
                    ext: ionosubMeta.ext,
                    tags: [
                        nosubMeta.ext,
                        ionosubMeta.ext,
                    ],
                    labels: [
                        nosubMeta.name,
                        ionosubMeta.name,
                    ],
                ]
                [chunk, deepcopy(nosubMeta) + newMeta, grids.flatten()]
            }

        def diffs2D_sub_ionosub = chipsLssa.out
            // group grids from vis name together
            .groupTuple(by: 0)
            .filter { _, metas, __ -> metas.size() >= 2 }
            .map { chunk, metas, grids ->
                def ionosubMeta = metas.find { it.sub == 'ionosub'} ?: [:]
                def subMeta = metas.find { it.sub == 'sub'} ?: [:]
                def newMeta = [
                    ptype: '2D_diff',
                    pol: 'both',
                    title: "crosspower diff (ionosub-sub)\\n${chunk}",
                    plot_name: 'chips2d_diff_ionosub_sub',
                    name: ionosubMeta.name,
                    ext: ionosubMeta.ext,
                    tags: [
                        ionosubMeta.ext,
                        subMeta.ext,
                    ],
                    labels: [
                        ionosubMeta.name,
                        subMeta.name,
                    ],
                ]
                [chunk, deepcopy(ionosubMeta) + newMeta, grids.flatten()]
            }

        def ratios2D_sub = chipsLssa.out
            // group grids from vis name together
            .groupTuple(by: 0)
            .filter { _, metas, __ -> metas.size() >= 2 }
            .map { chunk, metas, grids ->
                def nosubMeta = metas.find { it.sub == null} ?: [:]
                def subMeta = metas.find { it.sub == 'sub'} ?: [:]
                def newMeta = [
                    ptype: '2D_ratio',
                    pol: 'both',
                    title: "crosspower ratio (nosub:sub)\\n${chunk}",
                    plot_name: 'chips2d_ratio_nosub_sub',
                    name: subMeta.name,
                    ext: subMeta.ext,
                    tags: [
                        nosubMeta.ext,
                        subMeta.ext,
                    ],
                    labels: [
                        nosubMeta.name,
                        subMeta.name,
                    ],
                    max_power: 1e2,
                    min_power: 1e-2,
                ]
                [chunk, deepcopy(nosubMeta) + newMeta, grids.flatten()]
            }

        def ratios2D_ionosub = chipsLssa.out
            // group grids from vis name together
            .groupTuple(by: 0)
            .filter { _, metas, __ -> metas.size() >= 2 }
            .map { chunk, metas, grids ->
                def nosubMeta = metas.find { it.sub == null} ?: [:]
                def ionosubMeta = metas.find { it.sub == 'ionosub'} ?: [:]
                def newMeta = [
                    ptype: '2D_ratio',
                    pol: 'both',
                    title: "crosspower ratio (nosub:ionosub)\\n${chunk}",
                    plot_name: 'chips2d_ratio_nosub_ionosub',
                    name: ionosubMeta.name,
                    ext: ionosubMeta.ext,
                    tags: [
                        nosubMeta.ext,
                        ionosubMeta.ext,
                    ],
                    labels: [
                        nosubMeta.name,
                        ionosubMeta.name,
                    ],
                    max_power: 1e2,
                    min_power: 1e-2,
                ]
                [chunk, deepcopy(nosubMeta) + newMeta, grids.flatten()]
            }

        singles1D
            .mix(singles2D)
            .mix(comps1D)
            .mix(diffs2D_sub)
            .mix(diffs2D_iono)
            .mix(diffs2D_sub_ionosub)
            .mix(ratios2D_sub)
            .mix(ratios2D_ionosub)
            .filter { chunk, meta, __ ->
                if (meta.name == null) {
                    print("\n -> meta.name is null in ${chunk}: ${meta}")
                    return false
                }
                if (meta.tags == null || meta.tags.contains(null)) {
                    print("\n -> meta.tags contains null in ${chunk}: ${meta}")
                    return false
                }
                return true
            }
            | chipsPlot

    emit:
        frame = chipsPlot.out
            .flatMap { _, meta, pngs ->
                def suffix = ""
                if (meta.nobs > 1) {
                    suffix = "_x" + sprintf("%04d", meta.nobs)
                }
                coerceList(pngs).collect { png ->
                    ["${meta.plot_name}_${meta.name}${suffix}", png]
                }
            }.groupTuple()
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
    obsCSV = channel.of(obsids_file)
        .splitCsv()
        .filter { line -> !line[0].startsWith('#') }
        .map { line -> def (obsid, cluster) = line
            def meta = [:]
            if (cluster != null) {
                meta.cluster = cluster
            }
            [obsid, meta]
        }
        .view()

    obsids = obsCSV.map { obsid, meta -> obsid }

    // analyse obsids with web services
    obsids | ws

    ws.out.obsMetafits
        .map { obsid, _ -> obsid }
        .collectFile(
            name: "ws_obs_pass.csv", newLine: true, sort: true,
            storeDir: "${results_dir}"
        )
        | view { [it, it.readLines().size()] }

    // download preprocessed, unless noprep is set
    if (params.noprep) {
        prep(channel.empty(), channel.empty())
    } else {
        prep(ws.out.obsMetafits, ws.out.obsMeta)
    }

    // channel of obsids that pass the flag gate
    prep.out.obsMetaVis
        .map { obsid, _, __ -> obsid }
        .collectFile(
            name: "prep_obs_pass.csv", newLine: true, sort: true,
            storeDir: "${results_dir}"
        )
        | view { [it, it.readLines().size()] }

    prep.out.obsFlags
        .filter { obsid, flagAnts, prepAnts, manualAnts, newAnts, _ -> newAnts.size() > 0 }
        .map { obsid, flagAnts, prepAnts, manualAnts, newAnts, _ ->
            ([
                obsid,
                displayInts(flagAnts),
                displayInts(prepAnts),
                displayInts(manualAnts),
                displayInts(newAnts),
            ]).join("\t")
        }
        .collectFile(
            name: "prep_flags.tsv", newLine: true, sort: true,
            seed: ([ "OBS", "ORIGINAL ANTS", "PREP ANTS", "MANUAL ANTS", "NEW ANTS" ]).join("\t"),
            storeDir: "${results_dir}"
        )
        // display output path and number of lines
        | view { [it, it.readLines().size()] }

    // combine new flags from prep.out.obsFlags with ws.out.obsMeta
    obsMeta = prep.out.obsFlags.join(ws.out.obsMeta)
        .map {obsid, flagAnts, prepAnts, manualAnts, newAnts, fineChanFlags, meta ->
            def newMeta = [prepFlags: newAnts]
            if (fineChanFlags) {
                newMeta.fineChanFlags = fineChanFlags
            }
            [obsid, deepcopy(meta) + newMeta]
        }
    qaPrep(prep.out.obsMetaVis, obsMeta)
    if (!params.novideo) {
        ws.out.frame
            .mix(prep.out.frame)
            .mix(qaPrep.out.frame)
            | makeVideos
    }
    // make zips
    if (params.tarchive) {
        prep.out.zip.mix(qaPrep.out.zip) | makeTarchives
    }
}

// given a channel of tuple (obs, metafits, vis), calibrate and analyse
workflow qaPrep {
    take:
        obsMetafitsVis
        obsMeta
    main:
        // get sourcelists for each obs (only currently used in subtraction, not calibration)
        obsMetafitsVis.map { obsid, metafits, _ -> [obsid, metafits ] }
            | (hypSrclistAO & hypSrclistYaml)

        // obsMetaSrclist: channel of tuple(obsid, metafits, srclist)
        // cluster unless --nocluster
        if (params.nocluster) {
            channel.empty() | rexCluster
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
            channel.empty() | cal
        } else {
            obsMetafitsVis.join(obsMeta)
                .filter { _, metafits, vis, __ -> metafits != null && vis != null}
                .map { obs, metafits, vis, meta ->
                    [obs, meta, metafits, vis]
                }
                | cal
        }

        cal.out.obsMetaCalPass
            .map { it -> def (obsid, meta, soln) = it;
                def calflags = deepcopy(meta.calflags?:[])
                ([obsid, meta.name, displayInts(meta.prepflags?:[]), displayInts(calflags)]).join("\t")
            }
            .collectFile(
                name: "passcal.tsv", newLine: true, sort: true,
                seed: ([ "OBS", "NAME", "PREPFLAGS", "CALFLAGS" ]).join("\t"),
                storeDir: "${results_dir}${params.cal_suffix}"
            )
            | view { [it, it.readLines().size()] }

        // channel of arguments for hypApply{UV,MS}
        // - take tuple(obsid, meta, soln) from cal.out.obsMetaCalPass
        // - match with tuple(obsid, metafits, prepUVFits) by obsid
        if (params.noapply) {
            allApply = channel.empty()
        } else {
            allApply = obsMetafitsVis
                .join(cal.out.obsMetaCalPass)
                .map { obsid, metafits, prepUVFits, meta, soln ->
                    def newMeta = [
                        time_res: params.apply_time_res,
                        freq_res: params.apply_freq_res,
                        nodut1: params.nodut1,
                        apply_args: params.apply_args,
                    ]
                    def calflags = deepcopy(meta.calflags?:[])
                    if (calflags) {
                        newMeta.apply_args = "${params.apply_args} --tile-flags ${calflags.join(' ')}"
                    }
                    // calflags = [8,9,10,11,12,13,15,27,30,39,40,42,44,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,79,88,89,90,91,92,93,94,95,117]
                    // newMeta.apply_args = "${newMeta.apply_args} --tile-flags ${calflags.join(' ')}"
                    newMeta.apply_name = "${newMeta.time_res}s_${newMeta.freq_res}kHz"
                    [obsid, deepcopy(meta) + deepcopy(newMeta), metafits, prepUVFits, soln]
                }
        }

        // apply calibration solutions to uvfits and ms unless --nouv or --noms
        if (params.nouv) {
            channel.empty() | hypApplyUV
        } else {
            allApply | hypApplyUV
        }
        if (params.noms) {
            channel.empty() | hypApplyMS
        } else {
            allApply | hypApplyMS
        }
        // get uvfits subtraction arguments from apply output.
        subArgsUV = obsMetaSrclist.cross(hypApplyUV.out)
            .map { obsMetaSrclist_, hypApplyUV_ ->
                def (obsid, metafits, srclist) = obsMetaSrclist_
                def (__, meta, vis) = hypApplyUV_;
                def newMeta = [sub_nsrcs: params.sub_nsrcs]
                [obsid, deepcopy(meta) + newMeta, metafits, vis, srclist]
            }
        subArgsMS = obsMetaSrclist.cross(hypApplyMS.out)
            .map { obsMetaSrclist_, hypApplyUV_ ->
                def (obsid, metafits, srclist) = obsMetaSrclist_
                def (__, meta, vis) = hypApplyUV_;
                def newMeta = [sub_nsrcs: params.sub_nsrcs]
                [obsid, deepcopy(meta) + newMeta, metafits, vis, srclist]
            }
        if (params.nosub) {
            channel.empty() | (hypSubUV & hypSubMS)
        } else {
            subArgsUV | hypSubUV
            subArgsMS | hypSubMS
        }
        if (params.noionosub) {
            channel.empty() | (hypIonoSubUV & hypIonoSubMS)
        } else {
            subArgsUV.map {obsid, meta, metafits, vis, srclist ->
                def newMeta = [ionosub_nsrcs: params.ionosub_nsrcs]
                [obsid, deepcopy(meta) + newMeta, metafits, vis, srclist]}
                | hypIonoSubUV
            subArgsMS.map {obsid, meta, metafits, vis, srclist ->
                def newMeta = [ionosub_nsrcs: params.ionosub_nsrcs]
                [obsid, deepcopy(meta) + newMeta, metafits, vis, srclist]}
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

        // collect ionoqa results as .tsv
        cthulhuPlot.out
            // form row of tsv from json fields we care about
            .map { obsid, _, __, ___, json ->
                def stats = parseJson(json);
                [
                    obsid,
                    stats.PCA_EIGENVALUE?:'',
                    stats.MED_ABS_OFFSET?:'',
                    stats.QA?:'',
                    stats.MAX_SHIFT_SRC?:'',
                    stats.MAX_SHIFT_RA?:'',
                    stats.MAX_SHIFT_DEC?:'',
                    stats.N_TIMES?:'',
                    stats.N_SRCS?:'',
                ].join("\t")
            }
            .collectFile(
                name: "ionoqa.tsv", newLine: true, sort: true,
                seed: [
                    "OBS",
                    "HYP_IONO_PCA",
                    "HYP_IONO_MAG",
                    "HYP_IONO_QA",
                    "MAX_SHIFT_SRC",
                    "MAX_SHIFT_RA",
                    "MAX_SHIFT_DEC",
                    "N_TIMES",
                    "N_SRCS",
                ].join("\t"),
                storeDir: "${results_dir}"
            )
            // display output path and number of lines
            | view { [it, it.readLines().size()] }

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
                def newMeta = [
                    lowfreq:jsonMeta.freqs[0],
                    first_lst: jsonMeta.times[0]['lst_rad'],
                    first_jd: jsonMeta.times[0]['jd1'],
                    nchans:jsonMeta.freqs.size(),
                    ntimes:jsonMeta.times.size()
                ];
                ['eorband', 'eorfield', 'config'].each { key ->
                    if (jsonMeta[key] != null) {
                        newMeta[key] = jsonMeta[key]
                    }
                }
                [obsid, deepcopy(meta) + newMeta, vis] }

        // QA uvfits visibilities
        obsMetaUVSmart | uvfits

        // channel of calibrated, subtracted and ionosubtracted ms: tuple(obsid, name, vis)
        obsMetaMS = hypApplyMS.out.map { obsid, meta, vis, _ -> [obsid, meta, vis] }
            .mix(hypIonoSubMS.out.map { obsid, meta, vis, _, __ -> [obsid, meta, vis] })
            .mix(hypSubMS.out.map { obsid, meta, vis, _ -> [obsid, meta, vis] })

        // image and qa measurementsets or uvfits unless --noimage
        if (params.noimg) {
            channel.empty() | img
            obsMetaImgPass = uvfits.out.obsMetaUVPass.map { obsid, meta, _ -> [obsid, meta, [], []]}
        } else {
            // visibilities for imaging and power spectrum
            // filter by visibilities that pass uvfits QA if --nouv is not set
            if (params.nouv) {
                obsMetaVisPass = obsMetaMS
            } else {
                obsMetaVisPass = uvfits.out.obsMetaUVPass
            }

            obsMetaVisPass
                .map { obsid, meta, vis -> [obsid, meta, [vis]]}
                | img

            obsMetaImgPass = img.out.obsMetaImgPass
        }

        // tuple of (obsid, meta) for all visibilities that pass imgQA
        obsMetaPass = obsMetaImgPass.map { it -> def (obsid, meta) = it; [obsid, meta] }

        obsMetaPass
            .map { it -> def (obsid, meta) = it;
                ([obsid, meta.name]).join("\t")
            }
            .collectFile(
                name: "img_names_pass.tsv", newLine: true, sort: true,
                seed: ([ "OBS", "NAME"]).join("\t"),
                storeDir: "${results_dir}${params.img_suffix}${params.cal_suffix}"
            )
            | view { [it, it.readLines().size()] }

        // Apologies in advance for anyone who has to debug this.
        //
        // We want a tuple of (chunk, chunkMeta obsids[G]) for imaging and power spectrum.
        //
        // The pipeline can produce multiple visibilities for a given obsid. `meta` keys can be used to differentiate these. e.g.
        // - `sub` - type of subtraction (`null` means no subtraction)
        // - `cal_prog` - the program used to calibrate the visibility
        // - `poly` - whether a polyfit was applied
        //
        // We want to compare results between visibilities of different types within a single obsid, and within compatible groups of obsids.
        // Within the set of visibilities for each obsid, one is the "primary" visibility that others are compared to, e.g. the one with no subtractions.
        // A group of obsids is compatible if they have the same field, band and telescope configuration.
        // The metadata for each group of obsids comes from primary meta
        // Within each group of obsids, we sort by some metric in the primary meta  (e.g. window : wedge power ratio)
        // Then we split that group into chunks of size G.
        chunkMetaPass = obsMetaPass
            // TODO: chips only works with even timesteps ?
            // .filter { obsid, meta, vis ->
            //     def isEven = meta.ntimes % 2 == 0;
            //     if (!isEven) {
            //         println "can't run chips on ${obsid} ${meta.name} because it has an odd number of timesteps"
            //     }
            //     isEven
            // }
            // group metas by obsid
            .groupTuple(by: 0)
            // ensure metas are a list, even when there is only one meta
            .map { obsid, metas -> [ obsid, coerceList(metas) ]}
            // determine how to group each obsid, and how to sort that obsid within each group.
            .map { obsid, metas ->
                // find the meta of the primary subobservation (in this case, the unsubtracted visibilities)
                def metaPrime = metas.find { meta -> meta.sub == null }
                // determine the value to sort obsids by
                def sort = Float.NaN
                if (metaPrime.p_window?:Float.NaN != Float.NaN && metaPrime.p_wedge?:Float.NaN != Float.NaN) {
                    sort = metaPrime.p_window / metaPrime.p_wedge
                }
                // group by field, band, config
                def group_tokens = []
                def first_token = ""
                if (metaPrime.eorfield != null) {
                    first_token += "eor${metaPrime.eorfield}"
                }
                if (metaPrime.eorband != null) {
                    first_token += (metaPrime.eorband == 0 ? "low" : "high")
                }
                if (first_token.size() > 0) {
                    group_tokens << first_token
                }
                if (metaPrime.config != null) {
                    group_tokens << metaPrime.config
                }
                // if (metaPrime.pointing != null) {
                //     group_tokens << "p${metaPrime.pointing}"
                // }
                // if (metaPrime.lst != null) {
                //     def nearest_lst = ((metaPrime.lst.round().intValue() + 180) % 360 - 180)
                //     group_tokens << sprintf("lst%+03d", nearest_lst)
                // }
                def group = group_tokens.join('_')

                if (metas.size() > 1) {
                    [
                        "cal_prog", "time_res", "freq_res", "lowfreq", "nchans",
                        "eorband", "eorfield", "config"
                    ].each { key ->
                        assert metas[1..-1].every { meta -> meta[key] == metas[0][key] }, \
                            "meta key ${key} not consistent across subtractions for obsid ${obsid}: ${metas}"
                    }
                }
                [group, sort, obsid, metas]
            }
            .groupTuple(by: 0)
            // only keep groups with more than one obsid
            .filter { it -> it[1] instanceof List }
            // .view { it -> "\n -> grouped by obsid: ${it}"}
            // chunk obsid groups into subgroups of 20
            .flatMap { group, all_sorts, all_obsids, all_metas ->
                def groupMeta = [group: group]

                [all_sorts, all_obsids, all_metas].transpose()
                    .sort { it -> it[0] }
                    .collate(20, true)
                    .collect { chunk ->
                        def (sorts, obsids, metas) = chunk.transpose()
                        def obs_list = obsids.sort()
                        def hash = obs_list.join(' ').md5()[0..7]
                        def chunkMeta = [
                            // ntimes: metas.collect { meta -> meta[0].ntimes?:0 }
                            sort_bounds: [sorts[0], sorts[-1]],
                            hash: hash,
                        ]
                        ["${group}_${hash}", deepcopy(groupMeta) + chunkMeta, obs_list, metas]
                    }
            }
            // .view { it -> "\n -> chunked by obsid: ${it}"}
            .tap { groupChunkMetaPass }
            // flatten obsids out of each chunk
            .flatMap { chunk, _, obsids, all_metas ->
                [obsids, all_metas].transpose().collect { obsid, metas ->
                    [chunk, obsid, metas]
                }
            }
            // .view { it -> "\n -> flattened obsid: ${it}"}
            // flatten visibilities out of each obsid
            .flatMap { chunk, obsid, metas ->
                metas.collect { meta ->
                    [chunk, meta.name, obsid, meta]
                }
            }
            // .view { it -> "\n -> flattened vis type: ${it}"}
            // group by both obs group and visibility type name
            .groupTuple(by: 0..1)
            .map { chunk, name, obsids, metas ->
                def groupMeta = [name: name]
                [
                    "cal_prog", "time_res", "freq_res", "lowfreq", "nchans",
                    "eorband", "eorfield", "config", "sub"
                ].each { key ->
                    if (metas[0][key]) {
                        groupMeta[key] = metas[0][key]
                    }
                }
                [chunk, groupMeta, obsids]
            }
            // finally, we have a channel of visibilities from the same obs group and vis type
            // .view { it -> "\n -> chunkMetaPass: ${it}"}

        groupChunkMetaPass
            .map { chunk, chunkMeta, obsids, all_metas ->
                def sort_bounds = chunkMeta.sort_bounds?:[Float.NaN, Float.NaN]
                (
                    sort_bounds.collect { sprintf("%5.4f", it) }
                    + [chunk, obsids.size(), obsids.join(' ')]
                ).join("\t")
            }
            .collectFile(
                name: "obs_group_chunks.tsv", newLine: true, sort: true,
                seed: ([ "SORT LEFT", "SORT RIGHT", "GROUP CHUNK", "NOBS", "OBSIDS" ]).join("\t"),
                storeDir: "${results_dir}${params.img_suffix}${params.cal_suffix}"
            )
            | view { [it, it.readLines().size()] }

        chunkMetaPass.map { group, _, obsids -> [group, obsids] }
            | storeManifest

        chunkMetaPass
            .map { group, meta, obsids ->
                [group, meta.name, obsids.size(), obsids.join(' ')].join("\t")
            }
            .collectFile(
                name: "obs_groups.tsv", newLine: true, sort: true,
                seed: ([ "GROUP", "NAME", "NVIS", "VISS" ]).join("\t"),
                storeDir: "${results_dir}${params.img_suffix}${params.cal_suffix}"
            )
            | view { [it, it.readLines().size()] }

        // TODO: filter uvfits.out.obsMetaUVPass to only include obsids that pass imgQA
        obsMetaUVPass = uvfits.out.obsMetaUVPass.cross(obsMetaPass) { def (obsid, meta) = it; [obsid, meta.name] }
            .map { obsMetaUVPass_, obsMetaPass_ ->
                def (obsid, _, uvfits) = obsMetaUVPass_;
                def (__, meta) = obsMetaPass_;
                [obsid, meta, uvfits]
            }

        if (params.noimgcombine) {
            imgCombine(channel.empty(), channel.empty())
        } else {
            imgCombine(obsMetaUVPass, chunkMetaPass)
        }

        if (params.nochipscombine) {
            chips(obsMetaUVPass, channel.empty())
        } else {
            chips(obsMetaUVPass, chunkMetaPass)
        }

        // stack images from chunks
        if (!params.nostackthumbs) {
            obsMetaImgPass
                .cross(
                    chunkMetaPass.flatMap { group, groupMeta, obsids ->
                        obsids.collect { obsid -> [obsid, groupMeta, group] }
                    }
                ) { def (obsid, meta) = it; [obsid, meta.name] }
                .flatMap { obsMetaImgPass_, obsMetaPass_ ->
                    def (obsid, _, imgMetas, imgs) = obsMetaImgPass_;
                    def (__, groupMeta, group) = obsMetaPass_;
                    [imgMetas, imgs].transpose().collect { imgMeta, img ->
                        [group, groupMeta.name, imgMeta.suffix, imgMeta, img]
                    }
                }
                .groupTuple(by: 0..2)
                .map { group, _, __, imgMetas, imgs ->
                    [group, imgMetas[0], imgs]
                }
                | stackImgs
                | stackThumbnail
        }

        // archive data to object store
        if (params.archive) {
            if (params.archive_prep) {
                prep_archive = obsMetafitsVis.map { _, __, vis -> ["prep", vis] }
            } else {
                prep_archive = channel.empty()
            }
            if (params.archive_uvfits) {
                vis_archive = obsMetaUV.map { _, __, vis -> ["uvfits", vis]}
            } else {
                vis_archive = channel.empty()
            }
            prep_archive.mix(vis_archive)
                .mix(cal.out.obsMetaCalPass.map { _, __, soln -> ["soln", soln] })
                .mix(cal.out.archive)
                .mix(uvfits.out.archive)
                .mix(img.out.archive)
                | archive
        }

    emit:
        frame = cal.out.frame
            .mix(uvfits.out.frame)
            .mix(img.out.frame)
            .mix(chips.out.frame)
            .mix(cthulhuPlot.out.flatMap {
                def (_, meta, pngs) = it;
                pngs.collect { png -> ["cthulhuplot_${meta.name}", png] }
            }.groupTuple())
        zip = cal.out.zip
            .mix(uvfits.out.zip)
            .mix(img.out.zip)
            .mix(
                hypIonoSubUV.out.map { it -> def (obsid, meta, _, offsets) = it; ["offsets_${meta.name}", offsets]}
                .groupTuple()
            )

}

workflow makeVideos {
    take:
        frame
    emit:
        videos = frame.map { name, frames ->
                def latest = frames.collect { it.lastModified() }.max()
                def cachebust = "${latest}_x" + sprintf("%04d", frames.size())
                def sorted = frames.collect { path -> file(deepcopy(path.toString())) }.sort(false)
                [name, sorted, cachebust]
            }
            | ffmpeg
            | view { video, _ -> [video, video.size()] }
}

workflow makeTarchives {
    take:
        zip
    emit:
        zips = zip.map { name, files ->
                def latest = files.collect { file -> file.lastModified() }.max()
                def cachebust = "${latest}_x" + sprintf("%04d", files.size())
                // def sorted = files.collect { path -> file(deepcopy(path.toString())) }.sort(false)
                [name, files, cachebust]
            }
            | tarchive
            | view { zip_, cachebust -> [zip_, zip_.size()] }
}

// entrypoint: get externally preprocessed uvfits files from csv file and run qa
workflow extPrep {
    obsVis = channel.of(file("external${params.external_suffix}.csv"))
        .splitCsv(header: ["obsid", "vis"])
        .filter { !it.obsid.startsWith('#') }
        .map { [it.obsid, file(it.vis)] }

    obsids = obsVis.map { obsid, _ -> obsid }
    obsids | ws

    qaPrep(ws.out.obsMetafits.join(obsVis), channel.empty())

    if (!params.novideo) {
        ws.out.frame
            .mix(qaPrep.out.frame)
            | makeVideos
    }
    // make zips
    if (params.tarchive) {
        qaPrep.out.zip | makeTarchives
    }
}

