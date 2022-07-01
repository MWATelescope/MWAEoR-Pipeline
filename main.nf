#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Ensure the raw visibility files are present, or download via ASVO
process ensureRaw {
    // persist results in outdir, process will be skipped if files already present.
    storeDir "$params.outdir/$obsid/raw"
    // allow multiple retries
    maxRetries 5
    // exponential backoff: sleep for 2^attempt minutes after each fail
    errorStrategy { sleep(Math.pow(2, task.attempt) * 60000 as long); return 'retry' }

    input:
    val obsid

    output:
    tuple val(obsid), path("${obsid}.metafits"), path("${obsid}_2*.fits")

    script:
    """
    # echo commands, exit on any failures
    set -ex
    export MWA_ASVO_API_KEY="${params.asvo_api_key}"

    function ensure_disk_space {
        local needed=\$1
        while read -r avail; do
            if [[ \$avail -lt \$needed ]]; then
                echo "Not enough disk space available in \$PWD, need \$needed, have \$avail"
                exit 1
            fi
        done < <(df --output=avail . | tail -n 1)
    }

    # download any ASVO jobs for the obsid that are ready
    function get_first_ready_job {
        # extract id url and size from ready download vis jobs
        ${params.giant_squid} list -j --types download_visibilities --states ready -- $obsid \
            | tee /dev/stderr \
            | ${params.jq} -r '.[]|[.jobId,.files[0].fileUrl//"",.files[0].fileSize//""]|@tsv' \
            | tee ready.tsv
        while read jobid url size; do
            ensure_disk_space \$size || exit \$?
            if [ -n "\$url" ]; then
                echo "[obs:$obsid]: Downloading job \$jobid from \$url (\$size bytes)"
                curl \$url | tar -x
                exit 0
            fi
        done <ready.tsv
        return 1
    }

    # download any ready jobs, exit if success, else suppress errors
    get_first_ready_job && exit 0 || true

    # try and submit a job, if it's already there this will fail
    ${params.giant_squid} submit-vis --delivery "acacia" $obsid -w || true

    # extract id and state from pending download vis jobs
    ${params.giant_squid} list -j --types download_visibilities --states queued processing -- $obsid \
        | ${params.jq} -r '.[]|[.jobId,.jobState]|@tsv' \
        | tee pending.tsv

    # for each pending job for the obsid, wait for completion and try downloading again
    while read jobid state; do
        echo "[obs:$obsid]: waiting until \$jobid with state \$state is Ready"
        ${params.giant_squid} wait -j -- \$jobid
        # a new job is ready, so try and get it
        get_first_ready_job && exit 0 || true
    done <pending.tsv
    exit 1
    """
}

// Ensure preprocessed observation is present, or preprocess raw with Birli
process ensurePrep {
    storeDir "$params.outdir/${obsid}/prep"

    // label jobs that need a bigger cpu allocation
    label "cpu"

    input:
    tuple val(obsid), path("${obsid}.metafits"), path("${obsid}_2*.fits")

    output:
    tuple val(obsid), path("${obsid}.uvfits"), path("${obsid}*.mwaf")

    script:
    """
    ${params.birli} \
        -u "${obsid}.uvfits" \
        -f "${obsid}_%%.mwaf" \
        -m "${params.outdir}/${obsid}/raw/${obsid}.metafits" \
        "${params.outdir}/${obsid}/raw/${obsid}_2"*.fits
    """
}

// ensure calibration solutions and vis are present, or calibrate prep with hyperdrive
process ensureCal {
    storeDir "$params.outdir/${obsid}/cal"

    // label jobs that need a bigger gpu allocation
    label "gpu"

    input:
    tuple val(obsid), path("${obsid}.metafits"), path("${obsid}.uvfits"), \
        val(name), val(cal_args), val(apply_args)

    output:
    tuple val(obsid), val(name), path("soln_${obsid}_${name}.fits"), path("${obsid}_${name}.uvfits")

    script:
    """
    export HYPERDRIVE_CUDA_COMPUTE=${params.cuda_compute}
    ${params.hyperdrive} di-calibrate ${cal_args} \
        --data "${params.outdir}/$obsid/raw/${obsid}.metafits" "${obsid}.uvfits" \
        --beam "${params.beam_path}" \
        --source-list "${params.sourcelist}" \
        --outputs "soln_${obsid}_${name}.fits"
    ${params.hyperdrive} solutions-apply ${apply_args} \
        --data "${obsid}.metafits" "${obsid}.uvfits" \
        --solutions "soln_${obsid}_${name}.fits" \
        --outputs "${obsid}_${name}.uvfits"
    """
}

// QA tasks that can be run on each preprocessed obs or its flags
process prepQA {
    storeDir "$params.outdir/${obsid}/QA"
    input:
    tuple val(obsid), path("${obsid}.metafits"), path("${obsid}.uvfits"), path("${obsid}*.mwaf")

    output:
    // TODO: change this to whatever paths QA outputs
    path("qa_prep_results_${obsid}.txt")

    script:
    """
    echo "TODO: QA on preprocessed obs"
    ls -al ${obsid}.metafits ${obsid}.uvfits ${obsid}*.mwaf > qa_prep_results_${obsid}.txt
    """
}

// QA tasks that can be run on each calibration parameter set
process calQA {
    storeDir "$params.outdir/${obsid}/QA"
    input:
    tuple val(obsid), val(name), path("soln_${obsid}_${name}.fits"), path("${obsid}_${name}.uvfits")

    output:
    // TODO: change this to whatever paths QA outputs
    path("qa_cal_results_${obsid}_${name}.txt")

    script:
    """
    echo "TODO: QA on calibration set ${name} for ${obsid}"
    ls -al soln_${obsid}_${name}.fits ${obsid}_${name}.uvfits > qa_cal_results_${obsid}_${name}.txt
    """
}

// ensure calibrated visiblities are present, or apply solutions with hyperdrive

workflow {
    obsids = channel.fromPath(params.obsids_path).splitCsv().flatten()

    // calibration parameter sets:
    // - short name
    // - hyperdrive di-cal args
    //   - `-n` = `--num-sources`
    // - hyperdrive soln apply args
    // this is sourced from a CSV file with a parameter set on each line.
    // we want to run a different calibration job for each parameter set
    cal_param_sets = channel \
        .fromPath(params.cal_params_path) \
        .splitCsv(header: false, skip: 1, strip: true)

    ensureRaw(obsids)
    ensurePrep(ensureRaw.out)
    // prep QA on ensureRaw[obsid, metafits] cross ensurePrep[uvfits, mwaf]
    ensureRaw.out \
        .cross(ensurePrep.out) \
        .map( it -> [it[0][0], it[0][1], it[1][1], it[1][2]]) \
        | prepQA

    // calibrate on the cartesian product of:
    //  - ensureRaw[obsid, metafits] cross ensurePrep[uvfits]
    //  - cal_param_set[0:3]
    ensureRaw.out \
        .cross(ensurePrep.out) \
        .map( it -> [it[0][0], it[0][1], it[1][1]])\
        .combine(cal_param_sets) \
        | ensureCal \
        | calQA

    // TODO: gather QA results
}