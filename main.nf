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

    # save ASVO job info for this obsid to tsv
    # also save all jobs to json for debugging
    function cache_jobs {
        giant-squid list -j |
            tee jobs.json |
            jq -r '.[] \
                |select(.obsid==$obsid) \
                |[.jobId,.jobState,.files[0].fileUrl//""] \
                |@tsv' |
            tee jobs.tsv
    }

    # download any ASVO jobs for the obsid that are ready
    function get_ready_jobs {
        cache_jobs
        while read jobid state url size; do
            if [ "\$state" = "Ready" ] && [ -n "\$url" ]; then
                echo "[obs:$obsid]: Downloading from \$url"
                curl \$url | tar -x
                exit 0
            fi
        done <jobs.tsv
        return 1
    }

    # download any ready jobs, exit if success, else suppress errors
    get_ready_jobs && exit 0 || true

    # try and submit a job, if it's already there this will fail
    giant-squid submit-vis --delivery "acacia" $obsid -w || true

    cache_jobs

    # wait for any "Queued" or "Processing" jobs
    while read jobid state url size; do
        if [ "\$state" = "Queued" ] || [ "\$state" = "Processing" ]; then
            echo "[obs:$obsid]: waiting until \$jobid with state \$state is Ready"
            giant-squid wait \$jobid
            break
        fi
    done <jobs.tsv

    # download any ready jobs, don't suppress errors
    get_ready_jobs
    exit \$?
    """
}

// Ensure preprocessed observation is present, or preprocess raw with Birli
process ensurePrep {
    storeDir "$params.outdir/${obsid}/prep"

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