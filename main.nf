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

    env | sort
    df \$TMPDIR

    export http_proxy="http://proxy.per.dug.com:3128"
    export https_proxy="http://proxy.per.dug.com:3128"
    export all_proxy="proxy.per.dug.com:3128"
    export ftp_proxy="http://proxy.per.dug.com:3128"

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

    # try and submit a job, if it's already there this will fail, and we can suppress the
    # warning, otherwise we can download and exit
    ${params.giant_squid} submit-vis --delivery "acacia" $obsid -w \
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

    // this is necessary for singularity
    stageInMode "copy"
    module "singularity"

    maxRetries 1

    input:
    tuple val(obsid), path("${obsid}.metafits"), path("*") // <-preserves names of fits files

    output:
    tuple val(obsid), path("${obsid}.uvfits"), path("${obsid}*.mwaf")

    script:
    """
    export flag_template="${obsid}_%%.mwaf"
    if [ ${obsid} -gt 1300000000 ]; then
        flag_template="${obsid}_ch%%%.mwaf"
    fi
    ${params.birli} \
        -u "${obsid}.uvfits" \
        -f \$flag_template \
        -m "${obsid}.metafits" \
        ${obsid}*.fits
    """
}

// QA tasks that can be run on each preprocessed obs or its flags. The exit status of this job determines
// if this obs is good to use for calibration.
process prepQA {
    storeDir "$params.outdir/${obsid}/QA"
    input:
    tuple val(obsid), path("${obsid}.metafits"), path("${obsid}.uvfits"), path("*") // <-preserves names of mwaf files

    output:
    // TODO: change this to whatever paths QA outputs
    tuple val(obsid), path("qa_prep_results_${obsid}.txt") optional true

    // if this script returns failure, don't progress this obs to calibration, just ignore
    errorStrategy 'ignore'

    script:
    """
    echo "TODO: QA on preprocessed obs, possibly AOFlagger Occupancy, dead dipole fraction"
    # example: deliberately fail one obsid
    # if [${obsid} -eq "1322308000"]; then
    if false; then
        exit 1
    else
        ls -al ${obsid}.metafits ${obsid}.uvfits ${obsid}*.mwaf | tee qa_prep_results_${obsid}.txt
        ls -al qa_prep_results_${obsid}.txt
    fi
    """
}

// ensure calibration solutions and vis are present, or calibrate prep with hyperdrive
process ensureCal {
    storeDir "$params.outdir/${obsid}/cal"

    // label jobs that need a bigger gpu allocation
    label "gpu"

    module "cuda/11.3.1:gcc-rt/9.2.0"
    stageInMode "copy"

    input:
    tuple val(obsid), path("${obsid}.metafits"), path("${obsid}.uvfits"), \
        val(name), val(cal_args), val(apply_args)

    output:
    tuple val(obsid), val(name), path("soln_${obsid}_${name}.fits"), path("${obsid}_${name}.uvfits")

    script:
    """
    export HYPERDRIVE_CUDA_COMPUTE=${params.cuda_compute}
    ${params.hyperdrive} di-calibrate ${cal_args} \
        --data "${obsid}.metafits" "${obsid}.uvfits" \
        --beam "${params.beam_path}" \
        --source-list "${params.sourcelist}" \
        --outputs "soln_${obsid}_${name}.fits"
    ${params.hyperdrive} solutions-apply ${apply_args} \
        --data "${obsid}.metafits" "${obsid}.uvfits" \
        --solutions "soln_${obsid}_${name}.fits" \
        --outputs "${obsid}_${name}.uvfits"
    """
}

// QA tasks that can be run on each calibration parameter set
process calQA {
    storeDir "$params.outdir/${obsid}/QA"
    errorStrategy 'ignore'

    input:
    tuple val(obsid), val(name), path("soln_${obsid}_${name}.fits"), path("${obsid}_${name}.uvfits")

    output:
    // TODO: change this to whatever paths QA outputs
    tuple val(obsid), val(name), path("qa_cal_results_${obsid}_${name}.txt") optional true

    script:
    """
    echo "TODO: QA on calibration set ${name} for ${obsid}"

    # example: deliberately fail one obsid, name
    // if [${obsid} -eq "1322308000"] && [${name} == "50l_src4k_8s_80kHz"]; then
    if false; then
        exit 1
    else
        ls -al soln_${obsid}_${name}.fits ${obsid}_${name}.uvfits | tee qa_cal_results_${obsid}_${name}.txt
        ls -al qa_cal_results_${obsid}_${name}.txt
    fi
    """
}

// ensure calibrated visiblities are present, or apply solutions with hyperdrive

workflow {
    obsids = channel.fromPath(params.obsids_path).splitCsv().flatten()

    // ensure raw files are downloaded
    ensureRaw(obsids)

    // ensure preprocessed files have been generated from raw files
    ensurePrep(ensureRaw.out)

    // QA preprocessed files to prevent a bad obs from reaching calibration
    // - get the obsid and metafits from ensureRaw 
    // - cross with the uvfits and mwaf from ensurePrep
    ensureRaw.out \
        .cross(ensurePrep.out) \
        .map( it -> [it[0][0], it[0][1], it[1][1], it[1][2]]) \
        | prepQA

    // create a channel of calibration parameter sets from a csv with columns:
    // - short name that appears in the output filename
    // - args for hyperdrive di-cal, where `-n` = `--num-sources`
    // - args for hyperdrive soln apply args
    cal_param_sets = channel \
        .fromPath(params.cal_params_path) \
        .splitCsv(header: false, skip: 1, strip: true)

    // calibration for each obsid that passes QA, and each parameter set
    // - get obsid, metafits from ensureRaw cross with uvfits from ensurePrep
    // - cross with prepQA to reject jobs that failed QA
    // - combine with each calibration parameter set
    ensureRaw.out \
        .cross(ensurePrep.out) \
        .map( it -> [it[0][0], it[0][1], it[1][1]]) \
        .cross(prepQA.out) \
        .map( it -> [it[0][0], it[0][1], it[0][2]]) \
        .combine(cal_param_sets) \
        | ensureCal \
        | calQA

    // TODO: gather QA results
    // pick "best" param set? 
}