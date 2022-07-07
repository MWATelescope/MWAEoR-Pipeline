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
    tuple val(obsid), path("birli_${obsid}.uvfits"), path("birli_${obsid}.ms.zip"), path("${obsid}*.mwaf"), path("birli_prep.log")

    script:
    """
    export flag_template="${obsid}_%%.mwaf"
    if [ ${obsid} -gt 1300000000 ]; then
        flag_template="${obsid}_ch%%%.mwaf"
    fi
    ${params.birli} \
        -u "birli_${obsid}.uvfits" \
        -M "birli_${obsid}.ms" \
        -f \$flag_template \
        -m "${obsid}.metafits" \
        ${obsid}*.fits | tee birli_prep.log
    zip -r "birli_${obsid}.ms.zip" "birli_${obsid}.ms"
    """
}

// QA tasks for flags.
process flagQA {
    storeDir "$params.outdir/${obsid}/flag_qa"
    input:
    tuple val(obsid), path("${obsid}.metafits"), path("birli_${obsid}.uvfits"), path("*") // <-preserves names of mwaf files

    output:
    // TODO: change this to whatever paths QA outputs
    tuple val(obsid), path("*.json")

    // if this script returns failure, don't progress this obs to calibration, just ignore
    errorStrategy 'ignore'

    script:
    """
    echo "TODO: QA AOFlagger Occupancy, dead dipole fraction"
    ls -al ${obsid}.metafits *.uvfits *.mwaf
    touch flag_qa_${obsid}.json
    exit 0
    """
}

// ensure calibration solutions and vis are present, or calibrate prep with hyperdrive
process ensureCal {
    storeDir "$params.outdir/${obsid}/cal$params.cal_suffix"

    // label jobs that need a bigger gpu allocation
    label "gpu"

    module "cuda/11.3.1:gcc-rt/9.2.0"
    stageInMode "copy"

    input:
    tuple val(obsid), path("${obsid}.metafits"), path("${obsid}.uvfits"), path("cal_params.csv")

    output:
    tuple val(obsid), path("hyp_soln_${obsid}_*.fits"), path("hyp_${obsid}_*.uvfits"), path("hyp_${obsid}_*.ms.zip"), \
        path("hyp_di-cal_${obsid}_*.log"), path("hyp_apply_${obsid}_*.log")

    script:
    """
    export HYPERDRIVE_CUDA_COMPUTE=${params.cuda_compute}
    # use this control which GPU we're sending the job to
    export CUDA_VISIBLE_DEVICES=0

    set -ex
    export names=""
    while IFS=, read -r name cal_args apply_args; do
        # on first iteration, this does nothing. on nth, wait for all background jobs to finish
        if [[ \$CUDA_VISIBLE_DEVICES -eq 0 ]]; then
            wait \$(jobs -rp)
        fi
        # hyperdrive di-cal then solutions apply, backgrounded in a subshell on the same gpu
        (
            ${params.hyperdrive} di-calibrate \${cal_args} \
                --data "${obsid}.metafits" "${obsid}.uvfits" \
                --beam "${params.beam_path}" \
                --source-list "${params.sourcelist}" \
                --outputs "hyp_soln_${obsid}_\${name}.fits" \
                > hyp_di-cal_${obsid}_\${name}.log
            ${params.hyperdrive} solutions-apply \${apply_args} \
                --data "${obsid}.metafits" "${obsid}.uvfits" \
                --solutions "hyp_soln_${obsid}_\${name}.fits" \
                --outputs "hyp_${obsid}_\${name}.uvfits" "hyp_${obsid}_\${name}.ms" \
                > hyp_apply_${obsid}_\${name}.log
            zip -r "hyp_${obsid}_\${name}.ms.zip" "hyp_${obsid}_\${name}.ms"
        ) &
        names="\${names} \${name}"
        # increment the target device mod num_gpus
        CUDA_VISIBLE_DEVICES=\$((CUDA_VISIBLE_DEVICES+1%${params.num_gpus}))
    done < cal_params.csv

    # wait for all the background jobs to finish
    wait \$(jobs -rp)

    # check all the files we expect are present.
    for name in \$names; do
        # print out important info from log
        for log in hyp_di-cal_${obsid}_\${name}.log hyp_apply_${obsid}_\${name}.log; do
            if [ ! -f \$log ]; then
                echo "ERROR: \$log not found"
                exit 1
            else
                grep -iE "err|warn|hyperdrive|chanblocks|reading|writing|flagged" \$log
            fi
        done
        for file in "hyp_soln_${obsid}_\${name}.fits" "hyp_${obsid}_\${name}.uvfits" "hyp_${obsid}_\${name}.ms.zip"; do
            if [ ! -f \$file ]; then
                echo "Missing file \$file"
                exit 1
            fi
        done
    done
    """
}

// QA tasks that can be run on each calibration parameter set
process calQA {
    storeDir "$params.outdir/${obsid}/cal_qa"
    errorStrategy 'ignore'

    input:
    tuple val(obsid), path("${obsid}.metafits"), path("*") // <- hyp_soln_*.fits

    output:
    // TODO: change this to whatever paths QA outputs
    tuple val(obsid), path("*.json")

    stageInMode "copy"
    module "singularity"

    script:
    """
    export name="\$(basename *.fits)"
    name="\${name%.*}"
    singularity exec --bind \$PWD:/tmp --writable-tmpfs --pwd /tmp --no-home ${params.mwa_qa_sif} \
        cal_qa *.fits *.metafits X --out "\${name}_X.json"
    singularity exec --bind \$PWD:/tmp --writable-tmpfs --pwd /tmp --no-home ${params.mwa_qa_sif} \
        cal_qa *.fits *.metafits Y --out "\${name}_Y.json"
    """ 
}

// create dirty iamges of xx,yy,v
process wscleanDirty {
    storeDir "$params.outdir/${obsid}/img${params.img_suffix}"
    input:
    tuple val(obsid), val(name), path("*") // *.ms.zip

    output:
    tuple val(obsid), val(name), path("wsclean_${name}*-{XX,YY,V}-dirty.fits")

    label "cpu"
    label "img"

    script:
    """
    unzip *.ms.zip
    ${params.wsclean} \
        -weight briggs -1.0 \
        -name wsclean_${name} \
        -size ${params.img_size} ${params.img_size} \
        -scale ${params.img_scale} \
        -pol xx,yy,v \
        -abs-mem ${params.img_mem} \
        -channels-out ${params.img_channels_out} \
        *.ms
    """
}

// create source list from deconvolved image of stokes I
process wscleanSources {
    // - shallow deconvolution, 1000 iter just enough for removing sidelobes
    // - channels: need mfs
    storeDir "$params.outdir/${obsid}/img${params.img_suffix}"
    input:
    tuple val(obsid), val(name), path("*") // *.ms.zip

    output:
    tuple val(obsid), val(name), path("wscleanSources_${name}*-I-image.fits"), path("wscleanSources_${name}*-sources.txt")

    // label jobs that need a bigger cpu allocation
    label "cpu"
    label "img"

    script:
    """
    unzip *.ms.zip
    export name="\$(basename *.ms)"
    name="\${name%.*}"
    # todo: idg?
    ${params.wsclean} \
        -mgain 0.95 -weight briggs -1.0 -multiscale \
        -name wscleanSources_\${name} \
        -size ${params.img_size} ${params.img_size} \
        -scale ${params.img_scale} \
        -niter ${params.img_niter} \
        -pol i \
        -auto-threshold ${params.img_auto_threshold} \
        -auto-mask ${params.img_auto_mask} \
        -abs-mem ${params.img_mem} \
        -channels-out ${params.img_channels_out} \
        -save-source-list \
        *.ms
    """
}

process imgQA {
    // will need dirty V and clean XX, YY for a single vis file
    storeDir "$params.outdir/${obsid}/img_qa"
    errorStrategy 'ignore'

    input:
    tuple val(obsid), val(name), path("*") // <- "*-image.fits"

    output:
    tuple val(obsid), path("${name}.json")

    stageInMode "copy"
    module "singularity"

    script:
    """
    singularity exec --bind \$PWD:/tmp --writable-tmpfs --pwd /tmp --no-home ${params.mwa_qa_sif} \
        img_qa *.fits --out ${name}.json
    ls -al
    """
}

// todo: process visQa

// ensure calibrated visiblities are present, or apply solutions with hyperdrive

workflow {
    obsids = channel.fromPath(params.obsids_path).splitCsv().flatten()

    // ensure raw files are downloaded
    ensureRaw(obsids)

    // ensure preprocessed files have been generated from raw files
    ensurePrep(ensureRaw.out)

    // QA preprocessed flags
    // - get the obsid and metafits from ensureRaw 
    // - cross with the uvfits and mwaf from ensurePrep
    ensureRaw.out \
        .map(it -> it[0..1]) \
        .join(ensurePrep.out) \
        .map(it -> it[0..2] + [it[4]]) \
        | flagQA

    // create a channel of calibration parameter sets from a csv with columns:
    // - short name that appears in the output filename
    // - args for hyperdrive di-cal, where `-n` = `--num-sources`
    // - args for hyperdrive soln apply args
    cal_param_sets = channel \
        .fromPath(params.cal_params_path) \
        // .splitCsv(header: false, skip: 1, strip: true)
        .collect()

    // calibration for each obsid that passes QA, and each parameter set
    // - get obsid, metafits from ensureRaw cross with uvfits from ensurePrep
    // - combine with calibration parameter sets
    ensureRaw.out \
        .map(it -> it[0..1]) \
        .join(ensurePrep.out)
        .map(it -> it[0..2]) \
        .combine(cal_param_sets) \
        | ensureCal

    // calibration QA
    // - get obsid, metafits from ensureRaw cross with hyp_soln from ensureCal
    ensureRaw.out \
        .map(it -> it[0..1]) \
        .join(ensureCal.out) \
        .map( it -> it[0..2]) \
        .transpose() \
        | view 
        // TODO: | calQA

    // qa images of all preprocessed or calibrated measurement sets
    // - get ms from prep
    // - mix with ms from cal
    // - give each image a name (everything before the .ms.zip)
    ensurePrep.out \
        .map(it -> [it[0], it[2]]) \
        .mix(ensureCal.out.map( it -> [it[0], it[3]] )) \
        .transpose() \
        .map(it -> [\
            it[0], \
            (it[1].getName() =~ /^(.*)\.ms\.zip$/)[0][1], \
            it[1]\
        ]) \
        | wscleanDirty
        // TODO: | (wscleanDirty & wscleanSources)

    // for every group of dirty images, group by channel number, and imgQA for all pols
    wscleanDirty.out \
        .flatMap(it -> \
            it[2] \
                // .collect(path -> [path.getName().split('-')[..-3].join('-'), path]) \
                .collect(path -> [(path.getName() =~ /^(.*)-\w+-\w+.fits$/)[0][1], path]) \
                .groupBy(tup -> tup[0]) \
                .entrySet() \
                .collect(entry -> [it[0], entry.key, entry.value.collect(item -> item[1])]) \
        ) \
        | imgQA

    // what about analysis of di cal residuals?

    // TODO: gather QA results
}