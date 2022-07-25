#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Ensure the raw visibility files are present, or download via ASVO
process ensureRaw {
    // persist results in outdir, process will be skipped if files already present.
    storeDir "$params.outdir/$obsid/raw"
    // allow multiple retries
    maxRetries 5
    // exponential backoff: sleep for 2^attempt hours after each fail
    errorStrategy { sleep(Math.pow(2, task.attempt) * 60*60*1000 as long); return 'retry' }

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
            | ${params.jq} -r '.[]|[.jobId,.files[0].fileUrl//"",.files[0].fileSize//""]|@tsv' \
            | tee ready.tsv
        while read jobid url size; do
            ensure_disk_space \#((2 * \$size)) || exit \$?
            if [ -n "\$url" ]; then
                echo "[obs:$obsid]: Downloading job \$jobid from \$url (\$size bytes)"
                curl \$url | tar -x
                # TODO download and use giant squid to check hash match
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

    input:
    tuple val(obsid), path("${obsid}.metafits"), path("*") // <-preserves names of fits files

    output:
    tuple val(obsid), path("birli_${obsid}.uvfits"), path("${obsid}*.mwaf"), path("birli_prep.log")

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
    echo workflow.workDir=${workflow.workDir} >> birli_prep.log
    """
    // ms:
    // ${params.birli} \
    //     -M "birli_${obsid}.ms" \
    // zip -r "birli_${obsid}.ms.zip" "birli_${obsid}.ms"
}

// QA tasks for flags.
process flagQA {

    input:
    tuple val(obsid), path("${obsid}.metafits"), path("*") // <-preserves names of mwaf files

    output:
    // TODO: change this to whatever paths QA outputs
    tuple val(obsid), path("total_occupancy.csv")

    storeDir "$params.outdir/${obsid}/flag_qa"

    // if this script returns failure, don't progress this obs to calibration, just ignore
    errorStrategy 'ignore'
    // without this, nextflow checks for the output before fs has had time
    // to sync, causing it to think the job has failed.
    afterScript 'ls -al; sleep 5s'

    module 'python/3.9.7'

    script:
    """
    #!/usr/bin/env python
    # workflow.workDir=${workflow.workDir}

    from astropy.io import fits
    from glob import glob
    from time import sleep
    from sys import stderr

    paths = glob("*.mwaf")
    total_occupancy = 0
    num_coarse_chans = len(paths)
    for path in paths:
        with fits.open(path) as hdus:
            # hdus.info()
            flag_data = hdus[1].data['FLAGS']
            path_occupancy = flag_data.sum() / flag_data.size
            print(f"{path}\t{100 * path_occupancy:6.2f}%")
            total_occupancy += path_occupancy
    total_occupancy /= num_coarse_chans
    with open('total_occupancy.csv', 'w') as f:
        f.write(f"{total_occupancy}")
    print(f"--- total ---\t{100 * total_occupancy:6.2f}%")
    """
}

// ensure calibration solutions and vis are present, or calibrate prep with hyperdrive
process ensureCal {
    // TODO: figure out why this keeps failing
    errorStrategy 'ignore'

    storeDir "$params.outdir/${obsid}/cal$params.cal_suffix"

    // label jobs that need a bigger gpu allocation
    label "gpu"

    module "cuda/11.3.1:gcc-rt/9.2.0"
    stageInMode "copy"

    input:
    tuple val(obsid), path("${obsid}.metafits"), path("${obsid}.uvfits"), path("cal_params.csv")

    output:
    tuple val(obsid), path("hyp_soln_${obsid}_*.fits"), path("hyp_${obsid}_*.uvfits"), \
        path("hyp_di-cal_${obsid}_*.log"), path("hyp_apply_${obsid}_*.log")

    script:
    """
    export HYPERDRIVE_CUDA_COMPUTE=${params.cuda_compute}
    # use this control which GPU we're sending the job to
    export CUDA_VISIBLE_DEVICES=0

    set -ex
    ls -al
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
                --outputs "hyp_${obsid}_\${name}.uvfits" \
                > hyp_apply_${obsid}_\${name}.log
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
        for file in "hyp_soln_${obsid}_\${name}.fits" "hyp_${obsid}_\${name}.uvfits"; do
            if [ ! -f \$file ]; then
                echo "Missing file \$file"
                exit 1
            fi
        done
    done
    """
}

// // really simple gate for vis, to prevent "past end of file" errors
// process visGate {
//     input:
//     tuple val(obsid), path("*") // <- .uvfits
//     output:
//     tuple val(obsid), val(obsid)

//     // if this script returns failure, don't progress this vis
//     errorStrategy 'ignore'

//     module 'python/3.9.7'

//     script:
//     """
//     #!/usr/bin/env python
//     # workflow.workDir=${workflow.workDir}

//     from astropy.io import fits
//     from glob import glob
//     from time import sleep
//     from sys import stderr

//     paths = glob("*.uvfits")
//     for path in paths:
//         with fits.open(path) as hdus:
//             print(hdus[1].data[0].shape)
//     """
// }

// QA tasks that can be run on each visibility file.
process visQA {
    input:
    tuple val(obsid), val(name), path("*") // <- *.uvfits

    output:
    tuple val(obsid), path("${name}_vis_metrics.json")

    storeDir "$params.outdir/${obsid}/vis_qa"

    // TODO: figure out why this fails sometimes
    errorStrategy 'ignore'
    // maxRetries 1

    // needed for singularity
    stageInMode "copy"
    module "singularity"

    // without this, nextflow checks for the output before fs has had time
    // to sync, causing it to think the job has failed.
    beforeScript 'sleep 1s; ls -al'
    afterScript 'ls -al; sleep 1s'

    script:
    """
    set -ex
    singularity exec --bind \$PWD:/tmp --writable-tmpfs --pwd /tmp --no-home ${params.mwa_qa_sif} \
        vis_qa *.uvfits --out "${name}_vis_metrics.json"
    """
}

// QA tasks that can be run on each calibration parameter set
process calQA {

    input:
    tuple val(obsid), path("${obsid}.metafits"), val(name), path("*") // <- hyp_soln_*.fits

    output:
    tuple val(obsid), path("${name}_{X,Y}.json") // TODO: , path("${name}_{phases,amps}.png")

    storeDir "$params.outdir/${obsid}/cal_qa"

    // TODO: figure out why this fails sometimes
    errorStrategy 'ignore'
    // maxRetries 1

    stageInMode "copy"
    module "singularity"
    // TODO: figure out why this keeps failing
    errorStrategy 'ignore'

    // without this, nextflow checks for the output before fs has had time
    // to sync, causing it to think the job has failed.
    beforeScript 'sleep 1s; ls -al'
    afterScript 'ls -al; sleep 1s'

    script:
    """
    set -ex
    ls -al
    singularity exec --bind \$PWD:/tmp --writable-tmpfs --pwd /tmp --no-home ${params.mwa_qa_sif} \
        cal_qa *.fits *.metafits X --out "${name}_X.json"
    singularity exec --bind \$PWD:/tmp --writable-tmpfs --pwd /tmp --no-home ${params.mwa_qa_sif} \
        cal_qa *.fits *.metafits Y --out "${name}_Y.json"
    # TODO: recompile hyperdrive with plotting
    # hyperdrive solutions-plot -m "${obsid}.metafits" *.fits
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

process imgQA {
    // will need dirty V and clean XX, YY for a single vis file

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
        img_qa *.fits --out ${name}.json
    """
}

workflow {
    obsids = channel.fromPath(params.obsids_path).splitCsv().flatten()

    // ensure raw files are downloaded
    ensureRaw(obsids)

    // check raw file size / count
    ensureRaw.out \
        .map(it -> String.format( \
            "%s\t%s\t%s\t%s", \
            it[0], \
            it[1].size(), \
            it[2].size(), \
            it[2].collect(path -> path.size()).sum() \
        )) \
        .collectFile(name: "${projectDir}/results/prep_stats.tsv", newLine: true) \
        | view

    // ensure preprocessed files have been generated from raw files
    ensurePrep(ensureRaw.out)

    // quick gate to prevent past end of file
    // ensurePrep.out \
    //     .map(it -> [it[0], it[2]]) \
    //     | visGate

    // flag QA
    // - get the obsid, metafits from ensureRaw
    // - cross with the mwaf from ensurePrep
    ensureRaw.out \
        .map(it -> it[0..1]) \
        .join(ensurePrep.out) \
        .map(it -> it[0..1] + [it[3]]) \
        | flagQA

    // export flagQA results
    flagQA.out \
        .map(it -> String.format( \
            "%s\t%s", \
            it[0], \
            Float.parseFloat(file(it[1]).getText()), \
        )) \
        .collectFile(name: "${projectDir}/results/occupancy.tsv", newLine: true) \
        | view

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
        .join(ensurePrep.out) \
        .map(it -> it[0..2]) \
        // .join( \
        //     visGate.out \
        //     .map(it -> it[0]) \
        // )
        .join( \
            flagQA.out \
            .map(it -> [it[0], Float.parseFloat(file(it[1]).getText())])
            .filter(it -> it[1] < params.flag_occupancy_threshold) \
            .map(it -> [it[0]])
        ) \
        .combine(cal_param_sets) \
        | ensureCal

    // vis QA
    // - get obsid, uvfits from ensurePrep
    // - mix with obsid, uvfits from ensureCal
    // - give each vis a name from basename of uvfits

    // prep seems to be too big without averaging
    // ensurePrep.out \
    //     .map(it -> it[0..1]) \
    //     .mix( \
    //     ) \

    ensureCal.out \
        .map(it -> [it[0], it[2]]) \
        .transpose() \
        .map(it -> [\
            it[0],
            it[1].getBaseName(), \
            it[1] \
        ])\
        | (visQA & psMetrics)

    // calibration QA
    // - get obsid, metafits from ensureRaw cross with hyp_soln from ensureCal
    // - give each calibration a name from basename of uvfits
    ensureRaw.out \
        .map(it -> it[0..1]) \
        .join(ensureCal.out) \
        .map( it -> it[0..2]) \
        .transpose() \
        .map(it -> [\
            it[0],
            it[1],
            it[2].getBaseName(), \
            it[2] \
        ])\
        | calQA

    // img QA
    // - get obsid, metafits cross with uvfits from ensureCal
    // - give each image a name from basename of uvfits
    // - run wsclean
    ensureRaw.out \
        .map(it -> it[0..1]) \
        .join( \
            ensureCal.out \
            .map( it -> [it[0], it[2]] ) \
        )
        .transpose() \
        .map(it -> [\
            it[0], \
            it[2].getBaseName(), \
            it[1],
            it[2]\
        ]) \
        | wscleanDirty
    //     // TODO: | (wscleanDirty & wscleanSources)

    // for every group of dirty images, group by channel number, and imgQA for all pols
    wscleanDirty.out \
        .flatMap(it -> \
            it[2] \
                .collect(path -> [(path.getName() =~ /^(.*)-\w+-dirty.fits$/)[0][1], path]) \
                .groupBy(tup -> tup[0]) \
                .entrySet() \
                .collect(entry -> [it[0], entry.key, entry.value.collect(item -> item[1])]) \
        ) \
        | imgQA

    psMetrics.out \
        .map(it -> it[1].getText()) \
        .collectFile(name: "${projectDir}/results/ps_metrics.dat") \
        | view

    // what about analysis of di cal residuals?

    // TODO: gather QA results
}