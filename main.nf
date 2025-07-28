#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include {
    birli_argstr_suffix;
    calqa_pass;
    cmt_imgqa_pass_sub;
    cmt_imgqa_pass;
    cmt_ps_metrics_pass_sub;
    cmt_ps_metrics_pass;
    coerceList;
    contigRanges;
    decomposeImg;
    deepcopy;
    displayInts;
    displayRange;
    exitCodes;
    firstFail;
    get_seconds;
    getFailReason;
    groovy2bashAssocArray;
    groupMeta;
    hyp_apply_name;
    is_multichannel;
    is_multiinterval;
    isNaN;
    mapMerge;
    obsids_file;
    openWithDelay;
    parseCsv2;
    parseFloatOrNaN;
    parseJson;
    prepqa_pass;
    results_dir;
    wrap_angle;
    wscleanDConvParams;
    wscleanParams;
    wsSummarize;
    getFreqResHz;
    getFreqReskHz;
} from './modules/utils.nf'


// default entrypoint: get preprocessed vis from asvo and run qa
workflow {
    // get obsids from csv
    obsCSV = channel.of(obsids_file())
        .splitCsv()
        .filter { line -> !line[0].startsWith('#') }
        .map { line ->
            def (obsid, cluster) = line
            def meta = [:]
            if (cluster != null) {
                meta.cluster = cluster
            }
            [obsid, deepcopy(meta)]
        }

    obsids = obsCSV.map { obsid, _m -> obsid }

    // analyse obsids with web services
    obsids | ws

    ws.out.obsMetafits
        .map { obsid, _metafits -> obsid }
        .collectFile(
            name: "ws_obs_pass.csv", newLine: true, sort: true,
            storeDir: "${results_dir()}"
        )
        | view { [it, it.readLines().size()] }

    // download preprocessed, unless noprep is set
    if (params.noprep) {
        prep(channel.empty())
    } else {
        prep(ws.out.obsMeta.join(ws.out.obsMetafits))
    }

    // channel of obsids that pass the flag gate
    prep.out.subobsMetaVisPass
        .map { obsid, meta, vis -> [obsid, meta.subobs?:'', meta.name?:'', vis].join('\t') }
        .collectFile(
            name: "prep_subobs_name_pass.csv", newLine: true, sort: true,
            storeDir: "${results_dir()}"
        )
        | view { [it, it.readLines().size()] }

    qaPrep( prep.out.subobsMetaVisPass, ws.out.obsMetafits )
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

    all_fail_codes = ws.out.fail_codes
        .join(prep.out.fail_codes, remainder: true)
        .join(qaPrep.out.fail_codes, remainder: true)
        .map { obsid, ws_code, prep_code, prepqa_code ->
            def codes = [ws_code, prep_code, prepqa_code].findAll { it != null }
            def aFailCode = coerceList(firstFail([codes]))[0]
            [obsid, aFailCode]
        }

    // do lst counts, scatterplot
    all_fail_codes
        .filter { _o, fail_code -> fail_code == getFailReason(0x00) }
        .join(ws.out.obsMeta)
        .map { obsid, _fail_code, meta ->
            [obsid, meta.lst.round().intValue()]
        }
        .groupTuple(by: 1)
        .map { obsids_, lst ->
            [
                String.format("%+03d", lst),
                obsids_.size(),
            ].join("\t")
        }
        .collectFile(
            name: "lst_counts_all.tsv", newLine: true, sort: true,
            seed: ([ "LST", "COUNT" ]).join("\t"),
            storeDir: "${results_dir()}${params.img_suffix}${params.cal_suffix}"
        )
        .view { it.readLines().size() }
        .map { tsv ->
            def meta = [
                name: "lst_counts_all", title: "lst counts after qa",
                x: "LST", y: "COUNT"
            ]
            [meta, tsv]
        }
        | tsvScatterPlot

    // do ewp counts
    all_fail_codes
        .filter { _o, fail_code -> fail_code == getFailReason(0x00) }
        .join(ws.out.obsMeta)
        .map { obsid, _fail_code, meta ->
            [obsid, meta.ew_pointing]
        }
        .groupTuple(by: 1)
        .map { obsids_, ewp ->
            [
                ewp == null ? "" : String.format("%+1d", ewp),
                obsids_.size(),
            ].join("\t")
        }
        .collectFile(
            name: "ewp_counts_all.tsv", newLine: true, sort: true,
            seed: ([ "LST", "COUNT" ]).join("\t"),
            storeDir: "${results_dir()}${params.img_suffix}${params.cal_suffix}"
        )

    all_fail_codes.groupTuple(by: 1)
        .map { obsids_, fail_code ->
            [
                fail_code,
                obsids_.size(),
            ].join("\t")
        }
        .collectFile(
            name: "fail_counts_all.tsv", newLine: true, sort: true,
            seed: ([ "FAIL CODE", "COUNT" ]).join("\t"),
            storeDir: "${results_dir()}${params.img_suffix}${params.cal_suffix}"
        )
}

// download observation metadata from webservices in json format
process wsMeta {
    // persist results in outdir, process will be skipped if files already present.
    storeDir "${params.outdir}/meta"
    // tag to identify job in squeue and nf logs
    tag "${obsid}"
    time {5.minute * task.attempt * params.scratchFactor}

    // allow multiple retries
    maxRetries 2
    errorStrategy {
        if (task.attempt > 2) {
            return 'ignore'
        }
        return (task.exitStatus == 8 ? 'retry' : 'ignore')
    }

    input:
    val(obsid)
    output:
    tuple val(obsid), path(wsmeta), path(wsfiles)

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
    storeDir "${params.outdir}/meta"
    label "tap"
    tag "${obsid}"

    input:
    val(obsid)
    output:
    tuple val(obsid), path(tapmeta)

    script:
    tapmeta = "${obsid}_tapmeta.json"
    template "tapmeta.py"
}

// download observation metadata from webservices in metafits format
process wsMetafits {
    storeDir "${params.outdir}/${obsid}/raw"
    label "rate_limit"
    tag "${obsid}"

    input:
    val(obsid)
    output:
    tuple val(obsid), path(metafits)

    script:
    metafits = "${obsid}.metafits"
    """
    #!/bin/bash -eux
    ${params.proxy_prelude} # ensure proxy is set if needed
    wget -O "${metafits}" "http://ws.mwatelescope.org/metadata/fits?obs_id=${obsid}&include_ppds=${params.metafits_incl_ppds}"
    """
}

process wsSkyMap {
    storeDir "${params.outdir}/${obsid}/meta"
    label "rate_limit"
    tag "${obsid}"

    input:
    val(obsid)
    output:
    tuple val(obsid), path(skymap)

    script:
    skymap = "${obsid}_skymap.png"
    """
    #!/bin/bash -eux
    ${params.proxy_prelude} # ensure proxy is set if needed
    wget -O "${skymap}" "http://ws.mwatelescope.org/observation/skymap/?obs_id=${obsid}"
    """
}

process wsPPDs {
    storeDir "${params.outdir}/${obsid}/meta"
    tag "${obsid}"
    label "rate_limit"

    input:
    val(obsid)
    output:
    tuple val(obsid), path(ppds)

    when: !params.noppds

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
    storeDir "${params.outdir}/${obsid}/meta"
    tag "${obsid}"

    errorStrategy 'ignore'

    label "python"

    input:
        tuple val(obsid), path(metafits)
    output:
        tuple val(obsid), path(json), path(tsv)

    script:
    // metrics = "${obsid}_occupancy.json"
    json = "${obsid}_meta.json"
    tsv = "${obsid}_inputs.tsv"
    txt = "${obsid}_inputs.txt"
    template "metajson.py"
}

// temporary ASVO workaround, raw files from mwacache
process cacheBoxRaw {
    storeDir "${params.outdir}/${obsid}/raw"
    tag "${obsid}"

    input:
    tuple val(obsid), val(meta)
    output:
    tuple val(obsid), val(meta), path(metafits), path(raw)

    script:
    metafits = "${obsid}_metafits*.fits"
    raw = "${obsid}_2*.fits"
    """
    #!/bin/bash -eux
    for i in {01..08}; do
        for j in {1..3}; do
            rsync -auz --info=progress2 \
                mwacache\${i}:/volume\${j}/incoming/${obsid}'*'.fits . \
                || true
        done
    done
    ls -al
    """
}

process birliPrepUV {
    storeDir "${params.outdir}/${obsid}/prep"
    label "birli"
    label "cpu_half"
    time { params.scratchFactor * 2.hour * Math.pow(task.attempt, 2) }
    disk { 50.GB * Math.pow(task.attempt, 2) }
    stageInMode "symlink"
    memory { 200.GB * task.attempt }
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'ignore' }
    tag "${obsid}${birli_argstr_suffix()[2]?:''}"

    input:
    tuple val(obsid), val(meta_), path(metafits), path(raw)
    output:
    tuple val(obsid), val(meta), path(uvfits), path(metrics)
        // , path("${obsid}${spw}*.mwaf"), path("birli_prep.log")

    when: !params.noprep

    script:
    freq_res_khz = getFreqReskHz(meta_)
    if (freq_res_khz != null && params.prep_freq_res_khz != null && freq_res_khz > params.prep_freq_res_khz) {
        throw new Exception("error: target freq_res ${params.prep_freq_res_khz} < obs freq res ${freq_res_khz}. \nmeta=${meta}")
    }
    meta = deepcopy(meta_)
    prefix = "birli_"
    def (_argstr_asvo, argstr_cli, suffix) = birli_argstr_suffix()
    meta.prep_freq_res = params.prep_freq_res_khz?:freq_res_khz
    meta.prep_time_res = params.prep_time_res_s?:meta.time_res
    def coarse_chans = (meta.coarse_chans?:[]).sort(false)
    def not_contiguous = (coarse_chans.size > 1 && coarse_chans.size < coarse_chans[-1] - coarse_chans[0])
    def channel_glob = ""
    if (not_contiguous) {
        channel_glob = "_ch{??,???,??-??,??-???,???-???}"
    } else if (meta['subobs']) {
        subobs = meta['subobs']
        if (subobs =~ /ch[\d-]+/) {
            suffix += "_${subobs}"
        }
    }
    meta['channel_glob'] = channel_glob
    meta['birli_suffix'] = suffix
    uvfits = ''+"${prefix}${obsid}${suffix}${channel_glob}.uvfits"
    metrics = ''+"${prefix}${obsid}${suffix}${channel_glob}_metrics.fits"
    """
    set -eux
    ${params.birli} \
        ${argstr_cli} \
        -u "${prefix}${obsid}${suffix}.uvfits" \
        -m "${metafits}" \
        --metrics-out "${prefix}${obsid}${suffix}_metrics.fits" \
        -- ${raw}
    """
}

process hypPrepVisConvert {
    storeDir "${params.outdir}/${obsid}/prep"
    label "hyperdrive_cpu"
    time { params.scratchFactor * 1.hour * Math.pow(task.attempt, 2) }
    disk { 50.GB * Math.pow(task.attempt, 2) }
    stageInMode "symlink"
    memory { 200.GB * task.attempt }
    // errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'ignore' }
    errorStrategy 'terminate'
    tag "${obsid}${suffix}"

    input:
    tuple val(obsid), val(meta_), path(metafits), path(uvfits)
    output:
    tuple val(obsid), val(meta), path(vis)

    when: (params.prep_export_time_res_s != null || params.prep_export_freq_res_khz != null)

    script:
    meta = deepcopy(meta_)
    prefix = "birli_"
    suffix = '' + (params.prep_suffix ?: '')
    args = [:]
    if (meta['subobs'] != null) {
        suffix += "_${subobs}"
    }
    if (params.prep_export_time_res_s != null) {
        if (meta.prep_time_res != null && meta.prep_time_res > params.prep_export_time_res_s) {
            throw new Exception("error: target time_res ${params.prep_export_time_res_s} < prep time res ${meta.prep_time_res}. \nmeta=${meta}")
        }
        suffix += "_${params.prep_export_time_res_s}s"
        args['time-average'] = params.prep_export_time_res_s
    }
    if (params.prep_export_freq_res_khz != null) {
        if (meta.prep_freq_res != null && meta.prep_freq_res > params.prep_export_freq_res_khz) {
            throw new Exception("error: target freq_res ${params.prep_export_freq_res_khz} < prep freq res ${meta.prep_freq_res}. \nmeta=${meta}")
        }
        suffix += "_${params.prep_export_freq_res_khz}kHz"
        args['freq-average'] = params.prep_export_freq_res_khz
    }
    if (params.ssins_apply) {
        suffix += ".ssins"
    }
    prepFlags = meta.prepFlags?:[]
    fineChanFlags = meta.fineChanFlags?:[]
    if (prepFlags.size() > 0) {
        args["tile-flags"] = "${prepFlags.join(' ')}"
    }
    if (prepFlags.size() > 0) {
        args["tile-flags"] = "${prepFlags.join(' ')}"
    }
    if (fineChanFlags.size() > 0) {
        args["fine-chan-flags-per-coarse-chan"] = "${fineChanFlags.join(' ')}"
    }
    argstr = args
        .collect { k, v ->
            if (v == null) {
                ['--' + k]
            }
            else {
                ['--' + k] + v
            }
        }
        .flatten()
        .join(' ')
    vis = ''+"${prefix}${obsid}${suffix}.${params.prep_export_ext}"
    """
    set -eux
    ${params.hyperdrive_cpu} vis-convert \
        ${argstr} \
        -d $metafits $uvfits \
        -o "${vis}"
    """
}

process demo03_mwalib {
    stageInMode "symlink"
    storeDir "${params.outdir}/${obsid}/raw_qa"
    label "mwa_demo"
    tag "${metafits.baseName}"

    input:
    tuple val(obsid), val(meta), path(metafits)
    output:
    tuple val(obsid), val(meta), path("${metafits.baseName}-antennas.tsv"), path("${metafits.baseName}-channels.tsv")

    when: !params.nodemo

    script:
    """
    ${params.demo_prelude?:''}
    /demo/03_mwalib.py ${metafits}
    """
}

process demo04_ssins {
    stageInMode "copy"
    storeDir "${params.outdir}/${obsid}/${qa}_qa"
    label "ssins"
    label "mwa_demo"
    label "cpu_quarter"
    memory { MemoryUnit.of(7 * vis_bytes) } // TODO: make ssins memory footprint suck less
    time 6.hour
    tag "${base}${meta.plot_base?:''}"

    input:
    tuple val(obsid), val(meta), path(metafits), path(vis_)
    output:
    tuple val(obsid), val(meta), path("${base}${meta.plot_base?:''}.png"), path("${base}${meta.mask_base?:''}*_SSINS_mask.h5", optional: true)

    when: !params.nodemo

    script:
    vis = coerceList(vis_)
    vis_bytes = vis.collect { it.size() }.sum()
    firstVis = vis[0]
    if (firstVis.extension == "uvfits") {
        base = firstVis.baseName
        qa = "prep"
    } else {
        base = obsid
        qa = "raw"
    }
    """
    ${params.demo_prelude?:''}
    /demo/04_ssins.py ${meta.argstr?:''} ${metafits} ${vis.join(' ')}
    """
}

process demo10_phase {
    stageInMode "copy"
    storeDir "${params.outdir}/${obsid}/bf"
    label "mwa_demo"
    label "mem_super"
    tag "${base}${meta.plot_base?:''}"
    time 8.hours

    input:
    tuple val(obsid), val(meta), path(vis)
    output:
    tuple val(obsid), val(meta), path("${base}.${meta.plot_base?:''}*.fits")

    when: !(params.nodemo || params.noimg || params.nodemoallsky)

    script:
    base = coerceList(vis)[0].baseName
    suffix = ''
    argstr = ''
    if (meta.phase_centre != null) {
        suffix += ".${meta.phase_centre}"
    }
    if (meta.avg_time != null) {
        suffix += ".${meta.avg_time}"
    }
    if (meta.avg_time != null) {
        suffix += ".${meta.avg_time}"
    }
    """
    ${params.demo_prelude?:''}
    /demo/10_phase.py ${meta.argstr?:''} ${coerceList(vis).join(' ')}
    """
}

process demo11_allsky {
    stageInMode "copy"
    storeDir "${params.outdir}/${obsid}/img"
    label "mwa_demo"
    label "mem_super"
    tag "${base}${meta.plot_base?:''}"
    time 8.hours

    input:
    tuple val(obsid), val(meta), path(vis)
    output:
    tuple val(obsid), val(meta), path("${base}${meta.plot_base?:''}*.fits")

    when: !(params.nodemo || params.noimg || params.nodemoallsky)

    script:
    base = coerceList(vis)[0].baseName
    """
    ${params.demo_prelude?:''}
    /demo/11_allsky.py ${meta.argstr?:''} ${coerceList(vis).join(' ')}
    """
}


workflow extRaw {
    channel.of(obsids_file())
        .splitCsv()
        .filter { line -> !line[0].startsWith('#') }
        .map { line ->
            def (obsid, _comment) = line
            obsid
        }
        .unique()
        .map { obsid_ ->
            def obsid = coerceList(obsid_)[0]
            def meta = [obsid: obsid]
            // def raw = file("${params.outdir}/../raw/${obsid}_2?????????????_ch???_???.fits")
            def raw = file("${params.outdir}/${obsid}/raw/${obsid}_2*.fits")
            [ obsid, deepcopy(meta), raw ]
        }
        .filter { _o, _m, raw_ ->
            def raw = coerceList(raw_)
            raw.size() && raw.every{
                if (!it.exists()) { print("raw does not exist: ${it}") }
                it.exists()
            }
        }
        .tap { obsMetaRaw }
        .map { obsid, _m, _raw -> obsid }
        .tap { obsids }

    obsids | ws

    obsMetaRaw.join(ws.out.obsMeta).join(ws.out.obsMetafits)
        .map { obsid, meta, vis, wsMeta, metafits ->
            [ obsid,  mapMerge(meta, wsMeta), metafits, vis ]
        }
        | birliPrepUV

    birliPrepUV.out.flatMap { obsid, meta, uvfits_, _metrics ->
            def uvfits = coerceList(uvfits_)
            if (uvfits.size > 1) {
                uvfits.collect { f ->
                    def newMeta = [:]
                    def last_token = '' + f.baseName.split('_')[-1]
                    if (last_token =~ /ch[\d-]+/) {
                        newMeta.subobs = last_token
                    }
                    [obsid, mapMerge(meta, newMeta), f]
                }
            } else {
                [[obsid, meta, uvfits_]]
            }
        }
        .tap { subobsVis }
        | uvMeta

    ws.out.obsMeta.cross(uvMeta.out) { it[0] }
        .map { obsMeta_, uvMeta_ ->
            def (obsid, wsMeta) = obsMeta_
            def (_o, meta_, uvJson) = uvMeta_
            def uvmeta = parseJson(uvJson)
            def newMeta = [
                lowfreq: uvmeta.freqs[0],
                freq_res_hz: (uvmeta.freqs[1] - uvmeta.freqs[0]),
                nchans: (uvmeta.freqs?:[]).size(),
                ntimes: (uvmeta.times?:[]).size(),
            ]
            ['config', 'eorband', 'num_ants', 'total_weight'].each { key ->
                if (uvmeta[key] != null) {
                    newMeta[key] = uvmeta[key]
                }
            }
            if (newMeta.ntimes > 0) {
                newMeta.lst = Math.toDegrees(uvmeta.times[0].lst_rad)
            }
            [obsid, mapMerge(mapMerge(meta_, wsMeta), newMeta)]
        }
        .tap { subobsMeta }
        .cross(subobsVis) { def (obsid, meta) = it; [obsid, meta.subobs?:''] }
        .map { subobsMeta_, subobsVis_ ->
            def (obsid, meta) = subobsMeta_
            def (_o, _m, uvfits) = subobsVis_
            [obsid, meta, uvfits]
        }
        .tap { subobsMetaVis }

    flag(subobsMetaVis, ws.out.obsMetafits)

    // for all obs that pass flag:
    flag.out.subobsFlagmetaPass.map { obsid, meta, flagMeta ->
            // update meta with flagMeta
            [[obsid, meta.subobs?:''], mapMerge(meta, flagMeta)]
        }
        // join with vis from subobsMetaVis
        .join(subobsMetaVis.map { obsid, meta, uvfits ->
            [[obsid, meta.subobs?:''], uvfits]
        })
        .map { obsSubobs, meta, uvfits ->
            def (obsid, _subobs) = obsSubobs
            [obsid, meta, uvfits]
        }
        | ssinsQA

    qaPrep( ssinsQA.out.subobsMetaVisSSINs, ws.out.obsMetafits )

    ws.out.frame.mix(flag.out.frame)
        .mix(ssinsQA.out.frame)
        .mix(qaPrep.out.frame)
        .map { n, l -> [n, l as ArrayList] }
        | makeVideos
}

workflow asvoRawFlow {
    channel.of(obsids_file())
        .splitCsv()
        .filter { line -> !line[0].startsWith('#') }
        .map { line ->
            def (obsid, _comment) = line
            obsid
        }
        .unique()
        | asvoRaw

    asvoRaw.out
        .filter { _o, _metafits, raw_ ->
            def raw = coerceList(raw_)
            raw.size() && raw.every{
                if (!it.exists()) { print("raw does not exist: ${it}") }
                it.exists()
            }
        }
        .map { obsid, _m, raw ->
            def meta = [:]
            [obsid, meta, raw]
        }
        .tap { obsMetaRaw }
        .map { obsid, _m, _raw -> obsid }
        .tap { obsids }

// todo combine with extRaw
// }
// workflow raw {
    // take:
    //     obsMetaRaw
    // main:

    obsids | ws

    obsMetaRaw.join(ws.out.obsMeta).join(ws.out.obsMetafits)
        .map { obsid, meta, vis_, wsMeta, metafits ->
            [ obsid, mapMerge(meta, wsMeta), metafits, coerceList(vis_) ]
        }
        .tap { obsMetaMetafitsRaw }

    if (params.prepByCh) {
        rawByCh = obsMetaMetafitsRaw.join(ws.out.mwalibMeta)
            .flatMap { obsid, _m, _metafits, vis, mwalibMeta ->
                vis.collect { f ->
                    def channels = (mwalibMeta.channels?:[]).withIndex().collect { chan, idx ->
                            chan['idx'] = idx;
                            chan
                        }
                    def (obsid_, _datestamp, chan, _batch) = f.baseName.split('_')
                    if ((obsid_ as Integer) != (obsid as Integer)) {
                        throw new Exception("obsid mismatch in ${f}: ${obsid_} != ${obsid}")
                    }
                    def gpuboxMatch = chan =~ /gpubox([\d]+)/
                    def chMatch = chan =~ /ch([\d]+)/
                    def channelInfo = null
                    if (gpuboxMatch) {
                        channelInfo = channels.findAll { it['gpubox_number'] as Integer == gpuboxMatch[0][1] as Integer }[0]
                    } else if (chMatch) {
                        channelInfo = channels.findAll { it['rec_chan_number'] as Integer == chMatch[0][1] as Integer }[0]
                    } else {
                        throw new Exception("unknown channel in ${f}")
                    }
                    [obsid, channelInfo.rec_chan_number, channelInfo.idx, f]
                }
            }
            .groupTuple(by: 0..2)

        obsWsmetaVis = obsMetaMetafitsRaw.cross(rawByCh).map { obsMetaMetafitsRaw_, rawByCh_ ->
                def (obsid, meta, metafits, _raw) = obsMetaMetafitsRaw_
                def (_o, chan, chanIdx, vis) = rawByCh_
                def newMeta = [subobs:"ch${chan}", coarse_chans:[chan], birli_chan_ranges:["${chanIdx}"]]
                [obsid, mapMerge(meta, newMeta), metafits, vis]
            }
            .tap { obsMetafitsRawByCh }
            | birliPrepUV
    } else {
        obsMetafitsRawByCh = obsMetaMetafitsRaw
        obsWsmetaVis = obsMetaMetafitsRaw | birliPrepUV
    }

    obsMetafitsRawByCh.flatMap { obsid, meta, metafits, vis ->
            [
                ['--autos',                        '.diff.auto',  '.spectrum'],
                ['--no-diff --autos',              '.auto',       '.spectrum'],
                ['--crosses',                      '.diff.cross', '.spectrum'],
                ['--crosses --no-diff',            '.cross',      '.spectrum'],
                // ['--sigchain --autos',             '.diff.auto',  '.sigchain'],
                // ['--sigchain --no-diff --autos',   '.auto',       '.sigchain'],
                // ['--sigchain --crosses',           '.diff.cross', '.sigchain'],
                // ['--sigchain --no-diff --crosses', '.cross',      '.sigchain'],
            ].collect { argstr, plot_prefix, plot_suffix ->
                def newMeta = [argstr:argstr, plot_base:"${plot_prefix}${plot_suffix}"]
                if (meta.subobs != null) {
                    newMeta.argstr += " --suffix=${meta.subobs}"
                    newMeta.mask_base = "${plot_prefix}${meta.subobs}"
                    newMeta.plot_base = "${plot_prefix}${meta.subobs}${plot_suffix}"
                }
                [obsid, mapMerge(meta, newMeta), metafits, vis]
            }
        }
        .tap { obsMetaMetafitsRawSsins }

    birliPrepUV.out.flatMap { obsid, meta, uvfits_, _metrics ->
            def uvfits = coerceList(uvfits_)
            if (uvfits.size > 1) {
                uvfits.collect { f ->
                    def newMeta = [:]
                    def last_token = '' + f.baseName.split('_')[-1]
                    if (last_token =~ /ch[\d-]+/) {
                        newMeta.subobs = last_token
                    }
                    [obsid, mapMerge(meta, newMeta), f]
                }
            } else {
                [[obsid, meta, uvfits_]]
            }
        }
        .tap { subobsVis }
        | uvMeta


    subobsVis.join(ws.out.obsMetafits)
        .flatMap { obsid, meta, vis, metafits ->
            [
                // ['--autos',                        '.diff.auto',  '.spectrum'],
                // ['--autos --no-diff',              '.auto',       '.spectrum'],
                ['--crosses',                      '.diff.cross', '.spectrum'],
                // ['--crosses --no-diff',            '.cross',      '.spectrum'],
                // ['--sigchain --autos',             '.diff.auto',  '.sigchain'],
                // ['--sigchain --autos --no-diff',   '.auto',       '.sigchain'],
                // ['--sigchain --crosses',           '.diff.cross', '.sigchain'],
                // ['--sigchain --crosses --no-diff', '.cross',      '.sigchain'],
                // ['--flags --autos --no-diff',      '.auto',       '.flags'],
                // ['--flags --crosses --no-diff',    '.cross',      '.flags'],
            ].collect { argstr, plot_prefix, plot_suffix ->
                def newMeta = [argstr:argstr, plot_base:"${plot_prefix}${plot_suffix}"]
                if (meta.subobs != null) {
                    newMeta.argstr += " --suffix=${meta.subobs}"
                    newMeta.mask_base = "${plot_prefix}${meta.subobs}"
                    newMeta.plot_base = "${plot_prefix}${meta.subobs}${plot_suffix}"
                }
                [obsid, mapMerge(meta, newMeta), metafits, vis]
            }
        }
        .tap { obsMetaMetafitsPrepSsins }

    obsMetaMetafitsRawSsins.mix(obsMetaMetafitsPrepSsins)
        | demo04_ssins

    ws.out.obsMeta.cross(uvMeta.out) { it[0] }
        .map { obsMeta_, uvMeta_ ->
            def (obsid, wsMeta) = obsMeta_
            def (_o, meta_, uvJson) = uvMeta_
            def uvmeta = parseJson(uvJson)
            def newMeta = [
                lowfreq: uvmeta.freqs[0],
                freq_res: (uvmeta.freqs[1] - uvmeta.freqs[0]),
                nchans: (uvmeta.freqs?:[]).size(),
                ntimes: (uvmeta.times?:[]).size(),
            ]
            ['config', 'eorband', 'num_ants', 'total_weight'].each { key ->
                if (uvmeta[key] != null) {
                    newMeta[key] = uvmeta[key]
                }
            }
            if (newMeta.ntimes > 0) {
                newMeta.lst = Math.toDegrees(uvmeta.times[0].lst_rad)
            }
            [obsid, mapMerge(mapMerge(meta_, wsMeta), newMeta)]
        }
        .tap { subobsMeta }
        .cross(subobsVis) { def (obsid, meta) = it; [obsid, meta.subobs?:''] }
        .map { subobsMeta_, subobsVis_ ->
            def (obsid, meta) = subobsMeta_
            def (_o, _m, uvfits) = subobsVis_
            [obsid, meta, uvfits]
        }
        .tap { subobsMetaVis }

    flag(subobsMetaVis, ws.out.obsMetafits)

    // for all obs that pass flag:
    flag.out.subobsFlagmetaPass.map { obsid, meta, flagMeta ->
            // update meta with flagMeta
            [[obsid, meta.subobs?:''], mapMerge(meta, flagMeta)]
        }
        // join with vis from subobsMetaVis
        .join(subobsMetaVis.map { obsid, meta, uvfits ->
            [[obsid, meta.subobs?:''], uvfits]
        })
        .map { obsSubobs, meta, uvfits ->
            def (obsid, _subobs) = obsSubobs
            [obsid, meta, uvfits]
        }
        | ssinsQA

    qaPrep( ssinsQA.out.subobsMetaVisSSINs, ws.out.obsMetafits )

    ws.out.frame.mix(flag.out.frame)
        .mix(ssinsQA.out.frame)
        .mix(qaPrep.out.frame)
        .map { n, l -> [n, l as ArrayList] }
        | makeVideos
}

// workflow extCache {
//     def name = params.visName ?: "ssins"
//     channel.of(obsids_file())
//         .splitCsv()
//         .filter { line -> !line[0].startsWith('#') }
//         .map { line ->
//             def (obsid, _comment) = line
//             obsid
//         }
//         | ws
//     ws.out.obsMetafits.map { obsid, _metafits ->
//             def meta = [:]
//             [ obsid, deepcopy(meta) ]
//         }
//         | cacheBoxRaw
//         .tap { obsMetaMetafitsRaw }
//         | birliPrepUV
//     obsMetaMetafitsRaw | demo04_ssins
//     birliPrepUV.out.flatMap { obsid, meta, uvfits_ ->
//             def uvfits = coerceList(uvfits_)
//             if (uvfits.size > 1) {
//                 uvfits.collect { f ->
//                     def newMeta = [:]
//                     def last_token = '' + f.baseName.split('_')[-1]
//                     if (last_token =~ /ch[\d-]+/) {
//                         newMeta.subobs = last_token
//                     }
//                     [obsid, mapMerge(meta, newMeta), f]
//                 }
//             } else {
//                 [[obsid, meta, metafits, uvfits_]]
//             }
//         }
//         .tap { subobsVis }
//         | uvMeta
//     subobsMetaVis = ws.out.obsMeta.cross(uvMeta.out) { it[0] }
//         .map { obsMeta_, uvMeta_ ->
//             def (obsid, wsMeta) = obsMeta_
//             def (_, meta_, uvJson) = uvMeta_
//             def uvmeta = parseJson(uvJson)
//             def newMeta = [
//                 lowfreq: uvmeta.freqs[0],
//                 freq_res: (uvmeta.freqs[1] - uvmeta.freqs[0]),
//                 nchans: (uvmeta.freqs?:[]).size(),
//                 ntimes: (uvmeta.times?:[]).size(),
//             ]
//             ['config', 'eorband', 'num_ants', 'total_weight'].each { key ->
//                 if (uvmeta[key] != null) {
//                     newMeta[key] = uvmeta[key]
//                 }
//             }
//             if (newMeta.ntimes > 0) {
//                 newMeta.lst = Math.toDegrees(uvmeta.times[0].lst_rad)
//             }
//             [obsid, mapMerge(mapMerge(meta_, wsMeta), newMeta)]
//         }
//         .tap { subobsMeta }
//         .cross(subobsVis) { def (obsid, meta) = it; [obsid, meta.subobs?:''] }
//         .map { subobsMeta_, subobsVis_ ->
//             def (obsid, meta) = subobsMeta_
//             def (_, __, uvfits) = subobsVis_
//             [obsid, meta, uvfits]
//         }
//     qaPrep( subobsMetaVis, ws.out.obsMetafits )
//     qaPrep.out.frame | makeVideos
// }

// Ensure the raw visibility files are present, or download via ASVO
process asvoRaw {
    storeDir "${params.outdir}/${obsid}/raw"
    tag "${obsid}"
    time { 1.hour * Math.pow(task.attempt, 4) * params.scratchFactor }
    disk { 100.GB * Math.pow(task.attempt, 4) }
    memory { 60.GB * Math.pow(task.attempt, 4) }
    label "giant_squid"
    label "rate_limit"
    maxRetries 2
    errorStrategy {
        return 'ignore'
        // TODO: copy from asvoPrep, exponential backoff: sleep for 2^attempt hours after each fail
        // failure_reason = [
        //     5: "I/O error or hash mitch",
        //     28: "No space left on device",
        // ][task.exitStatus]
        // if (failure_reason) {
        //     println "task ${task.hash} failed with code ${task.exitStatus}: ${failure_reason}"
        //     return 'ignore'
        // }
        // retry_reason = [
        //     1: "general or permission",
        //     11: "Resource temporarily unavailable",
        //     75: "Temporary failure, try again"
        // ][task.exitStatus] ?: "unknown"
        // wait_hours = Math.pow(2, task.attempt)
        // println "sleeping for ${wait_hours} hours and retrying task ${task.hash}, which failed with code ${task.exitStatus}: ${retry_reason}"
        // sleep(wait_hours * 60*60*1000 as long)
        // return 'retry'
    }

    input:
    val obsid
    output:
    tuple val(obsid), path("${obsid}.metafits"), path("${obsid}_2*.fits")

    script:
    """
    # echo commands, exit on any failures
    set -eux
    export MWA_ASVO_API_KEY="${params.asvo_api_key}"
    ${params.proxy_prelude} # ensure proxy is set if needed
    # submit a job to ASVO, suppress failure if a job already exists.
    ${params.giant_squid} submit-vis --delivery "acacia" $obsid || true
    # extract id and state from pending download vis jobs
    ${params.giant_squid} list -j --types download_visibilities --states queued,error -- $obsid \
        | ${params.jq} -r '.[]|[.jobId,.jobState]|@tsv' \
        | tee pending.tsv
    # extract id url size hash from any ready download vis jobs for this obsid
    ${params.giant_squid} list -j --types download_visibilities --states ready -- $obsid \
        | tee /dev/stderr \
        | ${params.jq} -r '.[]|[.jobId,.files[0].fileUrl//"",.files[0].fileSize//"",.files[0].fileHash//""]|@tsv' \
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
        exit 0 # success
    fi
    echo "no ready jobs"
    exit 75 # temporary
    """
}

// asvoRaw was here: https://github.com/MWATelescope/MWAEoR-Pipeline/commit/e49ac2765c0fcb6bada9d3245d59a65b702f5707

// Download preprocessed files from asvo
//
// If the preprocessed files are not present, the `asvoPrep` process schedules
// a conversion job on ASVO using [Giant Squid](github.com/mwaTelescope/giant-squid).
// An exponential backoff is used to wait until the conversion job is ready to
// download from [Acacia](https://pawsey.org.au/systems/acacia/), then the
// archive is downloaded with `wget`, and hash-validated as it is inflated with
// `tar`.
process asvoPrep {
    storeDir "${params.outdir}/${obsid}/prep"
    tag "${obsid}"
    label "giant_squid"
    label "rate_limit"
    label "nvme"
    time { 1.hour * Math.pow(task.attempt, 4) * params.scratchFactor }
    disk { 50.GB * Math.pow(task.attempt, 4) }
    memory { 50.GB * Math.pow(task.attempt, 4) }
    maxRetries 2
    errorStrategy {
        def reason = exitCodes()[task.exitStatus] ?: "unknown"
        def result = 'ignore'
        if (task.exitStatus == 75) { // may take up to a day for asvo to finish
            // result = 'retry'
            // wait_minutes = Math.pow(5, task.attempt)
            // println "${result}: sleeping for ${wait_minutes} minutes. task ${task.hash}, failed with code ${task.exitStatus}: ${reason}"
            // sleep(wait_minutes * 60*1000 as long)
        }
        if (task.exitStatus in [11, 140]) { // temporary failure, or out of time
            result = 'retry'
        }
        println "WARN ${result} task ${task.getClass()} ${task.hash} ${reason}"
        return result
    }

    input:
    tuple val(obsid), val(meta)
    output:
    tuple val(obsid), val(meta), path(uvfits)

    script:
    prefix = "birli_"
    def (argstr, _argstr_cli, suffix) = birli_argstr_suffix();
    uvfits = ''+"${prefix}${obsid}*${suffix}.uvfits"
    print("workDir=${workflow.workDir}");
    // echo \$'meta=${meta}'
    """
    #!/bin/bash -eux
    export MWA_ASVO_API_KEY="${params.asvo_api_key}"
    export SINGULARITY_CACHEDIR="${workflow.workDir?(workflow.workDir+'/.singularity'):'/tmp'}"
    ${params.proxy_prelude} # ensure proxy is set if needed

    """ + ( params.pullPrep ? """
    # download if available in accacia
    ${params.rclone} ls "${params.bucket_prefix}.prep" --include '*${obsid}*'
    ${params.rclone} copy "${params.bucket_prefix}.prep/" . --include '${uvfits}' --progress --stats-one-line || true
    ls -al
    if [ \$(ls -1 ${uvfits} 2> /dev/null | wc -l) -gt 0 ]; then
        exit 0
    else
        echo "failed to download '${uvfits}' from bucket '${params.bucket_prefix}.prep'"
        echo "trying asvo"
    fi
    """ : "" ) + """

    # submit a job to ASVO, suppress failure if a job already exists.
    ${params.giant_squid} submit-conv -v -p ${argstr} ${obsid} || true

    # sleep 1s to prevent DoS
    sleep 1

    # list pending conversion jobs
    # ${params.giant_squid} list -j --types conversion --states queued,processing,error -- ${obsid}

    # sleep 1s to prevent DoS
    # sleep 1

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
        [ -f birli_manifest.sha1sum ] && cat birli_manifest.sha1sum
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
    storeDir "${params.outdir}/${obsid}/prep_qa"
    tag "${base}"
    label "python"
    label "nvme"
    memory {
        // max 30G Virtual for 60ts, 144T, 768chans
        [
            (50.GB * Math.pow(2,task.attempt) * (meta.ntimes?:60) * (meta.num_ants?:144) * (meta.num_chans?:768) \
                / (60 * 144 * 768)),
            360.GB
        ].min()
    }
    time {params.scratchFactor * 10.minute * Math.pow(2,task.attempt)}
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'ignore' }
    maxRetries 2

    input:
    tuple val(obsid), val(meta), path(metafits), path(uvfits)
    output:
    tuple val(obsid), val(meta), path(metrics)


    when: !(params.noprepqa || params.noqa)

    script:
    base = uvfits.baseName
    metrics = ''+"occupancy_${base}.json"
    template "flagqa.py"
}

process ssins {
    storeDir "${params.outdir}/${obsid}/prep_qa"
    tag "${base}"
    label "ssins"
    label "cpu_quarter"
    label "nvme"
    // TODO: make ssins memory footprint suck less
    // set memory limit from uvfits size. e.g. need 137GB for a 31GB file
    // should be 4.5×, but apparently it needs 6× (+1 if shm)
    memory { MemoryUnit.of(7 * uvfits.size()) }
    // set
    time { params.scratchFactor * 30.minute }
    // can't do this because no ntimes in prep meta
    // time { 30.min * ((meta.ntimes?:15) * (meta.nchans?:384) / (15 * 384) ) }

    input:
    tuple val(obsid), val(meta), path(uvfits)
    output:
    tuple val(obsid), val(meta),
        path("${output_prefix}_SSINS_mask.h5"),
        path("${output_prefix}{autos,cross,flagged}_SSINS*.png"),
        path("${output_prefix}ssins_occ.json")
        // path(ssins_uvfits)
        // todo: path("ssins_VDH.png"),
        // todo: path("match_events.json"),

    when: !params.nossins

    script:
    base = uvfits.baseName
    output_prefix = meta.output_prefix ?: "${base}${meta.subobs?:''}_"
    args = [ //noqa
        plot_title: meta.plot_title ?: "${base}${meta.subobs?:''}",
        sel_ants: meta.sel_ants ? meta.sel_ants.join(' ') : null,
        output_prefix: output_prefix,
        guard_width: meta.guard_width ?: (params.prep_freq_res_khz?:10) * 500,
        uvfits: uvfits
    ].findAll { _k, v -> v != null }
        .collect { k, v -> ''+"""--${k} ${v} """ }
        .join(" ")
    template "ssins.py"
}

process absolve {
    tag "${base}"
    label "ssins"
    label "nvme"
    label "mem_full"
    time { params.scratchFactor * 15.minute }
    storeDir "${params.outdir}/${obsid}/prep"

    input:
    tuple val(obsid), val(meta_), path(uvfits), path(mask)
    output:
    tuple val(obsid), val(meta), path(flagged)

    when: params.ssins_apply

    script:
    base = uvfits.baseName
    meta = deepcopy(meta_)
    if (meta_.name?:'' != '') {
        meta['name'] = ''+"${meta_.name}.ssins"
    } else {
        meta['name'] = "ssins"
    }
    flagged = "${base}.ssins.uvfits"
    args = "--src ${mask} --dst ${uvfits}"
    template "ssins_apply.py"
}

process autoplot {
    storeDir "${params.outdir}/${obsid}/prep_qa"
    tag "${base}${suffix}"
    label "python"
    label "nvme"
    label "mem_half"
    time { params.scratchFactor * 15.minute }

    input:
    tuple val(obsid), val(meta), path(metafits), path(uvfits)
    output:
    tuple val(obsid), val(meta), path(autoplot)


    when: !params.noautoplot

    script:
    base = uvfits.baseName
    suffix = meta.suffix?:""
    autoplot = ''+"autoplot_${base}${suffix}.png"
    args = meta.autoplot_args
    template "autoplot.py"
}

// make a reduced ao-format sourcelist with hyperdrive
process hypSrclistAO {
    storeDir "${params.outdir}/${obsid}/cal${params.cal_suffix}"
    tag "${obsid}"
    label "hyperdrive"
    label "cpu_full"
    time 15.minute

    input:
    tuple val(obsid), path(metafits)
    output:
    tuple val(obsid), path(reduced)

    when: !params.nosrclist

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
    storeDir "${params.outdir}/${obsid}/cal${params.cal_suffix}"
    tag "${obsid}"
    label "hyperdrive"
    label "cpu_full"
    time 15.minute

    input:
    tuple val(obsid), path(metafits)
    output:
    tuple val(obsid), path(reduced)

    when: !params.nosrclist

    script:
    reduced = ''+"${obsid}_reduced_n${params.sub_nsrcs}.yaml"
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
// process rexCluster {
//     storeDir "${params.outdir}/${obsid}/cal${params.cal_suffix}"
//     tag "${obsid}"
//     label "mwa_reduce"

//     input:
//     tuple val(obsid), path(reduced)
//     output:
//     tuple val(obsid), path(cluster)

//     script:
//     cluster = ''+"${obsid}_cluster_n${params.sub_nsrcs}_i${params.ionosub_nsrcs}.txt"
//     nclusters = params.sub_nsrcs / params.ionosub_nsrcs
//     """
//     #!/bin/bash -eux

//     ${params.mwa_reduce} cluster \
//         "${reduced}" \
//         "${cluster}" \
//         ${nclusters}
//     """
// }

// calibrate with mwa reduce
// process rexCalSol {
//     storeDir "${params.outdir}/${obsid}/cal${params.cal_suffix}"
//     tag "${obsid}${meta.subobs?:''}"
//     label "mwa_reduce"
//     label "nvme"
//     label "mem_full"
//     label "cpu_full"

//     input:
//     tuple val(obsid), val(dical_args), path(metafits), path(vis), val(tile_flags)
//     output:
//     tuple val(obsid),
//         path("rex_soln_${obsid}${meta.subobs?:''}_${name_glob}.bin"),
//         path("rex_di-cal_${obsid}${meta.subobs?:''}_${name_glob}.log")

//     script:
//     dical_names = dical_args.keySet().collect()
//     para = dical_names.size() > 1
//     name_glob = para ? "{" + dical_names.join(',') + "}" : dical_names[0]
//     flag_args = tile_flags.size() > 0 ? (''+"--tile-flags ${tile_flags.join(' ')}") : ""
//     vis_ms = "${uvfits.baseName}.ms"
//     """
//     #!/bin/bash -eux
//     """ + groovy2bashAssocArray(dical_args, "dical_args") + """
//     ${params.casa} -c "importuvfits('${vis}', '${vis_ms}')"
//     singularity exec ${params.cotter_sif} fixmwams vis.ms ${metafits}

//     for name in \${!dical_args[@]}; do
//         export soln_name="rex_soln_${obsid}${meta.subobs?:''}_\${name}.bin"
//         export log_name="rex_di-cal_${obsid}${meta.subobs?:''}_\${name}.log"
//         export args=\${dical_args[\$name]}
//         ${params.mwa_reduce} calibrate \${args} \
//             -applybeam -mwa-path ${params.img_mwa_path} \
//             -m "${params.sourcelist}" \
//             -i 50 \
//             -a 1e-4 1e-8 \
//             ${vis_ms} \
//             \${soln_name} | tee \${log_name}
//         ps=("\${PIPESTATUS[@]}")
//         if [ \${ps[0]} -ne 0 ]; then
//             echo "mwa_reduce failed. status=\${ps[0]}"
//             exit \${ps[0]}
//         fi
//     fi
//     """
// }

// ensure calibration solutions are present, or calibrate prep with hyperdrive
// do multiple calibration solutions for each obs, depending on dical_args
// meta is a hashmap containing extra arguments, think `**kwargs` in Python
process hypCalSol {
    storeDir "${params.outdir}/${obsid}/cal${params.cal_suffix}"
    tag "${obsid}${meta.subobs?:''}"

    // label jobs that need a bigger gpu allocation
    label "hyperdrive"
    memory { MemoryUnit.of(2 * uvfits.size()) }
    label "cpu_quarter"
    label "gpu_nvme"
    label "rate_limit_20"
    time { params.scratchFactor * 2.5.hour * Math.pow(task.attempt, 4) }
    errorStrategy {
        task.exitStatus in 137..140 ? 'retry' : 'ignore'
    }
    maxRetries 2
    stageInMode "copy"

    input:
    tuple val(obsid), val(meta), path(metafits), path(uvfits), val(dical_args)
    output:
    tuple val(obsid), val(meta),
        path("hyp_soln_${obsid}${meta.subobs?:''}_${name_glob}.fits"),
        path("hyp_di-cal_${obsid}${meta.subobs?:''}_${name_glob}.log", optional: true)
    // todo: model subtract: path("hyp_model_${obsid}_${name_glob}.uvfits")

    script:
    // stupid sanity check
    if (meta.obsid != null && obsid != meta.obsid) {
        throw new Exception("hypCalSol: obsid ${obsid} does not match meta.obsid ${meta.obsid}")
    }

    dical_names = dical_args.keySet().collect()
    para = dical_names.size() > 1
    name_glob = para ? "{" + dical_names.join(',') + "}" : dical_names[0]
    name_prefix = ''
    if (meta.name?:'' != '') {
        name_prefix = ''+"${meta.name}_"
    }
    name_glob = ''+"${name_prefix}${name_glob}"
    flag_args = ""
    prepFlags = meta.prepFlags?:[]
    fineChanFlags = meta.fineChanFlags?:[]
    if (prepFlags.size() > 0) {
        flag_args += " --tile-flags ${prepFlags.join(' ')}"
    }
    if (fineChanFlags.size() > 0) {
        flag_args += " --fine-chan-flags-per-coarse-chan ${fineChanFlags.join(' ')}"
    }

    // echo \$'meta=${meta}'
    """
    #!/bin/bash -eux
    """ + (para ? "export CUDA_VISIBLE_DEVICES=0" : "") + """
    export num_gpus="\$(${params.cmd_gpu_count})"
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
        export soln_name="hyp_soln_${obsid}${meta.subobs?:''}_${name_prefix}\${name}.fits"
        export log_name="hyp_di-cal_${obsid}${meta.subobs?:''}_${name_prefix}\${name}.log"
    """ + (params.pullCalSol ? """
        # download if available in accacia
        ${params.rclone} copy "${params.bucket_prefix}.soln/\${soln_name}" .
        if [ -f "\${soln_name}" ]; then
            touch \${log_name}
            continue
        fi
    """ : "") + """

        args=\${dical_args[\$name]}
        # hyperdrive di-cal backgrounded in a subshell
        (
            ${params.hyperdrive} di-calibrate \${args} \
                ${params.hyp_srclist_args} \
                ${params.hyp_dical_args} \
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
        # TODO: model subtract: --model-filenames "hyp_model_${obsid}${meta.subobs?:''}_\${name}.uvfits"
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
// process polyFit {
//     input:
//     tuple val(obsid), val(meta), path(soln)
//     output:
//     tuple val(obsid), val(meta), path(poly_soln), path(logs)

//     storeDir "${params.outdir}/${obsid}/cal${params.cal_suffix}"
//     tag "${obsid}${meta.subobs?:''}.${meta.dical_name}"

//     label "python"

//     script:
//     meta = mapMerge(meta, [name: ''+"poly_${meta.name}"])
//     poly_soln = "${meta.cal_prog}_soln_${obsid}${meta.subobs?:''}_${meta.name}.fits"
//     logs = "polyfit_${obsid}${meta.subobs?:''}_${meta.name}.log"
//     """
//     #!/bin/bash -eux
//     run_polyfit.py "${soln}" --outfile "${poly_soln}" | tee "${logs}"
//     ps=("\${PIPESTATUS[@]}")
//     if [ \${ps[0]} -ne 0 ]; then
//         echo "run_polyfit.py failed. status=\${ps[0]}"
//         exit \${ps[0]}
//     fi
//     """
// }

// do fourier fit of calibration solution phases
// meta is a map containing calibration solution metadata, fields:
// - name: name of the calibration (dical_name)
// - cal_prog: name of the calibration program (hyp or rex)
process phaseFit {
    storeDir "${params.outdir}/${obsid}/cal${params.cal_suffix}"
    tag "${obsid}${meta.subobs?:''}.${meta.dical_name}"
    label "fitcal"

    input:
    tuple val(obsid), val(meta_), path(soln)
    output:
    tuple val(obsid), val(meta), path(fit_soln) // , path(logs)

    when: params.fitphase

    script:
    meta = mapMerge(meta_, [name: ''+"${meta_.name}_fitted", fitted: true])
    base = soln.baseName
    fit_soln = ''+"${base}_fitted.fits"
    """
    #!/bin/bash -eux
    run_fitting.py "${soln}" --phase
    ps=("\${PIPESTATUS[@]}")
    if [ \${ps[0]} -ne 0 ]; then
        echo "run_fitting.py failed. status=\${ps[0]}"
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
    storeDir "${params.outdir}/${obsid}/cal${params.cal_suffix}"
    tag "${obsid}${meta.subobs?:''}.${meta.name}"
    label "hyperdrive_cpu"
    label "cpu_half"
    label "nvme"
    memory { MemoryUnit.of(2 * vis.size()) }

    input:
    tuple val(obsid), val(meta_), path(metafits), path(vis), path(soln)
    output:
    tuple val(obsid), val(meta), path(cal_vis), path(logs, optional: true)

    when: !params.nouv && !params.noapply

    script:
    meta = mapMerge(meta_, [name: ''+"${meta_.name}_${meta_.apply_name}"])
    cal_vis = ''+"hyp_${obsid}${meta.subobs?:''}_${meta.name}.uvfits"
    logs = ''+"hyp_apply_${meta.name}.log"
    args = meta.apply_args?:""
    if (meta.time_res != null) {
        args += " --time-average=${meta.time_res}s"
    }
    if (meta.freq_res_khz != null) {
        args += " --freq-average=${meta.freq_res_khz}kHz"
    }
    if (meta.nodut1) {
        args += " --ignore-dut1"
    }

    // echo \$'meta=${meta}'
    """
    #!/bin/bash -eux
    ${params.hyperdrive_cpu} solutions-apply ${args} \
        --data "${metafits}" "${vis}" \
        --solutions "${soln}" \
        --outputs "${cal_vis}" \
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
    storeDir "${params.outdir}/${obsid}/cal${params.cal_suffix}"
    // storeDir "/data/curtin_mwaeor/FRB_hopper/"

    tag "${obsid}${meta.subobs?:''}.${meta.name}"
    label "hyperdrive_cpu"
    label "cpu_half"
    label "gpu"

    input:
    tuple val(obsid), val(meta_), path(metafits), path(vis), path(soln)
    output:
    tuple val(obsid), val(meta), path(cal_vis), path(logs, optional: true)


    when: !params.noms && !params.noapply

    script:
    meta = mapMerge(meta_, [name: ''+"${meta_.name}_${meta_.apply_name}"])
    cal_vis = ''+"hyp_${obsid}${meta.subobs?:''}_${meta.name}.ms"
    logs = ''+"hyp_apply_${meta.name}_ms.log"
    // echo \$'meta=${meta}'
    """
    #!/bin/bash -eux
    ${params.hyperdrive_cpu} solutions-apply ${meta.apply_args} \
        --time-average=${meta.time_res}s \
        --freq-average=${meta.freq_res}kHz \
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
    storeDir "${params.outdir}/${obsid}/cal${params.cal_suffix}"
    tag "${obsid}${meta.subobs?:''}.${meta_.name}"
    label "hyperdrive"
    label "cpu_half"
    label "gpu"

    input:
    tuple val(obsid), val(meta_), path(metafits), path(vis), path(srclist)
    output:
    tuple val(obsid), val(meta), path(sub_vis), path(logs, optional: true)

    when: !params.nosub

    script:
    meta = mapMerge(meta_, [sub: "sub", name: ''+"sub_${meta_.name}"])
    sub_vis = ''+"hyp_${obsid}${meta.subobs?:''}_${meta.name}.uvfits"
    logs = ''+"hyp_vis-${meta.name}_uv.log"
    """
    #!/bin/bash -eux
    ${params.hyperdrive} vis-sub \
        --data "${metafits}" "${vis}" \
        --beam-file "${params.beam_path}" \
        --source-list "${srclist}" \
        --invert --num-sources ${meta.sub_nsrcs} \
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
    storeDir "${params.outdir}/${obsid}/cal${params.cal_suffix}"
    tag "${obsid}${meta.subobs?:''}.${meta_.name}"
    label "hyperdrive"
    label "cpu_half"
    label "gpu"

    input:
    tuple val(obsid), val(meta_), path(metafits), path(vis), path(srclist)
    output:
    tuple val(obsid), val(meta), path(sub_vis), path(logs, optional: true)


    when: !params.nosub

    script:
    meta = mapMerge(meta_, [sub: "sub", name: ''+"sub_${meta_.name}"])
    sub_vis = ''+"hyp_${obsid}${meta.subobs?:''}_${meta.name}.ms"
    logs = ''+"hyp_vis-${meta.name}_ms.log"
    """
    #!/bin/bash -eux
    ${params.hyperdrive} vis-sub \
        --data "${metafits}" "${vis}" \
        --beam-file "${params.beam_path}" \
        --source-list "${srclist}" \
        --invert --num-sources ${meta.sub_nsrcs} \
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
    storeDir "${params.outdir}/${obsid}/cal${params.cal_suffix}"
    tag "${obsid}${meta_.subobs?:''}.${meta_.name}_i${meta_.ionosub_nsrcs}"
    label "hyperdrive"
    label "cpu_half"
    label "mem_full" // need a full node otherwise gpu runs out of memory
    label "gpu"
    label "rate_limit_50"
    // 8h per 10GB file
    stageInMode "copy"
    time { params.scratchFactor * (coerceList(vis).collect { it.size() }.sum() / 10.GB.toBytes()) * 24.hour }

    input:
    tuple val(obsid), val(meta_), path(metafits), path(vis), path(srclist)
    output:
    tuple val(obsid), val(meta), path(sub_vis), path(json), path(logs, optional: true)

    when: !params.noionosub

    script:
    meta = mapMerge(meta_, [sub: "ionosub", name: ''+"ionosub_${meta_.name}${params.peel_suffix?:''}_i${meta_.ionosub_nsrcs}"])
    sub_vis = ''+"hyp_${obsid}${meta.subobs?:''}_${meta.name}.uvfits"
    logs = ''+"hyp_vis-${meta.name}_uv.log"
    json = ''+"hyp_peel_${obsid}${meta.subobs?:''}_${meta.name}_uv.json"
    """
    #!/bin/bash -eux
    export RUST_BACKTRACE=1
    ${params.hyperdrive} peel ${params.hyp_peel_args} \
        ${params.hyp_srclist_args} \
        -n ${params.sub_nsrcs} \
        --data "${metafits}" "${vis}" \
        --beam-file "${params.beam_path}" \
        --source-list "${srclist}" \
        --iono-sub ${meta.ionosub_nsrcs} \
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
    storeDir "${params.outdir}/${obsid}/cal${params.cal_suffix}"
    tag "${obsid}${meta_.subobs?:''}.${meta_.name}_i${meta_.ionosub_nsrcs}"
    label "hyperdrive"
    label "cpu_half"
    label "gpu"
    label "rate_limit_50"

    input:
    tuple val(obsid), val(meta_), path(metafits), path(vis), path(srclist)
    output:
    tuple val(obsid), val(meta), path(sub_vis), path(json), path(logs, optional: true)

    when: !params.noionosub

    script:
    meta = mapMerge(meta_, [sub: "ionosub", name: ''+"ionosub_${meta_.name}_i${meta_.ionosub_nsrcs}"])
    sub_vis = ''+"hyp_${obsid}${meta.subobs?:''}_${meta.name}.ms"
    logs = ''+"hyp_vis-${meta.name}_ms.log"
    json = ''+"hyp_peel_${obsid}${meta.subobs?:''}_${meta.name}_ms.json"
    """
    #!/bin/bash -eux
    export RUST_BACKTRACE=1
    ${params.hyperdrive}peel ${params.hyp_peel_args} \
        ${params.hyp_srclist_args} \
        --data "${metafits}" "${vis}" \
        --beam-file "${params.beam_path}" \
        --source-list "${srclist}" \
        --iono-sub ${meta.ionosub_nsrcs} \
        --sub ${meta.sub_nsrcs} \
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
    storeDir "${params.outdir}/${obsid}/iono_qa${params.cal_suffix}"
    tag "${obsid}${meta.subobs?:''}.${meta.name}"
    label "python"
    label "mem_half"
    time 1.hour

    input:
    tuple val(obsid), val(meta), path(srclist), path(offsets)
    output:
    tuple val(obsid), val(meta), path("cthuluplot_${title}*.png"), path(csv), path(json), path("tec_${title}*.png")

    script:
    title = ''+"${obsid}${meta.subobs?:''}_${meta.name}_i${meta.ionosub_nsrcs}"
    plot = ''+"cthuluplot_${title}.png"
    tec = ''+"tec_${title}.png"
    csv = ''+"ionoqa_${title}.csv"
    json = ''+"ionoqa_${title}.json"
    extra = "--plot-altaz 90 0 --offset-type arrow --show-axes --resolution 1024 --dpi 100"
    if (meta.time_res) {
        extra += " --time-res=${meta.time_res}"
    }
    if (meta.centre_freq) {
        extra += " --obs-frequency=${meta.centre_freq as int}"
    }
    if (meta.ra_phase_center != null && meta.dec_phase_center!= null) {
        extra += " --obs-radec ${meta.ra_phase_center} ${meta.dec_phase_center}"
    }
    template "cthulhuplot.py"
}

// process ionoqa {
//     input:
//     tuple val(obsid), val(meta), path(srclist), path(offsets)
//     output:
//     tuple val(obsid), val(meta), path(ionoqa)

//     storeDir "${params.outdir}/${obsid}/iono_qa${params.cal_suffix}"

//     tag "${obsid}.${meta.name}"

//     label "python"

//     time 10.minute

//     script:
//     title = "${obsid}_${meta.name}"
//     ionoqa = "ionoqa_${title}.json"
//     extra = ""
//     template "ionoqa.py"
// }

// QA tasks that can be run on preprocessed visibility files.
process prepVisQA {
    storeDir "${params.outdir}/${obsid}/prep_qa"
    tag "${base}"
    label "python"
    label "nvme"
    memory { MemoryUnit.of(uvfits.size() * 3) }
    time {params.scratchFactor * 20.minute * task.attempt}
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'ignore' }
    maxRetries 2

    input:
    tuple val(obsid), val(meta), path(metafits), path(uvfits)
    output:
    tuple val(obsid), val(meta), path(metrics)


    when: !(params.noprepqa || params.noqa)

    script:
    base = uvfits.baseName
    metrics = ''+"${base}_prepvis_metrics.json"
    """
    #!/bin/bash -eux
    run_prepvisqa.py ${uvfits} "${metafits}" --out "${metrics}"
    """
}

// QA tasks that can be run on calibrated visibility files.
process visQA {
    storeDir "${params.outdir}/${obsid}/vis_qa${params.cal_suffix}"
    tag "${obsid}${meta.subobs?:''}.${meta.name}"
    label "python"

    input:
    tuple val(obsid), val(meta), path(uvfits)
    output:
    tuple val(obsid), val(meta), path(json)

    script:
    base = uvfits.baseName
    json = ''+"${base}_vis_metrics.json"
    """
    #!/bin/bash -eux
    run_visqa.py ${uvfits} --out "${json}"
    """
}

process uvMeta {
    storeDir "${params.outdir}/${obsid}/vis_qa${params.cal_suffix}"
    tag "${base}"
    label "python"
    label "nvme"
    // needs 20min, 60GB of memory per 20GB of uvfits (at least 60GB)
    memory = { (task.attempt * volume * 60).GB }
    time {
        [
            (
                params.scratchFactor * task.attempt * 20.minute * volume
            ),
            task.attempt * 4.hour
        ].min()
    }
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'ignore' }
    maxRetries 2
    // stageInMode "symlink" // can't symlink any more, scratch too slow
    stageInMode "copy"

    input:
        tuple val(obsid), val(meta), path(vis)
    output:
        tuple val(obsid), val(meta), path(uvmeta)
        // tuple val(obsid), val(meta), path(vis), path(uvmeta), emit: obsMetaVisJson

    script:
    vis_=coerceList(vis)
    totalSizeGiga = vis_.collect { MemoryUnit.of((it.size())) }.sum()
    totalSizeGigaFloat = (float)(totalSizeGiga.toGiga())
    volume = (float)([totalSizeGigaFloat/20.0f,1.0f].max())
    base = vis_[0].baseName
    uvmeta = ''+"uvmeta_${base}.json"
    template "uvmeta.py"
}

// QA tasks that can be run on each calibration parameter set
process calQA {
    storeDir "${params.outdir}/${obsid}/cal_qa${params.cal_suffix}"
    tag "${obsid}${meta.subobs?:''}.${meta.name}"
    label "python"

    input:
    tuple val(obsid), val(meta), path(metafits), path(soln)
    output:
    tuple val(obsid), val(meta), path(metrics)


    when: !(params.nocalqa || params.noqa)

    script:
    metrics = ''+"metrics_${soln.baseName}_X.json"
    """
    #!/bin/bash -eux
    run_calqa.py "${soln}" "${metafits}" --pol X --out "${metrics}"
    """
}

process phaseFits {
    storeDir "${params.outdir}/${obsid}/cal_qa${params.cal_suffix}"
    tag "${obsid}.${meta.name}"
    label "mwax_mover"

    input:
    tuple val(obsid), val(meta), path(metafits), path(soln)
    output:
    tuple val(obsid), val(meta), path("${obsid}* ${name} phase_fits.tsv"), path("${obsid}*${name}*.png")

    when: !params.nophasefits

    script:
    name = ''+"${meta.cal_prog}_${meta.name}"
    """
    #!/bin/bash -eux
    python /app/scripts/cal_analysis.py \
        --name "${name}" \
        --metafits "${metafits}" --solns ${soln} \
        --phase-diff-path=/app/phase_diff.txt \
        --plot-residual --residual-vmax=0.5
    """
}

// write info from solutions to json
process solJson {
    storeDir "${params.outdir}/${obsid}/cal_qa${params.cal_suffix}"
    tag "${obsid}${meta.subobs?:''}.${meta.name}"
    label "python"

    input:
    tuple val(obsid), val(meta), path(soln)
    output:
    tuple val(obsid), val(meta), path(metrics)

    script:
    metrics = ''+"${meta.cal_prog}_soln_${obsid}${meta.subobs?:''}_${meta.name}.fits.json"
    template "soljson.py"
}

process plotPrepVisQA {
    storeDir "${params.outdir}/${obsid}/prep_qa"
    tag "${obsid}${meta.subobs?:''}"
    label "python"
    time {15.minute * task.attempt}

    input:
    tuple val(obsid), val(meta), path(metrics)
    output:
    tuple val(obsid), val(meta), path("${base}_{rms,modz}.png")

    when: !params.noplotprepqa

    script:
    base = ''+"prepvis_metrics_${metrics.baseName}"
    """
    #!/bin/bash -eux
    plot_prepvisqa.py "${metrics}" --out "${base}.png" --save
    """
}

process plotSols {
    storeDir "${params.outdir}/${obsid}/cal_qa${params.cal_suffix}"
    tag "${obsid}${meta.subobs?:''}.${meta.name}"
    label "hyperdrive_cpu"

    input:
    tuple val(obsid), val(meta), path(metafits), path(soln)
    output:
    tuple val(obsid), val(meta), path(plots_glob)

    script:
    plots_glob = ''+"${meta.cal_prog}_soln_${obsid}${meta.subobs?:''}*_${meta.name}_{phases,amps}*.png"
    """
    ${params.hyperdrive_cpu} solutions-plot ${params.hyp_sols_plot_args} -m "${metafits}" ${soln}
    """
}

process plotCalQA {
    storeDir "${params.outdir}/${obsid}/cal_qa${params.cal_suffix}"
    tag "${obsid}${meta.subobs?:''}.${meta.name}"
    label "python"

    input:
    tuple val(obsid), val(meta), path(metrics)
    output:
    tuple val(obsid), val(meta), path("${plot_base}_{rms,fft}.png")

    when: !params.noplotcalqa

    script:
    plot_base = ''+"calmetrics_${obsid}${meta.subobs?:''}_${meta.name}"
    """
    #!/bin/bash -eux
    plot_calqa.py "${metrics}" --out "${plot_base}" --save
    """
}

process plotVisQA {
    storeDir "${params.outdir}/${obsid}/vis_qa${params.cal_suffix}"
    tag "${obsid}"
    label "python"

    input:
    tuple val(obsid), val(meta), path(metrics)
    output:
    tuple val(obsid), val(meta), path("${meta.cal_prog}_${obsid}_${meta.name}_vis_metrics_*.png")

    when: !params.noplotvisqa

    script:
    """
    #!/bin/bash -eux
    plot_visqa.py "${metrics}" --out "${meta.cal_prog}_${obsid}_${meta.name}_vis_metrics.png" --save
    """
}

process plotImgQA {
    storeDir "${results_dir()}${params.img_suffix}${params.cal_suffix}"
    stageInMode "symlink"
    tag "${name}"
    label "python"

    input:
    tuple val(name), path("??????????.json")
    output:
    tuple val(name), path("${base}_*.png")

    script:
    base = ''+"wsclean_hyp_${name}-MFS"
    """
    #!/bin/bash
    set -ex
    plot_imgqa.py --out "${base}" --save ??????????.json
    """
}

process delaySpec {
    storeDir "${params.outdir}/${obsid}/vis_qa${params.cal_suffix}"
    tag "${obsid}${meta.subobs?:''}.${meta.name}"
    label "python"

    input:
    tuple val(obsid), val(meta), path(uvfits)
    output:
    tuple val(obsid), val(meta), path(dlyspec)

    when: !params.nodelayspec

    script:
    title = ''+"${obsid}${meta.subobs?:''}_${meta.name}"
    dlyspec = ''+"dlyspec_${title}.png"
    vmin = 3e13 // noqa
    vmax = 1e15 // noqa
    template "jline_delay_spec_from_uvfits.py"
}

// create dirty iamges of xx,yy,v
// TODO: merge img_params into meta?
process wscleanDirty {
    storeDir "${params.outdir}/${obsid}/img${params.img_suffix}${params.cal_suffix}"
    tag "${obsid}${meta.subobs?:''}${meta.inter_tok?:''}.${meta.name}"
    label "wsclean"
    label "cpu_half"
    label "mem_half"
    label "nvme"
    time { 1.8.minute * (1 + (multiplier * pix_mult * chan_mult * inter_mult)) }

    input:
    tuple val(obsid), val(meta), path(vis), val(img_params)
    output:
    tuple val(obsid), val(meta), path(img_glob)

    script:
    multiplier = vis.collect().size()
    // mult_suffix = multiplier > 1 ? (''+"_x${multiplier}") : ""
    // name_mult = "${name}${mult_suffix}"
    img_name = ''+"${meta.cal_prog}_${obsid}${meta.subobs?:''}${meta.subobs?:''}_${meta.name}${meta.inter_tok?:''}"

    // multipliers for determining compute resources
    pix_mult = 1 + (img_params.size / 1024) ** 2
    chan_mult = 1 + ("${img_params.channels_out}".split(' ')[0] as Double) / 25
    if (meta.interval) {
        inter_mult = meta.interval[1] - meta.interval[0]
    } else {
        inter_mult = 1 + ("${img_params.intervals_out}".split(' ')[0] as int) / 3
    }

    vis_ms = vis.collect {''+"${it.baseName}.ms"}
    vis = vis.collect()

    img_glob = ''+"wsclean_${img_name}"
    if (is_multiinterval() && !meta.inter_tok) {
        img_glob += "-t????"
    }
    if (is_multichannel()) {
        img_glob += "-{MFS,????}"
        // img_glob += "-*"
    }
    img_glob += "-{XX,YY,XY,XYi,I,Q,U,V}-{dirty,uv-real,uv-imag}.fits"
    img_args = ''+img_params.args
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
            vis_ms.collect {''+"${params.chgcentre} ${params.chgcentre_args} ${it}"}.join("\n") : \
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
// TODO: use IDG from Jack
// wsclean -name ./images/test_compact_idg -niter 10000 \
//    -scale 0.02 -size 2048 2048 \
//    -auto-threshold 0.5 -auto-mask 3 \
//    -pol I -multiscale -weight briggs 0  -j 10 -mgain 0.85 \
//    -no-update-model-required -abs-mem 30 \
//    -grid-with-beam -use-idg -idg-mode cpu -pb-undersampling 4 -mwa-path /astro/mwaeor/jline/software/ \
//    data/1088285720/compact_list-sky-model_8s_80kHz_fee-interp_band*.ms

// TODO: apply primary beam, save psf,
// NCHANS=768
// CHREM=$((NCHANS - CHAN))
// scale=0.25
// imsize=4096
// wsclean -mwa-path /pawsey/mwa/ \
//  -size ${imsize} ${imsize}  -scale ${scale}arcmin -channels-out ${CHREM} -channel-range ${CHAN} ${NCHANS} \
//  -apply-primary-beam -save-psf-pb -make-psf -reuse-primary-beam -weight natural -pol i \
//  /scratch/mwaeor/thimu/phase2_pawsey/moon_pawsey/cal_ms/${OBSID}.ms


process wscleanDConv {
    storeDir "${params.outdir}/${obsid}/img${img_params.suffix}${params.cal_suffix}"

    tag "${obsid}${meta.subobs?:''}${meta.inter_tok?:''}.${meta.name}.${(img_params.pol?:'').split(' ')[0]}"
    label "wsclean"
    label "cpu_quarter"
    label "nvme_full"
    // if (coerceList(vis).collect().size() > 1) {
    // } else {
    //     label "nvme"
    // }
    stageInMode "symlink"

    // cpus { max(4, min(36,1 + (Math.pow(task.attempt, 2) * multiplier))) as int }
    time { [23.hour, 40.minute * (1 + (Math.pow(task.attempt, 2) * multiplier))].min() }
    memory { [350.GB, 20.GB * (1 + (Math.pow(task.attempt, 2) * multiplier))].min() }
    maxRetries 3

    input:
    tuple val(obsid), val(meta), path(vis), val(img_params), path(dirtyImgs)
    output:
    tuple val(obsid), val(meta), path(img_glob)
        // path("wsclean_${img_name}-sources.txt") <- only works for stokes I

    when: !params.noimg

    script:
    multiplier = Math.sqrt(vis.collect().size())
    img_name = img_params.name?:''
    if (img_name == '') {
        img_name = ''+"${meta.cal_prog}_${obsid}${meta.subobs?:''}_${meta.name}${meta.inter_tok?:''}"
    }

    // multipliers for determining compute resources
    multiplier *= 1 + (img_params.size / 1024)
    multiplier *= 1 + ("${img_params.channels_out}".split(' ')[0] as int) / 25
    if (meta.interval) {
        inter_mult = meta.interval[1] - meta.interval[0]
    } else {
        inter_mult = 1 + ("${img_params.intervals_out}".split(' ')[0] as int) / 3
    }
    iter_mult = 1 + Math.sqrt(img_params.niter as Double) / 1000
    multiplier *= iter_mult
    vis_ms = vis.collect {''+"${it.baseName}.ms"}
    vis = vis.collect()
    img_glob = ''+(img_params.glob?:'')
    if (img_glob == '') {
        img_glob = "wsclean_${img_name}"
        if (is_multiinterval() && !meta.inter_tok) {
            img_glob += "-t????"
        }
        if (is_multichannel()) {
            img_glob += "-MFS"
        }
        if (img_params.pol) {
            img_glob += "-{XX,YY,XY,XYi,I,Q,U,V}"
        }
        img_glob += "-{image,psf,residual}.fits"
    }
    img_args = img_params.args
    if (img_params.weight) {
        img_args += " -weight ${img_params.weight}"
    }
    if (img_params.pol) {
        img_args += " -pol ${img_params.pol}"
    }
    // echo \$'meta=${meta}, img_params=${img_params}, img_glob=${img_glob}'
    """
    #!/bin/bash -eux
    echo \$'img_params=${img_params}, img_glob=${img_glob}, multiplier=${multiplier}'
    df -h .
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
        """ + ((meta.chan) ? "-channel-range ${meta.chan[0]} ${meta.chan[1]}" : "") + """ \
    """ + vis_ms.join(' ')
}

process imgQuantiles {
    storeDir "${params.outdir}/${obsid}/img_qa${params.img_suffix}${params.cal_suffix}"
    tag "${obsid}${meta.subobs?:''}.${meta.name}.${meta.suffix}"
    label "python"
    time {5.min * task.attempt}

    input:
    tuple val(obsid), val(meta), path(fits)
    output:
    tuple val(obsid), val(meta), path(csv)

    script:
    csv = "quantile_${fits.baseName}.csv"
    template "img_meta.py"
}

// power spectrum metrics via chips
process psMetrics {
    storeDir "${params.outdir}/${obsid}/ps_metrics${params.cal_suffix}"
    tag "${obsid}${meta.subobs?:''}.${meta.name}"
    label "chips"
    label "cpu_half"
    time 40.minute

    input:
    tuple val(obsid), val(meta), path("${meta.cal_prog}_${obsid}${meta.subobs?:''}_${meta.name}.uvfits")
    output:
    tuple val(obsid), val(meta), path(out_metrics)

    script:
    uvbase = ''+"${meta.cal_prog}_${obsid}${meta.subobs?:''}_${meta.name}"
    nchans = ''+meta.nchans
    eorband = ''+meta.eorband
    out_metrics = ''+"output_metrics_${uvbase}.dat"
    // echo \$'${meta}'
    """
    #!/bin/bash -eux
    export DATADIR="\$PWD"
    export OUTPUTDIR="\$PWD/"
    export OMP_NUM_THREADS=${task.cpus}
    which ${params.ps_metrics}
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
    tag "${group}"
    time {5.minute * task.attempt * params.scratchFactor}
    storeDir "${params.outdir}/${group}"

    input:
    tuple val(group), val(obsids)
    output:
    path(manifest)


    script:
    manifest = ''+"manifest_${group}.csv"
    """
    cat <<EOF > ${manifest}
${obsids.join('\n')}
EOF"""
}

// power spectrum with chips
process chipsGrid {
    storeDir "${params.outdir}/${group}/ps_metrics${params.cal_suffix}/${meta.name}"
    tag "${group}.${meta.name}.${pol}"
    label "chips"
    label "cpu_half"
    label "mem_full"
    label "nvme"
    maxRetries 3
    errorStrategy { return task.exitStatus > 1 ? "retry" : "ignore" }
    time { 50.minute * obsids.size() * params.scratchFactor }

    input:
    tuple val(group), val(meta), val(obsids), path(viss)
    output:
    tuple val(group), val(meta), path(grid)
        // path("vis_tot_${pol}.${ext}.dat"), \
        // path("vis_diff_${pol}.${ext}.dat"), \
        // path("noise_tot_${pol}.${ext}.dat"), \
        // path("noise_diff_${pol}.${ext}.dat"), \
        // path("weights_${pol}.${ext}.dat")
        // path("{vis_tot,vis_diff,noise_tot,noise_diff,weights}_${pol}.${ext}.dat")

    when: (!params.nopowerspec && !params.nochips)

    script:
    pol = meta.pol?:"xx"
    grid = "grid_${group}.${meta.name}.${pol}.tar.gz"
    lowfreq = meta.lowfreq
    ext = meta.ext
    nchans = meta.nchans
    eorband = meta.eorband
    eorfield = meta.eorfield
    period = meta.period?:"8.0"
    freq_res_hz = getFreqResHz(meta)
    freq_idx_start = meta.freq_idx_start?:0
    details = "details_${group}.${meta.name}.${pol}.txt"
    prepare_diff_cmd = """prepare_diff "${ext}" "${nchans}" "${freq_idx_start}" "${pol}" "${ext}" "${eorband}" -p "${period}" -c "${freq_res_hz}" -n "${lowfreq}" """

    // echo "meta=${meta}"
    """
    #!/bin/bash -ux
    """ + (eorfield == null || eorband == null || !nchans || !ext || lowfreq == null || !freq_res_hz || freq_idx_start == null || period == null ? "exit 2" : "" ) + """
    export DATADIR="\$PWD" INPUTDIR="\$PWD/" OUTPUTDIR="\$PWD/" OBSDIR="\$PWD/" OMP_NUM_THREADS=${task.cpus}

    echo eorfield=${eorfield} | tee -a ${details}
    echo eorband=${eorband} | tee -a ${details}
    echo nchans=${nchans} | tee -a ${details}
    echo ext=${ext} | tee -a ${details}
    echo pol=${pol} | tee -a ${details}
    echo lowfreq=${lowfreq} | tee -a ${details}
    echo freq_res=${freq_res_hz} | tee -a ${details}
    echo freq_idx_start=${freq_idx_start} | tee -a ${details}
    echo period=${period} | tee -a ${details}
    echo date=\$(date -Is) | tee -a ${details}
    echo host=\$(hostname) | tee -a ${details}
    echo gridvisdiff_bin=\$(which gridvisdiff) | tee -a ${details}

    """ + (
        [coerceList(obsids), coerceList(viss)].transpose().collect { obsid, vis ->
            def cmd = """gridvisdiff "${vis}" "${obsid}" "${ext}" "${eorband}" -f "${eorfield}" """
            """
    ${cmd}
    echo '${cmd}' >> ${details}
    export ps=\$?
    if [ \${ps} -ne 0 ] && [ \${ps} -ne 42 ]; then
        echo "gridvisdiff failed. status=\${ps}"
        exit \${ps}
    fi
    """
        }.join("\n")
    ) + """
    echo prepare_diff_bin=\$(which prepare_diff) | tee -a ${details}

    echo '${prepare_diff_cmd}' >> ${details}
    ${prepare_diff_cmd}
    export ps=\$?
    if [ \${ps} -ne 0 ] && [ \${ps} -ne 42 ]; then
        echo "prepare_diff failed. pol=${pol}, status=\${ps}"
        exit \${ps}
    fi
    # check for nulls
    for dat in *_${pol}.${ext}.dat; do
        if ! grep -m1 -qP "[^\\0]" \$dat; then
            echo "oops, all nulls \$dat"
            exit 3
        fi
    done

    echo "hashes" >> ${details}
    sha1sum {vis_tot,vis_diff,noise_tot,noise_diff,weights}_${pol}.${ext}.dat | tee -a ${details}
    tar -zcvf "${grid}" ${details} {vis_tot,vis_diff,noise_tot,noise_diff,weights}_${pol}.${ext}.dat
    exit 0
    """
}

process chipsCombine {
    storeDir "${params.outdir}/${group}/ps_metrics${params.cal_suffix}/${meta.name}"
    tag "${group}.${pol}.${meta.name}"
    label "chips"
    label "cpu_full"
    label "mem_super"
    label "nvme_full"
    stageInMode "symlink"
    // takes about 6 hours for 300, 16 hours for 1000,
    // double it for safety
    // clamp 1h < x < 24h
    time {
        def t = 3.hour
        if (exts.size() < 1000) {
            t = [2.hour, 6.hour * (exts.size() - 0.5) / 200].max()
        } else {
            t = [32.hour, 16.hour * (exts.size() - 0.5) / 1000].min()
        }
        [[params.scratchFactor * t, 1.hour].max(), 23.95.hour].min()
    }

    input:
    tuple val(group), val(meta), val(exts), path(grids)

    output:
    tuple val(group), val(meta), path("combine_${pol}.${combined_ext}.txt"), path(grid)

    script:
    // ext = "${group}_${meta.name}"
    combined_ext = meta.ext
    pol = meta.pol?:"xx"
    nchans = meta.nchans
    // echo "meta=${meta}"
    grid = "grid_${group}.${meta.name}.${pol}.tar.gz"
    details = "details_${group}.${meta.name}.${pol}.txt"
    """
    #!/bin/bash -eux
    """ + (!nchans || !combined_ext || exts.size() == 0 ? "exit 2" : "" ) + """

    export DATADIR="\$PWD" INPUTDIR="\$PWD/" OUTPUTDIR="\$PWD/" OBSDIR="\$PWD/" OMP_NUM_THREADS=${task.cpus}

    echo nchans=${nchans} | tee -a ${details}
    echo exts="${exts.join(' ')}" | tee -a ${details}
    echo combined_ext=${combined_ext} | tee -a ${details}
    echo pol=${pol} | tee -a ${details}
    echo date=\$(date -Is) | tee -a ${details}
    echo host=\$(hostname) | tee -a ${details}
    echo combine_data_bin=\$(which combine_data) | tee -a ${details}

    echo -n "" > combine_${pol}.${combined_ext}.txt
    cat > exts.txt <<EOF
""" + exts.join('\n') + """
EOF
    cat > grids.txt <<EOF
""" + grids.join('\n') + """
EOF
    set -euxo pipefail
    pwd
    df -h .
    while read -r ext grid; do
        tar -zxvf "\$grid" &
        echo "${pol}.\${ext}" | tee -a combine_${pol}.${combined_ext}.txt
    done < <(paste exts.txt grids.txt)
    # wait for all to finish.
    # the while ... done < <(paste ...) process substitution works with wait but
    # paste ... | while read ... done won't work with wait!
    wait

    # check error code, allow 0 or 42
    set +e
    combine_data "combine_${pol}.${combined_ext}.txt" ${nchans} "${pol}.${combined_ext}" 1
    export chips_return=\$?
    [ \$chips_return -eq 0 ] || [ \$chips_return -eq 42 ] || exit \$chips_return

    # check for nulls
    for dat in *_${pol}.${combined_ext}.dat; do
        if ! grep -m1 -qP "[^\\0]" \$dat; then
            echo "oops, all nulls \$dat"
            exit 3
        fi
    done

    echo "hashes" >> ${details}
    sha1sum {vis_tot,vis_diff,noise_tot,noise_diff,weights}_${pol}.${combined_ext}.dat | tee -a ${details}
    tar -zcvf "${grid}" ${details} {vis_tot,vis_diff,noise_tot,noise_diff,weights}_${pol}.${combined_ext}.dat
    exit 0
    """
}

process chipsLssa {
    storeDir "${params.outdir}/${group}/${params.lssa_bin}${params.cal_suffix}_b${bias_mode}/${meta.name}"
    tag "${group}.${params.lssa_bin}.${meta.name}.${pol}_${bias_mode}"
    label "chips"
    label "cpu_half"
    label "mem_half"
    label "nvme"
    time 1.hour

    input:
    tuple val(group), val(meta_), path(grid)

    output:
    tuple val(group), val(meta),
        path("{crosspower,residpower,residpowerimag,totpower,flagpower,fg_num,outputweights}_${pol}_${bias_number}.iter.${ext}.dat")

    script:
    maxu = params.lssa_maxu
    nbins = params.lssa_nbins
    bias_mode = (meta_.bias_mode ?: params.lssa_bias_mode ?:0 )
    // TODO: pretty sure meta_.lowfreq is midpoint but chips wants start frequency
    lowfreq = meta_.lowfreq
    nchan = meta_.nchans
    freq_res_hz = getFreqResHz(meta_)
    period = meta_.int_time

    // lssa = "lssa_${group}.${meta.name}.${pol}.tar.gz"

    if (bias_mode==0) {
        nchan_out = nchan
        start_chan = 0
    } else if (bias_mode==10) {
        nchan_out = Math.round(nchan / 2)
        start_chan = Math.round(nchan / 2)
    } else if (bias_mode==12) {
        nchan_out = Math.round(nchan / 2)
        start_chan = Math.round(nchan / 4)
    } else if (bias_mode==13) {
        nchan_out = Math.round(nchan / 2)
        start_chan = Math.round(nchan / 8)
    } else {
        throw new Exception("unknown bias_mode ${bias_mode}")
    }
    lssa_bin = params.lssa_bin
    // simple fft takes bias_mode
    bias_number = bias_mode // general fft takes  uses channel start in filename instead of bias_mode
    syntax_specific_args = "${bias_mode}"
    if (lssa_bin =~ /.*general.*/) {
        bias_number = start_chan
        syntax_specific_args = "${nchan_out} ${start_chan}"
    }
    meta = mapMerge(meta_, [
        nbins: nbins,
        maxu: maxu,
        bias_mode: bias_mode,
        bias_number: bias_number,
        nchan_out: nchan_out,
        start_chan: start_chan,
        lssa_bin: lssa_bin,
    ])

    ext = meta.ext
    pol = meta.pol?:"xx"
    eorband = meta.eorband
    details = "details_${group}.${meta.name}.${pol}.${params.lssa_bin}_b${bias_mode}.txt"

    if (eorband!=1) {
        throw new Exception("eorband=${eorband} currently hardcoded as 1 (high) in fft")
    }

    // echo "meta=${meta}"
    """
    #!/bin/bash -eux
    """ + (eorband == null || !nchan || nbins == null || !ext || maxu == null || bias_mode == null ? "exit 2" : "" ) + """

    tar -zxvf ${grid}

    export DATADIR="\$PWD" INPUTDIR="\$PWD/" OUTPUTDIR="\$PWD/" OBSDIR="\$PWD/" OMP_NUM_THREADS=${task.cpus}

    echo ext=${ext} | tee -a ${details}
    echo nchan=${nchan} | tee -a ${details}
    echo nchan_out=${nchan_out} | tee -a ${details}
    echo nbins=${nbins} | tee -a ${details}
    echo pol=${pol} | tee -a ${details}
    echo maxu=${maxu} | tee -a ${details}
    echo bias_mode=${bias_mode} | tee -a ${details}
    echo eorband=${eorband} | tee -a ${details}
    echo lowfreq=${lowfreq} | tee -a ${details}
    echo chanwidth=${freq_res_hz} | tee -a ${details}
    echo period=${period} | tee -a ${details}
    echo date=\$(date -Is) | tee -a ${details}
    echo host=\$(hostname) | tee -a ${details}
    echo lssa_bin=\$(which ${lssa_bin}) | tee -a ${details}

    ${lssa_bin} "${ext}" "${nchan}" "${nbins}" "${pol}" "${maxu}" "${ext}" ${syntax_specific_args} -p "${period}" -c "${freq_res_hz}" \
        2>&1 | tee ${lssa_bin}.txt

    export cross_size="\$(stat -c%s crosspower_${pol}_${bias_number}.iter.${ext}.dat)"
    if (( cross_size < 4097 )); then
        echo "crosspower_${pol}_${bias_number}.iter.${ext}.dat is too small (\$cross_size), exiting"
        exit 3
    fi
    """
    // tar -zcvf "${lssa}" ${details} {crosspower,residpower,residpowerimag,totpower,flagpower,fg_num,outputweights}_${pol}_${bias_mode}.iter.${ext}.dat
}

process chipsPlot {
    storeDir "${params.outdir}/${group}/${params.lssa_bin}${params.cal_suffix}_b${bias_mode}/${meta.name}"
    tag "${group}.${meta.name}.${ptype}.${pol}_${bias_mode}"
    label "chips_wrappers"
    time 15.minute

    input:
    tuple val(group), val(meta_), path(grid)
    output:
    tuple val(group), val(meta), path(plot)

    when: !params.noplotchips

    script:
    bias_mode = (meta_.bias_mode ?: params.lssa_bias_mode ?: 0 )
    pols_present = coerceList(grid).collect { (it.name =~ /[\w_]+_([xy]{2})_[\d]+.*/)[0][1] }.unique()
    if (pols_present.size() == 1) {
        pol = pols_present[0]
    } else if (pols_present.size() == 2) {
        pol = "both"
    } else {
        throw new Exception("unknown polarisations in ${grid}")
    }
    // meta.pol is the argument passed to plot, pol determines the filename
    meta = mapMerge(meta_, [pol: pol, bias_mode: bias_mode])
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
    freq_res_hz = getFreqResHz(meta)
    args = [
        // title: meta.title,
        // file group
        basedir: "./",
        bias_mode: bias_mode,
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
        kperp_max: (meta.kperp_max ?: params.kperp_max),
        kperp_min: (meta.kperp_min ?: params.kperp_min),
        kparra_min: (meta.kparra_min ?: params.kparra_min),
        plot_wedge_cut_2D: meta.plot_wedge_cut_2D,
        chips_tag_one_label: (meta.labels?:[])[0],
        chips_tag_two_label: (meta.labels?:[])[1],
        chips_tag_three_label: (meta.labels?:[])[2],
        // chips group
        lowerfreq: meta.lowfreq,
        chan_width: freq_res_hz,
        umax: meta.maxu,
        // density_correction: meta.density_correction,
        // omega_matter: meta.omega_matter,
        // omega_baryon: meta.omega_baryon,
        // omega_lambda: meta.omega_lambda,
        // hubble: meta.hubble,

    ].findAll { __, v -> v != null }
        .collect { k, v -> """--${k} "${v}" """ }
        .join(" ")
    if (meta.plot_delta) {
        args += " --plot_delta"
        suffix += "_delta"
    }
    if (meta.plot_wedge_cut) {
        args += " --plot_wedge_cut_2D"
    }
    if (meta.nchans) {
        args += " --N_chan_orig ${meta.nchans}"
    }
    if (meta.kperp) {
        args += " --K_perp ${meta.kperp}"
    }
    plot = "chips${dims}D_${pol}_${suffix}.png"

    // template "jline_plotchips.py"
    """
    plotchips_all.py ${args}
    mv *.png ${plot} || true # in case plot is not named correctly
    """
}

process chips1d_tsv {
    storeDir "${params.outdir}/${group}/${params.lssa_bin}${params.cal_suffix}_b${bias_mode}/${meta.name}"
    tag "${group}.${meta.name}.${pol}_${bias_mode}"
    label "chips_wrappers"
    label "nvme"
    time 15.minute

    input:
    tuple val(group), val(meta), path(grid)
    output:
    tuple val(group), val(meta), path(tsv)

    script:
    pol = meta.pol?:"null"

    suffix = "${meta.ext}"
    bias_mode = (meta.bias_mode ?: params.lssa_bias_mode ?: 0 )
    tsv = "1D_power_${pol}_${bias_mode}.${suffix}.tsv"
    lowerfreq = meta.lowfreq
    nchan = meta.nchans
    eorband = meta.eorband?:1

    freq_res_hz = getFreqResHz(meta)

    args = [
        title: meta.title,
        // file group
        basedir: "./",
        chips_tag: meta.ext,
        polarisation: meta.pol,

        wedge_factor: meta.wedge_factor,
        low_k_edge: meta.low_k_edge,
        high_k_edge: meta.high_k_edge,
        num_k_edges: meta.nbins,
        kperp_max: (meta.kperp_max ?: params.kperp_max),
        kperp_min: (meta.kperp_min ?: params.kperp_min),
        kparra_min: (meta.kparra_min ?: params.kparra_min),
        kparra_max: (meta.kparra_max ?: params.kparra_max),

        // chips group
        lowerfreq: lowerfreq,
        chan_width: freq_res_hz,
        umax: meta.maxu,
        bias_mode: bias_mode,
        density_correction: (meta.density_correction ?: params.density_correction),
        omega_matter: (meta.omega_matter ?: params.omega_matter),
        omega_baryon: (meta.omega_baryon ?: params.omega_baryon),
        omega_lambda: (meta.omega_lambda ?: params.omega_lambda),
        hubble: (meta.hubble ?: params.hubble),
        // num obs

    ].findAll { _k, v -> v != null }
        .collect { k, v -> """--${k} "${v}" """ }
        .join(" ")

    if (nchan) {
        args += " --N_chan_orig ${nchan}"
    }
    if (meta.kperp) {
        args += " --K_perp ${meta.kperp}"
    }
    if (params.k_edges != null) {
        args += " --ktot_bin_edges <(echo '${params.k_edges.join('\n')}')"
    }

    """
    chips1D_tsv.py ${args}
    cat *.tsv
    """
}

// analyse images of V,XX,YY
process imgQA {
    storeDir "${params.outdir}/${obsid}/img_qa${params.img_suffix}${params.cal_suffix}"
    tag "${obsid}${meta.subobs?:''}.${meta.name}"
    label "python"
    time {5.minute * task.attempt * params.scratchFactor}
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'ignore' }
    maxRetries 2

    input:
    tuple val(obsid), val(meta), path(fits)
    output:
    tuple val(obsid), val(meta), path(json)

    script:
    json = "wsclean_hyp_${obsid}${meta.subobs?:''}_${meta.name}-MFS.json"
    """
    #!/bin/bash -eux
    run_imgqa.py ${fits.join(' ')} --out ${json}
    """
}

// analyse images of V,XX,YY
process krVis {
    storeDir "${params.outdir}/${obsid}/img_qa${params.img_suffix}${params.cal_suffix}"
    tag "${obsid}${meta.subobs?:''}${meta.inter_tok?:''}.${meta.name}"
    label "python"

    input:
    tuple val(obsid), val(meta), path(fits), path(roma0)
    output:
    tuple val(obsid), val(meta), path("${out_prefix}{I,L,V,IV}-*.png")

    script:
    out_prefix = "${obsid}${meta.subobs?:''}${meta.inter_tok?"-${meta.inter_tok}":''}_${meta.name}-"
    template "krvis.py"
}

process uvPlot {
    storeDir "${params.outdir}/${obsid}/vis_qa${params.cal_suffix}"
    tag "${obsid}${meta.subobs?:''}.${meta.name}"
    label "python"

    input:
    tuple val(obsid), val(meta), path(uvfits)
    output:
    tuple val(obsid), val(meta), path("${uvplot}_{XX,YY,XY,YX}.png")

    when: !params.nouvplot

    script:
    title = "${obsid}${meta.subobs?:''}_${meta.name}"
    uvplot = "uvplot_${uvfits.baseName}"
    template "uvplot_2d.py"
}

// make a thumbnail png from a fits image
process thumbnail {
    storeDir "${params.outdir}/${obsid}/img_qa${params.img_suffix}${params.cal_suffix}"
    tag "${obsid}${meta.subobs?:''}.${meta.name}.${meta.suffix}"
    label "python"
    time {5.min * task.attempt * params.scratchFactor}
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'ignore' }
    maxRetries 2

    input:
    tuple val(obsid), val(meta), path(img)
    output:
    tuple val(obsid), val(meta), path(thumb)

    script:
    thumb = "${obsid}${meta.subobs?:''}_${meta.name}_${meta.suffix}.png"
    args = [
        fits: img,
        title: meta.title?:"${obsid}${meta.subobs?:''} ${meta.name} ${meta.suffix}",
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
        .collect { k, _v -> "--${k}" }

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
    storeDir "${params.outdir}/${obsid}/img_qa${params.img_suffix}${params.cal_suffix}"
    tag "${obsid}${meta.subobs?:''}${meta.subobs?:''}.${meta.name}.${meta.prod}.${meta.orderName}"
    label "python"

    input:
    tuple val(obsid), val(meta_), path(fits)
    output:
    tuple val(obsid), val(meta), path(polcomp)

    script:
    meta = mapMerge(meta_, ['orderName': meta_.order.join('')])
    polcomp = "${obsid}${meta.subobs?:''}${meta.subobs?:''}_${meta.name}_polcomp_${meta.prod}_${meta.orderName}.png"
    title = "${obsid}${meta.subobs?:''}${meta.subobs?:''} ${meta.name} ${meta.prod} ${meta.orderName}"
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
    storeDir "${params.outdir}/${obsid}/img_qa${params.img_suffix}${params.cal_suffix}"
    tag "${obsid}${meta.subobs?:''}.${meta.name}"
    label "imagemagick"

    input:
    tuple val(obsid), val(meta), path(thumbs)
    output:
    tuple val(obsid), val(meta), path(montage)

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
    storeDir "${results_dir()}${params.cal_suffix}"
    stageInMode "symlink"
    label "python"
    tag "${name}"

    input:
        tuple val(name), path("??????????.json")
    output:
        path("cal_qa*.png")

    when: !params.noplotcalqa

    script:
    """
    #!/bin/bash -eux
    plot_caljson.py --out cal_qa_${name} --save ??????????.json
    """
}

process stackImgs {
    storeDir "${params.outdir}/${chunk}/img${params.img_suffix}${params.cal_suffix}"
    tag "${chunk}${meta.subobs?:''}.${meta.name}"
    label 'imagemagick'

    input:
        tuple val(chunk), val(meta), path(imgs)
    output:
        tuple val(chunk), val(meta), path(stack)

    script:
    stack = "stack_${meta.name}_${meta.pol}.fits"
    """
    #!/bin/bash -eux
    convert ${imgs.join(' ')} -average -auto-gamma -auto-level ${stack}
    """
}

process stackThumbnail {
    storeDir "${params.outdir}/${chunk}/img_qa${params.img_suffix}${params.cal_suffix}"
    tag "${chunk}.${meta.name}.${meta.suffix}"
    label "python"
    time {10.min * task.attempt * params.scratchFactor}
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'ignore' }
    maxRetries 2

    input:
    tuple val(chunk), val(meta), path(img)
    output:
    tuple val(chunk), val(meta), path("${chunk}_${meta.name}_${meta.suffix}.png")

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

// process updateQuality {
//     input:
//         tuple val(obsid), val(name), val(params)
//     output:
//         tuple val(obsid), val(name), path("quality_${name}.json")

//     storeDir "${params.outdir}/${obsid}/meta"

//     tag "${obsid}.${name}"
// }

process tsvScatterPlot {
    storeDir "${results_dir()}${params.img_suffix}${params.cal_suffix}"
    stageInMode "copy"
    label "python"
    tag "${meta.name}"

    input:
        tuple val(meta), path(tsv)
    output:
        path(plot)

    script:
    plot = "${meta.name}.png"
    argstr = ([plot: plot, tsv: tsv] + meta).findAll { k, v ->
            v != null && [
                'tsv', 'x', 'y', 'c', 'plot', 'title', 'dpi', 'palette',
                'figwidth', 'figheight'
            ].contains(k)
        }
        .collect { k, v -> ["--${k}"] + coerceList(v) }
        .flatten()
        .collect { token -> "\"${token}\"" }
        .join(' ')
    template "plot_tsv.py"
}

process tarchive {
    storeDir "${results_dir()}${params.img_suffix}${params.cal_suffix}"
    stageInMode "copy"
    tag "${name}"

    input:
        tuple val(name), path(files), val(cachebust)
    output:
        tuple path(out), path("${dot_cachebust}")

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
    storeDir "${results_dir()}${params.img_suffix}${params.cal_suffix}"
    stageInMode "symlink"
    tag "${name}"
    label "ffmpeg"
    time { [8.hour, 1.hour * (task.attempt ** 2) * params.scratchFactor ].min() }
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'ignore' }
    maxRetries 2

    input:
        tuple val(name), path("??????????.png"), val(cachebust)
    output:
        tuple path("${name}.mp4"), path("${dot_cachebust}")

    when: !params.novideo
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
    tag "$x"
    storeDir "${params.outdir}/.fakebuckets/${bucket}"
    // label "rclone"
    label "datamover"

    // TODO: use size of File x to determine memory
    // memory = (x.size() * 1.5)

    input:
        tuple val(bucket_suffix), path(x)
    output:
        path("${x}.shadow")

    script:
    bucket = "${params.bucket_prefix}.${bucket_suffix}"
    """
    #!/bin/bash -eux
    touch "${x}.shadow"
    ${params.proxy_prelude} # ensure proxy is set if needed
    ${params.rclone} mkdir "${bucket}"
    ${params.rclone} copy --copy-links "$x" "${bucket}"
    """
}

process unArchive {
    tag "${storeDirSuffix}/${filename}"
    storeDir "${params.outdir}/${storeDirSuffix}"
    label "datamover_dl"

    input:
        tuple val(storeDirSuffix), val(bucket_suffix), val(filename)
    output:
        path("${filename}")

    script:
    if (bucket_suffix == null) {
        """
        touch "${filename}"
        """
    } else {
        bucket = "${params.bucket_prefix}.${bucket_suffix}"
        """
        #!/bin/bash -eux
        ${params.proxy_prelude} # ensure proxy is set if needed
        ${params.rclone} copy "${bucket}/${filename}" .
        """
    }
}

// collect ao quality metrics
process aoQuality {
    storeDir "${params.outdir}/${obsid}/vis_qa${params.cal_suffix}"

    input:
    tuple val(obsid), val(name), path("vis.ms")
    output:
    tuple val(obsid), val(name), path("${obsid}_${name}_aoquality_{sum,rfi,b,t,f,q}.tsv")

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

// use mwa webservices to gate obsids by pointing and faults
workflow ws {
    take:
        // channel of obs ids
        obsids

    main:
        obsids | wsMeta & tapMeta

        def quality_updates = channel.empty()
        if( file(params.quality_updates_path).exists() ) {
            quality_updates = file(params.quality_updates_path)
                .readLines()
                .findAll { !it.startsWith('#') && it.length() > 13 }
                .collectEntries { line ->
                    def (obsid, quality, comment) = line.split(',')
                    [obsid, [dataquality: quality, dataqualitycomment: comment]]
                }
        }

        def tile_updates = channel.empty()
        if( file(params.tile_updates_path).exists() ) {
            tile_updates = file(params.tile_updates_path)
                .readLines()
                .collect { line ->
                    def (firstObsid, lastObsid, tileIdxs) = line.split(',') + ['', '']
                    [firstObsid as int, lastObsid as int, (tileIdxs as String).split("\\|").collect {it as Integer} ]
                }
        }

        // print(params)

        def wsSummary = wsMeta.out.join(tapMeta.out).map { obsid, wsJson, filesJson, tapJson ->
                try {
                    def quality_update = quality_updates[obsid]?:[:]
                    def manualAnts = ([]) as Set
                    tile_updates.each {
                        def (firstObsid, lastObsid, tileIdxs, _comment) = it
                        if (obsid as int >= firstObsid && obsid as int <= lastObsid) {
                            manualAnts.addAll(tileIdxs)
                        }
                    }
                    def summary = wsSummarize(obsid, wsJson, filesJson, tapJson, quality_update, manualAnts)
                    [ obsid, deepcopy(summary) ]
                } catch (Exception e) {
                    println "error summarizing ${obsid}"
                    org.codehaus.groovy.runtime.StackTraceUtils.sanitize(e).printStackTrace()
                    throw e
                }
            }

        fail_codes = wsSummary
            // .filter { obsid, meta, uvfits, flagMeta -> (meta.prepFlags?:[]).size() > 0 }
            .map { obsid, summary ->
                [
                    obsid,
                    summary.fail_code
                ]
            }
        fail_codes.groupTuple(by: 1)
            .map { obsids_, fail_code ->
                [
                    fail_code,
                    obsids_.size(),
                ].join("\t")
            }
            .collectFile(
                name: "fail_counts_ws.tsv", newLine: true, sort: true,
                seed: ([ "FAIL CODE", "COUNT" ]).join("\t"),
                storeDir: "${results_dir()}"
            )
            | view { it.readLines().size() }

        wsSummary.map { obsid, summary ->
                [obsid, summary.lst.round().intValue()]
            }
            .groupTuple(by: 1)
            .map { obsids_, lst ->
                [
                    String.format("%+03d", lst),
                    obsids_.size(),
                ].join("\t")
            }
            .collectFile(
                name: "lst_counts_ws.tsv", newLine: true, sort: true,
                seed: ([ "LST", "COUNT" ]).join("\t"),
                storeDir: "${results_dir()}"
            )
            | view { it.readLines().size() }

        // display wsSummary
        wsStats = wsSummary.map { obsid, summary ->
                [
                    obsid,
                    (summary.fail_code==null||summary.fail_code==getFailReason(0x00))?'':summary.fail_code,
                    summary.fail_reasons.join('|'),
                    isNaN(summary.starttime_mjd)?'':String.format("%.5f", summary.starttime_mjd),
                    summary.starttime_utc?:'',

                    // pointing
                    String.format("% 5.2f", summary.ra_pointing),
                    String.format("% 5.2f", summary.dec_pointing),
                    String.format("% 5.2f", summary.az_pointing),
                    String.format("% 5.2f", summary.el_pointing),
                    isNaN(summary.ra_phase_center)?'':String.format("% 5.2f", summary.ra_phase_center),
                    isNaN(summary.dec_phase_center)?'':String.format("% 5.2f", summary.dec_phase_center),
                    isNaN(summary.ew_pointing)?'':String.format("% 2d", summary.ew_pointing),
                    isNaN(summary.sweet_pointing)?'':String.format("% 2d", summary.sweet_pointing),
                    String.format("% 5.2f", summary.lst),
                    summary.obs_name,
                    isNaN(summary.eorfield)?'':summary.eorfield,
                    isNaN(summary.sun_elevation)?'':summary.sun_elevation,
                    isNaN(summary.sun_pointing_distance)?'':summary.sun_pointing_distance,

                    // channels
                    summary.freq_res,
                    displayInts(summary.coarse_chans, ' '),
                    isNaN(summary.eorband)?'':summary.eorband,

                    // times
                    summary.int_time,
                    summary.nscans,

                    // config
                    summary.config,
                    summary.n_dead_dipoles,
                    summary.dead_dipole_frac,
                    summary.n_bad_tiles,
                    summary.n_good_tiles,
                    summary.bad_tile_frac,
                    displayInts(summary.tile_nums, ' '),
                    displayInts(summary.bad_tiles, ' '),
                    displayInts(summary.tile_rxs, ' '),

                    // iono quality
                    isNaN(summary.iono_magnitude)?'':summary.iono_magnitude,
                    isNaN(summary.iono_pca)?'':summary.iono_pca,
                    isNaN(summary.iono_qa)?'':summary.iono_qa,

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
                    "OBS", "FAIL CODE", "FAIL REASON", "START MJD", "START UTC",
                    "RA POINT", "DEC POINT", "AZ POINT", "EL POINT", "RA PHASE", "DEC PHASE",
                    "EW POINT", "SWEET POINT", "LST DEG", "OBS NAME", "EOR FIELD", "SUN ELEV", "SUN POINT",
                    "FREQ RES", "COARSE CHANS", "EOR BAND",
                    "TIME RES", "N SCANS",
                    "CONFIG", "N DEAD DIPOLES", "DEAD DIPOLE FRAC", "N BAD TILES", "N GOOD TILES", "BAD TILES FRAC",
                    "TILE NUMS", "FLAG TILES", "TILE RXS",
                    "IONO MAG", "IONO PCA", "IONO QA",
                    "N FILES", "N ARCHIVED",
                    "QUALITY", "QUALITY COMMENT",
                    "STATE FAULTS", "POINTING FAULTS", "FREQ FAULTS", "GAIN FAULTS", "BEAM FAULTS", "FAULT STR"
                ].join("\t"),
                storeDir: "${results_dir()}",
            )
            // display output path and number of lines
            | view { [it, it.readLines().size()] }

        wsStats.flatMap { file ->
                [
                    ["FAIL CODE", "tab10"],
                    ["IONO QA", "viridis"],
                    ["N BAD TILES", "viridis"],
                ].collect { col, palette ->
                    def short_name = col.replaceAll(/\s+/, '_').toLowerCase()
                    def meta = [
                        name: "ws_${short_name}", title: col, palette: palette,
                        x: 'OBS', y: 'LST DEG', c: col,
                    ]
                    [ meta, file ]
                }
            }
            | tsvScatterPlot

        wsSummary
            .filter { _o, summary -> summary.fail_reasons == [] }
            .map { obsid, _summary -> obsid }
            .tap { pass }
            | (wsMetafits & wsSkyMap & wsPPDs)

        wsMetafits.out | metaJson

        wsMetafits.out.map { obsid, metafits ->
                [ obsid, [:], metafits]
            }
            | demo03_mwalib

        demo03_mwalib.out.map { obsid, meta, antennasTsv_, channelsTsv_ ->
                def antennasTsv = coerceList(antennasTsv_)[0]
                def antennas = parseCsv2(antennasTsv, true, 0, '\t')
                def channelsTsv = coerceList(channelsTsv_)[0]
                def channels = parseCsv2(channelsTsv, true, 0, '\t')
                def newMeta = [
                    channels: channels,
                    antennas: antennas,
                ]
                [obsid, mapMerge(meta, newMeta)]
            }
            .tap { mwalibMeta }

    emit:
        // channel of good obsids with their metafits: tuple(obsid, metafits)
        obsMetafits = wsMetafits.out

        // channel of video name and frames to convert
        frame = pass.join(wsSkyMap.out)
            .map { _o, png -> ["skymap", png] }
            .mix( pass.join(wsPPDs.out.map { _o, png -> ["ppd", png] }) )
            .groupTuple()

        // channel of (obsid, metadata hashmap)

        obsMeta = wsSummary.map { obsid, summary ->
            def meta = [:]
            [
                // "groupid", "starttime_utc", "starttime_mjd", "obs_name"
                "ew_pointing", "centre_freq",
                "n_tiles", "bad_ants", "manual_ants", "tile_nums",
                "eorband", "eorfield", "lst", "int_time",
                "ra_phase_center", "dec_phase_center",
                "coarse_chans"
            ].each { key ->
                if (summary[key] != null) {
                    meta[key] = summary[key]
                }
            }
            if (summary['freq_res'] != null) {
                meta['freq_res_khz'] = summary['freq_res']
                meta['freq_res_hz'] = summary['freq_res'] * 1e3
            }
            if (summary["config"] != null) {
                meta["longconfig"] = summary["config"]
            }
            [obsid, deepcopy(meta)]
        }

        mwalibMeta = mwalibMeta

        fail_codes = fail_codes
}

workflow ssinsQA {
    take:
        subobsMetaVis

    main:
        // run ssins
        subobsMetaVis
            // - if ssins_by_rx, run ssins for each rx separately
            .flatMap { obsid, meta, uvfits ->
                def rxAnts = meta.unflaggedRxAnts?:[:]
                def allAnts = rxAnts.collect { _rx, ants -> ants }.flatten()
                def subobsAnts = [[meta.subobs, allAnts]]
                if (params.ssins_by_rx) {
                    subobsAnts += (rxAnts.collect {rx, ants ->
                        [(meta.subobs?:'') + String.format("_rx%02d", rx), ants]
                    })
                }
                def base = uvfits.baseName
                subobsAnts.findAll { _rx, ants -> ants.size() > 1 }
                    .collect { subobs, ants ->
                        // tile_idxs = ants.collect { it['Tile'] } as ArrayList
                        // ant_nums = ants.collect { it['Antenna'] + 1 } as ArrayList
                        def ssinsMeta = [
                            subobs: subobs,
                            plot_title: "\"${base}${subobs?:''}\\nn=${ants.size()}\"",
                            output_prefix: "${base}${subobs?:''}_",
                            // sel_ants: ant_nums
                            sel_ants: ants
                        ]
                        [obsid, mapMerge(meta, ssinsMeta), uvfits]
                    }
            }
            | ssins

        // analyse ssins occupancy
        ssinsOcc = ssins.out.map { def (obsid, meta, _mask, _plots, occ_json) = it;
                // def stats = [:]
                // parseJson(occ_json).each { item ->
                //     stats[item.key] = item.value
                // }
                def stats = parseJson(occ_json)
                stats.dab_total = stats.findAll {item ->
                        item.key.contains('DAB') && item.value != null
                    }
                    .collect { item -> item.value }
                    .sum()
                stats.narrow_total = stats.findAll {item ->
                        item.key.contains('narrow') && item.value != null
                    }
                    .collect { item -> item.value }
                    .sum()
                [obsid, meta, deepcopy(stats)]
            }

        ssinsOcc.map { def (obsid, meta, occ) = it;
                    ([
                        obsid,
                        meta.subobs?:'',
                        isNaN(occ.total)?'':String.format("%.2f", occ.total),
                        isNaN(occ.streak)?'':String.format("%.2f", occ.streak),
                        isNaN(occ.dab_total)?'':String.format("%.2f", occ.dab_total),
                        isNaN(occ.narrow_total)?'':String.format("%.2f", occ.narrow_total),
                    ]).join("\t")
                }
                .collectFile(
                    name: "ssins_occupancy.tsv", newLine: true, sort: true,
                    seed: (["OBS", "SUBOBS", "TOTAL", "STREAK", "DAB TOTAL", "NARROW TOTAL"]).join("\t"),
                    storeDir: "${results_dir()}"
                )
                | view { [it, it.readLines().size()] }

        allDABs = ssinsOcc.flatMap { _o, _m, occ ->
                occ.findAll { it.key.startsWith("DAB") }.collect { it.key }
            }
            .unique()
            .toSortedList()
            .map{it -> [it]}

        channel.of("OBS", "SUBOBS", "DAB TOTAL").concat(allDABs.flatten())
            .toList()
            .map { it.join("\t") }
            .concat(
                ssinsOcc.filter { _o, _m, occ ->
                        occ.dab_total > 0
                    }
                    .combine(allDABs)
                    .map { obs, meta, occ, dabs ->
                        ([
                            obs,
                            meta.subobs?:'',
                            occ.dab_total?:''
                        ] + dabs.collect { dab ->
                            occ[dab]?:''
                        }).join("\t")
                    }
                    .toSortedList()
                    .flatten()
            ).collectFile(
                name: "ssins_dab.tsv", newLine: true, sort: false,
                storeDir: "${results_dir()}"
            )
            .view { [it, it.readLines().size()] }

        allNarrows = ssinsOcc.flatMap { _o, _m, occ ->
                occ.findAll { it.key.contains('narrow') }.collect { it.key }
            }
            .unique()
            .toSortedList()
            .map{it -> [it]}

        channel.of("OBS", "SUBOBS", "NARROW TOTAL").concat(allNarrows.flatten())
            .toList()
            .map { it.join("\t") }
            .concat(
                ssinsOcc.filter { _o, _m, occ ->
                        occ.narrow_total > 0
                    }
                    .combine(allNarrows)
                    .map { obs, meta, occ, narrows ->
                        ([
                            obs,
                            meta.subobs?:'',
                            occ.narrow_total?:''
                        ] + narrows.collect { narrow ->
                            occ[narrow]?:''
                        }).join("\t")
                    }
                    .toSortedList()
                    .flatten()
            ).collectFile(
                name: "ssins_narrow.tsv", newLine: true, sort: false,
                storeDir: "${results_dir()}"
            )
            .view { [it, it.readLines().size()] }

        if (params.noprepqa || params.noqa) {
            subobsMetaVisPass = subobsMetaVis
            subobsMetaVisFlags = channel.empty()
            fail_codes = channel.empty()
        } else {
            // TODO: do we really need vis here?
            subobsMetaVisFlags = (subobsMetaVis.map { obsid, meta, uvfits -> [[obsid, meta.subobs?:''], meta, uvfits ]})
                .join(ssinsOcc.map { obsid, meta, occ -> [[obsid, meta.subobs?:''], occ ] }, remainder: false)
                .map { obsSubobs, meta, uvfits, ssinsStats ->
                    def (obsid, _subobs) = obsSubobs
                    def manualAnts = (meta.manual_ants?:[]) as Set
                    def flagMeta = [manualAnts: (manualAnts as ArrayList).sort(false)]
                    def newFlags = ([]) as Set
                    if (!params.noManualFlags && manualAnts.size() > 0) {
                        // println("${obsid} newFlags += manualAnts:${displayInts(manualAnts)}")
                        newFlags.addAll(manualAnts) // as Set?
                    }
                    if (ssinsStats != null) {
                        [ "dab_total", "narrow_total", "streak", "total" ].each { key ->
                            if (ssinsStats[key]) {
                                flagMeta[ "ssins_${key}" ] = ssinsStats[key]
                            }
                        }
                    }

                    def (fail_code, reason) = prepqa_pass(flagMeta)
                    // flag fail codes

                    flagMeta.fail_code = getFailReason(fail_code)
                    if (reason) {
                        flagMeta.reasons = reason
                    }

                    [obsid, meta, uvfits, deepcopy(flagMeta)]
                }

            fail_codes = subobsMetaVisFlags
                // .filter { obsid, meta, _vis, flagMeta -> (meta.prepFlags?:[]).size() > 0 }
                .map { obsid, _m, _vis, flagMeta ->
                    [ obsid, flagMeta.fail_code, flagMeta.reasons ]
                }

            fail_codes
                // .filter { _, fail_code, __ -> fail_code != getFailCode(0x00) }
                .groupTuple(by: 0)
                .map { obsid, obsCodes_, reasons ->
                    def (failCode, reason) = firstFail([
                        coerceList(obsCodes_),
                        coerceList(reasons),
                    ])
                    [ obsid, failCode, reason?:'' ].join("\t")
                }
                .collectFile(
                    name: "reasons_ssins.tsv", newLine: true, sort: true,
                    seed: ([ "OBSID", "FAIL CODE", "REASON" ]).join("\t"),
                    storeDir: "${results_dir()}"
                )
                | view { it.readLines().size() }

            fail_codes.groupTuple(by: 1)
                .map { obsids, fail_code, _reasons ->
                    [ fail_code, obsids.size() ].join("\t")
                }
                .collectFile(
                    name: "fail_counts_ssins.tsv", newLine: true, sort: true,
                    seed: ([ "FAIL CODE", "COUNT" ]).join("\t"),
                    storeDir: "${results_dir()}"
                )
                | view { it.readLines().size() }

            subobsMetaVisFlags
                // .filter { obsid, meta, uvfits, flagMeta -> (meta.ssinsFlagsTsv?:[]).size() > 0 }
                .map { obsid, meta, _vis, flagMeta ->
                    ([
                        obsid,
                        meta.subobs?:'',
                        meta.lst?:'',
                        isNaN(meta.ew_pointing)?'':meta.ew_pointing,
                        flagMeta.fail_code?:'',
                        flagMeta.reasons?:'',
                        flagMeta.ssins_dab_total==null?'':flagMeta.ssins_dab_total,
                        flagMeta.ssins_narrow_total==null?'':flagMeta.ssins_narrow_total,
                        flagMeta.ssins_streak==null?'':flagMeta.ssins_streak,
                        flagMeta.ssins_total==null?'':flagMeta.ssins_total,
                    ]).join("\t")
                }
                .collectFile(
                    name: "ssins_flags.tsv", newLine: true, sort: true,
                    seed: ([ "OBS", "SUBOBS", "LST DEG", "EW POINTING", "FAIL CODE", "REASONS", "SSINS DAB", "SSINS NARROW", "SSINS STREAK", "SSINS TOTAL" ]).join("\t"),
                    storeDir: "${results_dir()}"
                )
                // display output path and number of lines
                .view { [it, it.readLines().size()] }
                .flatMap { file ->
                    [
                        ["FAIL CODE", "tab10"],
                        ["SSINS TOTAL", "viridis"],
                        ["SSINS STREAK", "viridis"],
                        ["SSINS NARROW", "viridis"],
                        ["SSINS DAB", "viridis"],
                    ].collect { col, palette ->
                        def short_name = col.replaceAll(/\s+/, '_').toLowerCase()
                        def meta = [
                            name: "prep_${short_name}", title: col, palette: palette,
                            x: 'OBS', y: 'LST DEG', c: col,
                        ]
                        [ meta, file ]
                    }
                }
                | tsvScatterPlot

            subobsMetaVisPass = subobsMetaVisFlags.filter { _o, _m, _vis, flagMeta ->
                    (flagMeta.reasons?:'') == ''
                }
                .map { obsid, meta, uvfits, _flagMeta -> [obsid, meta, uvfits]}

            subobsMetaVisPass.map { obsid, meta, vis ->
                    ([
                        obsid,
                        meta.subobs?:'',
                        meta.name?:'',
                        meta.lst?:'',
                        isNaN(meta.ew_pointing)?'':meta.ew_pointing,
                        displayInts(meta.prepFlags?:[], ' '),
                        vis
                    ]).join("\t")
                }
                .collectFile(
                    name: "pass_ssins.tsv", newLine: true, sort: true,
                    seed: ([ "OBS", "SUBOBS", "NAME", "LST DEG", "EW POINTING", "NEW ANTS", "VIS" ]).join("\t"),
                    storeDir: "${results_dir()}"
                )
        }

        if (params.ssins_apply) {
            // TODO: dubious join might need cross
            subobsMetaVisSSINs = subobsMetaVisPass.join(
                    ssins.out.map { obsid, meta, mask, _plots, _json ->
                            [obsid, meta, mask]
                        }
                        .filter { _o, meta, _mask ->
                            !(((meta.subobs?:'').split('_')?:[''])[-1] =~ /rx\d+/)
                        }
                )
                .map { obsid, _meta, uvfits, meta, mask ->
                    [obsid, meta, uvfits, mask]
                }
                | absolve
            // TODO: flagqa
        } else {
            subobsMetaVisSSINs = subobsMetaVisPass
        }

    emit:
        subobsMetaVisSSINs
        subobsReasons = subobsMetaVisFlags.filter { _o, _m, _v, flagMeta ->
                (flagMeta.reasons?:'') != ''
            }
            .map { obsid, meta, _v, flagMeta -> [obsid, meta, flagMeta.reasons]}
        fail_codes
        frame = ssins.out.flatMap { _o, _m, _mask, imgs, _j ->
                imgs.collect { img ->
                    def tokens = coerceList(img.baseName.split('_') as ArrayList)
                    def suffix = tokens[-1]
                    def prefix = null
                    if (suffix == "SSINS") {
                        try {
                            prefix = tokens[-2]
                        } catch (Exception e) {
                            prefix = ''
                            println(tokens)
                            throw e
                        }
                        [''+"ssins_${prefix}", img]
                    } else {
                        try {
                            prefix = tokens[-3]
                        } catch (Exception e) {
                            prefix = ''
                            println(tokens)
                            throw e
                        }
                        [''+"ssins_${prefix}_${suffix}", img]
                    }
                }
            }
            .groupTuple()
        archive = ssins.out.map { _o, _m, mask, _imgs, _j -> ["ssins", mask]}
}

// ensure preprocessed uvfits are downloaded
workflow prep {
    take:
        // channel of obsids with their metafits: tuple(obsid, metafits)
        obsMetaMetafits
    main:
        // download preprocessed uvfits
        obsMetaMetafits.map { obsid, meta, __ -> [obsid, meta] }
            | asvoPrep

        asvoPrep.out.flatMap { obsid, meta, uvfits ->
                coerceList(uvfits).collect { f ->
                    def newMeta = [:]
                    def base_tokens = f.baseName.split('_') as ArrayList
                    if (base_tokens.size() > 2 && base_tokens[2] =~ /ch\d+/ ) {
                        newMeta['subobs'] = base_tokens[2]
                    }
                    [obsid, mapMerge(meta, newMeta), f]
                }
            }
            | uvMeta

        // TODO: dubious join might need cross
        obsMetaMetafits.join(asvoPrep.out).flatMap { obsid, _meta, metafits, meta, uvfits_ ->
                def uvfits = coerceList(uvfits_)
                if (uvfits.size > 1) {
                    uvfits.collect { f ->
                        def newMeta = [:]
                        def last_token = '' + f.baseName.split('_')[-1]
                        if (last_token =~ /ch\d+/) {
                            newMeta.subobs = last_token
                        }
                        [obsid, mapMerge(meta, newMeta), metafits, f]
                    }
                } else {
                    [[obsid, meta, metafits, uvfits_]]
                }
            }
            .tap { subobsMetaMetafitsPrep }

        subobsMetaMetafitsPrep.map { obsid, meta, _mf, uvfits -> [ obsid, meta, uvfits] }.tap { subobsMetaVis }
        obsMetaMetafits.map { obsid, _meta, metafits -> [obsid, metafits] }.tap { obsMetafits }

        flag(subobsMetaVis, obsMetafits)

        // for all obs that pass flag:
        flag.out.subobsFlagmetaPass.map { obsid, meta, flagMeta ->
                // update meta with flagMeta
                [[obsid, meta.subobs?:''], mapMerge(meta, flagMeta)]
            }
            // join with vis from subobsMetaVis
            .join(subobsMetaVis.map { obsid, meta, uvfits ->
                [[obsid, meta.subobs?:''], uvfits]
            })
            .map { obsSubobs, meta, uvfits ->
                def (obsid, _subobs) = obsSubobs
                [obsid, meta, uvfits]
            }
            | ssinsQA

    emit:
        // channel of obsids which pass the flag gate: tuple(obsid, meta, metafits, uvfits)
        subobsMetaVisPass = ssinsQA.out.subobsMetaVisSSINs
        subobsReasons = ssinsQA.out.subobsReasons
        // channel of video name and frames to convert
        frame = flag.out.frame
            .mix(ssinsQA.out.frame)

        // channel of files to archive, and their buckets
        archive = flag.out.archive
            .mix(ssinsQA.out.archive)

        zip = flag.out.zip

        fail_codes =
            flag.out.fail_codes.join(
                    ssinsQA.out.fail_codes.map { obsid, fail_code, _reasons -> [obsid, fail_code] },
                    remainder: true
                )
                .map { obsid, flag_code, prep_code ->
                    def codes = [flag_code, prep_code].findAll { it != null }
                    [obsid, coerceList(firstFail([codes]))[0]]
                }
}

// analyse (absolved) subobs prep uvfits for antenna flags, occupancy
workflow flag {
    take:
        // channel of subobservation visibilities and metadata: tuple(obsid, meta, uvfits)
        subobsMetaVis
        // channel of metafits for each observation
        obsMetafits

    main:
        // analyse preprocessed vis qa and aoflagger occupancy
        obsMetafits.cross(subobsMetaVis)
            .map { obsMetafits_, subobsMetaVis_ ->
                def (obsid, metafits) = obsMetafits_
                def (_o, meta, uvfits) = subobsMetaVis_
                [obsid, meta, metafits, uvfits]
            }
            .tap { subobsMetaMetafitsPrep }
            // TODO: flagQA first, then exclude ants from prepVisQA
            | (prepVisQA & flagQA)

        // collect flagQA results
        // TODO: collect sky_chans from flagQA
        sky_chans = (params.sky_chans).collect { ch -> "$ch".toString() }
        flagQA.out
            // form row of tsv from json fields we care about
            .map { obsid, meta, json ->
                def stats = parseJson(json)
                def flagged_sky_chans = stats.flagged_sky_chans?:[]
                def chan_occupancy = sky_chans.collect { ch ->
                    def occ = ((stats.channels?:[:])[ch]?:[:]).rfi_occupancy
                    isNaN(occ) ? '' : occ
                }
                def rfi_occ = stats.total_rfi_occupancy
                def flagged_sky_chan_idxs = displayInts(flagged_sky_chans, ' ')
                def flagged_fchan_idxs = displayInts(stats.flagged_fchan_idxs?:[], ' ')
                def flagged_timestep_idxs = displayInts(stats.flagged_timestep_idxs?:[], ' ')
                def preflagged_ants = displayInts(stats.preflagged_ants?:[], ' ')
                ([
                    obsid,
                    meta.subobs?:'',
                    stats.num_chans?:'',
                    flagged_sky_chan_idxs,
                    flagged_fchan_idxs,
                    flagged_sky_chans.size/chan_occupancy.size,
                    stats.num_times?:'',
                    stats.num_baselines?:'',
                    flagged_timestep_idxs,
                    preflagged_ants,
                    stats.total_occupancy,
                    isNaN(rfi_occ)?'':rfi_occ,
                ] + chan_occupancy).join("\t")
            }
            .collectFile(
                name: "occupancy.tsv", newLine: true, sort: true,
                seed: ([
                    "OBS", "SUBOBS",
                    "N CHAN", "FLAGGED SKY CHANS", "FLAGGED FINE CHANS", "FLAGGED SKY CHAN FRAC",
                    "N TIME","N BL", "FLAGGED TIMESTEPS", "FLAGGED ANTS",
                    "TOTAL OCCUPANCY", "NON PREFLAGGED"
                ] + sky_chans).join("\t"),
                storeDir: "${results_dir()}"
            )
            // display output path and number of lines
            | view { [it, it.readLines().size()] }

        // collect prepVisQA results as .tsv
        prepVisQA.out
            // form row of tsv from json fields we care about
            .map { obsid, meta, json ->
                def stats = parseJson(json);
                def bad_ants = stats.BAD_ANTS?:[]
                [
                    obsid,
                    isNaN(meta.ew_pointing)?'':meta.ew_pointing,
                    meta.subobs?:'',
                    stats.STATUS?:'',
                    stats.NANTS?:'',
                    stats.NTIMES?:'',
                    stats.NCHAN?:'',
                    stats.NPOLS?:'',
                    bad_ants.size(),
                    displayInts(bad_ants, ' '),
                ].join("\t")
            }
            .collectFile(
                name: "prepvis_metrics.tsv", newLine: true, sort: true,
                seed: [
                    "OBS",
                    "EWP",
                    "SUBOBS",
                    "STATUS",
                    "NANTS",
                    "NTIMES",
                    "NCHAN",
                    "NPOLS",
                    "N_BAD_ANTS",
                    "BAD_ANTS",
                ].join("\t"),
                storeDir: "${results_dir()}"
            )
            // display output path and number of lines
            | view { [it, it.readLines().size()] }

        [
            ["XX_RMS", { it -> (it.XX?:[:]).RMS?:[] }],
            ["YY_RMS", { it -> (it.YY?:[:]).RMS?:[] }],
            ["XX_MODZ0", { it -> ((it.XX?:[:]).MODZ_SCORE?:[:])["0"]?:[] }],
            ["XX_MODZ1", { it -> ((it.XX?:[:]).MODZ_SCORE?:[:])["1"]?:[] }],
            ["YY_MODZ0", { it -> ((it.YY?:[:]).MODZ_SCORE?:[:])["0"]?:[] }],
            ["YY_MODZ1", { it -> ((it.YY?:[:]).MODZ_SCORE?:[:])["1"]?:[] }],
        ].each { metric, getMetric ->
            prepVisQA.out
                // form row of tsv from json fields we care about
                .map { obsid, meta, json ->
                    def stats = parseJson(json);
                    ([ obsid, meta.subobs?:'' ] + getMetric(stats)).join("\t")
                }
                .collectFile(
                    name: "prepvis_${metric}.tsv", newLine: true, sort: true,
                    seed: [
                        "OBS",
                        "SUBOBS",
                        metric,
                    ].join("\t"),
                    storeDir: "${results_dir()}"
                )
                // display output path and number of lines
                | view { [it, it.readLines().size()] }
        }

        // plot prepvisQA
        prepVisQA.out | plotPrepVisQA

        if (params.noprepqa || params.noqa) {
            channel.empty().tap { subobsMetaFlags }
            subobsFlagmetaPass = subobsMetaVis.map { obsid, meta, _uvfits -> [obsid, meta, [:]] }
            subobsMetaPass = subobsFlagmetaPass.map { obsid, meta, _flagMeta ->
                    [obsid, meta]
                }

            fail_codes = channel.empty()

        } else {
            (subobsMetaMetafitsPrep.map { obsid, meta, _mf, _prep -> [[obsid, meta.subobs?:''], meta] })
                .join(flagQA.out.map { obsid, meta, flagJson -> [[obsid, meta.subobs?:''], parseJson(flagJson)]})
                .join(prepVisQA.out.map { obsid, meta, prepJson -> [[obsid, meta.subobs?:''], parseJson(prepJson)]})
                .map { obsidSubobs, meta, flagStats, prepStats ->
                    def (obsid, _subobs) = obsidSubobs
                    def manualAnts = (meta.manual_ants?:[]) as Set
                    def newMeta = [:]
                    def flagMeta = [manualAnts: (manualAnts as ArrayList).sort(false)]
                    if (meta.n_tiles != null) {
                        flagMeta.n_tiles = meta.n_tiles
                    }
                    def newFlags = ([]) as Set
                    if (!params.noManualFlags && manualAnts.size() > 0) {
                        newFlags.addAll(manualAnts) // as Set?
                    }
                    def prepAnts = []
                    if (prepStats != null) {
                        prepAnts = prepStats.BAD_ANTS?:[]
                        flagMeta += [ prepAnts: prepAnts ]
                        if (prepStats.STATUS) {
                            flagMeta += [prep_status: prepStats.STATUS]
                        }
                        if (!params.noPrepFlags && prepAnts.size() > 0) {
                            newFlags.addAll(prepAnts) // as Set?
                        }
                    }
                    def flagAnts = []
                    if (flagStats != null) {
                        flagAnts = flagStats.preflagged_ants?:[]
                        flagMeta += [
                            total_occ: flagStats.total_occupancy,
                            total_non_preflagged_bl_occ: flagStats.total_rfi_occupancy,
                        ]
                        if (flagStats.preflagged_ants) {
                            flagMeta += [ flagAnts: flagAnts ]
                        }
                        if (flagStats.flagged_fchan_idxs) {
                            newMeta += [ fineChanFlags: flagStats.flagged_fchan_idxs?:[] ]
                        }
                        if (flagStats.times != null) {
                            newMeta += [ ntimes: flagStats.num_times ]
                        }
                        if (flagStats.freqs != null) {
                            newMeta += [ nchans: flagStats.num_chans ]
                        }
                        if (flagAnts.size() > 0) {
                            newFlags.removeAll(flagAnts)
                        }
                        def inputs = (flagStats.INPUTS?:[])
                        if (inputs.size() > 0) {
                            def unflaggedAnts = inputs.findAll {
                                it['Pol'] == "X" && it['Flag'] == 0 && !newFlags.contains(it['Antenna'])
                            }
                            flagMeta += [
                                unflaggedAnts: unflaggedAnts,
                                unflaggedRxAnts: unflaggedAnts.groupBy { it['Rx'] }
                                    .collectEntries { k, v ->
                                        [(k): v.collect{ it_ -> it_['Antenna'] }]
                                    },
                                unflaggedFlavourAnts: unflaggedAnts.groupBy { it['Flavors'] }
                                    .collectEntries { k, v ->
                                        [(k): v.collect{ it_ -> it_['Antenna'] }]
                                    }
                            ]
                        }
                    }
                    if (newFlags.size() > 0) {
                        // println("${obsid} newFlags:${displayInts(newFlags)}")
                        newMeta += [prepFlags: (newFlags as ArrayList).sort(false)]
                        flagMeta += [prepFlags: (newFlags as ArrayList).sort(false)]
                    }

                    def (fail_code, reason) = prepqa_pass(flagMeta)

                    flagMeta.fail_code = getFailReason(fail_code)
                    if (reason) {
                        flagMeta.reasons = reason
                    }

                    [obsid, mapMerge(meta, newMeta), deepcopy(flagMeta)]
                }
                .tap { subobsMetaFlags }
                .map { obsid, meta, flagMeta ->
                    [ obsid, meta, flagMeta.fail_code, flagMeta.reasons ]
                }
                .tap { subobsMetaReasons }
                .groupTuple(by: 0)
                .map { obsid, _meta, obsCodes_, _reasons ->
                    def failCode = coerceList(firstFail([coerceList(obsCodes_)]))[0]
                    [ obsid, failCode ]
                }
                .tap { fail_codes }

            subobsFlagmetaPass = subobsMetaFlags.filter { _o, _m, flagMeta ->
                    (flagMeta.reasons?:'') == ''
                }
                .map { obsid, meta, flagMeta ->
                    [obsid, meta, flagMeta]
                }
            subobsMetaPass = subobsFlagmetaPass.map { obsid, meta, _flagMeta ->
                    [obsid, meta]
                }

            subobsMetaReasons
                .groupTuple(by: 0)
                .map { obsid, metas, obsCodes_, reasons ->
                    def ff = [null, null]
                    try {
                        ff = firstFail([
                            coerceList(obsCodes_),
                            coerceList(reasons),
                        ])
                        def (failCode, reason) = ff
                        return [ obsid, failCode, reason?:'' ].join("\t")
                    } catch (Exception e) {
                        println("obsid=${obsid} metas=${metas} obsCodes_=${obsCodes_} reasons=${reasons}")
                        println("ff=${ff}")
                        org.codehaus.groovy.runtime.StackTraceUtils.sanitize(e).printStackTrace()
                        throw e
                    }
                }
                .collectFile(
                    name: "reasons_flag.tsv", newLine: true, sort: true,
                    seed: ([ "OBSID", "FAIL CODE", "REASON" ]).join("\t"),
                    storeDir: "${results_dir()}"
                )
                | view { it.readLines().size() }

            fail_codes.groupTuple(by: 1)
                .map { obsids, fail_code ->
                    [ fail_code, obsids.size() ].join("\t")
                }
                .collectFile(
                    name: "fail_counts_flag.tsv", newLine: true, sort: true,
                    seed: ([ "FAIL CODE", "COUNT" ]).join("\t"),
                    storeDir: "${results_dir()}"
                )
                | view { it.readLines().size() }

            subobsMetaFlags
                .map { obsid, meta, flagMeta ->
                    ([
                        obsid,
                        meta.subobs?:'',
                        meta.lst?:'',
                        isNaN(meta.ew_pointing)?'':meta.ew_pointing,
                        flagMeta.fail_code?:'',
                        flagMeta.reasons?:'',
                        flagMeta.total_occ==null?'':flagMeta.total_occ,
                        flagMeta.total_non_preflagged_bl_occ==null?'':flagMeta.total_non_preflagged_bl_occ,
                        displayInts(flagMeta.flagAnts?:[], ' '),
                        displayInts(flagMeta.prepAnts?:[], ' '),
                        displayInts(flagMeta.manualAnts?:[], ' '),
                        displayInts(meta.prepFlags?:[], ' '),
                    ]).join("\t")
                }
                .collectFile(
                    name: "prep_flags.tsv", newLine: true, sort: true,
                    seed: ([
                        "OBS", "SUBOBS", "LST DEG", "EW POINTING", "FAIL CODE", "REASONS",
                        "PREP OCC", "RFI OCC", "ORIGINAL ANTS", "PREP ANTS",
                        "MANUAL ANTS", "NEW ANTS"
                    ]).join("\t"),
                    storeDir: "${results_dir()}"
                )
                // display output path and number of lines
                | view { [it, it.readLines().size()] }

            subobsMetaPass.map { obsid, meta ->
                    ([
                        obsid,
                        meta.subobs?:'',
                        displayInts(meta.prepFlags?:[], ' '),
                    ]).join("\t")
                }
                .collectFile(
                    name: "pass_flag.tsv", newLine: true, sort: true,
                    seed: ([ "OBS", "SUBOBS", "NEW ANTS" ]).join("\t"),
                    storeDir: "${results_dir()}"
                )
        }

        autoplotByRx = subobsFlagmetaPass.map {obsid, meta, flagMeta ->
                [[obsid, meta.subobs?:''], mapMerge(meta, flagMeta)]
            }.join(subobsMetaMetafitsPrep.map{obsid, meta, metafits, vis ->
                [[obsid, meta.subobs?:''], metafits, vis]
            }).flatMap { obsSubobs, meta, metafits, uvfits ->
                def (obsid, _subobs) = obsSubobs
                (meta.unflaggedRxAnts?:[:]).collect { rx, ants ->
                    def suffix = String.format("_rx%02d", rx);
                    def plotMeta = [
                        suffix: suffix,
                        "autoplot_args": "${params.autoplot_args} --sel_ants ${ants.join(' ')} --plot_title '${obsid}${suffix}'"
                    ]
                    [obsid, mapMerge(meta, plotMeta), metafits, uvfits]
                }
            }
        // autoplotByFlavor = subobsMetaFlags.map {obsid, meta, flagMeta ->
        //         [[obsid, meta.subobs?:''], flagMeta.unflaggedAnts]
        //     }.join(subobsMetaMetafitsPrep.map{obsid, meta, metafits, vis ->
        //         [[obsid, meta.subobs?:''], metafits, vis]
        //     }).flatMap { obsSubobs, meta, metafits, uvfits, flavorAnts ->
        //         def (obsid, _subobs) = obsSubobs
        //         flavorAnts.collect { flavor, ants ->
        //             suffix = String.format("_flavor-%s", flavor);
        //             sel_ants = ants.collect{it['Antenna']}
        //             plotMeta = [
        //                 suffix: suffix,
        //                 "autoplot_args": "${params.autoplot_args} --sel_ants ${sel_ants.join(' ')} --plot_title '${obsid}${suffix}'"
        //             ]
        //             [obsid, mapMerge(meta, plotMeta), metafits, uvfits]
        //         }
        //     }
        // autoplotByRx.mix(autoplotByFlavor)
        //     | autoplot
        // channel.empty() | autoplot
        autoplotByRx | autoplot

    emit:
        // channel of good subobs with their metafits: tuple(obsid, meta, uvfits)
        subobsMetaPass
        subobsFlagmetaPass
        // channel of video name and frames to convert
        frame = channel.empty()
            // TODO: fix
            // [f4/db1b29] NOTE: Process `extPrep:makeVideos:ffmpeg (prepvisqa_1322653536)` terminated with an error exit status (1) -- Error is ignored
            // [3b/4e3aa8] NOTE: Process `extPrep:makeVideos:ffmpeg (prepvisqa_1321964256)` terminated with an error exit status (1) -- Error is ignored
            // [5b/b0f82a] NOTE: Process `extPrep:makeVideos:ffmpeg (prepvisqa_1321447264)` terminated with an error exit status (1) -- Error is ignored
            // [76/f2946e] NOTE: Process `extPrep:makeVideos:ffmpeg (prepvisqa_1322653896)` terminated with an error exit status (1) -- Error is ignored
            // [59/5eb9f4] NOTE: Process `extPrep:makeVideos:ffmpeg (prepvisqa_1321791968)` terminated with an error exit status (1) -- Error is ignored
            // [4c/ee5983] NOTE: Process `extPrep:makeVideos:ffmpeg (prepvisqa_1321792688)` terminated with an error exit status (1) -- Error is ignored
            // plotPrepVisQA.out.flatMap { _, __, imgs ->
            //     imgs.collect { img ->
            //         def suffix = img.baseName.split('_')[-1]
            //         [''+"prepvisqa_${suffix}", img]
            //     }
            // }
            .mix(autoplot.out.map {_o, meta, img -> [''+"prepvisqa_autoplot${meta.suffix?:''}", img]})
            .groupTuple()
        archive = channel.empty() // TODO: archive flag jsons
        zip = prepVisQA.out.map { _o, _m, json -> ["prepvisqa", json]}
            .mix(flagQA.out.map { _o, _m, json -> ["flagqa", json]})
            .groupTuple()
        fail_codes
}

workflow cal {
    take:
        // channel of metafits and preprocessed uvfits: tuple(obsid, meta, metafits, uvfits)
        subobsMetaMetafitsVis
        // TODO: obsMetafits and subobsMetaVis instead
    main:
        // channel of metafits for each obsid: tuple(obsid, metafits)
        obsMetafits = subobsMetaMetafitsVis.map { obsid, _meta, metafits, _vis -> [obsid, metafits] }
        // hyperdrive di-calibrate on each obs
        subobsMetaMetafitsVis
            .map { def (obsids, meta, metafits, uvfits) = it
                [obsids, meta, metafits, uvfits, ["${params.dical_name}": params.dical_args]]
            }
            | hypCalSol

        // - hypCalSol gives multiple solutions, flatMap gives 1 tuple per solution.
        hypCalSol.out.flatMap { obsid, meta, solns, _log ->
                coerceList(solns).unique().collect { soln ->
                    // give each calibration a name from basename of solution fits.
                    // this is everything after the obsid
                    def dical_name = soln.baseName.split('_')[3..-1].join('_');
                    // print("dical_name(${soln.baseName}) => ${dical_name}")
                    def newMeta = [
                        dical_name: dical_name,
                        name: dical_name,
                        cal_prog: "hyp"
                    ]
                    [obsid, mapMerge(meta, newMeta), soln]
                }
            }
            // channel of individual dical solutions: tuple(obsid, meta, soln)
            .tap { eachCal }


        // hyperdrive dical log analysis
        // hypCalSol.out
        //     .flatMap { obsid, meta, solns, logs ->
        //         [
        //             (solns instanceof List ? solns : [solns]),
        //             (logs instanceof List ? logs : [logs])
        //         ].transpose().collect { soln, log ->
        //             name = soln.getBaseName().split('_')[3..-1].join('_')
        //             def (convergedDurationSec, convergedNumerator, convergedDenominator) = ['', '', '']
        //             if (!diCalLog.isEmpty()) {
        //                 startMatch = diCalLog.getText() =~ /([\d: -]+) INFO  hyperdrive di-calibrate/
        //                 startTime = null
        //                 if (startMatch) {
        //                     startTime = logDateFmt().parse(startMatch[0][1])
        //                 }
        //                 convergedMatch = (diCalLog.getText() =~ /([\d: -]+) INFO  All timesteps: (\d+)\/(\d+)/)
        //                 if (startTime && convergedMatch) {
        //                     convergedTime = logDateFmt().parse(convergedMatch[0][1])
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
        //         storeDir: "${results_dir()}${params.cal_suffix}"
        //     )
        //     | view { [it, it.readLines().size()] }

        // phase fit on calibration solutions
        obsMetafits.cross(
                eachCal.map { obsid, meta, soln -> [[obsid, meta.name], meta, soln] }
                    .groupTuple(by: 0)
                    .map{ obsName, metas, solns ->
                        def (obsid, _n) = obsName
                        [obsid, coerceList(metas)[0], solns]
                    }
            )
            .map { obsMetafits_, hypCalSol_ ->
                def (obsid, metafits) = obsMetafits_
                def (_o, meta, solns) = hypCalSol_
                [obsid, meta, metafits, solns]
            }
            | phaseFits

        // generate json from solns and fit phases
        eachCal | (solJson & phaseFit)
        // channel of all dical (and polyfit) solutions: tuple(obsid, meta, soln)
        allCal = eachCal.mix(
            phaseFit.out.map { obsid, meta, soln -> [obsid, meta, soln] }
        )

        // collect solJson results as .tsv
        solJson.out
            // form row of tsv from json fields we care about
            .map { obsid, meta, json ->
                def stats = parseJson(json)
                def results = (stats.RESULTS?:[])[0]?:[]
                def (nans, _convergences) = results.split { it == "NaN" }
                [
                    obsid,
                    meta.subobs?:'',
                    meta.lst,
                    meta.name,
                    (results.size() > 0 ? (nans.size() / results.size()) : ''),
                    results.join('\t')
                ].join("\t")
            }
            .collectFile(
                name: "cal_results${params.cal_suffix}.tsv", newLine: true, sort: true,
                seed: [
                    "OBS", "SUBOBS", "LST", "CAL NAME", "NAN FRAC", "RESULTS BY CH"
                ].join("\t"),
                storeDir: "${results_dir()}${params.cal_suffix}"
            )
            // display output path and number of lines
            | view { [it, it.readLines().size()] }

        // calibration QA and plot solutions
        obsMetafits.cross(allCal)
            // cross with solutions from eachCal and eachPolyCal
            // TODO: eliminate the need for metafits in plotSols and calQA
            .map { obsMetafits_, allCal_ ->
                def (obsid, metafits) = obsMetafits_
                def (_o, meta, soln) = allCal_
                [obsid, meta, metafits, soln]
            }
            | (plotSols & calQA)

        // plot each calQA result
        calQA.out | plotCalQA
        calQA.out.map { _o, meta, json ->
                def name = meta.name
                [name, json] }
            .groupTuple(by: 0)
            .map { name, jsons ->
                [name, jsons.sort(false)]
            }
            | plotCalJsons

        [
            ["XX_RMS", { stats -> (stats.XX?:[:]).RMS?:[] }],
            ["YY_RMS", { stats -> (stats.YY?:[:]).RMS?:[] }],
            ["XX_MODZ", { stats -> (stats.XX?:[:]).RMS_MODZ?:[] }],
            ["YY_MODZ", { stats -> (stats.YY?:[:]).RMS_MODZ?:[] }],
        ].each { metric, getMetric ->
            calQA.out
                // form row of tsv from json fields we care about
                .map { obsid, meta, json -> [obsid, meta, parseJson(json)] }
                .filter { _o, _m, stats -> stats != null }
                .map { obsid, meta, stats ->
                    ([ obsid, meta.subobs?:'', meta.name ] + getMetric(stats)).join("\t")
                }
                .collectFile(
                    name: "calqa_${metric}.tsv", newLine: true, sort: true,
                    seed: [
                        "OBS", "SUBOBS", "CAL NAME", metric,
                    ].join("\t"),
                    storeDir: "${results_dir()}${params.cal_suffix}"
                )
                // display output path and number of lines
                | view { [it, it.readLines().size()] }
        }

        allTiles = subobsMetaMetafitsVis.flatMap { _o, meta, _mf, _v ->
                meta.tile_nums
            }
            .unique()
            .toSortedList()
            .map{it -> [it]}
            // .view { it -> "allTiles ${it}" }

        // can't collectFile in a channel, so this is not possible?
        // keys = channel.of("intercept", "length", "chi2dof", "sigma_resid")
        // pols = channel.of("xx", "yy")
        // keys.combine(pols)
        [
            "intercept_xx", "intercept_yy",
            "length_xx", "length_yy",
            "chi2dof_xx", "chi2dof_yy",
            "sigma_resid_xx", "sigma_resid_yy",
        ].collect { key ->
            channel.of("OBS", "LST").concat(allTiles.flatten())
                .toList()
                .map { it.join("\t") }
                .concat(
                    phaseFits.out.combine(allTiles).map { obsid, meta, tsv, _pngs, tiles ->
                            phaseFits = parseCsv2(coerceList(tsv)[0], true, 0, '\t')
                            ([obsid, meta.lst] + tiles.collect { tile ->
                                def tilePhaseFits = phaseFits.find { "${it['tile_id']}" == "${tile}" }?:[:]
                                tilePhaseFits[key]?:''
                            }).join("\t")
                        }
                        .toSortedList()
                        .flatten()
                ).collectFile(
                    name: "phase_fits_${params.dical_name}_${key}.tsv", newLine: true, sort: false,
                    storeDir: "${results_dir()}${params.cal_suffix}"
                )
                .view { [it, it.readLines().size()] }
        }

        if (params.noqa) {
            obsMetaStats = allCal.map { obsid, meta, _soln ->
                    [obsid, meta]
                }
        } else {
            obsMetaStats = calQA.out
                .map { obsid, meta, json ->
                    def stats = parseJson(json)
                    def newMeta = [
                        unused_bls: parseFloatOrNaN(stats.PERCENT_UNUSED_BLS) / 100,
                        unconvg_chs: parseFloatOrNaN(stats.PERCENT_NONCONVERGED_CHS) / 100,
                        rms_convg: parseFloatOrNaN(stats.RMS_CONVERGENCE),
                        skewness: parseFloatOrNaN(stats.SKEWNESS),
                        rx_var: parseFloatOrNaN(stats.RECEIVER_VAR),
                        dfft_pow: parseFloatOrNaN(stats.DFFT_POWER),
                        calqa_bad_ants: stats.BAD_ANTS?:[],
                        calqa_status: stats.STATUS?:'',
                        calqa_fail_reason: stats.FAILURE_REASON?:'',
                    ]
                    def (fail_code, reason) = calqa_pass(newMeta)
                    newMeta.fail_code = fail_code
                    newMeta.reason = reason
                    if (meta.calqa_bad_ants) {
                        def prepFlags = (meta.prepFlags?:[]) as Set
                        def calFlags = ([]) as Set
                        if (!params.noCalFlags) {
                            calFlags.addAll(meta.calqa_bad_ants?:[])
                        }
                        def newflags = (calFlags - prepFlags) as ArrayList
                        if (newflags) {
                            newMeta.calFlags = newflags.sort(false)
                        }
                    }
                    [obsid, mapMerge(meta, newMeta)]
                }

            obsMetaStats.map { obsid, meta ->
                    def bad_ants = (meta.calqa_bad_ants?:[])
                    [
                        obsid,
                        meta.subobs?:'',
                        isNaN(meta.lst)?'':meta.lst,
                        isNaN(meta.ew_pointing)?'':meta.ew_pointing,
                        meta.name,
                        meta.calqa_status?:'',
                        bad_ants.size,
                        displayInts(bad_ants, ' '),
                        isNaN(meta.unused_bls)?'':meta.unused_bls,
                        isNaN(meta.unconvg_chs)?'':meta.unconvg_chs,
                        isNaN(meta.rms_convg)?'':meta.rms_convg,
                        isNaN(meta.skewness)?'':meta.skewness,
                        isNaN(meta.rx_var)?'':meta.rx_var,
                        isNaN(meta.dfft_pow)?'':meta.dfft_pow,
                        meta.calqa_fail_reason?:'',
                    ].join("\t")
                }
                .collectFile(
                    name: "cal_metrics_${params.dical_name}${params.cal_suffix}.tsv", newLine: true, sort: true,
                    seed: [
                        "OBS", "SUBOBS", "LST", "EWP", "CAL NAME", "STATUS",
                        "N BAD ANTS", "BAD ANTS",
                        "BL UNUSED FRAC", "NON CONVG CHS FRAC",
                        "RMS CONVG",
                        "SKEW",
                        "RECV VAR",
                        "DFFT POW",
                        "FAILURE_REASON"
                    ].join("\t"),
                    storeDir: "${results_dir()}${params.cal_suffix}"
                )
                // display output path and number of lines
                | view { [it, it.readLines().size()] }
        }
        // channel of obsids and names that pass qa. tuple(obsid, name)
        // - take tuple(obsid, cal_name, json) from calQA.out
        // - filter on json.STATUS == "PASS"
        // - take obsid and name
        obsMetaPass = obsMetaStats
            // TODO: reintroduce status filter
            // .filter { _, meta, stats ->
            //     stats.STATUS == null || stats.STATUS == "PASS"
            // }
            .filter { _o, meta -> meta.fail_code == 0x00 }
            .map { obsid, meta -> [obsid, meta] }

        fail_codes = obsMetaStats.map{ obs, meta ->
                [obs, getFailReason(meta.fail_code?:0x00), meta.reason?:""]
            }

        fail_codes.groupTuple(by: 0)
                .map { obsid, obsCodes_, reasons ->
                    def (failCode, reason) = firstFail([
                        coerceList(obsCodes_),
                        coerceList(reasons),
                    ])
                    [ obsid, failCode, reason?:'' ].join("\t")
                }
                .collectFile(
                    name: "reasons_calqa_${params.dical_name}${params.cal_suffix}.tsv", newLine: true, sort: true,
                    seed: ([ "OBSID", "FAIL CODE", "REASON" ]).join("\t"),
                    storeDir: "${results_dir()}${params.cal_suffix}"
                )
                | view { it.readLines().size() }

    emit:
        fail_codes = fail_codes.filter { _o, fail_code, _reason -> fail_code != getFailReason(0x00) }
        // channel of calibration solutions that pass qa. tuple(obsid, name, cal)
        // - take tuple(obsid, meta, soln) from allCal
        // - match with obsMetaPass on (obsid, cal_name)
        obsMetaCalPass = allCal
            .cross(obsMetaPass) {def (obsid, meta) = it; [obsid, meta.name]}
            .map{ allCal_, obsMetaPass_ ->
                def (obsid, _m, soln) = allCal_
                def (_o, meta) = obsMetaPass_
                [obsid, meta, soln]
            }
        // channel of files to archive, and their buckets
        archive =
            hypCalSol.out
                .flatMap { _o, _m, solns, _log ->
                    coerceList(solns).collect { soln ->
                        [ "soln", soln ]
                    }
                }
                .mix(
                    calQA.out.map { _o, _m, json -> ["calqa", json]}
                )

        // channel of video name and frames to convert
        frame = plotCalQA.out.mix(plotSols.out)
            .flatMap { _o, meta, pngs ->
                coerceList(pngs).collect { png ->
                    def tokens = png.baseName.split('_').collect()
                    def suffix = tokens.removeLast()
                    if(suffix =~ /\d+/) {
                        // def timeblock = suffix
                        suffix = tokens.removeLast()
                    }
                    [''+"calqa${params.cal_suffix?:''}_${meta.name}_${suffix}", png]
                }
            }
            .mix(phaseFits.out.flatMap { _o, _m, _tsv, pngs ->
                pngs.collect { png ->
                    def suffix = png.baseName.split(' ')[-1]
                    ["phasefits_${suffix}", png]
                }
            })
            .groupTuple()
        // channel of files to zip
        zip = calQA.out.map { _o, meta, json -> ["calqa_${meta.name}", json] }
            .mix(solJson.out.map { _o, meta, json -> ["soljson_${meta.name}", json] })
            .mix(phaseFits.out.map { _o, meta, tsv, _pngs -> ["phasefits_${meta.name}", tsv] })
            .mix(allCal.map { _o, meta, soln -> ["${meta.cal_prog}_soln_${meta.name}", soln] })
            .groupTuple()
}

// process uvfits visibilities
workflow uvfits {
    take:
        obsMetaUV
    main:
        // vis QA on each compact obsid, no subtractions
        obsMetaUV
            .filter { _o, meta, _vis -> meta.sub == null && meta.config =~ /.*Compact.*/ }
            | visQA

        obsMetaUV | uvPlot

        // collect visQA results as .tsv
        visQA.out
            // form row of tsv from json fields we care about
            .map { obsid, meta, json ->
                // todo: ignore manual and bad ants from meta
                // + echo 'meta=[name:ssins_30l_src4k_300it_8s_80kHz, cal_prog:hyp, obsid:1094404168, ew_pointing:0, obs_name:high_season2_2456911, starttime_utc:2014-09-10T17:09:12.000Z, starttime_mjd:56910.715, centre_freq:1.8240E+8, n_tiles:128, bad_ants:[85], manual_ants:[], tile_nums:[88, 111, 112, ...], eorband:1, eorfield:0, lst:3.6285414695739746, int_time:2.0, freq_res:80, ra_phase_center:0.0, dec_phase_center:-27.0, fineChanFlags:[0, 1, 16, 30, 31], unflaggedRxAnts:[1:[0, 1, 2, 3, 4, 5, 6, 7], ...], unflaggedFlavourAnts:[RG6_90:[0, 1, 5, 7, 8, 9, 10, 15, 16, 17, 18, 29, 31, 32, 34, 39, 54, 63, 120], ...], dical_name:ssins_30l_src4k_300it, fail_code:0, reason:[], calFlags:[35], time_res:8, nodut1:false, apply_args: --tile-flags 35, apply_name:8s_80kHz]'
                def stats = parseJson(json)
                def vis_ants = stats.POOR_ANTS?:[]
                def new_ants = vis_ants as Set
                def db_ants = meta.bad_ants?:[]
                new_ants -= db_ants as Set
                def manual_ants = meta.manual_ants?:[]
                new_ants -= manual_ants as Set
                def prep_ants = meta.prepFlags?:[]
                new_ants -= prep_ants as Set
                def cal_ants = meta.calFlags?:[]
                [
                    obsid,
                    meta.name,
                    displayInts(db_ants, ' '),
                    displayInts(manual_ants, ' '),
                    displayInts(prep_ants, ' '),
                    displayInts(cal_ants, ' '),
                    stats.NPOOR_ANTS?:'',
                    displayInts(vis_ants, ' '),
                    stats.NPOOR_BLS?:'',
                    displayInts(new_ants as ArrayList, ' '),
                ].join("\t")
            }
            .collectFile(
                name: "vis_metrics.tsv", newLine: true, sort: true,
                seed: [
                    "OBS", "VIS NAME",
                    "DB_ANTS", "MANUAL_ANTS", "PREP_ANTS", "CAL_ANTS",
                    "NPOOR_ANTS", "POOR ANTS", "NPOOR_BLS",
                    "NEW_ANTS",
                ].join("\t"),
                storeDir: "${results_dir()}${params.cal_suffix}"
            )
            // display output path and number of lines
            | view { [it, it.readLines().size()] }

        // write bands and fields to a file
        obsMetaUV.map { obsid, meta, vis ->
            [
                obsid, meta.name, meta.config?:'', meta.eorband?:'', meta.eorfield?:'',
                meta.nchans?:'', meta.lowfreq?:'', meta.freq_res_khz?:'', vis
            ].join('\t') }
            .collectFile(
                name: "eor_params.tsv", newLine: true, sort: true,
                seed: [
                    "OBSID", "NAME", "CONFIG", "BAND", "FIELD", "NCHAN", "LOWFREQ", "FREQRES", "VIS"
                ].join('\t'),
                storeDir: "${results_dir()}${params.cal_suffix}"
            )
            // display output path and number of lines
            | view { [it, it.readLines().size()] }

        visQA.out | plotVisQA

        if (params.noeor) {
            obsMetaUVEoR = obsMetaUV
        } else {
            obsMetaUVEoR = obsMetaUV
                .filter { _o, meta, _vis -> (meta.eorband != null && meta.eorfield != null) }
        }

        // ps_metrics
        if (params.nopsmetrics) {
            subobsMetaUVPass_ = obsMetaUVEoR
            fail_codes = channel.empty()
        } else {
            obsMetaUVEoR | psMetrics

            // collect psMetrics as a .dat
            psMetrics.out
                // read the content of each ps_metrics file including the trailing newline
                .map { _o, _meta, dat -> dat.getText() }
                .collectFile(
                    name: "ps_metrics.dat",
                    storeDir: "${results_dir()}${params.cal_suffix}"
                )
                // display output path and number of lines
                | view { [it, it.readLines().size()] }

            // collect psMetrics as a .tsv
            psMetrics.out
                // form each row of tsv
                .map { obsid, meta, dat ->
                    def dat_values = dat.getText().split('\n')[0].split(' ')[1..-1]
                    ([
                        obsid,
                        isNaN(meta.lst)?'':meta.lst,
                        isNaN(meta.ew_pointing)?'':meta.ew_pointing,
                        meta.longconfig?:'',
                        meta.name?:'',
                    ] + dat_values).join("\t")
                }
                .collectFile(
                    name: "ps_metrics.tsv", newLine: true, sort: true,
                    seed: [
                        "OBS", "LST", "EWP", "CONF", "CAL NAME", "P_WEDGE", "NUM_CELLS", "P_WINDOW", "NUM_CELLS",
                        "P_ALL", "D3"
                    ].join("\t"),
                    storeDir: "${results_dir()}${params.cal_suffix}"
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
                    [obsid, mapMerge(meta, newMeta)]
                }

            passReasons = psMeta.map { obsid, meta -> [obsid, meta.subobs?:'', meta]}
                // get cmt reduced metrics failures
                .groupTuple(by: 0..1)
                .flatMap { obsid, _subobs, metas ->
                    def nosubMeta = metas.find { it.sub == null} ?: [:]
                    def (nosub_fail_code, nosubReason) = cmt_ps_metrics_pass(nosubMeta)
                    def subMetas = metas.findAll { it.sub != null } ?: []
                    def subReasons = subMetas.collect { subMeta -> cmt_ps_metrics_pass_sub(nosubMeta, subMeta) }
                    def (asub_fail_code, asubReason) = subReasons.find { fail_code, _reason -> fail_code != 0x00 }?:[null, null]
                    if (nosub_fail_code == 0x00 && asubReason != null) {
                        nosub_fail_code = asub_fail_code
                        nosubReason = asubReason
                    }
                    return [
                        [obsid, nosubMeta, getFailReason(nosub_fail_code), nosubReason]
                    ] + [subMetas, subReasons].transpose().collect { subMeta, subReason ->
                        def (sub_fail_code, sub_reason) = subReason
                        [obsid, subMeta, getFailReason(sub_fail_code), sub_reason?:nosubReason]
                    }
                }

            passReasons
                .map { obsid, meta, fail_code, reason ->
                    [
                        obsid,
                        meta.name,
                        meta.subobs?:'',
                        meta.p_window,
                        meta.p_wedge,
                        meta.num_cells,
                        fail_code,
                        reason
                    ].join("\t")
                }
                .collectFile(
                    name: "reasons_ps.tsv", newLine: true, sort: true,
                    seed: [
                        "OBSID",
                        "NAME",
                        "SUBOBS",
                        "P_WINDOW", "P_WEDGE", "NUM_CELLS",
                        "FAIL CODE", "REASON"
                    ].join("\t"),
                    storeDir: "${results_dir()}${params.cal_suffix}"
                )
                // display output path and number of lines
                | view { [it, it.readLines().size()] }


            fail_codes = passReasons.map { obsid, _m, fail_code, _reason ->
                    [
                        obsid,
                        fail_code
                    ]
                }
                .groupTuple(by: 0)
                .map { obsid, obs_fail_codes ->
                    def aFailCode = coerceList(firstFail([ coerceList(obs_fail_codes) ]))[0]
                    [obsid, aFailCode]
                }

            fail_codes.groupTuple(by: 1)
                .map { obsids, fail_code ->
                    [
                        fail_code,
                        obsids.size(),
                    ].join("\t")
                }
                .collectFile(
                    name: "fail_counts_ps.tsv", newLine: true, sort: true,
                    seed: ([ "FAIL CODE", "COUNT" ]).join("\t"),
                    storeDir: "${results_dir()}${params.cal_suffix}"
                )
                | view { it.readLines().size() }

            obsMetaPass = passReasons
                .filter { _o, _m, _fail_code, reason -> reason == null }
                .map { obsid, meta, _fail_code, _reason -> [obsid, meta] }

            subobsMetaUVPass_ = obsMetaUV.cross(obsMetaPass) { def (obsid, meta) = it; [obsid, meta.name, meta.subobs?:''] }
                .map { obsMetaUV_ , obsMetaPass_ ->
                    def (obsid, _m, vis) = obsMetaUV_
                    def (_o, meta) = obsMetaPass_
                    [obsid, meta, vis]
                }
        }

        // delay spectrum
        obsMetaUVEoR | delaySpec

    emit:
        // channel of files to archive, and their buckets
        archive = visQA.out.map { _o, _m, json -> ["visqa", json]}
        // channel of video name and frames to convert
        frame = plotVisQA.out
            .flatMap { _o, meta, pngs ->
                coerceList(pngs).collect { png ->
                    def suffix = png.baseName.split('_')[-1]
                    [''+ "visqa_${meta.name}_${suffix}", png]
                }
            }
            .mix(uvPlot.out.flatMap {_o, name, pngs -> coerceList(pngs).collect { png ->
                def pol = png.baseName.split('_')[-2..-1].join('_')
                [''+"visqa_${name}_${pol}", png]
            }})
            .mix(delaySpec.out.map { _o, meta, png -> ["visqa_dlyspec_${meta.name}", png] })
            .groupTuple()
        // channel of files to zip
        zip = visQA.out.map { _o, meta, json -> [''+"visqa_${meta.name}", json] }
            .groupTuple()
        obsMetaUVPass = subobsMetaUVPass_
        fail_codes = fail_codes.filter { _o, fail_code -> fail_code != getFailReason(0x00) }
}


// img considerations
//
// ## pol
// Default: 'I'. Possible values: XX, XY, YX, YY, I, Q, U, V, RR, RL, LR or LL (case insensitive).
// It is allowed but not necessary to separate with commas, e.g.: 'xx,xy,yx,yy'.
// Two or four polarizations can be joinedly cleaned (see '-joinpolarizations'), but
// this is not the default. I, Q, U and V polarizations will be directly calculated from
// the visibilities, which might require correction to get to real IQUV values. The
// 'xy' polarization will output both a real and an imaginary image, which allows calculating
// true Stokes polarizations for those telescopes.
//
// ## join
//
//
// ## link
// requires wsclean 2.6
//
// ## continuing deconv <https://wsclean.readthedocs.io/en/latest/continue_deconvolution.html>
// The constructed model visibilities need to be stored in the measurement set.
// This implies that -mgain should be used (see the Selfcal instructions for more info on mgain) in the first run, or you have to manually predict the model image from the first run, before continuing.
// WSClean does not verify whether this assumption holds.
// As a result, it is not directly possible to use -no-update-model-required in the first run, because the second run requires the MODEL_DATA column to be filled. If -no-update-model-required was enabled, or only a model image without the corresponding predicted visibilities is available, it is still possible to continue the run by first predicting the model data (using wsclean -predict ...) from the model image before continuing.
//
// ## idg <https://wsclean.readthedocs.io/en/latest/image_domain_gridding.html>
//



// image visibilities and QA images
workflow img {
    take:
        // tuple of (obsid, meta, vis)
        obsMetaVis
    main:

        // def polModes = [
        //     "I": [glob: "-I", pol: "I", idg: true], // Total intensity science, i.e., not interested in QUV: only image I with -pol I.
        //     "Q": [glob: "-Q", pol: "Q", idg: true],
        //     "U": [glob: "-U", pol: "U", idg: true],
        //     "V": [glob: "-V", pol: "V", idg: true],
        //     "IQUV": [glob: "-{I,Q,U,V}", pol: "IQUV", idg: true], // Total intensity science, but “nice to have QUV” for e.g. sensivity analysis
        //     "IQUV_join": [glob: "-{I,Q,U,V}", pol: "IQUV", join_pols: true, idg: true], // Interested in all stokes parameter, cleaning each polarization in a joined way
        //     "IV_join": [glob: "-{I,V}", pol: "i,v", join_pols: true, ],
        //     "IQUV_lqu": [glob: "-{Q,U}", pol: "IQUV", link_pols: "q,u", idg: true ], // Interested in rotation measure synthesis
        //     "QU_sq": [glob: "-{Q,U}", pol: "q,u", join_pols: true, sq_chan_join: true ],
        //     "XXYY_join": [glob: "-{XX,YY}", pol: "xx,yy", join_pols: true, ],
        //     "XXXYXX_lxxyy": [glob: "-{XX,YY,XY,XYi}", pol: "xx,xy,yx,yy", link_pols: "xx,yy"],
        // ]

        obsMetaVis.map { obsid, meta, vis ->
                def newMeta = [
                    argstr:"--no-diff --crosses --combine-freq --pix 501",
                    plot_base:".cross"
                ]
                [obsid, mapMerge(meta, newMeta), vis]
            }
            | demo11_allsky

        // wsclean: make deconvolved images
        if (params.img_split_intervals) {
            // add -t???? suffix to name, to split wsclean over multiple jobs
            splitObsMetaVis = obsMetaVis.flatMap { obsid, meta, vis ->
                (0..(meta.ntimes-1)).collect { i ->
                    [obsid, mapMerge(meta, [interval: [i, i+1], inter_tok: String.format("-t%04d", i)]), vis]
                }
            }
        } else {
            splitObsMetaVis = obsMetaVis
        }

        // print("params.nodeconv: ${params.nodeconv}")
        if (params.nodeconv) {
            splitObsMetaVis.map {obsid, meta, vis ->
                    def imgParams = [ img_name: ''+"${meta.cal_prog}_${obsid}${meta.subobs?:''}_${meta.name}${meta.inter_tok?:''}" ]
                    [obsid, meta, vis, mapMerge(wscleanParams(), imgParams)]
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
            //         dirtyImgs = imgs.findAll {img -> img.baseName.split('-')[-1] == 'dirty'}
            //         newImgParams = mapMerge(imgParams, wscleanDConvParams())
            //         [obsid, meta, vis, newImgParams, dirtyImgs]
            //     }

            splitObsMetaVis.map { obsid, meta, vis ->
                    def imgParams = [ img_name: ''+"${meta.cal_prog}_${obsid}${meta.subobs?:''}_${meta.name}${meta.inter_tok?:''}" ]
                    [obsid, meta, vis, mapMerge(wscleanDConvParams(), imgParams), []]
                }
                .flatMap {obsid, meta, vis, imgParams_, dirtyImgs ->
                    [
                        // ["", "i"],
                        ["-{XX,YY}", "xx,yy -join-polarizations"],
                        // "IQUV_join": [glob: "-{I,Q,U,V}", pol: "IQUV", join_pols: true, idg: true], // Interested in all stokes parameter, cleaning each polarization in a joined way
                        // ["-{I,Q,U,V}", "IQUV -join-polarizations -grid-with-beam -use-idg -idg-mode cpu -pb-undersampling 4 -mwa-path /astro/mwaeor/jline/software/"], // Interested in all stokes parameter, cleaning each polarization in a joined way
                        // "QU_sq": [glob: "-{Q,U}", pol: "q,u", join_pols: true, sq_chan_join: true ],
                        // ["-{Q,U}", "q,u -join-polarizations -squared-channel-joining "],
                        ["-{I,V}", "i,v -join-polarizations"],
                        // ["-{XX,YY,XY,XYi}", "xx,xy,yx,yy -link-polarizations xx,yy"],
                        // todo: qu mostly just times out
                        // ["-{Q,U}", "q,u -join-polarizations -squared-channel-joining"],
                    ].collect { polGlob, polArg ->
                        def imgParams = [
                            pol: polArg,
                            glob: "wsclean_${imgParams_.img_name}",
                        ]
                        if (is_multiinterval() && !meta.inter_tok) {
                            imgParams.glob += "-t????"
                        }
                        if (is_multichannel()) {
                            imgParams.glob += "-MFS"
                        }
                        imgParams.glob += polGlob
                        imgParams.glob += "-image.fits"
                        [obsid, meta, vis, mapMerge(imgParams_, imgParams), dirtyImgs]
                    }
                }
                | wscleanDConv
        }

        obsMetaImg = wscleanDirty.out.mix(wscleanDConv.out)
            .flatMap { obsid, meta, imgs ->
                coerceList(imgs).collect { img ->
                    [obsid, mapMerge(meta, decomposeImg(img)), img] }}
            .branch { _o, meta, _img ->
                // target product is image unless deconv is disabled, separate MFS
                imgMfs: meta.prod == (params.nodeconv ? "dirty" : "image") && meta.chan == -1
                imgNoMfs: meta.prod == (params.nodeconv ? "dirty" : "image")
                // target product is uv grids
                gridMfs: meta.prod ==~ /uv-.*/ && meta.chan == -1
                gridNoMfs: meta.prod ==~ /uv-.*/
                // target product is psf
                psfMfs: meta.prod == "psf" && meta.chan == -1
                psfNoMfs: meta.prod == "psf"
                // other potential products: residual, model
            }

        if (!params.nopolcomp || params.thumbnail_limits) {
            // calculate quantiles (what values are at nth percentile)
            obsMetaImg.imgMfs \
                .mix(obsMetaImg.psfMfs)
                // .mix(obsMetaGrid) \
                | imgQuantiles

            // limits are used to set the color scale of each type of image.
            imgLimits = imgQuantiles.out
                .map { obsid, meta, hist ->
                    def high = parseCsv2(hist, true, 2, ',')
                        .find { row -> Float.compare(parseFloatOrNaN(row.quantile), params.thumbnail_quantile as Float) == 0 }
                    high = parseFloatOrNaN(high == null ? null : high.value)
                    def low = parseCsv2(hist, true, 2, ',')
                        .find { row -> Float.compare(parseFloatOrNaN(row.quantile), (1-params.thumbnail_quantile) as Float) == 0 }
                    low = parseFloatOrNaN(low == null ? null : low.value)
                    // subobs = meta.subobs?:''
                    [obsid, meta.inter_tok, meta.name, meta.inter_suffix, high, low]
                }
                .filter { _o, _interval, _name, _suff, high, low -> !(isNaN(high) || isNaN(low)) }
                // get max limit for each obs, interval, name, suffix (deals with multiple channels)
                .groupTuple(by: 0..3)
                .map { obsid, interval, name, suff, highs, lows ->
                    def high = coerceList(highs).max()
                    def low = coerceList(lows).max()
                    [obsid, interval, name, suff, high, low]
                }
                // for each obs, interval, name, produce a mapping from suff to limit
                .groupTuple(by: 0..2)
                .map { obsid, interval, name, suffs, highs, lows ->
                    def suffLimits_ = [
                            coerceList(suffs),
                            coerceList(highs),
                            coerceList(lows),
                        ]
                        .transpose()
                        .collect { suff, high, low -> [suff, [high, low]]}
                        .collectEntries()
                    [obsid, interval, name, suffLimits_]
                }
                // .view { "imgLimits ${it}" }

            obsMetaImgMfsPass = obsMetaImg.imgMfs.map { obsid, meta, img ->
                    [obsid, meta.inter_tok, meta.name, meta.inter_suffix, meta, img]
                }
                .cross(imgLimits.flatMap { obsid, interval, name, suffLimits_ ->
                    suffLimits_.collect { k, _v -> [ obsid, interval, name, k ]}}
                ) { it[0..3] }
                .map { obsMetaImgMfs_, _imgLimits ->
                    def (obsid, _intv, _name, _ints, meta, img) = obsMetaImgMfs_
                    [obsid, meta, img]
                }
                // .view { "obsMetaImgMfs ${it}" }

            // channel of all suffixes in
            suffs_ = obsMetaImg.imgMfs
                .map { _o, meta, _img -> meta.inter_suffix }
                .unique()
                .toSortedList()
                .map{it -> [it]}

            // write img limits to file
            imgLimits.combine(suffs_)
                .map { obsid, interval, name, limits, suffs ->
                        ([obsid, interval, name] + suffs.collect { suff ->
                                def hilo = ['', '']
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
                    storeDir: "${results_dir()}${params.img_suffix}${params.cal_suffix}"
                )
                // display output path and number of lines
                | view { [it, it.readLines().size()] }

            // value channel containing a map from img suffix to max limit excluding outliers
            // suffLimits = channel.from([:])
            suffLimits = imgLimits.combine(suffs_)
                .flatMap{ _o, _interval, name, limits, suffs ->
                    suffs.collect { suff ->
                        def hilo = limits[suff]?: [Float.NaN, Float.NaN]
                        [new Tuple(name, suff)] + hilo
                    }
                }
                .groupTuple(by: 0)
                .map { group, highs_, lows_ ->
                    def thumbnail_vmax = params.thumbnail_vmax?: {
                        def highs = coerceList(highs_)
                            .findAll { !isNaN(it) }
                            .sort(false)
                        def high_index = ((highs.size() - 1) * params.thumbnail_quantile).round() as Integer
                        highs[high_index]
                    }()

                    def thumbnail_vmin = params.thumbnail_vmin ?: {
                        def lows = coerceList(lows_)
                            .findAll { !isNaN(it) }
                            .sort(false)
                        def low_index = (lows.size() * (1-params.thumbnail_quantile)).round() as Integer
                        lows[low_index]
                    }()

                    [group, [thumbnail_vmax.toFloat(), thumbnail_vmin.toFloat()]]
                }
                .toList()
                .map { it.collectEntries() }
        } else {
            suffLimits = channel.from([:])
            obsMetaImgMfsPass = obsMetaImg.imgMfs
            obsMetaGrid = obsMetaImg.gridMfs
            obsMetaPsf = obsMetaImg.psfMfs
        }

        // suffLimits.view { "suffLimits \n${it.collect().join('\n')}"}

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
                .filter { _o, _m, _img, high, low ->
                    !(isNaN(high) || isNaN(low))
                }
                // group (pol, img, limit) by (obs, interval, name, prod).
                .map { obsid, meta, img, high, low ->
                    [obsid, meta.inter_tok, meta.name, meta.prod, meta.pol, img, high, low]
                }
                .groupTuple(by: 0..3)
                // make a hashmap of pol -> (img, limit) for each (obs, interval, name, prod)
                .map { obsid, interval, name, prod, pols, imgs, highs, lows ->
                    def polImgLimits = [
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
                            def (imgs, highs, _lows)  = order.collect { pol -> polImgLimits[pol] }.transpose()
                            def polcompMeta = [
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
                .mix(params.thumbnail_psfs ?obsMetaPsf : channel.empty())
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
                        def ratio = 0.1
                        if (newMeta.vmin as Float < 0 && newMeta.vmax as Float > 0) {
                            ratio = (-newMeta.vmin) / (newMeta.vmax - newMeta.vmin)
                        }
                        // newMeta.norm_args = '{\\\\"stretch\\\\":\\\\"asinh\\\\",\\\\"asinh_a\\\\":'+ "${ratio}" + '}'
                        newMeta.norm_args = '{\\\\"stretch\\\\":\\\\"power\\\\",\\\\"power\\\\":2,\\\\"clip\\\\":true}'
                    }
                    [obsid, mapMerge(meta, newMeta), img]
                }
                // .view { "thumbnail ${it}" }
                | thumbnail
        }

        // montage of polarizations
        if (!params.nomontage) {
            thumbnail.out.flatMap { obsid, meta, png ->
                // visibility name, without sub*
                def newMeta = mapMerge(meta, [vis_name: meta.name, sub: meta.sub?:'nosub'])
                def m = newMeta.vis_name =~ /(sub|ionosub)_(.*)/
                if (m) {
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
                // montage_name = [sub: suffix, pol: sub_name].get(params.montage_by)
                def subobs = "${newMeta.subobs?:''}${newMeta.inter_tok?:''}".toString()
                def montages = []
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

        // if (params.krvis) {
        //     obsMetaImgMfsPass
        //         .map { obsid, meta, img -> [obsid, meta.name, meta.inter_tok?:'', meta, img] }
        //         .groupTuple(by: 0..2)
        //         .filter { _o, _name, _intt, metas, _imgs -> metas.find { it.pol == 'I' } != null }
        //         .map { obsid, _name, _intt, metas, imgs ->
        //             def iMeta = metas.find { it.pol == 'I' }
        //             def (_iquvMetas, iquvImgs) = [metas, imgs].transpose().findAll { meta, _img ->
        //                 ['I', 'Q', 'U', 'V'].contains(meta.pol)
        //             }.transpose()
        //             [obsid, iMeta, iquvImgs, file('/pawsey/mwa/mwaeor/dev/telescope_data_visualisation/romaO.npy')]
        //         }
        //         | krVis
        // } else {
        //     channel.empty() | krVis
        // }

        // each vis name can have multiple images in obsMetaImgMfs
        // group by obsid, vis name using original meta from obsMetaVis
        obsMetaImgGroup = obsMetaVis.cross(obsMetaImgMfsPass) { def (obsid, meta) = it; [obsid, meta.name] }
            .map { obsMetaVis_, obsMetaImgMfs_ ->
                def (obsid, meta) = obsMetaVis_
                def (_o, imgMeta, img) = obsMetaImgMfs_
                [obsid, meta.name, meta, imgMeta, img]
            }
            .filter { _o, _name, _m, imgMeta, _img -> ['XX', 'YY', 'V'].contains(imgMeta.pol)}
            .groupTuple(by: 0..1)
            .map { obsid, _name, metas, imgMetas, imgs -> [obsid, metas[0], imgMetas, imgs] }
            // | view { def (obsid, meta) = it; "obsMetaImgGroup ${obsid}, ${meta.name}"}

        // imgQA for MFS images and groups of images
        //  - imgs need to be deconvolved for QA
        //  - can't handle multiple intervals
        if (params.noimgqa || params.nodeconv || params.img_split_intervals) {
            channel.empty() | imgQA
            obsMetaImgPass_ = obsMetaImgGroup
            fail_codes = channel.empty()
        } else {
            obsMetaImgGroup
                .map { obsid, meta, _imgMetas, imgs -> [obsid, meta, imgs] }
                .filter { _o, meta, imgs ->
                    // imgQA looks at pks flux, which needs to be deconvolved, only works with eor0
                    // meta.prod == 'image' &&
                    meta.eorfield == 0 &&
                    imgs.size() == 3
                }
                // .map { obsid, meta, img -> [obsid, meta.name, meta, img] }
                // .groupTuple(by: 0..1)
                // .map {obsid, name, metas, imgs -> [obsid, metas[0], imgs]}
                | imgQA

            imgQA.out.map { _o, meta, json -> [meta.name, json] }
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
                        meta.ew_pointing,
                        isNaN(meta.lst)?'':meta.lst,
                        meta.longconfig?:'',
                        meta.name,
                        xx.RMS_ALL?:"", xx.RMS_BOX?:"", xx_pks.PEAK_FLUX?:"", xx_pks.INT_FLUX?:"",
                        yy.RMS_ALL?:"", yy.RMS_BOX?:"", yy_pks.PEAK_FLUX?:"", yy_pks.INT_FLUX?:"",
                        v.RMS_ALL?:"", v.RMS_BOX?:"", v_pks.PEAK_FLUX?:"", v_pks.INT_FLUX?:"",
                        // stats.V_XX?.RMS_RATIO, stats.V_XX?.RMS_RATIO_BOX
                    ].join("\t")
                }
                .collectFile(
                    name: "img_metrics.tsv", newLine: true, sort: true,
                    seed: [
                        "OBS", "EWP", "LST", "CONF", "IMG NAME",
                        "XX RMS ALL", "XX RMS BOX", "XX PKS0023_026 PEAK", "XX PKS0023_026 INT",
                        "YY RMS ALL", "YY RMS BOX", "YY PKS0023_026 PEAK", "YY PKS0023_026 INT",
                        "V RMS ALL","V RMS BOX", "V PKS0023_026 PEAK", "V PKS0023_026 INT" ,
                        // "V:XX RMS RATIO", "V:XX RMS RATIO BOX"
                    ].join("\t"),
                    storeDir: "${results_dir()}${params.img_suffix}${params.cal_suffix}"
                )
                // display output path and number of lines
                | view { [it, it.readLines().size()] }

            imgQA.out
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
                    [obsid, mapMerge(meta, newMeta)]
                }
                .groupTuple(by: 0)
                .map { obsid, metas -> [obsid, coerceList(metas)] }
                .flatMap { obsid, metas ->
                    def nosubMeta = metas.find { meta -> meta.sub == null }
                    if (nosubMeta == null) {
                        return []
                    }
                    def (nosub_fail_code, nosubReason) = cmt_imgqa_pass(nosubMeta)
                    def subMetas = metas.findAll { it != null && it.sub?:"" != "" } ?: []
                    def subReasons = subMetas.collect { subMeta -> cmt_imgqa_pass_sub(nosubMeta, subMeta) }
                    def (asub_fail_code, asubReason) = subReasons.find { fail_code, _reason -> fail_code != 0x00 }?:[null, null]
                    if (nosub_fail_code == 0x00 && asubReason != null) {
                        nosub_fail_code = asub_fail_code
                        nosubReason = asubReason
                    }
                    return [
                        [obsid, nosubMeta, getFailReason(nosub_fail_code), nosubReason]
                    ] + [subMetas, subReasons].transpose().collect { subMeta, subReason ->
                        def (sub_fail_code, sub_reason) = subReason
                        [obsid, subMeta, getFailReason(sub_fail_code), sub_reason?:nosubReason]
                    }
                }
                .tap { obsMetaReasons }

            obsMetaReasons
                .map { obsid, _m, fail_code, _reason ->
                    [ obsid, fail_code ]
                }
                .groupTuple(by: 0)
                .map { obsid, obs_fail_codes ->
                    def aFailCode = coerceList(firstFail([ coerceList(obs_fail_codes) ]))[0]
                    [obsid, aFailCode]
                }
                .tap { fail_codes }

            fail_codes.groupTuple(by: 1)
                .map { obsids, fail_code ->
                    [
                        fail_code,
                        obsids.size(),
                    ].join("\t")
                }
                .collectFile(
                    name: "fail_counts_imgqa.tsv", newLine: true, sort: true,
                    seed: ([ "FAIL CODE", "COUNT" ]).join("\t"),
                    storeDir: "${results_dir()}${params.img_suffix}${params.cal_suffix}"
                )
                | view { it.readLines().size() }

            obsMetaReasons
                // .filter { _, __, fail_code, ___ -> fail_code != getFailCode(0x00) }
                .groupTuple( by: 0 )
                .map { obsid, metas, obsCodes_, reasons ->
                    def (failCode, reason, meta) = firstFail([
                        coerceList(obsCodes_),
                        coerceList(reasons),
                        coerceList(metas),
                    ])
                    [
                        obsid,
                        meta.name,
                        meta.xx_pks_int?:"",
                        meta.yy_pks_int?:"",
                        meta.v_pks_int?:"",
                        meta.v_rms_box?:"",
                        failCode,
                        reason?:''
                    ].join("\t")
                }
                .collectFile(
                    name: "reasons_imgqa.tsv", newLine: true, sort: true,
                    seed: [
                        "OBSID", "NAME",
                        "XX PKS INT", "YY PKS INT", "V PKS INT", "V RMS BOX",
                        "FAIL CODE", "REASON"
                    ].join("\t"),
                    storeDir: "${results_dir()}${params.img_suffix}${params.cal_suffix}"
                )
                // display output path and number of lines
                | view { [it, it.readLines().size()] }

            obsMetaPass = obsMetaReasons
                .filter { _o, _meta, _fail_code, reason -> reason == null }
                .map { obsid, meta, _fail_code, _reason -> [obsid, meta] }

            obsMetaImgPass_ = obsMetaImgGroup.cross(obsMetaPass) { def (obsid, meta) = it; [obsid, meta.name] }
                .map { obsMetaImgGroup_, obsMetaPass_ ->
                    def (obsid, _meta, imgMetas, imgs) = obsMetaImgGroup_
                    def (_o, meta) = obsMetaPass_
                    [obsid, meta, imgMetas, imgs]
                }
        }

    emit:
        // channel of files to archive, and their buckets
        archive = obsMetaImgMfsPass.filter { obsid, _meta, _img -> !obsid.startsWith("e") }
            .transpose()
            .map { _o, __, img -> ["img", img]}
            .mix( imgQA.out.map { _o, _m, json -> ["imgqa", json]} )
        // channel of video name and frames to convert
        frame = polComp.out.map { _o, meta, png ->
                [''+"imgqa_${meta.name}_polcomp_${meta.prod}_${meta.orderName}", png]
            }
            .mix(polMontage.out.map { _o, meta, png ->
                [''+"imgqa_${meta.name}_polmontage", png]
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
            // .mix(krVis.out.flatMap {_, meta, pngs ->
            //     pngs.collect { png ->
            //         def suffix = png.baseName.split('-')[-2..-1].join('-');
            //         ["krvis_${meta.name}_${suffix}", png]
            //     }
            // })
            .groupTuple()
        // channel of files to zip
        zip = imgQA.out.map { _o, meta, json -> ["imgqa_${meta.name}", json] }
            .mix(obsMetaImgMfsPass.map { _o, meta, img -> ["img_${meta.name}", img] })
            .groupTuple()
        // channel of tuple (obsid, imgMeta, img) that pass qa
        obsMetaImgPass = obsMetaImgPass_
        fail_codes = fail_codes
}

workflow imgCombine {
    take:
        // tuple of (obsid, meta, uvfits)
        obsMetaUV
        // tuple of (chunk, chunkMeta, obsids) for grouping
        chunkMetaObs

    main:

        // if (params.img_split_coarse_chans) {
        //     if (params.sky_chans == null || params.sky_chans.size == 0) {
        //         throw new Exception("img_split_coarse_chans is enabled but params.sky_chans=${params.sky_chans}")
        //     }
        //     nCoarseChans = params.sky_chans.size
        //     obsMetaUV.map { obsid, meta, uvfits ->
        //             if (meta.nchans == null || meta.nchans % nCoarseChans != 0) {
        //                 throw new Exception("nchans ${meta.nchans} is not a multiple of ${nCoarseChans}")
        //             }
        //             meta.nchans
        //         }
        //         .unique()
        //         .collect()
        //         .subscribe { it ->
        //             if (it.size() > 1) {
        //                 throw new Exception("obsMetaUV has multiple nchans ${it}")
        //             }
        //         }

        //     // splitObsMetaUV = obsMetaUV.flatMap { obsid, meta, uvfits ->
        //     //         params.sky_chans.withIndex().collect { ch, idx ->
        //     //             fineChansPerCoarse = meta.nchans / nCoarseChans
        //     //             def newMeta = deepcopy(meta)
        //     //             newMeta.chan = [idx * fineChansPerCoarse, (idx+1) * fineChansPerCoarse]
        //     //             newMeta.chan_tok = String.format("-ch%03d", ch)
        //     //             // TODO: newMeta.name = "${meta.name}${meta.chan_tok}"
        //     //             [obsid, newMeta, uvfits]
        //     //         }
        //     //     }
        // }
        // // else {
        // //     splitObsMetaUV = obsMetaUV
        // // }

        // combined images
        groupMetaVisImgParams = obsMetaUV.cross(
                chunkMetaObs.flatMap { chunk, chunkMeta, obsids ->
                    obsids.collect { obsid -> [obsid, chunkMeta, chunk] }
                }
            ) { def (obsid, meta) = it; [obsid, meta.name] }
            .map { obsMetaUV_, chunkMetaObs_ ->
                def (_o, _m, vis) = obsMetaUV_
                def (_o_, chunkMeta, chunk) = chunkMetaObs_
                [chunk, chunkMeta.name, chunkMeta, vis]
            }
            .groupTuple(by: 0..1)
            // .view { it -> "\n -> imgCombine obsMetaUV x chunkMetaObs ${it}\n" }
            .map { chunk, _name, chunkMetas, viss ->
                [chunk, coerceList(chunkMetas)[0], viss.flatten(), deepcopy(wscleanParams())]
            }
            // .view { it -> "\n -> imgCombine before filter ${it}\n" }
            .filter { _chunk, chunkMeta, _viss, _imgParams -> chunkMeta.name }

        if (params.img_split_coarse_chans && !params.noimgcombine) {
            if (params.sky_chans == null || params.sky_chans.size == 0) {
                throw new Exception("img_split_coarse_chans is enabled but params.sky_chans=${params.sky_chans}")
            }
            nCoarseChans = params.sky_chans.size
            print(" -> deleteme nCoarseChans=${nCoarseChans}")
            fineChansPerCoarse = obsMetaUV.map { _o, meta, _uvfits ->
                    if (meta.nchans == null || meta.nchans % nCoarseChans != 0) {
                        throw new Exception("nchans ${meta.nchans} is not a multiple of ${nCoarseChans}")
                    }
                    meta.nchans
                }
                .unique()
                // .collect()
                .toList()
                .map { it ->
                    if (it.size() > 1) {
                        throw new Exception("obsMetaUV has multiple nchans ${it}")
                    }
                    if (it[0] == null) {
                        throw new Exception("nchans is null ${it}")
                    }
                    it[0] / nCoarseChans
                }

            splitGroupMetaVisImgParams = groupMetaVisImgParams.combine(fineChansPerCoarse)
                .flatMap { chunk, chunkMeta, viss, imgParams, fchans ->
                    // println("  -> deleteme fchans=${fchans} chunk=${chunk} chunkMeta=${chunkMeta} imgParams=${imgParams}")
                    params.sky_chans.withIndex().collect { ch, idx ->
                        // println("   -> deleteme ch=${ch} idx=${idx}")
                        def chan_tok = String.format("-ch%03d", ch)
                        def newMeta = [
                            chan: [idx * fchans, (idx+1) * fchans],
                            chan_tok: chan_tok,
                            subobs: ''+"${chunkMeta.subobs?:''}${chan_tok}",
                        ]
                        [chunk, mapMerge(chunkMeta, newMeta), viss, imgParams]
                    }
                }
        } else {
            splitGroupMetaVisImgParams = groupMetaVisImgParams
        }

        if (params.nodeconv) {
            splitGroupMetaVisImgParams | wscleanDirty
            channel.empty() | wscleanDConv
        } else {
            channel.empty() | wscleanDirty
            splitGroupMetaVisImgParams.map { chunk, meta, vis, imgParams ->
                    [chunk, meta, vis, mapMerge(wscleanDConvParams(), imgParams), []]
                }
                | wscleanDConv
        }

        // make thumbnails
        // - only MFS images unless thumbnail_all_chans
        // - exclude dirty if deconv is enabled
        wscleanDConv.out.mix(wscleanDirty.out)
            .flatMap { group, meta, imgs ->
                imgs.collect { img ->
                    [group, mapMerge(meta, decomposeImg(img)), img] }
            }
            .filter { _group, meta, _img ->
                meta.prod !=~ /uv-.*/ \
                && (meta.chan?:-1 == -1 || params.thumbnail_all_chans) \
                && (meta.prod == (params.nodeconv ? "dirty" : "image"))
            }
            | thumbnail

        if (params.img_split_coarse_chans) {
            frame = thumbnail.out
                .filter { _group, meta, _png ->
                    meta.name =~ /ionosub/
                }
                .map { group, meta, png ->
                    [''+"imgCombine_${group}_${meta.name}", png]
                }
            }
        else {
            frame = channel.empty()
        }
    emit:
        frame = frame.groupTuple()
}

// WARN: joins assume single obs per vis
workflow extChips {
    // obsVis = channel.of(file("external.csv"))
    //     .splitCsv(header: ["obsid", "vis"])
    //     .filter { !it.obsid.startsWith('#') }
    //     .map { [it.obsid, file(it.vis)] }
    def baseMeta = [
        name: params.visName ?: "ionosub_30l_src4k_300it_8s_80kHz_i1000",
        cal_prog: params.cal_prog ?: "hyp",
        eorfield: params.eorfield ?: 0,
    ]

    channel.of(obsids_file())
        .splitCsv()
        .filter { line -> !line[0].startsWith('#') }
        .map { line ->
            def (obsid, _comment) = line
            def vis = file("${params.outdir}/${obsid}/cal${params.cal_suffix}/${baseMeta.cal_prog}_${obsid}_${baseMeta.name}.uvfits")
            [ obsid, deepcopy(baseMeta), vis ]
        }
        .filter { obs, __, vis -> obs.size() > 0 && vis.exists() }
        .tap { obsMetaVis }
        // .view { "before uvMeta ${it}"}
        | uvMeta

    obsVis = obsMetaVis.map { obs, _meta, vis -> [obs, vis] }
    obsVis.join(uvMeta.out)
        .map { obs, vis, meta_, metaJson ->
            def uvmeta = parseJson(metaJson)
            def newMeta = [
                lowfreq: uvmeta.freqs[0],
                freq_res_hz: (uvmeta.freqs[1] - uvmeta.freqs[0]),
                nchans: (uvmeta.freqs?:[]).size(),
                ntimes: (uvmeta.times?:[]).size(),
                int_time: (uvmeta.times[1].gps - uvmeta.times[0].gps),
            ]
            ['eorband', 'num_ants', 'total_weight'].each { key ->
                if (uvmeta[key] != null) {
                    newMeta[key] = uvmeta[key]
                }
            }
            if (uvmeta['config'] != null) {
                newMeta['shortconfig'] = uvmeta['config']
            }
            if (meta_.ntimes > 0) {
                newMeta.lst = Math.toDegrees(uvmeta.times[0].lst_rad)
            }
            [obs, mapMerge(meta_, newMeta), vis]
        }
        // .view { "before psMetrics ${it}" }
        | psMetrics

    obsVis.join(psMetrics.out)
        .map { obs, vis, meta_, dat ->
            def (p_wedge, num_cells, p_window) = dat.getText().split('\n')[0].split(' ')[1..-1]
            def newMeta = [
                p_wedge: p_wedge,
                num_cells: num_cells,
                p_window: p_window,
            ]
            [
                "nchans", "eorband", "eorfield", "lowfreq", "freq_res_hz", "int_time"
            ].each { key ->
                if (meta_[key] == null) {
                    throw new Exception("obs=${obs} meta.${key} is null")
                }
            }
            [obs, mapMerge(meta_, newMeta), vis]
        }
        .tap { obsMetaPSVis }
        .map { obs, meta, vis ->
            def gmeta = groupMeta(meta)
            // print("obs=${obs} group=${gmeta.group}")
            [''+"${gmeta.group}", gmeta.sort, obs, gmeta, vis]
        }
        // .view { "before grouping ${it}" }
        .groupTuple(by: 0)
        // only keep groups with more than one obsid
        .filter { it -> it[1] instanceof List }
        // .view { it -> "\n -> grouped: ${it}"}
        // TODO: assumes one vis per obsid
        .flatMap { group, sorts, obss, gmetas, _viss ->
            [sorts, obss, gmetas].transpose()
                .sort { it -> it[0] }
                .collate(params.chunkSize, params.chunkRemainder)
                .take(params.chunkCount)
                .collect { chunk ->
                    def (chunkSorts_, chunkObss, chunkMetas) = chunk.transpose()
                    def meta = coerceList(chunkMetas)[0]
                    def chunkSorts = coerceList(chunkSorts_)
                    // meta.obsids = chunkObss.sort(false)
                    def obsids = coerceList(chunkObss).sort(false)
                    def newMeta = mapMerge(meta, [
                        nobs: obsids.size(),
                        hash: obsids.join(' ').md5()[0..7],
                        sort_bounds: [chunkSorts[0], chunkSorts[-1]],
                    ])
                    ["${group}_${newMeta.hash}", newMeta, obsids]
                }
        }
        // .view { it -> "\n -> chunked by obsid: ${it}"}
        .set { chunkMetaObs }

    chips(obsMetaPSVis, chunkMetaObs)
    // if (params.nochipscombine) {
    //     chips(obsMetaPSVis, channel.empty())
    // } else {
    //     chips(obsMetaPSVis, chunkMetaObs)
    // }
}

// make combined grids and images using chips and wsclean
workflow chips {
    take:
        // tuple of (obsid, meta, uvfits)
        obsMetaUV
        // tuple of (chunk, chunkMeta, obsids) for grouping
        chunkMetaObs

    main:
        pols = channel.of('xx', 'yy')

        // power spectrum
        obsMetaUV
            .combine(pols)
            .map { obsid, meta, viss, pol ->
                [obsid, mapMerge(meta, [nobs: 1, ext: ''+"${obsid}_${meta.name}", pol: pol]), [obsid], viss]
            }
            | chipsGrid

        chipsGrid.out
            .map { obsid, meta, grid ->
                [obsid, meta.name?:'', meta, [grid]]
            }
            .groupTuple(by: 0..1)
            .map { obsid, _name, metas, grids ->
                [obsid, coerceList(metas)[0], grids.flatten()]
            }
            .cross(
                chunkMetaObs.flatMap { chunk, chunkMeta, obsids ->
                    obsids.collect { obsid -> [obsid, deepcopy(chunkMeta), chunk] }
                }
            ) { def (obsid, meta) = it; [obsid, meta.name] }
            .map { chipsGrid_, chunkMetaObs_ ->
                def (_o, obsMeta, grid) = chipsGrid_
                def (__, chunkMeta, chunk) = chunkMetaObs_
                [chunk, chunkMeta.name, chunkMeta, obsMeta.ext, grid]
            }
            .groupTuple(by: 0..1)
            // .view { it -> "\n -> chipsGrid x chunkMetaObs ${it}\n" }
            .map { chunk, _name, chunkMetas, exts, grids ->
                def chunkMeta = chunkMetas[0]
                def newMeta = [ ext: ''+"${chunk}_${chunkMeta.name}"]
                // print("\n -> chunkMetas 0: ${chunkMeta}\n")
                [chunk, mapMerge(chunkMeta, newMeta), exts, grids.flatten()]
            }
            // .view { it -> "\n -> before chipsCombine ${it}\n" }
            .filter { _chunk, chunkMeta, _exts, _grid -> chunkMeta.name }
            .combine(pols)
            .map { chunk, chunkMeta, exts, grid, pol ->
                [chunk, mapMerge(chunkMeta, [pol:pol]), exts, grid.findAll { it.baseName.contains(pol) }]
            }
            .filter { _chunk, _chunkMeta, _exts, grids -> grids.size() > 0 }
            | chipsCombine

        if (params.lssaObs) {
            lssaExtra = chipsGrid.out.map { obsid, meta, grid ->
                    [obsid, mapMerge(meta, [sort_bounds: [0, 0]]), [grid]]
                }
        } else {
            lssaExtra = channel.empty()
        }

        chipsCombine.out.map { group, meta, _txt, dats -> [ group, meta, dats ] }
            .mix(lssaExtra)
            | chipsLssa
            | chips1d_tsv

        chips1d_tsv.out
            // .map { println("chips1d_tsv out ${it[0]}") }

        all_k_modes = chips1d_tsv.out.flatMap { _chunk, _meta, tsv ->
                tsv.readLines().collect {
                    // try { it.split("\t")[0].toFloat() } catch(NumberFormatException e) { null }
                    it.split("\t")[0]
                }
                .findAll {
                    try {
                        it.split("\t")[0].toFloat()
                    } catch(NumberFormatException _e) {
                        return false
                    }
                    return it != null
                }
            }
            .unique()
            .collect()
            .toList()
            .view()

        // combine all tsv files into one
        channel.of("SORT_LEFT", "SORT_RIGHT", "GRID", "NAME", "NOBS", "POL").concat(all_k_modes.flatten())
            .toList()
            .map { it.join("\t") }
            .concat(
                chips1d_tsv.out.combine(all_k_modes).map { group, meta, tsv, k_modes ->
                        def table = parseCsv2(coerceList(tsv)[0], true, 0, '\t')
                        def (deltas, _noises) = k_modes.collect { k_mode ->
                            def row = table.find { it.k_modes == k_mode }
                            if (row == null) {
                                return [Float.NaN, Float.NaN]
                            }
                            [parseFloatOrNaN(row.delta), parseFloatOrNaN(row.noise)]
                        }.transpose()
                        def (sort_left, sort_right) = meta.sort_bounds?:[Float.NaN, Float.NaN]
                        ([
                            sort_left, sort_right, group, meta.name, meta.nobs, meta.pol
                        ] + deltas).join('\t')
                    }
                    .collect(sort: true)
                    .flatMap { it -> it }
            )
            .collectFile(
                name: ''+"chips1d_delta_${params.lssa_bin}${params.lssa_bias_mode?:0}${params.obsids_suffix}" + (params.visName?".${params.visName}":"") + "${params.result_suffix}.tsv",
                newLine: true, sort: false,
                storeDir: "${results_dir()}${params.img_suffix}${params.cal_suffix}"
            )
            .tap { chips1d_delta_tsv }
            | view { [it, it.readLines().size()] }

        // chips1d_delta_tsv.out.map { }

        singles1D = chipsLssa.out.map { chunk, meta, lssa ->
                def newMeta = [
                    ptype: '1D',
                    // pol: 'both',
                    title: ''+"crosspower\\n${chunk}\\n${meta.name}",
                    plot_name: 'chips1d',
                    max_power: params.chips_max_power,
                    min_power: params.chips_min_power,
                    tags: [],
                ]
                [chunk, mapMerge(meta, newMeta), lssa]
            }

        singles1DDelta = chipsLssa.out.map { chunk, meta, lssa ->
                def newMeta = [
                    ptype: '1D',
                    // pol: 'both',
                    plot_delta: true,
                    title: ''+"crosspower\\n${chunk}\\n${meta.name}",
                    plot_name: 'chips1d',
                    max_power: params.chips_max_power,
                    min_power: params.chips_min_power,
                    tags: [],
                ]
                [chunk, mapMerge(meta, newMeta), lssa]
            }

        singles2D = chipsLssa.out.map { chunk, meta, lssa ->
                def newMeta = [
                    ptype: '2D',
                    // pol: 'both',
                    title: ''+"crosspower\\n${chunk}\\n${meta.name}",
                    plot_name: 'chips2d',
                    max_power: params.chips_max_power,
                    min_power: params.chips_min_power,
                    tags: [],
                ]
                [chunk, mapMerge(meta, newMeta), lssa]
            }

        comps1D = chipsLssa.out
            // group lssas from vis name together
            .groupTuple(by: 0)
            // .view { chunk, metas, _ -> "\n -> comps1D metas: \n${metas.join('\n')}"}
            .filter { _chunk, metas, _lssas -> metas.size() >= 3 }
            .map { chunk, metas, lssas ->
                def nosubMeta = metas.find { it.sub == null} ?: [:]
                def subMeta = metas.find { it.sub == 'sub'} ?: [:]
                def ionosubMeta = metas.find { it.sub == 'ionosub'} ?: [:]
                def newMeta = [
                    ptype: '1D_comp',
                    // pol: 'both',
                    title: ''+"crosspower\\n${chunk}",
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
                    max_power: params.chips_max_power,
                    min_power: params.chips_min_power,
                ]
                [chunk, mapMerge(nosubMeta, newMeta), lssas.flatten()]
            }
        diffs2D_sub = chipsLssa.out
            // group lssas from vis name together
            .groupTuple(by: 0)
            .filter { _chunk, metas, _lssas -> metas.size() >= 2 }
            .map { chunk, metas, lssas ->
                def nosubMeta = metas.find { it.sub == null} ?: [:]
                def subMeta = metas.find { it.sub == 'sub'} ?: [:]
                def newMeta = [
                    ptype: '2D_diff',
                    // pol: 'both',
                    title: ''+"crosspower diff (nosub-sub)\\n${chunk}",
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
                [chunk, mapMerge(nosubMeta, newMeta), lssas.flatten()]
            }
        diffs2D_iono = chipsLssa.out
            // group lssas from vis name together
            .groupTuple(by: 0)
            .filter { _chunk, metas, _lssas ->
                def hasUnsubMeta = metas.find { it.sub == null} != null
                def hasIonoSubMeta = metas.find { it.sub == 'ionosub'} != null
                metas.size() >= 2 & hasUnsubMeta & hasIonoSubMeta
            }
            .map { chunk, metas, lssas ->
                def nosubMeta = metas.find { it.sub == null} ?: [:]
                def ionosubMeta = metas.find { it.sub == 'ionosub'} ?: [:]
                def newMeta = [
                    ptype: '2D_diff',
                    // pol: 'both',
                    title: ''+"crosspower diff (nosub-ionosub)\\n${chunk}",
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
                [chunk, mapMerge(nosubMeta, newMeta), lssas.flatten()]
            }

        diffs2D_sub_ionosub = chipsLssa.out
            // group lssas from vis name together
            .groupTuple(by: 0)
            .filter { _chunk, metas, _lssas ->
                def hasSubMeta = metas.find { it.sub == 'sub'} != null
                def hasIonoSubMeta = metas.find { it.sub == 'ionosub'} != null
                metas.size() >= 2 & hasSubMeta & hasIonoSubMeta
            }
            .map { chunk, metas, lssas ->
                def ionosubMeta = metas.find { it.sub == 'ionosub'} ?: [:]
                def subMeta = metas.find { it.sub == 'sub'} ?: [:]
                def newMeta = [
                    ptype: '2D_diff',
                    // pol: 'both',
                    title: ''+"crosspower diff (ionosub-sub)\\n${chunk}",
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
                [chunk, mapMerge(ionosubMeta, newMeta), lssas.flatten()]
            }

        ratios2D_sub = chipsLssa.out
            // group lssas from vis name together
            .groupTuple(by: 0)
            .filter { _chunk, metas, _lssas ->
                def hasSubMeta = metas.find { it.sub == 'sub'} != null
                def hasUnsubMeta = metas.find { it.sub == null} != null
                metas.size() >= 2 & hasSubMeta & hasUnsubMeta
            }
            .map { chunk, metas, lssas ->
                def nosubMeta = metas.find { it.sub == null} ?: [:]
                def subMeta = metas.find { it.sub == 'sub'} ?: [:]
                def newMeta = [
                    ptype: '2D_ratio',
                    // pol: 'both',
                    title: ''+"crosspower ratio (nosub:sub)\\n${chunk}",
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
                [chunk, mapMerge(nosubMeta, newMeta), lssas.flatten()]
            }

        ratios2D_ionosub = chipsLssa.out
            // group lssas from vis name together
            .groupTuple(by: 0)
            .filter { _chunk, metas, _lssas ->
                def hasUnsubMeta = metas.find { it.sub == null} != null
                def hasIonoSubMeta = metas.find { it.sub == 'ionosub'} != null
                metas.size() >= 2 & hasUnsubMeta & hasIonoSubMeta
            }
            .map { chunk, metas, lssas ->
                def nosubMeta = metas.find { it.sub == null} ?: [:]
                def ionosubMeta = metas.find { it.sub == 'ionosub'} ?: [:]
                def newMeta = [
                    ptype: '2D_ratio',
                    // pol: 'both',
                    title: ''+"crosspower ratio (nosub:ionosub)\\n${chunk}",
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
                [chunk, mapMerge(nosubMeta, newMeta), lssas.flatten()]
            }

        singles1D
            .mix(singles1DDelta)
            .mix(singles2D)
            .mix(comps1D)
            .mix(diffs2D_sub)
            .mix(diffs2D_iono)
            .mix(diffs2D_sub_ionosub)
            .mix(ratios2D_sub)
            .mix(ratios2D_ionosub)
            .filter { chunk, meta, __ ->
                if (meta.name == null) {
                    print("\n -> meta.name is null in ${chunk}") // ${meta}")
                    return false
                }
                if (meta.tags == null || meta.tags.contains(null)) {
                    print("\n -> meta.tags contains null in ${chunk}") //)" ${meta}")
                    return false
                }
                return true
            }
            | chipsPlot

    emit:
        frame = chipsPlot.out
            .flatMap { _group, meta, pngs ->
                def suffix = ""
                if (meta.nobs > 1) {
                    suffix = "_x" + String.format("%04d", meta.nobs)
                }
                coerceList(pngs).collect { png ->
                    [''+"${meta.plot_name}_${meta.name}${suffix}", png]
                }
            }.groupTuple()
}

// entrypoint: move unfiltered preprocessed uvfits from asvo accacia to mwaeor accacia
workflow archivePrep {
    obsids = channel.of(obsids_file())
        .splitCsv()
        .flatten()
        .filter { line -> !line.startsWith('#') }
        .unique()

    obsids | asvoPrep

    if (params.archive) {
        asvoPrep.out.map { _o, _m, vis -> ["prep", vis] }
            | archive
    }
}

// given a channel of tuple (obs, metafits, vis), calibrate and analyse
workflow qaPrep {
    take:
        subobsMetaVis
        obsMetafits

    main:
        obsMetafits.cross(subobsMetaVis).map{ obsMetafits_, subobsMetaVis_ ->
                def (obsid, metafits) = obsMetafits_
                def (_o, meta, vis) = subobsMetaVis_
                [obsid, meta, metafits, vis]
            }
            .filter { _o, _m, metafits, vis -> metafits != null && vis != null }
            // .view { it -> "subobsMetaMetafitsVis ${it[0]} ${it[1].subobs?:''}" }
            .tap { subobsMetaMetafitsVis }
            .map { obsid, meta, _metafits, _vis -> [obsid, meta] }
            .tap { obsMeta }

        // subobsMetaMetafitsVis | hypPrepVisConvert

        // get sourcelists for each obs (only currently used in subtraction, not calibration)
        subobsMetaMetafitsVis.map { obsid, _m, metafits, _vis -> [obsid, metafits ] }
            .unique()
            | (hypSrclistAO & hypSrclistYaml)

        obsMetafits.cross(hypSrclistAO.out)
            .map { obsMetafits_, hypSrclistAO_ ->
                def (obsid, metafits) = obsMetafits_
                def (_o, srclist) = hypSrclistAO_
                [obsid, metafits, srclist]
            }

        // cluster unless --nocluster
        // if (params.nocluster) {
        //     channel.empty() | rexCluster
        //     obsMetafitsSrclist = obsMetafits.cross(hypSrclistAO.out)
        //         .map { obsMetafits_, hypSrclistAO_ ->
        //             def (obsid, metafits) = obsMetafits_
        //             def (_, srclist) = hypSrclistAO_
        //             [obsid, metafits, srclist] }
        // } else {
        //     hypSrclistAO.out | rexCluster
        //     obsMetafitsSrclist = obsMetafits.cross(rexCluster.out)
        //         .map { obsMetafits_, rexCluster_ ->
        //             def (obsid, metafits) = obsMetafits_
        //             def (_, srclist) = rexCluster_
        //             [obsid, metafits, srclist] }
        // }

        // calibrate each obs that passes flag gate unless --nocal:
        if (params.nocal) {
            // empty channel disables a process
            channel.empty() | cal
        } else {
            subobsMetaMetafitsVis | cal
        }

        cal.out.obsMetaCalPass
            .map { def (obsid, meta, _soln) = it;
                def calFlags = meta.calFlags?:[]
                ([obsid, meta.name, meta.subobs?:'', displayInts(meta.prepFlags?:[]), displayInts(calFlags)]).join("\t")
            }
            .collectFile(
                name: "pass_cal.tsv", newLine: true, sort: true,
                seed: ([ "OBS", "NAME", "SUBOBS", "PREPFLAGS", "CALFLAGS" ]).join("\t"),
                storeDir: "${results_dir()}${params.cal_suffix}"
            )
            | view { [it, it.readLines().size()] }

        // channel of arguments for hypApply{UV,MS}
        // - take tuple(obsid, meta, soln) from cal.out.obsMetaCalPass
        // - match with tuple(obsid, metafits, prepUVFits) by obsid

        // apply calibration solutions to uvfits and ms unless --nouv or --noms
        subobsMetaMetafitsVis.cross(cal.out.obsMetaCalPass) { [it[0], it[1].subobs?:''] }
            .map { subobsMetaMetafitsVis_, obsMetaCalPass_ ->
                def (obsid, _m, metafits, prepUVFits) = subobsMetaMetafitsVis_
                def (_o, meta, soln) = obsMetaCalPass_
                def newMeta = [
                    // warning: freq_res is overloaded, could be in hz or kHz
                    time_res: params.apply_time_res,
                    freq_res_khz: params.apply_freq_res,
                    nodut1: params.nodut1,
                    apply_args: params.apply_args?:'',
                ]
                def calFlags = meta.calFlags?:[]
                if (calFlags) {
                    newMeta.apply_args = ''+"${params.apply_args?:''} --tile-flags ${calFlags.join(' ')}"
                }
                // calFlags = [8,9,10,11,12,13,15,27,30,39,40,42,44,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,79,88,89,90,91,92,93,94,95,117]
                // newMeta.apply_args = "${newMeta.apply_args} --tile-flags ${calFlags.join(' ')}"
                newMeta.apply_name = hyp_apply_name(newMeta.time_res, newMeta.freq_res_khz)
                [obsid, mapMerge(meta, newMeta), metafits, prepUVFits, soln]
            }
            | (hypApplyUV & hypApplyMS)
        // TODO: make hyperdrive not crash if there are not num srcs
        subArgsUV = obsMetafits.cross(hypApplyUV.out)
            .map { obsMetafits_, hypApplyUV_ ->
                def (obsid, metafits) = obsMetafits_
                def (__, meta, vis) = hypApplyUV_;
                def srclist = file(params.sourcelist)
                def newMeta = [sub_nsrcs: params.sub_nsrcs]
                [obsid, mapMerge(meta, newMeta), metafits, vis, srclist]
            }
        subArgsMS = obsMetafits.cross(hypApplyMS.out)
            .map { obsMetafits_, hypApplyUV_ ->
                def (obsid, metafits) = obsMetafits_
                def (__, meta, vis) = hypApplyUV_;
                def srclist = file(params.sourcelist)
                def newMeta = [sub_nsrcs: params.sub_nsrcs]
                [obsid, mapMerge(meta, newMeta), metafits, vis, srclist]
            }
        subArgsUV | hypSubUV
        subArgsMS | hypSubMS
        subArgsUV.map {obsid, meta, metafits, vis, srclist ->
            def newMeta = [ionosub_nsrcs: params.ionosub_nsrcs]
            [obsid, mapMerge(meta, newMeta), metafits, vis, srclist]}
            | hypIonoSubUV
        subArgsMS.map {obsid, meta, metafits, vis, srclist ->
            def newMeta = [ionosub_nsrcs: params.ionosub_nsrcs]
            [obsid, mapMerge(meta, newMeta), metafits, vis, srclist]}
            | hypIonoSubMS

        // make cthulhu plots
        hypIonoSubUV.out.cross(hypSrclistYaml.out) {it[0]}
            .map { hypIonoSubUV_, hypSrclistYaml_ ->
                def (obsid, meta, _vis, offsets, __) = hypIonoSubUV_;
                def (_o, srclist) = hypSrclistYaml_;
                [obsid, meta, srclist, offsets]
            }
            | cthulhuPlot

        // collect ionosub constants as tsvs
        // ["alphas", "betas", /**"gains"**/].collect {
        //     hypIonoSubUV.out.flatMap { obsid, meta, _, offsetsJson, __ ->
        //             offsets = parseJson(offsetsJson)
        //             offsets.collect { src, data ->
        //                 pos = data['weighted_catalogue_pos_j2000']
        //                 ([src, pos.ra?:'', pos.dec?:''] + data[it]).join("\t")
        //             }
        //         }
        //         .collectFile(
        //             name: "iono_${it}.tsv", newLine: true, sort: false,
        //             seed: ["SRC", "RA", "DEC", it.toUpperCase()].join("\t"),
        //             storeDir: "${results_dir()}${params.cal_suffix}"
        //         )
        //         | view { [it, it.readLines().size()] }
        // }

        // collect ionoqa results as .tsv
        cthulhuPlot.out
            // form row of tsv from json fields we care about
            .map { obsid, meta, __, ___, json, ____ ->
                def stats = parseJson(json);
                [
                    obsid,
                    isNaN(meta.lst)?'':meta.lst,
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
                    "LST",
                    "HYP_IONO_PCA",
                    "HYP_IONO_MAG",
                    "HYP_IONO_QA",
                    "MAX_SHIFT_SRC",
                    "MAX_SHIFT_RA",
                    "MAX_SHIFT_DEC",
                    "N_TIMES",
                    "N_SRCS",
                ].join("\t"),
                storeDir: "${results_dir()}${params.cal_suffix}"
            )
            // display output path and number of lines
            | view { [it, it.readLines().size()] }

        // TODO: use iono pass filter
        // def obsIonoPass = cthulhuPlot.out.map { obsid, _, __, ___, json, ____ ->
        //         stats = parseJson(json);
        //         [obsid, stats.QA?:null]
        //     }
        //     .filter { _, qa ->
        //         qa != null && (params.filter_max_hyp_ionoqa == null || qa < params.filter_max_hyp_ionoqa)
        //     }
        //     .map { obsid, qa -> [obsid] }

        // channel of calibrated, subtracted and ionosubtracted uvfits: tuple(obsid, name, uvfits)
        hypApplyUV.out.map { obsid, meta, vis, _logs -> [obsid, meta, vis] }
            .mix(hypSubUV.out.map { obsid, meta, vis, _logs -> [obsid, meta, vis] })
            .mix(hypIonoSubUV.out.map { obsid, meta, vis, _json, _logs -> [obsid, meta, vis] })
            .tap{ obsMetaUV }
            | uvMeta

        // improve uvfits meta
        obsMetaUVSmart = obsMetaUV.cross(uvMeta.out) {[it[0], it[1].name, it[1].subobs?:'']}
            .map { obsMetaUV_, uvMeta_ ->
                def (obsid, meta, vis) = obsMetaUV_;
                def (_o, _m, json) = uvMeta_;
                def jsonMeta = parseJson(json)
                def newMeta = [
                    lowfreq:jsonMeta.freqs[0],
                    first_lst: jsonMeta.times[0]['lst_rad'],
                    first_jd: jsonMeta.times[0]['jd1'],
                    nchans:jsonMeta.freqs.size(),
                    ntimes:jsonMeta.times.size(),
                ];
                ['eorband', 'eorfield', 'config', 'total_weight'].each { key ->
                    if (jsonMeta[key] != null) {
                        newMeta[key] = jsonMeta[key]
                    }
                }
                [obsid, mapMerge(meta, newMeta), vis] }

        // QA uvfits visibilities
        obsMetaUVSmart | uvfits

        // channel of calibrated, subtracted and ionosubtracted ms: tuple(obsid, name, vis)
        obsMetaMS = hypApplyMS.out.map { obsid, meta, vis, _log -> [obsid, meta, vis] }
            .mix(hypIonoSubMS.out.map { obsid, meta, vis, _json, _log -> [obsid, meta, vis] })
            .mix(hypSubMS.out.map { obsid, meta, vis, _log -> [obsid, meta, vis] })

        // image and qa measurementsets or uvfits unless --noimage
        if (params.noimg) {
            channel.empty() | img
            obsMetaImgPass = uvfits.out.obsMetaUVPass.map { obsid, meta, _vis -> [obsid, meta, [], []]}
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
                .filter { _o, meta, _vis ->
                    if (params.noimgunsub) {
                        return (meta.sub?:'') != ''
                    }
                    return true
                }
                | img

            obsMetaImgPass = img.out.obsMetaImgPass
        }

        // tuple of (obsid, meta) for all visibilities that pass imgQA
        obsMetaPass = obsMetaImgPass.map { def (obsid, meta) = it; [obsid, meta] }

        obsMetaPass
            .map { def (obsid, meta) = it;
                ([obsid, meta.name]).join("\t")
            }
            .collectFile(
                name: "img_names_pass.tsv", newLine: true, sort: true,
                seed: ([ "OBS", "NAME"]).join("\t"),
                storeDir: "${results_dir()}${params.img_suffix}${params.cal_suffix}"
            )
            | view { [it, it.readLines().size()] }

        // Apologies in advance for anyone who has to debug this.
        //
        // We want a tuple of (chunk, chunkMeta, obsids[G]) for imaging and power spectra.
        //
        // The pipeline can produce multiple visibilities for a given obsid. `meta` keys can be used to differentiate these. e.g.
        // - `sub` - type of subtraction (`null` means no subtraction)
        // - `cal_prog` - the program used to calibrate the visibility
        // - `poly` - whether a polyfit was applied
        //
        // We want to compare results between visibilities of different types within a single obsid,
        // and within compatible groups of obsids.
        // Within the set of visibilities for each obsid, one is the "primary" visibility that others
        // are compared to, e.g. the one with no subtractions.
        // A group of obsids is compatible if they have the same field, band and telescope configuration.
        // The metadata for each group of obsids comes from primary meta
        // Within each group of obsids, we sort by some metric in the primary meta  (e.g. window : wedge power ratio)
        // Then we split that group into chunks of size G.
        obsMetaPass
            // group metas by obsid
            .groupTuple(by: 0)
            // ensure metas are a list, even when there is only one meta
            .map { obsid, metas -> [ obsid, coerceList(metas) ]}
            // determine how to group each obsid, and how to sort that obsid within each group.
            .map { obsid, metas ->
                def gmetas = metas.collect { meta -> mapMerge(meta, groupMeta(meta)) }
                // find the meta of the primary subobservation (in this case, the unsubtracted visibilities)
                def gmetaPrime = gmetas.find { meta -> meta.sub == null }
                if (gmetaPrime == null) {
                    gmetaPrime = gmetas[0]
                }
                ["${gmetaPrime.group}", gmetaPrime.sort, obsid, gmetas]
            }
            .groupTuple(by: 0)
            // only keep groups with more than one obsid
            .filter { it -> it[1] instanceof List }
            // .view { it -> "\n -> grouped by obsid: ${it}"}
            // chunk obsid groups
            .flatMap { group, all_sorts, all_obsids, all_metas ->
                groupMeta = [group: group]

                [all_sorts, all_obsids, all_metas].transpose()
                    .sort { it -> it[0] }
                    .collate(params.chunkSize, params.chunkRemainder)
                    .take(params.chunkCount)
                    .collect { chunk ->
                        def (sorts, obsids, metas) = chunk.transpose()
                        def obs_list = obsids.sort(false)
                        def hash = obs_list.join(' ').md5()[0..7]
                        def chunkMeta = [
                            // ntimes: metas.collect { meta -> meta[0].ntimes?:0 }
                            sort_bounds: [sorts[0], sorts[-1]],
                            hash: hash,
                            nobs: obs_list.size(),
                        ]
                        // if $group already ends with $hash then don't add it to the group
                        def group_hash = group
                        if (!group_hash.endsWith(hash)) {
                            group_hash = "${group}_${hash}"
                        }
                        [group_hash, mapMerge(groupMeta, chunkMeta), obs_list, metas]
                    }
            }
            // .view { it -> "\n -> chunked by obsid: ${it}"}
            .tap { groupChunkMetaPass }
            // flatten obsids out of each chunk
            .flatMap { chunk, chunkMeta, obsids, all_metas ->
                [obsids, all_metas].transpose().collect { obsid, metas ->
                    [chunk, obsid, metas.collect { meta -> mapMerge(meta, chunkMeta)}]
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
            .map { chunk, _name, obsids, metas ->
                [chunk, metas[0], obsids]
            }
            // finally, we have a channel of visibilities from the same obs group and vis type
            // .view { it -> "\n -> chunkMetaPass: ${it}"}
            .tap { chunkMetaPass }

        groupChunkMetaPass
            .map { chunk, chunkMeta, obsids, _all_metas ->
                def sort_bounds = chunkMeta.sort_bounds?:[Float.NaN, Float.NaN]
                (
                    sort_bounds.collect { String.format("%9.9f", it) } \
                    + [chunk, obsids.size(), obsids.join(' ')]
                ).join("\t")
            }
            .collectFile(
                name: "obs_group_chunks.tsv", newLine: true, sort: true,
                seed: ([ "SORT LEFT", "SORT RIGHT", "GROUP CHUNK", "NOBS", "OBSIDS" ]).join("\t"),
                storeDir: "${results_dir()}${params.img_suffix}${params.cal_suffix}"
            )
            | view { [it, it.readLines().size()] }

        chunkMetaPass.map { group, _m, obsids -> [group, obsids] }
            .unique()
            | storeManifest

        chunkMetaPass
            .map { group, meta, obsids ->
                [group, meta.name, obsids.size(), obsids.join(' ')].join("\t")
            }
            .collectFile(
                name: "obs_groups.tsv", newLine: true, sort: true,
                seed: ([ "GROUP", "NAME", "NVIS", "VISS" ]).join("\t"),
                storeDir: "${results_dir()}${params.img_suffix}${params.cal_suffix}"
            )
            | view { [it, it.readLines().size()] }

        uvfits.out.obsMetaUVPass.cross(obsMetaPass) { def (obsid, meta) = it; [obsid, meta.name, meta.subobs?:''] }
            .map { subobsMetaUVPass_, obsMetaPass_ ->
                def (obsid, _m, uvfits) = subobsMetaUVPass_;
                def (_o, meta) = obsMetaPass_;
                if (meta.obsid != null && meta.obsid != obsid) {
                    throw new Exception("obsid mismatch: meta.obsid ${meta.obsid} != obsid ${obsid}")
                }
                [obsid, meta, uvfits]
            }
            .tap { obsMetaUVPass }

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
                    def (_o, _m, imgMetas, imgs) = obsMetaImgPass_;
                    def (_o_, groupMeta, group) = obsMetaPass_;
                    [imgMetas, imgs].transpose().collect { imgMeta, img ->
                        [group, groupMeta.name, imgMeta.suffix, imgMeta, img]
                    }
                }
                .groupTuple(by: 0..2)
                .map { group, _name, _suffix, imgMetas, imgs ->
                    [group, imgMetas[0], imgs]
                }
                | stackImgs
                | stackThumbnail
        }

        // archive data to object store
        if (params.archive) {
            if (params.archive_prep) {
                prep_archive = subobsMetaMetafitsVis.map { _o, _m, _metafits, vis -> ["prep", vis] }
            } else {
                prep_archive = channel.empty()
            }
            if (params.archive_uvfits) {
                vis_archive = obsMetaUV.map { _o, _m, vis -> ["uvfits", vis]}
            } else {
                vis_archive = channel.empty()
            }
            prep_archive.mix(vis_archive)
                .mix(cal.out.obsMetaCalPass.map { _o, _m, soln -> ["soln", soln] })
                .mix(cal.out.archive)
                .mix(uvfits.out.archive)
                .mix(img.out.archive)
                | archive
        }

    emit:
        obsMetaUVPass
        frame = cal.out.frame
            .mix(uvfits.out.frame)
            .mix(img.out.frame)
            .mix(chips.out.frame)
            .mix(cthulhuPlot.out.flatMap {
                def (_o, meta, pngs, _csv, _json, tecs) = it;
                pngs.collect { png ->
                    [''+"cthulhuplot_${meta.name}", png]
                } + tecs.collect { tec ->
                    [''+"tec_${meta.name}", tec]
                }
            }.groupTuple())
            .mix(imgCombine.out.frame)
        zip = cal.out.zip
            .mix(uvfits.out.zip)
            .mix(img.out.zip)
            .mix(
                hypIonoSubUV.out.map { _o, meta, _vis, offsets, _log ->
                    [''+"offsets_${meta.name}", offsets]
                }
                .groupTuple()
            )

        fail_codes = uvfits.out.fail_codes.join(img.out.fail_codes, remainder: true)
            .map { obsid, uv_code, img_code ->
                def codes = [uv_code, img_code].findAll { it != null }
                def aFailCode = coerceList(firstFail([codes]))[0]
                [obsid, aFailCode]
            }
}

workflow makeVideos {
    take:
        frame
    emit:
        videos = frame.map { name, frames_ ->
                def frames = coerceList(frames_).flatten()
                def latest = frames.collect { it.lastModified() }.max()
                def cachebust = "${latest}_x" + String.format("%04d", frames.size())
                def sorted = frames.collect { path -> file(path.toString()) }.sort(false)
                [name, sorted, cachebust]
            }
            | ffmpeg
            | view { video, _cachebust -> [video, video.size()] }
}

workflow makeTarchives {
    take:
        zip
    emit:
        zips = zip.map { name, files ->
                def latest = files.collect { file -> file.lastModified() }.max()
                def cachebust = ''+"${latest}_x" + String.format("%04d", files.size())
                // sorted = files.collect { path -> file(path.toString()) }.sort(false)
                [name, files, cachebust]
            }
            | tarchive
            | view { zip_, _cachebust -> [zip_, zip_.size()] }
}

// entrypoint: get externally preprocessed uvfits files from csv file and run qa
workflow extPrep {

    def name = params.visName ?: "ssins"
    def cal_prog = params.cal_prog ?: "hyp"

    channel.of(obsids_file())
        .splitCsv()
        .filter { line -> !line[0].startsWith('#') }
        .map { line ->
            def (obsid, _comment) = line
            obsid
        }
        .unique()
        .map { obsid_ ->
            def obsid = coerceList(obsid_)[0]
            def (_argstr_asvo, _argstr_cli, suffix) = birli_argstr_suffix()
            def meta = [name:name, cal_prog:cal_prog, obsid: obsid, birli_suffix: suffix]
            def vis = file("${params.outdir}/${obsid}/prep/birli_${obsid}${suffix}.${name}.uvfits")
            // def vis = file("${params.outdir}/${obsid}.uvfits")
            [ obsid, deepcopy(meta), vis ]
        }
        .filter { _o, _m, vis -> vis.exists() }
        .tap { obsMetaVis }
        .map { obsid, _m, _vis -> obsid }
        .tap { obsids }

    obsids | ws

    obsWsmetaVis = obsMetaVis.join(ws.out.obsMeta)
        .map { obsid, meta, vis, wsMeta ->
            [ obsid, mapMerge(meta, wsMeta), vis ]
        }

    flag(obsWsmetaVis, ws.out.obsMetafits)

    obsFlagmetaVis = flag.out.subobsMetaPass.map { obsid, meta -> [[obsid, meta.subobs?:''], meta] }
            .join(obsWsmetaVis.map { obsid, meta, uvfits -> [[obsid, meta.subobs?:''], uvfits] })
            // .join(flag.out.subobsMetaRxAnts.map { obsid, meta, rxAnts -> [[obsid, meta.subobs?:''], rxAnts] })
            .map { obsSubobs, meta, uvfits ->
                def (obsid, _subobs) = obsSubobs
                [obsid, meta, uvfits]
            }

    qaPrep(obsFlagmetaVis, ws.out.obsMetafits)

    ws.out.fail_codes
        .join(flag.out.fail_codes, remainder: true)
        .join(qaPrep.out.fail_codes, remainder: true)
        .map { obsid, ws_code, flag_code, qaprep_code ->
            def codes = [ws_code, flag_code, qaprep_code].findAll { it != null }
            [obsid, coerceList(firstFail([codes]))[0]]
        }
        .tap { all_fail_codes }
        .groupTuple(by: 1)
        .map { obsids_, fail_code ->
            [ fail_code, obsids_.size() ].join("\t")
        }
        .collectFile(
            name: "fail_counts_all.tsv", newLine: true, sort: true,
            seed: ([ "FAIL CODE", "COUNT" ]).join("\t"),
            storeDir: "${results_dir()}${params.img_suffix}${params.cal_suffix}"
        )

    ws.out.frame
        .mix(flag.out.frame)
        .mix(qaPrep.out.frame)
        | makeVideos
    // make zips
    if (params.tarchive) {
        qaPrep.out.zip | makeTarchives
    }
}

workflow extPSMetrics {
    obsVis = channel.of(file("external.csv"))
        .splitCsv(header: ["obsid", "vis"])
        .filter { !it.obsid.startsWith('#') }
        .map { [it.obsid, file(it.vis)] }

    obsVis.map { obs, vis -> [obs, [:], vis ]} | psMetrics
}

workflow bootstrap {
    def name = params.visName ?: "ssins"
    def cal_prog = params.cal_prog ?: "hyp"

    channel.of(obsids_file())
        .splitCsv()
        .filter { line -> !line[0].startsWith('#') }
        .map { line ->
            def (obsid, _comment) = line
            obsid
        }
        .unique()
        .flatMap { obsid_ ->
            def obsid = coerceList(obsid_)[0]
            def (_argstr_asvo, _argstr_cli, prep_suffix) = birli_argstr_suffix()
            def apply_name = hyp_apply_name(params.apply_time_res, params.apply_freq_res)
            ["${params.dical_name}": params.dical_args].keySet().collect { dical_name ->
                def meta = [
                    obsid: obsid, cal_prog:cal_prog,
                    prep_suffix:prep_suffix,
                    apply_suffix: (apply_name ? "_${apply_name}" : ""),
                    dical_suffix:"_${dical_name}",
                    peel_suffix: "${params.peel_suffix}" + (params.ionosub_nsrcs != null ? "_i${params.ionosub_nsrcs}" : ""),
                ]
                [ obsid, deepcopy(meta) ]
            }
        }
        .flatMap { obsid, meta ->
            def prep_suffix = meta.prep_suffix ?: ""
            def dical_suffix = meta.dical_suffix ?: ""
            def peel_suffix = meta.peel_suffix ?: ""
            def apply_suffix = meta.apply_suffix ?: ""
            if (params.ssins_apply) {
                prep_suffix = "${prep_suffix}.ssins"
                dical_suffix = "_ssins${dical_suffix}"
            }
            def result = [
                [obsid, "prep", file("${params.outdir}/${obsid}/prep/birli_${obsid}${prep_suffix}.uvfits")],
            ]
            if (params.nocal) { return result }
            result.add([obsid, "soln", file("${params.outdir}/${obsid}/cal/${meta.cal_prog}_soln_${obsid}${dical_suffix}.fits")])
            result.add([obsid, null, file("${params.outdir}/${obsid}/cal/${meta.cal_prog}_di-cal_${obsid}${dical_suffix}.log")])
            if (params.noapply) { return result }
            result.add([obsid, "uvfits", file("${params.outdir}/${obsid}/cal/${meta.cal_prog}_${obsid}${dical_suffix}${apply_suffix}.uvfits")])
            result.add([obsid, null, file("${params.outdir}/${obsid}/cal/${meta.cal_prog}_apply${dical_suffix}${apply_suffix}.log")])
            if (params.noionosub) { return result }
            result.add([obsid, "uvfits", file("${params.outdir}/${obsid}/cal/${meta.cal_prog}_${obsid}_ionosub${dical_suffix}${apply_suffix}${peel_suffix}.uvfits")])
            result.add([obsid, "peel", file("${params.outdir}/${obsid}/cal/${meta.cal_prog}_peel_${obsid}_ionosub${dical_suffix}${apply_suffix}${peel_suffix}_uv.json")])
            result.add([obsid, null, file("${params.outdir}/${obsid}/cal/${meta.cal_prog}_vis-ionosub${dical_suffix}${apply_suffix}${peel_suffix}_uv.log")])
            return result
        }
        .branch { _, _bucket, f ->
            exists: f.exists()
            not_exists: true
        }
        .set { files }

    // def now = new Date();
    // files.exists.map { obsid, bucket, f ->
    //         def modDate = new Date(f.lastModified())
    //         def age = (now.getTime() - modDate.getTime()) / 1000 / 60 / 60 / 24
    //         print("age: ${age} days: ${f}")
    //     }

    files.not_exists.map { obsid, bucket, f ->
            def store = "${obsid}/${f.parent.name}"
            [store, bucket, f.name]
        }
        | unArchive
}
