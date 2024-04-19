#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// def REMEMBER TO ALWAYS DEF https://stackoverflow.com/questions/184002/groovy-whats-the-purpose-of-def-in-def-x-0

// TODO: evaluate(new File('./lib/lib.groovy'))

def obsids_file = file(params.obsids_path)
if (params.obsids_suffix) {
    obsids_file = file("${obsids_file.parent}/${obsids_file.baseName}${params.obsids_suffix}.${obsids_file.extension}")
}
def results_dir = '' + "${params.resultsdir}/results${params.obsids_suffix}${params.result_suffix}"

// whether imaging is configured for multiple channels
img_channels_out = (params.img_channels_out instanceof String ? \
    params.img_channels_out.split(' ')[0] :\
    params.img_channels_out)
multichannel = (img_channels_out as int > 1)
// whether imaging is configured for multiple intervals
multiinterval = (params.img_intervals_out as int > 1) || params.img_split_intervals

def coerceList(x) {
    if (x instanceof nextflow.util.ArrayBag) {
        x.toList()
    } else if (x instanceof List) {
        x
    } else {
        [x]
    }
}

import java.text.SimpleDateFormat
import java.util.LinkedHashMap
import groovy.transform.Synchronized
import org.codehaus.groovy.runtime.StackTraceUtils

def deepcopy(orig) {
    def bos = new ByteArrayOutputStream()
    def oos = new ObjectOutputStream(bos)
    try {
        oos.writeObject(orig); oos.flush()
    } catch (Exception e) {
        println("error deepcopying ${orig} ${e}")
        StackTraceUtils.sanitize(e).printStackTrace()
        throw e
    }
    def bin = new ByteArrayInputStream(bos.toByteArray())
    def ois = new ObjectInputStream(bin)
    return ois.readObject()
}

def mapMerge(a, b) {
    // > a = [a:1]
    // > b = a + [a:2]
    // > println a
    // [a:1]
    return deepcopy( a + b )
}

// sources:
// - man 7 signal
// - https://www.gnu.org/software/wget/manual/html_node/Exit-Status.html
exitCodes = [
    1: "general or permission",
    2: "incorrect usage",
    3: "File I/O error",
    4: "Network failure",
    // 5: "I/O error or hash mismatch",
    5: "Hash mismatch or SSL verification",
    6: "Authentication failure",
    8: "Server issued an error response",
    11: "Resource temporarily unavailable",
    28: "No space left on device",
    75: "Temporary failure, try again",
    101: "Panic! at the kernel",
    127: "File or directory not found",
    135: "??? qaPrep:img:thumbnail",
    137: "SIGKILL - killed (OOM)",
    139: "SIGSEGV - segmentation fault",
    140: "SIGUSR2 - out of time",
]

// download observation metadata from webservices in json format
process wsMeta {
    input:
    val(obsid)
    output:
    tuple val(obsid), path(wsmeta), path(wsfiles)

    // persist results in outdir, process will be skipped if files already present.
    storeDir "${params.outdir}/meta"
    // tag to identify job in squeue and nf logs
    tag "${obsid}"

    time {5.minute * task.attempt}

    // allow multiple retries
    maxRetries 2
    errorStrategy {
        return (task.exitStatus == 8 ? 'retry' : 'ignore')
    }

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
    input:
        tuple val(obsid), path(metafits)
    output:
        tuple val(obsid), path(json), path(tsv)

    storeDir "${params.outdir}/${obsid}/meta"
    tag "${obsid}"

    errorStrategy 'ignore'

    label "python"

    script:
    // metrics = "${obsid}_occupancy.json"
    json = "${obsid}_meta.json"
    tsv = "${obsid}_inputs.tsv"
    txt = "${obsid}_inputs.txt"
    template "metajson.py"
}

// temporary ASVO workaround, raw files from mwacache
process cacheBoxRaw {
    input:
    tuple val(obsid), val(meta)
    output:
    tuple val(obsid), val(meta), path(metafits), path(raw)

    storeDir "${params.outdir}/${obsid}/raw"

    tag "${obsid}"

    // errorStrategy "terminate"

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
    input:
    tuple val(obsid), val(meta_), path(metafits), path(raw)
    output:
    tuple val(obsid), val(meta), path(uvfits)
        // , path("${obsid}${spw}*.mwaf"), path("birli_prep.log")

    storeDir "${params.outdir}/${obsid}/prep"

    tag "${obsid}"

    script:
    meta = deepcopy(meta_)
    prefix = "birli_"
    suffix = ''
    args = [:]
    if (params.prep_time_res_s != null) {
        args['avg-time-res'] = params.prep_time_res_s
        suffix += "_${params.prep_time_res_s}s"
    }
    if (params.prep_freq_res_khz != null) {
        args['avg-freq-res'] = params.prep_freq_res_khz
        suffix += "_${params.prep_freq_res_khz}kHz"
    }
    if (params.prep_rfi != null && !params.prep_rfi) {
        args['no-rfi'] = null
        suffix += '_norfi'
    }
    if (params.prep_pc != null) {
        // --phase-centre <RA> <DEC>
        args['phase-centre'] = params.prep_pc
    }
    if (params.prep_edge_width != null) {
        args['flag-edge-width'] = params.prep_edge_width
        suffix += "_edg${params.prep_edge_width}"
    }
    meta['birli_args'] = args
    meta['birli_suffix'] = suffix
    argstr = args.collect { k, v ->
            if (v == null) {
                ['--'+k]
            } else {
                ['--'+k] + v
            }
        }
        .flatten()
        .join(' ');
    uvfits = ''+"${prefix}${obsid}${suffix}*.uvfits"
    """
    set -eux
    ${params.birli} \
        ${argstr} \
        -u "${prefix}${obsid}${suffix}.uvfits" \
        -m "${metafits}" \
        ${raw}
    """
}

workflow extRaw {

    channel.of(obsids_file)
        .splitCsv()
        .filter { line -> !line[0].startsWith('#') }
        .map { line ->
            def (obsid, _) = line
            obsid
        }
        .unique()
        .map { obsid_ ->
            def obsid = coerceList(obsid_)[0]
            def meta = [obsid: obsid]
            // def raw = file("${params.outdir}/../raw/${obsid}_2?????????????_ch???_???.fits")
            def raw = file("${params.outdir}/../raw/${obsid}_2?????????????_ch???_???.fits")
            [ obsid, deepcopy(meta), raw ]
        }
        .filter { _, __, raw_ ->
            def raw = coerceList(raw_)
            raw.size() && raw.every{
                if (!it.exists()) { print("raw does not exist: ${it}") }
                it.exists()
            }
        }
        .tap { obsMetaRaw }
        .map { obsid, meta, raw -> obsid }
        .tap { obsids }

    obsids | ws

    obsWsmetaVis = obsMetaRaw.join(ws.out.obsMeta).join(ws.out.obsMetafits)
        .map { obsid, meta, vis, wsMeta, metafits ->
            [ obsid,  mapMerge(meta, wsMeta), metafits, vis ]
        }
        | birliPrepUV

    birliPrepUV.out.flatMap { obsid, meta, uvfits ->
            coerceList(uvfits).collect { f ->
                [obsid, mapMerge(meta, [subobs: f.baseName.split('_')[-1]]), f]
            }
        }
        .tap { subobsVis }
        | uvMeta

    ws.out.obsMeta.cross(uvMeta.out) { it[0] }
        .map { obsMeta_, uvMeta_ ->
            def (obsid, wsMeta) = obsMeta_
            def (_, meta_, uvJson) = uvMeta_
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
            def (_, __, uvfits) = subobsVis_
            [obsid, meta, uvfits]
        }
        .tap { subobsMetaVis }

    subobsMetaVis | ssins

    flag(obsWsmetaVis, ws.out.obsMetafits)

    flag.out.subobsMetaPass.map { obsid, meta -> [[obsid, meta.subobs?:''], meta] }
        .join(obsWsmetaVis.map { obsid, meta, uvfits -> [[obsid, meta.subobs?:''], uvfits] })
        // .join(flag.out.subobsMetaRxAnts.map { obsid, meta, rxAnts -> [[obsid, meta.subobs?:''], rxAnts] })
        .map { obsSubobs, meta, uvfits ->
            def (obsid, _) = obsSubobs
            [obsid, meta, uvfits]
        }
        .tap { subobsFlagmetaVis }

    qaPrep( subobsFlagmetaVis, ws.out.obsMetafits )
    // subobsMetaVis.mix(qaPrep.out.obsMetaUVPass) | ssins

    ws.out.frame.mix(flag.out.frame).mix(qaPrep.out.frame)
        .mix(
            ssins.out.flatMap { _, __, ___, imgs, ____ ->
                imgs.collect { img ->
                    def tokens = coerceList(img.baseName.split('_') as ArrayList)
                    def suffix = tokens[-1]
                    if (suffix == "SSINS") {
                        try {
                            prefix = tokens[-2]
                        } catch (Exception e) {
                            prefix = ''
                            println(tokens)
                            throw e
                        }
                        ["ssins_${prefix}", img]
                    } else {
                        try {
                            prefix = tokens[-3]
                        } catch (Exception e) {
                            prefix = ''
                            println(tokens)
                            throw e
                        }
                        ["ssins_${prefix}_${suffix}", img]
                    }
                }
            }.groupTuple()
        )
        .map { n, l -> [n, l as ArrayList] } | makeVideos
}


workflow extCache {
    def name = params.visName ?: "ssins"

    channel.of(obsids_file)
        .splitCsv()
        .filter { line -> !line[0].startsWith('#') }
        .map { line ->
            def (obsid, _) = line
            obsid
        }
        | ws

    ws.out.obsMetafits.map { obsid, _ ->
            def meta = [:]
            [ obsid, deepcopy(meta) ]
        }
        | cacheBoxRaw
        | birliPrepUV

    birliPrepUV.out.flatMap { obsid, meta, uvfits ->
            uvfits.collect { f ->
                [obsid, mapMerge(meta, [subobs: f.baseName.split('_')[-1]]), f]
            }
        }
        .tap { subobsVis }
        | uvMeta

    subobsMetaVis = ws.out.obsMeta.cross(uvMeta.out) { it[0] }
        .map { obsMeta_, uvMeta_ ->
            def (obsid, wsMeta) = obsMeta_
            def (_, meta_, uvJson) = uvMeta_
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
            def (_, __, uvfits) = subobsVis_
            [obsid, meta, uvfits]
        }

    qaPrep( subobsMetaVis, ws.out.obsMetafits )

    qaPrep.out.frame | makeVideos

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
    input:
    tuple val(obsid), val(meta)
    output:
    tuple val(obsid), val(meta), path(uvfits)

    storeDir "${params.outdir}/${obsid}/prep"
    // tag to identify job in squeue and nf logs
    tag "${obsid}"

    if (params.pullPrep) {
        label "rclone"
    }
    // label "datamover" # memory limit 8G not enough

    label "rate_limit_20"
    label "nvme"
    // label "mem_full"

    time { 1.hour * Math.pow(task.attempt, 4) }
    disk { 50.GB * Math.pow(task.attempt, 4) }
    memory { 50.GB * Math.pow(task.attempt, 4) }

    // allow multiple retries
    maxRetries 2
    // exponential backoff: sleep for 5^attempt minutes after each fail
    errorStrategy {
        def reason = exitCodes[task.exitStatus] ?: "unknown"
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
        println "WARN ${result} task ${task.getClass()} ${task.hash}"
        return result
    }

    script:
    prefix = "birli_"
    suffix = ''
    args = [
        output: "uvfits"
    ]
    if (params.prep_time_res_s != null) {
        args.avg_time_res = params.prep_time_res_s
        suffix += "_${params.prep_time_res_s}s"
    }
    if (params.prep_freq_res_khz != null) {
        args.avg_freq_res = params.prep_freq_res_khz
        suffix += "_${params.prep_freq_res_khz}kHz"
    }
    if (params.prep_rfi != null && !params.prep_rfi) {
        args.no_rfi = "true"
        suffix += '_norfi'
    }
    argstr = args.collect { k, v -> ''+"${k}=${v}" }.join(',');
    uvfits = ''+"${prefix}${obsid}*${suffix}.uvfits"

    // echo \$'meta=${meta}'
    """
    #!/bin/bash -eux
    export MWA_ASVO_API_KEY="${params.asvo_api_key}"
    export SINGULARITY_CACHEDIR="${workflow.workDir?(workflow.workDir+'/.singularity'):'/tmp'}"
    ${params.proxy_prelude} # ensure proxy is set if needed

    """ + ( params.pullPrep ? """
    # download if available in accacia
    ${params.rclone} ls "${params.bucket_prefix}.prep" --include '*${obsid}*'
    ${params.rclone} copy "${params.bucket_prefix}.prep/" . --include '${uvfits}' --progress --stats-one-line
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
    input:
    tuple val(obsid), val(meta), path(metafits), path(uvfits)
    output:
    tuple val(obsid), val(meta), path(metrics)

    storeDir "${params.outdir}/${obsid}/prep_qa"

    tag "${base}"

    label "python"
    label "nvme"
    // label "mem_half"
    memory = {
        // max 30G Virtual for 60ts, 144T, 768chans
        [
            (30.GB * Math.pow(2,task.attempt) * (meta.ntimes?:60) * (meta.num_ants?:144) * (meta.num_chans?:768)
                / (60 * 144 * 768)),
            360.GB
        ].min()
    }
    time {20.minute * Math.pow(2,task.attempt)}
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'ignore' }
    maxRetries 2

    when: !params.noprepqa

    script:
    base = uvfits.baseName
    metrics = ''+"occupancy_${base}.json"
    template "flagqa.py"
}

process ssins {
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

    storeDir "${params.outdir}/${obsid}/prep_qa"

    tag "${base}"

    label "ssins"
    label "nvme"

    label "mem_half" // TODO: set these from uvfits size // "mem_full"
    time { 15.minute } // 2.hour
    // can't do this because no ntimes in prep meta
    // time { 30.min * ((meta.ntimes?:15) * (meta.nchans?:384) / (15 * 384) ) }

    when: !params.nossins

    script:
    base = uvfits.baseName
    output_prefix = meta.output_prefix ?: "${base}${meta.subobs?:''}_"
    args = [
        plot_title: meta.plot_title ?: "${base}${meta.subobs?:''}",
        sel_ants: meta.sel_ants ? meta.sel_ants.join(' ') : null,
        output_prefix: output_prefix,
        guard_width: meta.guard_width ?: (params.prep_freq_res_khz?:10) * 500,
        uvfits: uvfits
    ].findAll { _, v -> v != null }
        .collect { k, v -> ''+"""--${k} ${v} """ }
        .join(" ")
    template "ssins.py"
}

process absolve {
    input:
    tuple val(obsid), val(meta_), path(uvfits), path(mask)
    output:
    tuple val(obsid), val(meta), path(flagged)

    tag "${base}"

    label "ssins"
    label "nvme"
    label "mem_full"

    time { 30.minute }

    // errorStrategy "terminate"

    storeDir "${params.outdir}/${obsid}/prep"

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
    input:
    tuple val(obsid), val(meta), path(metafits), path(uvfits)
    output:
    tuple val(obsid), val(meta), path(autoplot)

    storeDir "${params.outdir}/${obsid}/prep_qa"

    tag "${base}${suffix}"

    label "python"
    label "nvme"
    label "mem_half"

    when: !params.noautoplot

    script:
    base = uvfits.baseName
    suffix = meta.suffix?:""
    autoplot = ''+"autoplot_${obsid}${suffix}.png"
    args = meta.autoplot_args
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
process rexCluster {
    input:
    tuple val(obsid), path(reduced)
    output:
    tuple val(obsid), path(cluster)

    storeDir "${params.outdir}/${obsid}/cal${params.cal_suffix}"

    tag "${obsid}"
    label "mwa_reduce"

    script:
    cluster = ''+"${obsid}_cluster_n${params.sub_nsrcs}_i${params.ionosub_nsrcs}.txt"
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
    ''+"declare -A ${name}=(" + map.collect { k, v -> "[${k}]=\"${v}\"".toString() }.join(" ") + ")".toString()
}

// calibrate with mwa reduce
process rexCalSol {
    input:
    tuple val(obsid), val(dical_args), path(metafits), path(vis), val(tile_flags)
    output:
    tuple val(obsid),
        path("rex_soln_${obsid}${meta.subobs?:''}_${name_glob}.bin"),
        path("rex_di-cal_${obsid}${meta.subobs?:''}_${name_glob}.log")

    storeDir "${params.outdir}/${obsid}/cal${params.cal_suffix}"

    tag "${obsid}${meta.subobs?:''}"

    label "mwa_reduce"
    label "nvme"
    label "mem_full"
    label "cpu_full"

    script:
    dical_names = dical_args.keySet().collect()
    para = dical_names.size() > 1
    name_glob = para ? "{" + dical_names.join(',') + "}" : dical_names[0]
    flag_args = tile_flags.size() > 0 ? (''+"--tile-flags ${tile_flags.join(' ')}") : ""
    vis_ms = "${uvfits.baseName}.ms"
    """
    #!/bin/bash -eux
    """ + groovy2bashAssocArray(dical_args, "dical_args") + """
    ${params.casa} -c "importuvfits('${vis}', '${vis_ms}')"
    singularity exec ${params.cotter_sif} fixmwams vis.ms ${metafits}

    for name in \${!dical_args[@]}; do
        export soln_name="rex_soln_${obsid}${meta.subobs?:''}_\${name}.bin"
        export log_name="rex_di-cal_${obsid}${meta.subobs?:''}_\${name}.log"
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
    tuple val(obsid), val(meta),
        path("hyp_soln_${obsid}${meta.subobs?:''}_${name_glob}.fits"),
        path("hyp_di-cal_${obsid}${meta.subobs?:''}_${name_glob}.log", optional: true)
    // todo: model subtract: path("hyp_model_${obsid}_${name_glob}.uvfits")

    storeDir "${params.outdir}/${obsid}/cal${params.cal_suffix}"

    tag "${obsid}${meta.subobs?:''}"

    // label jobs that need a bigger gpu allocation
    label "hyperdrive"
    label "mem_quarter" // TODO: set from uvfits size
    label "cpu_quarter"
    label "gpu_nvme"
    label "rate_limit_50"
    // if (params.pullCalSol) {
    //     label "rclone"
    // }

    time { 1.hour * Math.pow(task.attempt, 4) }
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'ignore' }
    maxRetries 2

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
        export soln_name="hyp_soln_${obsid}${meta.subobs?:''}_${name_prefix}\${name}.fits"
        export log_name="hyp_di-cal_${obsid}${meta.subobs?:''}_${name_prefix}\${name}.log"
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
    input:
    tuple val(obsid), val(meta_), path(soln)
    output:
    tuple val(obsid), val(meta), path(fit_soln) // , path(logs)

    storeDir "${params.outdir}/${obsid}/cal${params.cal_suffix}"
    tag "${obsid}${meta.subobs?:''}.${meta.dical_name}"

    label "fitcal"

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
    input:
    tuple val(obsid), val(meta_), path(metafits), path(vis), path(soln)
    output:
    tuple val(obsid), val(meta), path(cal_vis), path(logs)

    storeDir "${params.outdir}/${obsid}/cal${params.cal_suffix}"

    tag "${obsid}${meta.subobs?:''}.${meta.name}"
    label "hyperdrive"
    label "cpu_half"
    label "nvme"
    label "mem_half"

    when: !params.nouv && !params.noapply

    script:
    meta = mapMerge(meta_, [name: ''+"${meta_.name}_${meta_.apply_name}"])
    cal_vis = ''+"hyp_${obsid}${meta.subobs?:''}_${meta.name}.uvfits"
    logs = ''+"hyp_apply_${meta.name}.log"
    // echo \$'meta=${meta}'
    """
    #!/bin/bash -eux
    ${params.hyperdrive} solutions-apply ${meta.apply_args} \
        --time-average=${meta.time_res}s \
        --freq-average=${meta.freq_res}kHz \
        --data "${metafits}" "${vis}" \
        --solutions "${soln}" \
        --outputs "${cal_vis}" \
        ${meta.nodut1 ? "--ignore-dut1" : ""} \
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
    tuple val(obsid), val(meta_), path(metafits), path(vis), path(soln)
    output:
    tuple val(obsid), val(meta), path(cal_vis), path(logs)

    storeDir "${params.outdir}/${obsid}/cal${params.cal_suffix}"
    // storeDir "/data/curtin_mwaeor/FRB_hopper/"

    tag "${obsid}${meta.subobs?:''}.${meta.name}"
    label "hyperdrive"
    label "cpu_half"
    label "gpu"

    when: !params.noms && !params.noapply

    script:
    meta = mapMerge(meta_, [name: ''+"${meta_.name}_${meta_.apply_name}"])
    cal_vis = ''+"hyp_${obsid}${meta.subobs?:''}_${meta.name}.ms"
    logs = ''+"hyp_apply_${meta.name}_ms.log"
    // echo \$'meta=${meta}'
    """
    #!/bin/bash -eux
    ${params.hyperdrive} solutions-apply ${meta.apply_args} \
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
    input:
    tuple val(obsid), val(meta_), path(metafits), path(vis), path(srclist)
    output:
    tuple val(obsid), val(meta), path(sub_vis), path(logs)

    storeDir "${params.outdir}/${obsid}/cal${params.cal_suffix}"

    tag "${obsid}${meta.subobs?:''}.${meta_.name}"
    label "hyperdrive"
    label "cpu_half"
    label "gpu"

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
    input:
    tuple val(obsid), val(meta_), path(metafits), path(vis), path(srclist)
    output:
    tuple val(obsid), val(meta), path(sub_vis), path(logs)

    storeDir "${params.outdir}/${obsid}/cal${params.cal_suffix}"

    tag "${obsid}${meta.subobs?:''}.${meta_.name}"
    label "hyperdrive"
    label "cpu_half"
    label "gpu"

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
    input:
    tuple val(obsid), val(meta_), path(metafits), path(vis), path(srclist)
    output:
    tuple val(obsid), val(meta), path(sub_vis), path(json), path(logs)

    storeDir "${params.outdir}/${obsid}/cal${params.cal_suffix}"

    tag "${obsid}${meta_.subobs?:''}.${meta_.name}_i${meta_.ionosub_nsrcs}"
    label "hyperdrive"
    label "cpu_half"
    label "mem_full" // need a full node otherwise gpu runs out of memory
    label "gpu"
    label "rate_limit_50"

    time 2.5.hour

    when: !params.noionosub

    script:
    meta = mapMerge(meta_, [sub: "ionosub", name: ''+"ionosub_${meta_.name}_i${meta_.ionosub_nsrcs}"])
    sub_vis = ''+"hyp_${obsid}${meta.subobs?:''}_${meta.name}.uvfits"
    logs = ''+"hyp_vis-${meta.name}_uv.log"
    json = ''+"hyp_peel_${obsid}${meta.subobs?:''}_${meta.name}_uv.json"
    """
    #!/bin/bash -eux
    export RUST_BACKTRACE=1
    ${params.hyperdrive} peel ${params.hyp_peel_args} \
        ${params.hyp_srclist_args} \
        --data "${metafits}" "${vis}" \
        --beam-file "${params.beam_path}" \
        --source-list "${srclist}" \
        --iono-sub ${meta.ionosub_nsrcs} \
        --sub ${meta.sub_nsrcs} \
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
    tuple val(obsid), val(meta_), path(metafits), path(vis), path(srclist)
    output:
    tuple val(obsid), val(meta), path(sub_vis), path(json), path(logs)

    storeDir "${params.outdir}/${obsid}/cal${params.cal_suffix}"

    tag "${obsid}${meta_.subobs?:''}.${meta_.name}_i${meta_.ionosub_nsrcs}"
    label "hyperdrive"
    label "cpu_half"
    label "gpu"
    label "rate_limit_50"

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
    input:
    tuple val(obsid), val(meta), path(srclist), path(offsets)
    output:
    tuple val(obsid), val(meta), path("cthuluplot_${title}*.png"), path(csv), path(json), path("tec_${title}*.png")

    storeDir "${params.outdir}/${obsid}/iono_qa${params.cal_suffix}"

    tag "${obsid}${meta.subobs?:''}.${meta.name}"

    label "python"
    label "mem_half"

    time 1.hour

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
    input:
    tuple val(obsid), val(meta), path(metafits), path(uvfits)
    output:
    tuple val(obsid), val(meta), path(metrics)

    storeDir "${params.outdir}/${obsid}/prep_qa"

    tag "${base}"

    label "python"
    label "nvme"
    // label "mem_half"
    memory = {
        // max 50G Virtual for 60ts, 144T, 768chans
        [
            (50.GB * (meta.ntimes?:60) * (meta.num_ants?:144) * (meta.num_chans?:768)
                / (60 * 144 * 768)),
            360.GB
        ].min()
    }
    time {20.minute * task.attempt}
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'ignore' }
    maxRetries 2


    when: !params.noprepqa

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
    input:
    tuple val(obsid), val(meta), path(uvfits)
    output:
    tuple val(obsid), val(meta), path(json)

    storeDir "${params.outdir}/${obsid}/vis_qa${params.cal_suffix}"

    tag "${obsid}${meta.subobs?:''}.${meta.name}"

    label "python"

    script:
    base = uvfits.baseName
    json = ''+"${base}_vis_metrics.json"
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
        // tuple val(obsid), val(meta), path(vis), path(uvmeta), emit: obsMetaVisJson

    storeDir "${params.outdir}/${obsid}/vis_qa${params.cal_suffix}"

    tag "${base}"

    label "python"
    label "nvme"
    memory = 40.GB
    // memory = (vis.size * 3)
    time = {
        // max 25 min for 60ts, 144T, 768chans
        [
            (25.minute * (meta.ntimes?:60) * (meta.num_ants?:144) * (meta.num_chans?:768)
                / (60 * 128 * 768)),
            4.hour
        ].min()
    }
    // stageInMode "symlink"
    // can't symlink any more, scratch too slow

    script:
    base = vis.baseName
    uvmeta = ''+"uvmeta_${base}.json"
    template "uvmeta.py"
}

// QA tasks that can be run on each calibration parameter set
process calQA {
    input:
    tuple val(obsid), val(meta), path(metafits), path(soln)
    output:
    tuple val(obsid), val(meta), path(metrics)

    storeDir "${params.outdir}/${obsid}/cal_qa${params.cal_suffix}"

    tag "${obsid}${meta.subobs?:''}.${meta.name}"

    label "python"

    script:
    metrics = ''+"metrics_${soln.baseName}_X.json"
    """
    #!/bin/bash -eux
    run_calqa.py "${soln}" "${metafits}" --pol X --out "${metrics}"
    """
}

process phaseFits {
    input:
    tuple val(obsid), val(meta), path(metafits), path(soln)
    output:
    tuple val(obsid), val(meta), path("${obsid}* ${name} phase_fits.tsv"), path("${obsid}*${name}*.png")

    storeDir "${params.outdir}/${obsid}/cal_qa${params.cal_suffix}"

    tag "${obsid}.${meta.name}"

    label "mwax_mover"

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
    input:
    tuple val(obsid), val(meta), path(soln)
    output:
    tuple val(obsid), val(meta), path(metrics)

    storeDir "${params.outdir}/${obsid}/cal_qa${params.cal_suffix}"
    tag "${obsid}${meta.subobs?:''}.${meta.name}"

    label "python"

    script:
    metrics = ''+"${meta.cal_prog}_soln_${obsid}${meta.subobs?:''}_${meta.name}.fits.json"
    template "soljson.py"
}

process plotPrepVisQA {
    input:
    tuple val(obsid), val(meta), path(metrics)
    output:
    tuple val(obsid), val(meta), path("${base}_{rms,modz}.png")

    storeDir "${params.outdir}/${obsid}/prep_qa"

    tag "${obsid}${meta.subobs?:''}"

    label "python"
    time {15.minute * task.attempt}
    when: !params.noplotprepqa

    script:
    base = ''+"prepvis_metrics_${metrics.baseName}"
    """
    #!/bin/bash -eux
    plot_prepvisqa.py "${metrics}" --out "${base}.png" --save
    """
}

process plotSols {
    input:
    tuple val(obsid), val(meta), path(metafits), path(soln)
    output:
    tuple val(obsid), val(meta), path(plots_glob)

    storeDir "${params.outdir}/${obsid}/cal_qa${params.cal_suffix}"

    tag "${obsid}${meta.subobs?:''}.${meta.name}"

    label "hyperdrive"

    script:
    plots_glob = ''+"${meta.cal_prog}_soln_${obsid}${meta.subobs?:''}*_${meta.name}_{phases,amps}*.png"
    """
    ${params.hyperdrive} solutions-plot ${params.hyp_sols_plot_args} -m "${metafits}" ${soln}
    """
}

process plotCalQA {
    input:
    tuple val(obsid), val(meta), path(metrics)
    output:
    tuple val(obsid), val(meta), path("${plot_base}_{rms,fft}.png")

    storeDir "${params.outdir}/${obsid}/cal_qa${params.cal_suffix}"

    tag "${obsid}${meta.subobs?:''}.${meta.name}"

    label "python"

    when: !params.noplotcalqa

    script:
    plot_base = ''+"calmetrics_${obsid}${meta.subobs?:''}_${meta.name}"
    """
    #!/bin/bash -eux
    plot_calqa.py "${metrics}" --out "${plot_base}" --save
    """
}

process plotVisQA {
    input:
    tuple val(obsid), val(meta), path(metrics)
    output:
    tuple val(obsid), val(meta), path("${meta.cal_prog}_${obsid}_${meta.name}_vis_metrics_*.png")

    storeDir "${params.outdir}/${obsid}/vis_qa${params.cal_suffix}"

    tag "${obsid}"

    label "python"

    when: !params.noplotvisqa

    script:
    """
    #!/bin/bash -eux
    plot_visqa.py "${metrics}" --out "${meta.cal_prog}_${obsid}_${meta.name}_vis_metrics.png" --save
    """
}

process plotImgQA {
    input:
    tuple val(name), path("??????????.json")
    output:
    tuple val(name), path("${base}_*.png")

    storeDir "${results_dir}${params.img_suffix}${params.cal_suffix}"
    stageInMode "symlink"

    tag "${name}"

    label "python"

    script:
    base = ''+"wsclean_hyp_${name}-MFS"
    """
    #!/bin/bash
    set -ex
    plot_imgqa.py --out "${base}" --save ??????????.json
    """
}

process delaySpec {
    input:
    tuple val(obsid), val(meta), path(uvfits)
    output:
    tuple val(obsid), val(meta), path(dlyspec)

    storeDir "${params.outdir}/${obsid}/vis_qa${params.cal_suffix}"

    tag "${obsid}${meta.subobs?:''}.${meta.name}"

    label "python"

    when: !params.nodelayspec

    script:
    title = ''+"${obsid}${meta.subobs?:''}_${meta.name}"
    dlyspec = ''+"dlyspec_${title}.png"
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

    tag "${obsid}${meta.subobs?:''}${meta.inter_tok?:''}.${meta.name}"
    label "wsclean"
    label "cpu_half"
    label "mem_half"
    label "nvme"

    time { 1.8.minute * (1 + (multiplier * pix_mult * chan_mult * inter_mult)) }

    script:
    multiplier = vis.collect().size()
    mult_suffix = multiplier > 1 ? (''+"_x${multiplier}") : ""
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
    if (multiinterval && !meta.inter_tok) {
        img_glob += "-t????"
    }
    if (multichannel) {
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
process wscleanDConv {
    input:
    tuple val(obsid), val(meta), path(vis), val(img_params), path(dirtyImgs)
    output:
    tuple val(obsid), val(meta), path(img_glob)
        // path("wsclean_${img_name}-sources.txt") <- only works for stokes I

    storeDir "${params.outdir}/${obsid}/img${img_params.suffix}${params.cal_suffix}"

    tag "${obsid}${meta.subobs?:''}${meta.inter_tok?:''}.${meta.name}.${(img_params.pol?:'').split(' ')[0]}"
    label "wsclean"
    label "cpu_quarter"
    label "nvme"
    stageInMode "symlink"

    // cpus { max(4, min(36,1 + (Math.pow(task.attempt, 2) * multiplier))) as int }
    time { [23.hour, 1.hour * (1 + (Math.pow(task.attempt, 2) * multiplier))].min() }
    memory { [350.GB, 35.GB * (1 + (Math.pow(task.attempt, 2) * multiplier))].min() }
    maxRetries 3
    maxForks 60

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
    multiplier *= iter_mult = 1 + Math.sqrt(img_params.niter as Double) / 1000
    vis_ms = vis.collect {''+"${it.baseName}.ms"}
    vis = vis.collect()
    img_glob = ''+(img_params.glob?:'')
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
    echo \$'img_params=${img_params}, img_glob=${img_glob}'
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
    input:
    tuple val(obsid), val(meta), path(fits)
    output:
    tuple val(obsid), val(meta), path(csv)

    storeDir "${params.outdir}/${obsid}/img_qa${params.img_suffix}${params.cal_suffix}"

    tag "${obsid}${meta.subobs?:''}.${meta.name}.${meta.suffix}"

    label "python"

    time {5.min * task.attempt}

    script:
    csv = "quantile_${fits.baseName}.csv"
    template "img_meta.py"
}

// power spectrum metrics via chips
process psMetrics {
    input:
    tuple val(obsid), val(meta), path("${meta.cal_prog}_${obsid}${meta.subobs?:''}_${meta.name}.uvfits")
    output:
    tuple val(obsid), val(meta), path(out_metrics)

    storeDir "${params.outdir}/${obsid}/ps_metrics${params.cal_suffix}"

    tag "${obsid}${meta.subobs?:''}.${meta.name}"

    label "chips"
    label "cpu_half"

    time 40.minute

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
    input:
    tuple val(group), val(obsids)
    output:
    path(manifest)

    tag "${group}"

    time {5.minute * task.attempt}

    storeDir "${params.outdir}/${group}"

    script:
    manifest = ''+"manifest_${group}.csv"
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
    tuple val(group), val(meta), \
        path("vis_tot_${pol}.${ext}.dat"), \
        path("vis_diff_${pol}.${ext}.dat"), \
        path("noise_tot_${pol}.${ext}.dat"), \
        path("noise_diff_${pol}.${ext}.dat"), \
        path("weights_${pol}.${ext}.dat")
        // path("{vis_tot,vis_diff,noise_tot,noise_diff,weights}_${pol}.${ext}.dat")

    storeDir "${params.outdir}/${group}/ps_metrics${params.cal_suffix}/${meta.name}"

    tag "${group}.${meta.name}.${pol}"

    label "chips"
    label "cpu_half"
    label "mem_half"
    label "nvme"

    maxRetries 3
    errorStrategy { return task.exitStatus > 1 ? "retry" : "ignore" }

    time { 30.minute * obsids.size() }

    when: (!params.nopowerspec && !params.nochips)

    script:

    // meta = deepcopy(meta_)
    lowfreq = meta.lowfreq
    ext = meta.ext
    nchans = meta.nchans
    eorband = meta.eorband
    eorfield = meta.eorfield
    period = meta.period?:"8.0"
    freq_res = "${(meta.freq_res?:80) * 1000}"
    pol = meta.pol?:"xx"
    freq_idx_start = meta.freq_idx_start?:0

    // echo "meta=${meta}"
    """
    #!/bin/bash -ux
    """ + (eorfield == null || eorband == null || !nchans || !ext || lowfreq == null || !freq_res || freq_idx_start == null || period == null ? "exit 2" : "" ) + """
    export DATADIR="\$PWD" INPUTDIR="\$PWD/" OUTPUTDIR="\$PWD/" OBSDIR="\$PWD/"
    export OMP_NUM_THREADS=${task.cpus}

    # copy data files
    for f in present_vals.dat missing_vals.dat krig_weights.dat; do
        [ -f \$f ] || cp "/astro/mwaeor/ctrott/output/\$f" .
    done

    which gridvisdiff

    """ + (
        [coerceList(obsids), coerceList(viss)].transpose().collect { obsid, vis ->
            """
    gridvisdiff "${vis}" "${obsid}" "${ext}" "${eorband}" -f "${eorfield}"
    export ps=\$?
    if [ \${ps} -ne 0 ] && [ \${ps} -ne 42 ]; then
        echo "gridvisdiff failed. status=\${ps}"
        exit \${ps}
    fi
    """
        }.join("\n")
    ) + """
    which prepare_diff
    prepare_diff "${ext}" "${nchans}" "${freq_idx_start}" "${pol}" "${ext}" "${eorband}" -p "${period}" -c "${freq_res}" -n "${lowfreq}"
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
    """
}

process chipsCombine {
    input:
    tuple val(group), val(meta), val(exts), path(grids)

    output:
    tuple val(group), val(meta), path("combine_${pol}.${combined_ext}.txt"),
        path("{vis_tot,vis_diff,noise_tot,noise_diff,weights}_${pol}.${combined_ext}.dat")

    storeDir "${params.outdir}/${group}/ps_metrics${params.cal_suffix}/${meta.name}"

    tag "${group}.${pol}.${meta.name}"

    label "chips"
    label "cpu_full"
    label "mem_full"
    label "nvme_full"
    stageInMode { (exts.size() > 80) ? "symlink" : "copy" }

    // takes about 6 hours for 300, 16 hours for 1000,
    // double it for safety
    // clamp 1h < x < 24h
    time {
        t = 1.hour
        if (exts.size() < 1000) {
            t = [1.hour, 6.hour * (exts.size() - 0.5) / 300].max()
        } else {
            t = [32.hour, 16.hour * (exts.size() - 0.5) / 1000].min()
        }
        [[t * 2, 1.hour].max(), 24.95.hour].min()
    }

    script:
    // ext = "${group}_${meta.name}"
    combined_ext = meta.ext
    pol = meta.pol?:"xx"
    nchans = meta.nchans
    // echo "meta=${meta}"
    """
    #!/bin/bash -eux
    """ + (!nchans || !combined_ext || exts.size() == 0 ? "exit 2" : "" ) + """

    export DATADIR="\$PWD"
    export INPUTDIR="\$PWD/"
    export OUTPUTDIR="\$PWD/"
    export OBSDIR="\$PWD/"
    export OMP_NUM_THREADS=${task.cpus}


    echo -n "" > combine_${pol}.${combined_ext}.txt
    for ext in ${exts.join(' ')}; do
        # check for nulls
        for dat in *_${pol}.\${ext}.dat; do
            ls -al \$dat
            if ! grep -m1 -qP "[^\\0]" \$dat; then
                echo "oops, all nulls \$dat"
                exit 3
            fi
        done
        echo "${pol}.\${ext}" >> combine_${pol}.${combined_ext}.txt
    done

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

    exit 0
    """
}

process chipsLssa {
    input:
    tuple val(group), val(meta_), path(grid)

    output:
    tuple val(group), val(meta),
        path("{crosspower,residpower,residpowerimag,totpower,flagpower,fg_num,outputweights}_${pol}_${bias_mode}.iter.${ext}.dat")

    storeDir "${params.outdir}/${group}/${params.lssa_bin}${params.cal_suffix}/${meta.name}"

    tag "${group}.${params.lssa_bin}.${meta.name}"

    label "chips"
    label "cpu_full"
    label "mem_full"
    label "nvme"

    time 1.hour

    script:
    maxu = params.lssa_maxu
    nbins = params.lssa_nbins
    meta = mapMerge(meta_, [nbins: nbins, maxu: maxu])
    ext = meta.ext
    nchans = meta.nchans
    eorband = meta.eorband
    pol = meta.pol?:"xx"
    freq_idx_start = 0
    bias_mode = params.lssa_bias_mode

    // echo "meta=${meta}"
    """
    #!/bin/bash -eux
    """ + (eorband == null || !nchans || nbins == null || !ext || maxu == null || bias_mode == null ? "exit 2" : "" ) + """

    export DATADIR="\$PWD"
    export INPUTDIR="\$PWD/"
    export OUTPUTDIR="\$PWD/"
    export OBSDIR="\$PWD/"
    export OMP_NUM_THREADS=${task.cpus}

    ${params.lssa_bin} "${ext}" "${nchans}" "${nbins}" "${pol}" "${maxu}" "${ext}" "${bias_mode}" "${eorband}" \
        2>&1 | tee syslog_lssa_simple_${pol}.txt

    export cross_size="\$(stat -c%s crosspower_${pol}_${bias_mode}.iter.${ext}.dat)"
    if (( cross_size < 4097 )); then
        echo "crosspower_${pol}_${bias_mode}.iter.${ext}.dat is too small (\$cross_size), exiting"
        exit 3
    fi
    """
}

process chipsPlot {
    input:
    tuple val(group), val(meta), path(grid)
    output:
    tuple val(group), val(meta), path("chips${dims}D_${pol}_${suffix}.png")

    storeDir "${params.outdir}/${group}/${params.lssa_bin}${params.cal_suffix}/${meta.name}"

    tag "${group}.${meta.name}.${ptype}"

    label "python"

    time 15.minute

    script:
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
        suffix += "_delta"
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


process chips1d_tsv {
    input:
    tuple val(group), val(meta), path(grid)
    output:
    tuple val(group), val(meta), path("1D_power_${pol}.${suffix}.tsv")

    storeDir "${params.outdir}/${group}/${params.lssa_bin}${params.cal_suffix}/${meta.name}"

    tag "${group}.${meta.name}.${pol}"

    label "chips_wrappers"
    label "nvme"

    time 15.minute

    script:
    pol = meta.pol?:"null"

    suffix = "${meta.ext}"
    basedir = "/astro/mwaeor/dev/nfdata/eor0high_phase2a-128T_lst+04_06a3107b/ps_metrics/ionosub_30l_src4k_300it_8s_80kHz_i1000./"
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
        kperp_max: meta.kperp_max,
        kperp_min: meta.kperp_min,
        kparra_min: meta.kparra_min,

        // chips group
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
    if (meta.nchans) {
        args += " --N_chan ${meta.nchans}"
    }
    if (meta.kperp) {
        args += " --K_perp ${meta.kperp}"
    }

    // template "chips1D_tsv.py"
    script:
    """
    chips1D_tsv.py ${args}
    cp "1D_power_${pol}.tsv" "1D_power_${pol}.${suffix}.tsv"
    """
}

// analyse images of V,XX,YY
process imgQA {
    input:
    tuple val(obsid), val(meta), path(fits)
    output:
    tuple val(obsid), val(meta), path(json)

    storeDir "${params.outdir}/${obsid}/img_qa${params.img_suffix}${params.cal_suffix}"

    tag "${obsid}${meta.subobs?:''}.${meta.name}"

    label "python"
    time {5.minute * task.attempt}

    script:
    json = "wsclean_hyp_${obsid}${meta.subobs?:''}_${meta.name}-MFS.json"
    """
    #!/bin/bash -eux
    run_imgqa.py ${fits.join(' ')} --out ${json}
    """
}

// analyse images of V,XX,YY
process krVis {
    input:
    tuple val(obsid), val(meta), path(fits), path(roma0)
    output:
    tuple val(obsid), val(meta), path("${out_prefix}{I,L,V,IV}-*.png")

    // errorStrategy "terminate"

    storeDir "${params.outdir}/${obsid}/img_qa${params.img_suffix}${params.cal_suffix}"

    tag "${obsid}${meta.subobs?:''}${meta.inter_tok?:''}.${meta.name}"

    label "python"

    script:
    out_prefix = "${obsid}${meta.subobs?:''}${meta.inter_tok?"-${meta.inter_tok}":''}_${meta.name}-"
    template "krvis.py"
}



process uvPlot {
    input:
    tuple val(obsid), val(meta), path(uvfits)
    output:
    tuple val(obsid), val(meta), path("${uvplot}_{XX,YY,XY,YX}.png")

    storeDir "${params.outdir}/${obsid}/vis_qa${params.cal_suffix}"

    tag "${obsid}${meta.subobs?:''}.${meta.name}"

    label "python"

    when: !params.nouvplot

    script:
    title = "${obsid}${meta.subobs?:''}_${meta.name}"
    uvplot = "uvplot_${uvfits.baseName}"
    template "uvplot_2d.py"
}

// make a thumbnail png from a fits image
process thumbnail {
    input:
    tuple val(obsid), val(meta), path(img)
    output:
    tuple val(obsid), val(meta), path(thumb)

    storeDir "${params.outdir}/${obsid}/img_qa${params.img_suffix}${params.cal_suffix}"

    tag "${obsid}${meta.subobs?:''}.${meta.name}.${meta.suffix}"

    label "python"
    time {5.min * task.attempt}

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
    tuple val(obsid), val(meta_), path(fits)
    output:
    tuple val(obsid), val(meta), path(polcomp)

    storeDir "${params.outdir}/${obsid}/img_qa${params.img_suffix}${params.cal_suffix}"

    tag "${obsid}${meta.subobs?:''}${meta.subobs?:''}.${meta.name}.${meta.prod}.${meta.orderName}"

    label "python"

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
        tuple val(name), path("??????????.json")
    output:
        path("cal_qa*.png")

    storeDir "${results_dir}${params.cal_suffix}"
    stageInMode "symlink"

    label "python"
    tag "${name}"

    when: !params.noplotcalqa

    script:
    """
    #!/bin/bash -eux
    plot_caljson.py --out cal_qa_${name} --save ??????????.json
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

// process updateQuality {
//     input:
//         tuple val(obsid), val(name), val(params)
//     output:
//         tuple val(obsid), val(name), path("quality_${name}.json")

//     storeDir "${params.outdir}/${obsid}/meta"

//     tag "${obsid}.${name}"
// }

process tsvScatterPlot {
    input:
        tuple val(meta), path(tsv)
    output:
        path(plot)

    storeDir "${results_dir}${params.img_suffix}${params.cal_suffix}"
    stageInMode "copy"

    label "python"

    tag "${meta.name}"

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
    input:
        tuple val(bucket_suffix), path(x)
    output:
        path("${x}.shadow")

    tag "$x"
    storeDir "${params.outdir}/.fakebuckets/${bucket}"

    // label "rclone"
    label "datamover"

    // TODO: use size of File x to determine memory
    // memory = (x.size() * 1.5)

    script:
    bucket = "${params.bucket_prefix}.${bucket_suffix}"
    """
    #!/bin/bash -eux
    touch "${x}.shadow"
    ${params.proxy_prelude} # ensure proxy is set if needed
    rclone mkdir "${bucket}"
    rclone copy --copy-links "$x" "${bucket}"
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

// import groovy.json.JsonSlurper
import groovy.json.JsonSlurperClassic
import groovy.json.StringEscapeUtils
jslurp = new JsonSlurperClassic()

def parseJson(path) {
    // TODO: fix nasty hack to deal with NaNs
    try {
        // try reading the file first, it may not be ready yet.
        File(path).size()
    } catch (Exception _) {
        sleep (5)
    }
    try {
        def text = path.getText().replaceAll(/(NaN|-?Infinity)/, '"$1"')
        def parsed = jslurp.parseText(text)
        return deepcopy(parsed)
    } catch (Exception e) {
        println("error parsing ${path} ${e}")
    }
}

def parseCsv(path, header = true, skipBytes = 0, delimeter=',') {
    def allLines = path.readLines()
    if (!header) {
        return allLines.collect { it.split(delimeter) }
    } else {
        def head = allLines[0][skipBytes..-1].split(delimeter)
        def body = allLines[1..-1]
        def parsed = body.collect { row ->
            [head, row.split(delimeter)].transpose().collectEntries()
        }
        return deepcopy(parsed)
    }
}

def parseFloatOrNaN(s) {
    // apparently there's a difference between `float` numbers and `Float` objects.
    // Float.valueOf(String) and Float.parseFloat(String) should be equivalent
    if (s == null || s == 'NaN') {
        return Float.NaN
    } else if (s instanceof Float) {
        return s
    } else if (s instanceof java.math.BigDecimal) {
        return s.floatValue()
    }
    try {
        return Float.parseFloat(s) // .toFloat()? .floatValue()? Float.valueOf? WHY IS THIS SO HARD
    } catch (NumberFormatException e) {
        return Float.NaN
    }
}

SimpleDateFormat logDateFmt = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss");

def isNaN(n) {
    if (n==null || n=="") {
        return true
    }
    // def cls = null
    try {
        // cls = n.getClass();
        return (n as float).isNaN()
    } catch (groovy.lang.MissingMethodException e) {
        // print("WARN isNaN(${n})<${cls}>: ${e}")
        StackTraceUtils.sanitize(e).printStackTrace()
        throw e
        return true
    }
}

// display a long list of ints, replace bursts of consecutive numbers with ranges
// start, end, delim
def displayRange(Integer s, Integer e, String d=',') {
    return s == e ? "${s}${d}" : s == e - 1 ? "${s}${d}${e}${d}" : "${s}-${e}${d}"
}

max_ints_display = 50
mid_ints_display = (max_ints_display-1).intdiv(2)

def displayInts(l_, delim=',') {
    def l = (l_ as ArrayList).sort(false).unique()
    switch (l) {
        case { l.size == 0 }: return "";
        case { l.size == 1 }: return "${l[0]}";
        default:
            def (sb, start, end) = [''<<'', l[0], l[0]]
            for (i in l[1..-1]) {
                (sb, start, end) = i == end + 1 ? [sb, start, i] : [sb << displayRange(start, end, delim), i, i]
            }
            def result = (sb << displayRange(start, end, delim))[0..-2].toString()
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
            def (sb, start, end) = [[], l[0], l[0]]
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

def prepqa_pass(flagMeta) {
    def reasons = []
    def fail_code = 0x00 // no error
    if (params.flag_occupancy_threshold != null && flagMeta.total_occ > params.flag_occupancy_threshold) {
        reasons += "total_occ(${String.format('%.2f', flagMeta.total_occ)})>${params.flag_occupancy_threshold}"
        fail_code = fail_code==0x00 ? 0x31 : fail_code
    }
    if (params.rfi_occupancy_threshold != null && flagMeta.total_non_preflagged_bl_occ > params.rfi_occupancy_threshold) {
        reasons += "rfi_occ(${String.format('%.2f', flagMeta.total_non_preflagged_bl_occ)})>${params.rfi_occupancy_threshold}"
        fail_code = fail_code==0x00 ? 0x32 : fail_code
    }
    if (params.ssins_occupancy_threshold != null && flagMeta.ssins_total > params.ssins_occupancy_threshold) {
        reasons += "ssins_total(${String.format('%.2f', flagMeta.ssins_total)})>${params.ssins_occupancy_threshold}"
        fail_code = fail_code==0x00 ? 0x33 : fail_code
    }
    if (params.ssins_streak_threshold != null && flagMeta.ssins_streak > params.ssins_streak_threshold) {
        reasons += "ssins_streak(${String.format('%.2f', flagMeta.ssins_streak)})>${params.ssins_streak_threshold}"
        fail_code = fail_code==0x00 ? 0x34 : fail_code
    }
    if (params.ssins_narrow_threshold != null && flagMeta.ssins_narrow_total > params.ssins_narrow_threshold) {
        reasons += "ssins_narrow_total(${String.format('%.2f', flagMeta.ssins_narrow_total)})>${params.ssins_narrow_threshold}"
        fail_code = fail_code==0x00 ? 0x35 : fail_code
    }
    if (params.ssins_dab_threshold != null && flagMeta.ssins_dab_total > params.ssins_dab_threshold) {
        reasons += "ssins_dab_total(${String.format('%.2f', flagMeta.ssins_dab_total)})>${params.ssins_dab_threshold}"
        fail_code = fail_code==0x00 ? 0x36 : fail_code
    }
    if (!params.noprepqafilter && (flagMeta.prep_status?:"GOOD") != "GOOD") {
        reasons += "prep_status=${flagMeta.prep_status}"
        fail_code = fail_code==0x00 ? 0x37 : fail_code
    }

    if (params.filter_bad_tile_frac != null && flagMeta.n_tiles != null && (flagMeta.prepFlags?:[]).size() > 0) {
        // print("prepqa_pass filter_bad_tile_frac=${params.filter_bad_tile_frac} n_tiles=${flagMeta.n_tiles} flagAnts=${flagMeta.flagAnts?:[]} prepFlags=${flagMeta.prepFlags?:[]}")
        def n_tiles = flagMeta.n_tiles
        def n_bad_tiles = ((((flagMeta.flagAnts?:[]) as Set) + ((flagMeta.prepFlags?:[]) as Set)) as ArrayList).size()
        if (n_bad_tiles > params.filter_bad_tile_frac * n_tiles) {
            reasons += "n_bad_tiles(${n_bad_tiles}) > ${params.filter_bad_tile_frac} * n_tiles(${n_tiles})"
            fail_code = fail_code==0x00 ? 0x3F : fail_code
        }
    }

    def reason = null
    if (reasons.size() > 0) {
        reason = reasons.join('|')
    }
    return [fail_code, reason]
}

def calqa_pass(meta) {
    def reasons = []
    def fail_code = 0x00
    if (params.filter_max_rms_convg != null && meta.rms_convg != null) {
        if (meta.rms_convg > params.filter_max_rms_convg) {
            reasons += "rms_convg(${meta.rms_convg}) > ${params.filter_max_rms_convg}"
            fail_code = fail_code==0x00 ? 0x40 : fail_code
        }
    }
    // TODO: this
    // if (params.filter_bad_tile_frac != null && flagMeta.n_tiles != null && (flagMeta.prepFlags?:[]).size() > 0) {
    //     n_tiles = flagMeta.n_tiles
    //     n_bad_tiles = ((Set(flagMeta.bad_ants?:[]) + Set(flagMeta.prepFlags?:[])) as ArrayList).size()
    //     if (n_bad_tiles > params.filter_bad_tile_frac * n_tiles) {
    //         reasons += "n_bad_tiles(${n_bad_tiles}) > ${params.filter_bad_tile_frac} * n_tiles(${n_tiles})"
    //         fail_code = fail_code==0x00 ? 0x3F : fail_code
    //     }
    // }
    def reason = null
    if (reasons.size() > 0) {
        reason = reasons.join('|')
    }
    return [fail_code, reasons]
}

// POWER	P_win_sub, P_win	< 20	Small window power overall
def cmt_ps_metrics_pass(meta) {
    def fail_code = 0x00
    if (params.filter_max_ps_window != null && meta.p_window != null) {
        if (meta.p_window > params.filter_max_ps_window) {
            fail_code = ((meta.sub?:"")=="") ? 0x60 : 0x61
            return [fail_code, "${meta.sub?:""}_p_win(${meta.p_window}) > max_ps_window(${params.filter_max_ps_window})"]
        }
    }
    return [fail_code, null]
}

// POWER	Normal P_win/P_wg	< 0.1	Window power a small fraction of wedge power
// POWER	P_wg_sub/P_wg	< 0.3	More than 70% wedge power subtracted
// POWER	P_win_sub/P_win	< 1.0, >0.1	Window power not crazy after subtraction
def cmt_ps_metrics_pass_sub(nosubMeta, subMeta) {
    def (fail_code, subReason) = cmt_ps_metrics_pass(subMeta)
    if (fail_code != 0x00) {
        return [fail_code, subReason]
    }
    if (params.filter_max_ps_ratio != null && nosubMeta.p_window != null && nosubMeta.p_wedge != null) {
        p_win_p_wg = nosubMeta.p_window / nosubMeta.p_wedge
        if (p_win_p_wg > params.filter_max_ps_ratio) {
            return [0x62, "p_win:p_wg(${p_win_p_wg}) > max_ps_ratio(${params.filter_max_ps_ratio})"]
        }
    }
    if (params.filter_max_ps_wedge_sub != null && subMeta.p_wedge != null && nosubMeta.p_wedge != null) {
        sub_p_wg = subMeta.p_wedge / nosubMeta.p_wedge
        if (sub_p_wg > params.filter_max_ps_wedge_sub) {
            return [0x64, "sub_p_wg:p_wg(${sub_p_wg}) > max_ps_wedge_sub{${params.filter_max_ps_wedge_sub}}"]
        }
    }
    if (subMeta.p_window != null && nosubMeta.p_window != null) {
        sub_p_win = subMeta.p_window / nosubMeta.p_window
        if (params.filter_min_ps_win_sub != null && sub_p_win < params.filter_min_ps_win_sub) {
            return [0x66, "sub_p_win:p_win(${sub_p_win}) < min_ps_win_sub(${params.filter_min_ps_win_sub})"]
        }
        if (params.filter_max_ps_win_sub != null && sub_p_win > params.filter_max_ps_win_sub) {
            return [0x67, "sub_p_win:p_win(${sub_p_win}) > max_ps_win_sub(${params.filter_max_ps_win_sub})"]
        }
    }
    return [fail_code, subReason]
}

def cmt_imgqa_pass(meta) {
    def fail_code = 0x00
    if (params.filter_max_vrms_box!= null && meta.v_rms_box != null) {
        if (meta.v_rms_box > params.filter_max_vrms_box) {
            fail_code = ((meta.sub?:"")=="") ? 0x70 : 0x71
            return [fail_code, "${meta.sub?:""}_v_rms_box(${meta.v_rms_box}) > max_vrms_box(${params.filter_max_vrms_box})"]
        }
    }
    if (params.filter_max_pks_int_v_ratio != null && meta.xx_pks_int != null && meta.yy_pks_int != null && meta.v_pks_int != null) {
        pks_int_v_ratio = meta.v_pks_int / (meta.xx_pks_int + meta.yy_pks_int)
        if (pks_int_v_ratio > params.filter_max_pks_int_v_ratio) {
            fail_code = ((meta.sub?:"")=="") ? 0x72 : 0x73
            return [fail_code, "${meta.sub?:""}_pks_int v/(xx+yy) (${pks_int_v_ratio}) > max_pks_int_v_ratio(${params.filter_max_pks_int_v_ratio})"]
        }
    }
    if (params.filter_max_pks_int_diff != null && meta.xx_pks_int != null && meta.yy_pks_int != null) {
        pks_int_diff = (meta.xx_pks_int - meta.yy_pks_int).abs()
        if (pks_int_diff > params.filter_max_pks_int_diff) {
            fail_code = ((meta.sub?:"")=="") ? 0x74 : 0x75
            return [fail_code, "${meta.sub?:""}_pks_int |xx-yy| (${pks_int_diff}) > max_pks_int_diff(${params.filter_max_pks_int_diff})"]
        }
    }
    return [fail_code, null]
}

// IMG	Vrms box	< 0.05	RMS V should be small
// IMG	V/(XX+YY) PKS int	< 0.001	V should be small compared with XX and YY
// IMG	PKS XX and YY	|XX-YY| < 10.0	XX and YY integrated should be similar
// IMG	XX_sub/XX integ	< 0.2	Most flux subtracted
// IMG	YY_sub/YY integ	< 0.2
// IMG	XX_sub integ	< 0.5	Integrated remaining flux after subtraction is small
// IMG	YY_sub integ	< 0.5	Integrated remaining flux after subtraction is small
def cmt_imgqa_pass_sub(nosubMeta, subMeta) {
    def fail_code = 0x00
    // if (params.filter_max_pks_int_v_ratio != null && nosubMeta.xx_pks_int != null && nosubMeta.yy_pks_int != null && nosubMeta.v_pks_int != null) {
    //     pks_int_v_ratio = nosubMeta.v_pks_int / (nosubMeta.xx_pks_int + nosubMeta.yy_pks_int)
    //     if (pks_int_v_ratio > params.filter_max_pks_int_v_ratio) {
    //         return [0x72, "pks_int v/(xx+yy) (${pks_int_v_ratio}) > max_pks_int_v_ratio(${params.filter_max_pks_int_v_ratio})"]
    //     }
    // }
    if (params.filter_max_pks_int_sub_ratio != null && nosubMeta.xx_pks_int != null && subMeta.xx_pks_int != null) {
        pks_int_sub_ratio = subMeta.xx_pks_int / nosubMeta.xx_pks_int
        if (pks_int_sub_ratio > params.filter_max_pks_int_sub_ratio) {
            return [0x76, "${subMeta.sub?:""}_xx_pks_int:xx_pks_int (${pks_int_sub_ratio}) > max_pks_int_sub_ratio(${params.filter_max_pks_int_sub_ratio})"]
        }
    }
    if (params.filter_max_pks_int_sub_ratio != null && nosubMeta.yy_pks_int != null && subMeta.yy_pks_int != null) {
        pks_int_sub_ratio = subMeta.yy_pks_int / nosubMeta.yy_pks_int
        if (pks_int_sub_ratio > params.filter_max_pks_int_sub_ratio) {
            return [0x77, "${subMeta.sub?:""}_yy_pks_int:yy_pks_int (${pks_int_sub_ratio}) > max_pks_int_sub_ratio(${params.filter_max_pks_int_sub_ratio})"]
        }
    }
    if (params.filter_max_pks_int_sub != null && subMeta.xx_pks_int != null) {
        if (subMeta.xx_pks_int > params.filter_max_pks_int_sub) {
            return [0x78, "${subMeta.sub?:""}_xx_pks_int(${subMeta.v_rms_box}) > max_pks_int_sub(${params.filter_max_pks_int_sub})"]
        }
    }
    if (params.filter_max_pks_int_sub != null && subMeta.yy_pks_int != null) {
        if (subMeta.yy_pks_int > params.filter_max_pks_int_sub) {
            return [0x79, "${subMeta.sub?:""}_yy_pks_int(${subMeta.v_rms_box}) > max_pks_int_sub(${params.filter_max_pks_int_sub})"]
        }
    }
    return [fail_code, null]
}

// decompose an image filename into:
// - interval if multiple intervals or -1 if combined
// - channel of multiple channels or -1 if combined
// - polarization
// - image product name (dirty, image, model, uv-{real,imag})
def decomposeImg(img) {
    def tokens = img.baseName.split('-').collect()
    // defaults:
    def meta = [chan: -1, inter: -1, chan_tok: "MFS"]
    // this handles the case where product is "uv-{real,imag}"
    if (tokens.size > 1 && tokens[-2] == "uv") {
        tokens = tokens[0..-3] + ["uv-${tokens[-1]}"]
    }
    meta.prod = tokens.removeLast()
    if (tokens.size > 0 && tokens[-1] != "MFS") {
        meta.pol = tokens.removeLast()
    }
    // channel is only present in multi-frequency imaging
    if (multichannel && (chan_tok = tokens.removeLast()) != "MFS") {
        meta.chan_tok = chan_tok
        try {
            meta.chan = (chan_tok as int)
        } catch (java.lang.NumberFormatException e) {
            print("error parsing channel ${chan_tok} from ${img}. \n${meta}")
            throw e
        }
    }
    // suffix without interval
    meta.inter_suffix = [meta.chan_tok, meta.pol, meta.prod].join('-')
    meta.suffix = meta.inter_suffix
    if (multiinterval && (inter_tok = tokens.removeLast()) =~ /t\d{4}/) {
        meta.inter_tok = inter_tok
        meta.inter = inter_tok[1..-1] as int
        meta.suffix = "${inter_tok}-${meta.inter_suffix}"
    }
    deepcopy(meta)
}

def groupMeta(meta) {
    def newMeta = [:]
    newMeta.sort = 1
    if (meta.p_window?:Float.NaN != Float.NaN && meta.p_wedge?:Float.NaN != Float.NaN) {
        newMeta.sort *= parseFloatOrNaN(meta.p_window) / parseFloatOrNaN(meta.p_wedge)
    }
    if (meta.total_weight?:Float.NaN != Float.NaN) {
        weight = parseFloatOrNaN(meta.total_weight) / 1e9
        newMeta.sort /= weight
    }
    // group by field, band, config
    def group_tokens = []
    def first_token = ""
    if (meta.eorfield != null) {
        first_token += "eor${meta.eorfield}"
    }
    if (meta.eorband != null) {
        first_token += (meta.eorband == 0 ? "low" : "high")
    }
    if (first_token.size() > 0) {
        group_tokens << first_token
    }
    if (meta.config != null) {
        group_tokens << meta.config
    }
    if (params.groupByPointing && meta.ew_pointing != null) {
        group_tokens << sprintf("ewp%+1d", meta.ew_pointing)
    }
    if (params.groupByLst && meta.lst != null) {
        nearest_lst = ((meta.lst.round().intValue() + 180) % 360 - 180)
        group_tokens << sprintf("lst%+03d", nearest_lst)
    }

    newMeta.group = group_tokens.join('_')

    [
        "cal_prog", "time_res", "freq_res", "lowfreq", "nchans",
        "eorband", "eorfield", "config", "name", "sub"
    ].each { key ->
        if (meta[key] != null) {
            newMeta[key] = meta[key]
        }
    }
    deepcopy(newMeta)
}

// mapping from eor1 pointing to sweet pointing
eor12sweet = [
    0: 0, // [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
    1: 2, // [0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3]
    2: 4, // [3,2,1,0,3,2,1,0,3,2,1,0,3,2,1,0]
    3: 10, // [0,2,4,6,0,2,4,6,0,2,4,6,0,2,4,6]
    4: 12, // [6,4,2,0,6,4,2,0,6,4,2,0,6,4,2,0]
    5: 26, // [0,3,6,9,0,3,6,9,0,3,6,9,0,3,6,9]
    6: 28, // [9,6,3,0,9,6,3,0,9,6,3,0,9,6,3,0]
    7: 51, // [0,4,8,12,1,5,9,13,2,6,10,14,3,7,11,15]
    8: 54, // [12,8,4,0,13,9,5,1,14,10,6,2,15,11,7,3]
    9: 83, // [0,5,10,15,1,6,11,16,2,7,12,17,3,8,13,18]
    10: 87, // [18,13,8,3,17,12,7,2,16,11,6,1,15,10,5,0]
]

cleanCode = 0x00;
failCodes = [
    0x00: "000 - clean",

    // 0x0X - observational
    0x01: "001 - phase centre",
    0x02: "002 - channel selection",
    0x03: "003 - pointing",
    0x04: "004 - sun elevation",
    0x05: "005 - capture mode",

    // 0x1X - runtime
    0x10: "016 - bad tiles",
    0x11: "017 - dead dipoles",

    // 0x2X - quality
    0x20: "032 - data quality",
    0x21: "033 - large iono qa",
    0x22: "034 - null iono qa",
    0x23: "035 - no files",

    // 0x3X - prep
    0x31: "049 - total occupancy",
    0x32: "050 - rfi occupancy",
    0x33: "051 - ssins occupancy",
    0x34: "052 - ssins streak",
    0x35: "053 - ssins narrow",
    0x36: "054 - ssins dab",
    0x37: "055 - prep status",
    0x3F: "063 - prepqa bad tiles",

    // 0x4X - calibration
    0x40: "064 - high rms convg",
    0x4F: "079 - calqa bad tiles",
    // 0x5X - visQA
    // 0x6X - ps_metrics
    0x60: "096 - large unsub p_win",
    0x61: "097 - large sub p_win",
    0x62: "098 - small unsub p_win:p_wg",
    // 0x63: "099 - small sub p_win:p_wg",
    0x64: "100 - large sub_p_win:p_wg",
    // 0x65
    0x66: "102 - small sub_p_win:p_win",
    0x67: "103 - large sub_p_win:p_win",

    // 0x7X - imgQA
    0x70: "112 - large unsub v_rms_box",
    0x71: "113 - large sub v_rms_box",
    0x72: "114 - large unsub pks_int v/(xx+yy)",
    0x73: "115 - large sub pks_int v/(xx+yy)",
    0x74: "116 - large unsub pks_int |xx-yy|",
    0x75: "117 - large sub pks_int |xx-yy|",
    0x76: "118 - large sub_xx_pks_int:xx_pks_int",
    0x77: "119 - large sub_yy_pks_int:yy_pks_int",
    0x78: "120 - large sub xx_pks_int",
    0x79: "121 - large sub yy_pks_int",
]

// pass in a list of codes and (optianally) other associated items,
// return the first non-clean code with its associated items
//
// example:
// codes = [0x00, 0x31]
// reasons = ["clean", "total_occ"]
// firstFail([codes, reasons]) => [0x31, "total_occ"]
def firstFail(pt) {
    if (pt == null || pt.size() == 0) {
        return [cleanCode, null]
    }
    def tp = pt.transpose()
    if (tp == null || tp.size() == 0) {
        return [cleanCode, null]
    }
    try {
        def failures = tp.findAll { it -> it[0] != failCodes[cleanCode] }
        if (failures == null || failures.size() == 0) {
            return tp[0]
        } else {
            return failures[0]
        }
    } catch (Exception e) {
        print("pt=${pt}")
        StackTraceUtils.sanitize(e).printStackTrace()
        throw e
    }
}

def wrap_angle(a) {
    return (Float.valueOf(a) + 180) % 360 - 180
}

def wsSummarize(obsid, wsJson, filesJson, tapJson, quality_update, manualAnts) {
    def wsStats = parseJson(wsJson)
    def fileStats = parseJson(filesJson)
    def tapStats = parseJson(tapJson)

    def obs_name = wsStats.obsname;
    def groupid = wsStats.groupid;
    def projectid = wsStats.projectid;

    def ra_phase_center = wsStats.ra_phase_center;
    if (ra_phase_center == null) {
        ra_phase_center = (wsStats.metadata?:[:]).ra_pointing?:Float.NaN
    }
    ra_phase_center = wrap_angle(ra_phase_center)

    def dec_phase_center = wsStats.dec_phase_center;
    if (dec_phase_center == null) {
        dec_phase_center = (wsStats.metadata?:[:]).dec_pointing?:Float.NaN
    }
    dec_phase_center = wrap_angle(dec_phase_center)

    def az_pointing = (wsStats.metadata?:[:]).azimuth_pointing?:Float.NaN
    az_pointing = wrap_angle(az_pointing)
    def el_pointing = (wsStats.metadata?:[:]).elevation_pointing
    def ew_pointing = null
    // pointings: east-west vs EOR1 vs sweet:
    // | EW | E1 | SW |  az |    el |        delays |    lsts  |
    // |----+----+----+-----+-------+-------------- | -------- |
    // | -3 |  5 | 26 |  90 | 69.16 | {0,3,6,9,...} | -28..-20 |
    // | -2 |  3 | 10 |  90 | 76.28 | {0,2,4,6,...} | -19..-12 |
    // | -1 |  1 |  2 |  90 | 83.19 | {0,1,2,3,...} | -12..-4  |
    // |  0 |  0 |  0 |   0 | 90.00 | {0,0,0,0,...} |  -4..+4  |
    // | +1 |  2 |  4 | 270 | 83.19 | {3,2,1,0,...} |  +4..+12 |
    // | +2 |  4 | 12 | 270 | 76.28 | {6,4,2,0,...} | +11..+19 |
    // | +3 |  6 | 28 | 270 | 69.16 | {9,6,3,0,...} | +19..+28 |
    if ([-90,0,90].contains(az_pointing.round() as int)) {
        ew_pointing = ((90-el_pointing) / 7).round() as int
        if (az_pointing > 0) ew_pointing *= -1
    } else {
        println "unknown azel pointing az=${az_pointing} el=${el_pointing}"
    }
    def gridpoint_number = (wsStats.metadata?:[:]).gridpoint_number
    def gridpoint_name = (wsStats.metadata?:[:]).gridpoint_name
    def sweet_pointing = null
    if (gridpoint_name == "sweet") {
        sweet_pointing = gridpoint_number
    } else if (gridpoint_name == "EOR1") {
        if (eor12sweet[gridpoint_number] != null) {
            sweet_pointing = eor12sweet[gridpoint_number]
        } else {
            println "unknown EOR1 gridpoint_number ${gridpoint_number}"
        }
    } else {
        println "unknown gridpoint_name ${gridpoint_name}"
    }
    def lst = wrap_angle((wsStats.metadata?:[:]).local_sidereal_time_deg.floatValue())

    def eorfield = null;
    def eorband = null;
    // eor fields
    // | field | ra h | ra d | dec |
    // | ----- | ---- | ---- | --- |
    // | EoR0  | 0    |    0 | -27 |
    // | EoR1  | 4    |   60 | -30 |
    // | EoR2  | 10.3 |  155 | -10 |
    // | EoR3  | 1    |   15 | -27 |
    def nearest_ra = ra_phase_center.round().intValue()
    def nearest_dec = dec_phase_center.round().intValue()
    if (nearest_ra == 0 && nearest_dec == -27) {
        eorfield = 0
    } else if (nearest_ra == 60 && nearest_dec == -30) {
        eorfield = 1
    } else if (nearest_ra == 155 && nearest_dec == -10) {
        eorfield = 2
    } else if (nearest_ra == 15 && nearest_dec == -27) {
        eorfield = 3
    } else {
        println "unknown eor field for ${nearest_ra} ${nearest_dec}"
    }

    // eor bands
    // |   band | cent          | range                |
    // | ------ | ------------- | -------------------- |
    // | 0 low  |  120 (154MHz) | 109-132 (139-170MHz) |
    // | 1 high |  142 (182MHz) | 131-154 (167-198MHz) |
    // | 2 ulow |   70  (90MHz) |  59-88   (75-113MHz) |

    def coarse_chans = ((wsStats.rfstreams?:[:])["0"]?:[:]).frequencies?:[]
    def center_chan = null
    if (tapStats.center_channel_number != null) {
        center_chan = tapStats.center_channel_number as int
    } else if(coarse_chans != null && coarse_chans.size() > 0) {
        center_chan = coarse_chans[Math.max(coarse_chans.size()/2 as int - 1, 0)]
    }

    if (center_chan == 142) {
        eorband = 1
    } else if (center_chan == 120) {
        eorband = 0
    } else if (center_chan == 70) {
        eorband = 2
        print("unknown eor band for ${center_chan}")
    }

    def nscans = ((wsStats.stoptime?:0) - (wsStats.starttime?:0)) / (wsStats.int_time?:1)
    def delays = (wsStats.alldelays?:[:]).values().flatten()
    def quality = wsStats.quality?:[:]
    def tiles = wsStats.tdict?:[:]
    def tile_nums = tiles.collect { k, _ -> k as int }
    def tile_names = tiles.collect { _, v -> v[0] }
    def tile_rxs = tiles.collect { _, v -> v[1] }
    def n_tiles = tile_names.size()
    def n_lb = tile_names.count { it =~ /(?i)lb/ }
    def n_hex = tile_names.count { it =~ /(?i)hex/ }

    def bad_tiles = wsStats.bad_tiles?:[:]
    def n_bad_tiles = bad_tiles.size()
    def n_good_tiles = n_tiles - n_bad_tiles
    def bad_tile_frac = n_bad_tiles / n_tiles
    def n_dead_dipoles = delays.count { it == 32 }
    def dead_dipole_frac = null
    if (delays.size() != null && delays.size() > 0) {
        dead_dipole_frac = n_dead_dipoles / delays.size()
    }
    def dataquality = Float.valueOf(wsStats.dataquality?:0)
    def dataqualitycomment = wsStats.dataqualitycomment?:''
    def manual_dataquality = Float.valueOf(quality_update.dataquality?:0)
    if (quality_update.dataquality != null && dataquality != manual_dataquality) {
        dataquality = manual_dataquality
        dataqualitycomment = "manual: ${quality_update.dataqualitycomment?:''}"
    }
    def faults = wsStats.faults?:[:]
    def badstates = (faults.badstates?:[:]).values().flatten()
    def badpointings = (faults.badpointings?:[:]).values().flatten()
    def badfreqs = (faults.badfreqs?:[:]).values().flatten()
    def badgains = (faults.badgains?:[:]).values().flatten()
    def badbeamshape = (faults.badbeamshape?:[:]).values().flatten()
    def fail_reasons = []
    def capture_mode = wsStats.mode

    def config = tapStats.mwa_array_configuration
    def sun_elevation = Float.valueOf(tapStats.sun_elevation?:'NaN')
    def sun_pointing_distance = Float.valueOf(tapStats.sun_pointing_distance?:'NaN')

    def bad_ants = bad_tiles.collect { tile_nums.indexOf(it) + 1 }

    // fail codes
    def fail_code = 0x00 // no error

    // 0x0X - observational
    if (params.filter_eorfield != null && eorfield != params.filter_eorfield) {
        fail_reasons += [sprintf("phase_radec(%+.1f,%+.1f)!=eor%d", ra_phase_center, dec_phase_center, params.filter_eorfield)]
        fail_code = fail_code==0x00 ? 0x01 : fail_code
    }
    if (params.filter_ra != null && (ra_phase_center - params.filter_ra).abs() > 0.1) {
        fail_reasons += [sprintf("phase_ra(%+.1f)!=ra(%+.1f)", ra_phase_center, params.filter_ra)]
        fail_code = fail_code==0x00 ? 0x01 : fail_code
    }
    if (params.filter_dec != null && (dec_phase_center - params.filter_dec).abs() > 0.1) {
        fail_reasons += [sprintf("phase_dec(%+.1f)!=dec(%+.1f)", dec_phase_center, params.filter_dec)]
        fail_code = fail_code==0x00 ? 0x01 : fail_code
    }
    if (params.filter_eorband != null && eorband != params.filter_eorband) {
        fail_reasons += ["center_chan=${center_chan}"]
        fail_code = fail_code==0x00 ? 0x02 : fail_code
    }
    if (params.filter_sweet_pointings != null && !params.filter_sweet_pointings.contains(sweet_pointing)) {
        fail_reasons += ["sweet_pointing=${sweet_pointing}"]
        fail_code = fail_code==0x00 ? 0x03 : fail_code
    }
    if (params.filter_ew_pointings != null && !params.filter_ew_pointings.contains(ew_pointing)) {
        fail_reasons += ["ew_pointing=${ew_pointing}"]
        fail_code = fail_code==0x00 ? 0x03 : fail_code
    }
    if (params.filter_sun_elevation != null && sun_elevation != null && sun_elevation > params.filter_sun_elevation) {
        fail_reasons += [sprintf("sun_elevation(%+.1f)>%+.1f", Float.valueOf(sun_elevation), Float.valueOf(params.filter_sun_elevation))]
        fail_code = fail_code==0x00 ? 0x04 : fail_code
    }
    if (params.filter_min_sun_pointing_distance != null && sun_pointing_distance != null && sun_pointing_distance > params.filter_min_sun_pointing_distance) {
        fail_reasons += [sprintf("sun_pointing_distance(%+.1f)>%+.1f", Float.valueOf(sun_pointing_distance), Float.valueOf(params.filter_min_sun_pointing_distance))]
        fail_code = fail_code==0x00 ? 0x04 : fail_code
    }
    if (capture_mode == "NO_CAPTURE") {
        fail_reasons += ["no cap fr fr"]
        fail_code = fail_code==0x00 ? 0x05 : fail_code
    }

    // 0x1X - runtime
    if (params.filter_bad_tile_frac != null && bad_tile_frac > params.filter_bad_tile_frac) {
        fail_reasons += ["bad_tiles(${bad_tiles.size()})=${displayInts(bad_tiles)}"]
        fail_code = fail_code==0x00 ? 0x10 : fail_code
    }
    if (params.filter_dead_dipole_frac != null && dead_dipole_frac > params.filter_dead_dipole_frac) {
        fail_reasons += ["dead_dipole_frac=${dead_dipole_frac}"]
        fail_code = fail_code==0x00 ? 0x11 : fail_code
    }

    // 0x2X - quality
    if (params.filter_quality != null && dataquality > params.filter_quality) {
        fail_reasons += ["dataquality=${dataquality} (${dataqualitycomment})"]
        fail_code = fail_code==0x00 ? 0x20 : fail_code
    }
    if (params.filter_ionoqa && quality.iono_qa != null && quality.iono_qa > params.filter_ionoqa) {
        fail_reasons += [sprintf("rts_iono_qa(%.1f)>%.1f", Float.valueOf(quality.iono_qa), Float.valueOf(params.filter_ionoqa))]
        fail_code = fail_code==0x00 ? 0x21 : fail_code
    }
    if (params.filter_ionoqa && quality.iono_qa == null) {
        fail_reasons += ["rts_iono_qa is null"]
        fail_code = fail_code==0x00 ? 0x22 : fail_code
    }
    if (fileStats.num_data_files < 2) {
        fail_reasons += ["no data files"]
        fail_code = fail_code==0x00 ? 0x23 : fail_code
    }

    fail_code = failCodes[fail_code]

    [
        fail_code: fail_code,
        fail_reasons: fail_reasons,
        // obs metadata
        obs_name: obs_name,
        groupid: groupid,
        corrmode: wsStats.mode,
        delaymode: wsStats.delaymode_name,
        starttime_mjd: parseFloatOrNaN(tapStats.starttime_mjd),
        starttime_utc: tapStats.starttime_utc,
        sun_elevation: sun_elevation,
        sun_pointing_distance: sun_pointing_distance,

        // pointing
        ra_pointing: wsStats.metadata.ra_pointing,
        dec_pointing: wsStats.metadata.dec_pointing,
        az_pointing: az_pointing,
        el_pointing: el_pointing,
        ra_phase_center: ra_phase_center,
        dec_phase_center: dec_phase_center,
        ew_pointing: ew_pointing,
        sweet_pointing: sweet_pointing,
        lst: lst,
        eorfield: eorfield,

        // channels
        freq_res: wsStats.freq_res,
        coarse_chans: coarse_chans,
        eorband: eorband,
        centre_freq: ((coarse_chans[0] + coarse_chans[-1]) * 1.28e6 / 2),

        // times
        int_time: wsStats.int_time,
        nscans: nscans,

        // tiles
        config: config,
        n_tiles: n_tiles,
        n_good_tiles: n_good_tiles,
        tile_nums: tile_nums,
        tile_rxs: tile_rxs,
        bad_tiles: bad_tiles,
        bad_ants: bad_ants,
        manual_ants: (manualAnts as ArrayList),
        // iono quality
        iono_magnitude: quality.iono_magnitude,
        iono_pca: quality.iono_pca,
        iono_qa: quality.iono_qa,
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
}

// use mwa webservices to gate obsids by pointing and faults
workflow ws {
    take:
        // channel of obs ids
        obsids

    main:
        obsids | wsMeta & tapMeta

        def quality_updates = file(params.quality_updates_path)
            .readLines()
            .findAll { !it.startsWith('#') && it.length() > 13 }
            .collectEntries { line ->
                def (obsid, quality, comment) = line.split(',')
                [obsid, [dataquality: quality, dataqualitycomment: comment]]
            }

        def tile_updates = file(params.tile_updates_path)
            .readLines()
            .collect { line ->
                def (firstObsid, lastObsid, tileIdxs) = line.split(',') + ['', '']
                [firstObsid as int, lastObsid as int, (tileIdxs as String).split("\\|").collect {it as Integer} ]
            }

        // print(params)

        def wsSummary = wsMeta.out.join(tapMeta.out).map { obsid, wsJson, filesJson, tapJson ->
                try {
                    def quality_update = quality_updates[obsid]?:[:]
                    def manualAnts = ([]) as Set
                    tile_updates.each {
                        def (firstObsid, lastObsid, tileIdxs, comment) = it
                        if (obsid as int >= firstObsid && obsid as int <= lastObsid) {
                            manualAnts.addAll(tileIdxs)
                        }
                    }
                    def summary = wsSummarize(obsid, wsJson, filesJson, tapJson, quality_update, manualAnts)
                    [ obsid, deepcopy(summary) ]
                } catch (Exception e) {
                    println "error summarizing ${obsid}"
                    StackTraceUtils.sanitize(e).printStackTrace()
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
            .map { obsids, fail_code ->
                [
                    fail_code,
                    obsids.size(),
                ].join("\t")
            }
            .collectFile(
                name: "fail_counts_ws.tsv", newLine: true, sort: true,
                seed: ([ "FAIL CODE", "COUNT" ]).join("\t"),
                storeDir: "${results_dir}"
            )
            | view { it.readLines().size() }

        wsSummary.map { obsid, summary ->
                [obsid, summary.lst.round().intValue()]
            }
            .groupTuple(by: 1)
            .map { obsids, lst ->
                [
                    sprintf("%+03d", lst),
                    obsids.size(),
                ].join("\t")
            }
            .collectFile(
                name: "lst_counts_ws.tsv", newLine: true, sort: true,
                seed: ([ "LST", "COUNT" ]).join("\t"),
                storeDir: "${results_dir}"
            )
            | view { it.readLines().size() }

        // display wsSummary
        wsStats = wsSummary.map { obsid, summary ->
                [
                    obsid,
                    (summary.fail_code==null||summary.fail_code==failCodes[0x00])?'':summary.fail_code,
                    summary.fail_reasons.join('|'),
                    isNaN(summary.starttime_mjd)?'':sprintf("%.5f", summary.starttime_mjd),
                    summary.starttime_utc?:'',

                    // pointing
                    sprintf("% 5.2f", summary.ra_pointing),
                    sprintf("% 5.2f", summary.dec_pointing),
                    sprintf("% 5.2f", summary.az_pointing),
                    sprintf("% 5.2f", summary.el_pointing),
                    isNaN(summary.ra_phase_center)?'':sprintf("% 5.2f", summary.ra_phase_center),
                    isNaN(summary.dec_phase_center)?'':sprintf("% 5.2f", summary.dec_phase_center),
                    isNaN(summary.ew_pointing)?'':sprintf("% 2d", summary.ew_pointing),
                    isNaN(summary.sweet_pointing)?'':sprintf("% 2d", summary.sweet_pointing),
                    sprintf("% 5.2f", summary.lst),
                    summary.obs_name,
                    isNaN(summary.eorfield)?'':summary.eorfield,
                    isNaN(summary.sun_elevation)?'':summary.sun_elevation,
                    isNaN(summary.sun_pointing_distance)?'':summary.sun_pointing_distance,

                    // channels
                    summary.freq_res,
                    displayInts(summary.coarse_chans),
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
                    displayInts(summary.tile_nums, delim=' '),
                    displayInts(summary.bad_tiles, delim=' '),
                    displayInts(summary.tile_rxs, delim=' '),

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
                storeDir: "${results_dir}",
            )
            // display output path and number of lines
            | view { [it, it.readLines().size()] }

        wsStats.flatMap { file ->
                [
                    ["FAIL CODE", "tab10"],
                    ["IONO QA", "viridis"],
                    ["N BAD TILES", "viridis"],
                ].collect { col, palette ->
                    short_name = col.replaceAll(/\s+/, '_').toLowerCase()
                    meta = [
                        name: "ws_${short_name}", title: col, palette: palette,
                        x: 'OBS', y: 'LST DEG', c: col,
                    ]
                    [ meta, file ]
                }
            }
            | tsvScatterPlot

        wsSummary
            .filter { _, summary -> summary.fail_reasons == [] }
            .map { obsid, _ -> obsid }
            .tap { pass }
            | (wsMetafits & wsSkyMap & wsPPDs)

        wsMetafits.out | metaJson

    emit:
        // channel of good obsids with their metafits: tuple(obsid, metafits)
        obsMetafits = wsMetafits.out

        // channel of video name and frames to convert
        frame = pass.join(wsSkyMap.out)
            .map { _, png -> ["skymap", png] }
            .mix( pass.join(wsPPDs.out.map { _, png -> ["ppd", png] }) )
            .groupTuple()

        // channel of (obsid, metadata hashmap)

        obsMeta = wsSummary.map { obsid, summary ->
            def meta = [:]
            [
                // "groupid", "starttime_utc", "starttime_mjd", "obs_name"
                "ew_pointing", "centre_freq",
                "n_tiles", "bad_ants", "manual_ants", "tile_nums",
                "eorband", "eorfield", "lst", "int_time", "freq_res",
                "ra_phase_center", "dec_phase_center",
                "coarse_chans"
            ].each { key ->
                if (summary[key] != null) {
                    meta[key] = summary[key]
                }
            }
            [obsid, deepcopy(meta)]
        }

        fail_codes = fail_codes
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

        // subobsMetaMetafitsPrep = obsMetaUV.cross(uvMeta.out) {[it[0], it[1].name]}
        //     .map { obsMetaUV_, uvMeta_ ->
        //         def (obsid, meta, vis) = obsMetaUV_; (_, __, json) = uvMeta_;
        //         jsonMeta = parseJson(json)
        //         newMeta = [
        //             lowfreq:jsonMeta.freqs[0],
        //             first_lst: jsonMeta.times[0]['lst_rad'],
        //             first_jd: jsonMeta.times[0]['jd1'],
        //             nchans:jsonMeta.freqs.size(),
        //             ntimes:jsonMeta.times.size(),
        //         ];
        //         ['eorband', 'eorfield', 'config', 'total_weight'].each { key ->
        //             if (jsonMeta[key] != null) {
        //                 newMeta[key] = jsonMeta[key]
        //             }
        //         }
        //         [obsid, mapMerge(meta, newMeta), vis] }

        obsMetaMetafits.join(asvoPrep.out).flatMap { obsid, _, metafits, meta, uvfits_ ->
                def uvfits = coerceList(uvfits_)
                if (uvfits.size > 1) {
                    uvfits.collect { f ->
                        [obsid, mapMerge(meta, [subobs: f.baseName.split('_')[2]]), metafits, f]
                    }
                } else {
                    [[obsid, meta, metafits, uvfits_]]
                }
            }
            .tap { subobsMetaMetafitsPrep }

        subobsMetaMetafitsPrep.map { obsid, meta, _, uvfits -> [ obsid, meta, uvfits] }.tap { subobsMetaVis }
        obsMetaMetafits.map { obsid, meta, metafits -> [obsid, metafits] }.tap { obsMetafits }

        flag(subobsMetaVis, obsMetafits)

        flag.out.subobsFlagmetaPass.map { obsid, meta, flagMeta ->
                [[obsid, meta.subobs?:''], mapMerge(meta, flagMeta)]
            }
            .join(subobsMetaVis.map { obsid, meta, uvfits ->
                [[obsid, meta.subobs?:''], uvfits]
            })
            // .join(flag.out.subobsMetaRxAnts.map { obsid, meta, rxAnts -> [[obsid, meta.subobs?:''], rxAnts] })
            .map { obsSubobs, meta, uvfits ->
                def (obsid, _) = obsSubobs
                [obsid, meta, uvfits]
            }
            .tap { subobsFlagmetaVis }
            // ssins
            .flatMap { obsid, meta, uvfits ->
                def rxAnts = meta.unflaggedRxAnts?:[:]
                def allAnts = rxAnts.collect { rx, ants -> ants }.flatten()
                def subobsAnts = [[meta.subobs, allAnts]]
                if (params.ssins_by_rx) {
                    subobsAnts += (rxAnts.collect {rx, ants ->
                        [(meta.subobs?:'') + sprintf("_rx%02d", rx), ants]
                    })
                }
                def base = uvfits.baseName
                subobsAnts.findAll { _, ants -> ants.size() > 1 }
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
        ssinsOcc = ssins.out.map { def (obsid, meta, _, __, occ_json) = it;
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
                        occ.total?:'',
                        occ.streak?:'',
                        occ.dab_total?:'',
                        occ.narrow_total?:'',
                    ]).join("\t")
                }
                .collectFile(
                    name: "ssins_occupancy.tsv", newLine: true, sort: true,
                    seed: (["OBS", "SUBOBS", "TOTAL", "STREAK", "DAB TOTAL", "NARROW TOTAL"]).join("\t"),
                    storeDir: "${results_dir}"
                )
                | view { [it, it.readLines().size()] }

        allDABs = ssinsOcc.flatMap { obs, meta, occ ->
                occ.findAll { it.key.startsWith("DAB") }.collect { it.key }
            }
            .unique()
            .toSortedList()
            .map{it -> [it]}

        channel.of("OBS", "SUBOBS", "DAB TOTAL").concat(allDABs.flatten())
            .toList()
            .map { it.join("\t") }
            .concat(
                ssinsOcc.filter { obs, meta, occ ->
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
                storeDir: "${results_dir}"
            )
            .view { [it, it.readLines().size()] }

        allNarrows = ssinsOcc.flatMap { obs, meta, occ ->
                occ.findAll { it.key.contains('narrow') }.collect { it.key }
            }
            .unique()
            .toSortedList()
            .map{it -> [it]}

        channel.of("OBS", "SUBOBS", "NARROW TOTAL").concat(allNarrows.flatten())
            .toList()
            .map { it.join("\t") }
            .concat(
                ssinsOcc.filter { obs, meta, occ ->
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
                storeDir: "${results_dir}"
            )
            .view { [it, it.readLines().size()] }

        if (params.noprepqa) {
            subobsMetaVisPass = subobsFlagmetaVis
            subobsMetaVisFlags = channel.empty()
            fail_codes = channel.empty()
        } else {
            // TODO: do we really need vis here?
            subobsMetaVisFlags = (subobsFlagmetaVis.map { obsid, meta, uvfits -> [[obsid, meta.subobs?:''], meta, uvfits ]})
                .join(ssinsOcc.map { obsid, meta, occ -> [[obsid, meta.subobs?:''], occ ] }, remainder: false)
                .map { obsSubobs, meta, uvfits, ssinsStats ->
                    def (obsid, _) = obsSubobs
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

                    flagMeta.fail_code = failCodes[fail_code]
                    if (reason) {
                        flagMeta.reasons = reason
                    }

                    [obsid, meta, uvfits, deepcopy(flagMeta)]
                }

            fail_codes = subobsMetaVisFlags
                // .filter { obsid, meta, uvfits, flagMeta -> (meta.prepFlags?:[]).size() > 0 }
                .map { obsid, meta, uvfits, flagMeta ->
                    [ obsid, flagMeta.fail_code, flagMeta.reasons ]
                }

            fail_codes
                // .filter { _, fail_code, __ -> fail_code != failCodes[0x00] }
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
                    storeDir: "${results_dir}"
                )
                | view { it.readLines().size() }

            fail_codes.groupTuple(by: 1)
                .map { obsids, fail_code, _ ->
                    [ fail_code, obsids.size() ].join("\t")
                }
                .collectFile(
                    name: "fail_counts_ssins.tsv", newLine: true, sort: true,
                    seed: ([ "FAIL CODE", "COUNT" ]).join("\t"),
                    storeDir: "${results_dir}"
                )
                | view { it.readLines().size() }

            subobsMetaVisFlags
                // .filter { obsid, meta, uvfits, flagMeta -> (meta.ssinsFlagsTsv?:[]).size() > 0 }
                .map { obsid, meta, _, flagMeta ->
                    ([
                        obsid,
                        meta.subobs?:'',
                        meta.lst?:'',
                        meta.ew_pointing?:'',
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
                    storeDir: "${results_dir}"
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
                        short_name = col.replaceAll(/\s+/, '_').toLowerCase()
                        meta = [
                            name: "prep_${short_name}", title: col, palette: palette,
                            x: 'OBS', y: 'LST DEG', c: col,
                        ]
                        [ meta, file ]
                    }
                }
                | tsvScatterPlot

            subobsMetaVisPass = subobsMetaVisFlags.filter { obsid, meta, uvfits, flagMeta ->
                    (flagMeta.reasons?:'') == ''
                }
                .map { obsid, meta, uvfits, flagMeta -> [obsid, meta, uvfits]}

            subobsMetaVisPass.map { obsid, meta, vis ->
                    ([
                        obsid,
                        meta.subobs?:'',
                        meta.name?:'',
                        meta.lst?:'',
                        meta.ew_pointing?:'',
                        displayInts(meta.prepFlags?:[], delim=' '),
                        vis
                    ]).join("\t")
                }
                .collectFile(
                    name: "pass_ssins.tsv", newLine: true, sort: true,
                    seed: ([ "OBS", "SUBOBS", "NAME", "LST DEG", "EW POINTING", "NEW ANTS", "VIS" ]).join("\t"),
                    storeDir: "${results_dir}"
                )
        }

        if (params.ssins_apply) {
            subobsMetaVisSSINs = subobsMetaVisPass.join(
                    ssins.out.map { obsid, meta, mask, _, __ ->
                            [obsid, meta, mask]
                        }
                )
                .map { obsid, _, uvfits, meta, mask ->
                    [obsid, meta, uvfits, mask]
                }
                | absolve
            // TODO: flagqa
        } else {
            subobsMetaVisSSINs = subobsMetaVis
        }

    emit:
        // channel of obsids which pass the flag gate: tuple(obsid, meta, metafits, uvfits)
        subobsMetaVisPass = subobsMetaVisSSINs
        subobsReasons = subobsMetaVisFlags.filter { obsid, meta, uvfits, flagMeta ->
                (flagMeta.reasons?:'') != ''
            }
            .map { obsid, meta, uvfits, flagMeta -> [obsid, meta, flagMeta.reasons]}
        // channel of video name and frames to convert
        frame = flag.out.frame
            .mix(ssins.out.flatMap { _, __, ___, imgs, ____ ->
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
            })
            .groupTuple()

        // channel of files to archive, and their buckets
        archive = flag.out.archive
            .mix(ssins.out.map { _, __, mask, ___, ____ -> ["ssins", mask]})

        zip = flag.out.zip

        fail_codes =
            flag.out.fail_codes.join(
                    fail_codes.map { obsid, fail_code, reasons -> [obsid, fail_code] },
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
                def (_, meta, uvfits) = subobsMetaVis_
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
                    occ = stats.channels?[ch]?.rfi_occupancy
                    isNaN(occ)?'':occ
                }
                def rfi_occ = stats.total_rfi_occupancy
                def flagged_sky_chan_idxs = displayInts(flagged_sky_chans, delim=' ')
                def flagged_fchan_idxs = displayInts(stats.flagged_fchan_idxs?:[], delim=' ')
                def flagged_timestep_idxs = displayInts(stats.flagged_timestep_idxs?:[], delim=' ')
                def preflagged_ants = displayInts(stats.preflagged_ants?:[], delim=' ')
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
                storeDir: "${results_dir}"
            )
            // display output path and number of lines
            | view { [it, it.readLines().size()] }

        // collect prepVisQA results as .tsv
        prepVisQA.out
            // form row of tsv from json fields we care about
            .map { obsid, meta, json ->
                stats = parseJson(json);
                bad_ants = stats.BAD_ANTS?:[]
                [
                    obsid,
                    meta.subobs?:'',
                    stats.STATUS?:'',
                    stats.NANTS?:'',
                    stats.NTIMES?:'',
                    stats.NCHAN?:'',
                    stats.NPOLS?:'',
                    bad_ants.size(),
                    displayInts(bad_ants, delim=' '),
                ].join("\t")
            }
            .collectFile(
                name: "prepvis_metrics.tsv", newLine: true, sort: true,
                seed: [
                    "OBS",
                    "SUBOBS",
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
                .map { obsid, meta, json ->
                    stats = parseJson(json);
                    ([ obsid, meta.subobs?:'' ] + getMetric(stats)).join("\t")
                }
                .collectFile(
                    name: "prepvis_${metric}.tsv", newLine: true, sort: true,
                    seed: [
                        "OBS",
                        "SUBOBS",
                        metric,
                    ].join("\t"),
                    storeDir: "${results_dir}"
                )
                // display output path and number of lines
                | view { [it, it.readLines().size()] }
        }

        // plot prepvisQA
        prepVisQA.out | plotPrepVisQA

        if (params.noprepqa) {
            subobsMetaFlags = channel.empty()
            subobsMetaPass = subobsMetaVis.map { obsid, meta, uvfits -> [obsid, meta] }
            fail_codes = channel.empty()

        } else {
            (subobsMetaMetafitsPrep.map { obsid, meta, _, __ -> [[obsid, meta.subobs?:''], meta] })
                .join(flagQA.out.map { obsid, meta, flagJson -> [[obsid, meta.subobs?:''], parseJson(flagJson)]})
                .join(prepVisQA.out.map { obsid, meta, prepJson -> [[obsid, meta.subobs?:''], parseJson(prepJson)]})
                .map { obsidSubobs, meta, flagStats, prepStats ->
                    def (obsid, subobs) = obsidSubobs
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
                        inputs = (flagStats.INPUTS?:[])
                        if (inputs.size() > 0) {
                            unflaggedAnts = inputs.findAll {
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

                    flagMeta.fail_code = failCodes[fail_code]
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
                .map { obsid, _, obsCodes_, reasons ->
                    failCode = coerceList(firstFail([coerceList(obsCodes_)]))[0]
                    [ obsid, failCode ]
                }
                .tap { fail_codes }

            subobsMetaFlags.filter { obsid, meta, flagMeta ->
                    (flagMeta.reasons?:'') == ''
                }
                .map { obsid, meta, flagMeta ->
                    [obsid, meta, flagMeta]
                }
                .tap { subobsFlagmetaPass }
                .map { obsid, meta, flagMeta ->
                    [obsid, meta]
                }
                .tap { subobsMetaPass }

            subobsMetaReasons
                .groupTuple(by: 0)
                .map { obsid, metas, obsCodes_, reasons ->
                    try {
                        def ff = firstFail([
                            coerceList(obsCodes_),
                            coerceList(reasons),
                        ])
                        def (failCode, reason) = ff
                        return [ obsid, failCode, reason?:'' ].join("\t")
                    } catch (Exception e) {
                        println("obsid=${obsid} metas=${metas} obsCodes_=${obsCodes_} reasons=${reasons}")
                        try {println("ff=${ff}")} catch (Exception f) {}
                        StackTraceUtils.sanitize(e).printStackTrace()
                        throw e
                    }
                }
                .collectFile(
                    name: "reasons_flag.tsv", newLine: true, sort: true,
                    seed: ([ "OBSID", "FAIL CODE", "REASON" ]).join("\t"),
                    storeDir: "${results_dir}"
                )
                | view { it.readLines().size() }

            fail_codes.groupTuple(by: 1)
                .map { obsids, fail_code ->
                    [ fail_code, obsids.size() ].join("\t")
                }
                .collectFile(
                    name: "fail_counts_flag.tsv", newLine: true, sort: true,
                    seed: ([ "FAIL CODE", "COUNT" ]).join("\t"),
                    storeDir: "${results_dir}"
                )
                | view { it.readLines().size() }

            subobsMetaFlags
                .map { obsid, meta, flagMeta ->
                    ([
                        obsid,
                        meta.subobs?:'',
                        meta.lst?:'',
                        meta.ew_pointing?:'',
                        flagMeta.fail_code?:'',
                        flagMeta.reasons?:'',
                        flagMeta.total_occ==null?'':flagMeta.total_occ,
                        flagMeta.total_non_preflagged_bl_occ==null?'':flagMeta.total_non_preflagged_bl_occ,
                        displayInts(flagMeta.flagAnts?:[], delim=' '),
                        displayInts(flagMeta.prepAnts?:[], delim=' '),
                        displayInts(flagMeta.manualAnts?:[], delim=' '),
                        displayInts(meta.prepFlags?:[], delim=' '),
                    ]).join("\t")
                }
                .collectFile(
                    name: "prep_flags.tsv", newLine: true, sort: true,
                    seed: ([
                        "OBS", "SUBOBS", "LST DEG", "EW POINTING", "FAIL CODE", "REASONS",
                        "PREP OCC", "RFI OCC", "ORIGINAL ANTS", "PREP ANTS",
                        "MANUAL ANTS", "NEW ANTS"
                    ]).join("\t"),
                    storeDir: "${results_dir}"
                )
                // display output path and number of lines
                | view { [it, it.readLines().size()] }

            subobsMetaPass.map { obsid, meta ->
                    ([
                        obsid,
                        meta.subobs?:'',
                        displayInts(meta.prepFlags?:[], delim=' '),
                    ]).join("\t")
                }
                .collectFile(
                    name: "pass_flag.tsv", newLine: true, sort: true,
                    seed: ([ "OBS", "SUBOBS", "NEW ANTS" ]).join("\t"),
                    storeDir: "${results_dir}"
                )
        }

        // TODO: auto plots
        autoplotByRx = subobsFlagmetaPass.map {obsid, meta, flagMeta ->
                [[obsid, meta.subobs?:''], mapMerge(meta, flagMeta)]
            }.join(subobsMetaMetafitsPrep.map{obsid, meta, metafits, vis ->
                [[obsid, meta.subobs?:''], metafits, vis]
            }).flatMap { obsSubobs, meta, metafits, uvfits ->
                def (obsid, _) = obsSubobs
                (meta.unflaggedRxAnts?:[:]).collect { rx, ants ->
                    def suffix = sprintf("_rx%02d", rx);
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
        //         def (obsid, _) = obsSubobs
        //         flavorAnts.collect { flavor, ants ->
        //             suffix = sprintf("_flavor-%s", flavor);
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
        frame = plotPrepVisQA.out.flatMap { _, __, imgs ->
                imgs.collect { img ->
                    def suffix = img.baseName.split('_')[-1]
                    [''+"prepvisqa_${suffix}", img]
                }
            }
            .mix(autoplot.out.map {_, meta, img -> [''+"prepvisqa_autoplot${meta.suffix?:''}", img]})
            .groupTuple()
        archive = channel.empty() // TODO: archive flag jsons
        zip = prepVisQA.out.map { _, __, json -> ["prepvisqa", json]}
            .mix(flagQA.out.map { _, __, json -> ["flagqa", json]})
            .groupTuple()
        fail_codes
}

workflow cal {
    take:
        // channel of metafits and preprocessed uvfits: tuple(obsid, meta, metafits, uvfits)
        obsMetaMetafitsVis
        // TODO: obsMetafits and subobsMetaVis instead
    main:
        // channel of metadata for each obsid: tuple(obsid, meta)
        obsMeta = obsMetaMetafitsVis.map { obsid, meta, _, __ -> [obsid, meta] }
        // channel of metafits for each obsid: tuple(obsid, metafits)
        obsMetafits = obsMetaMetafitsVis.map { obsid, _, metafits, __ -> [obsid, metafits] }
        // hyperdrive di-calibrate on each obs
        obsMetaMetafitsVis
            .map { def (obsids, meta, metafits, uvfits) = it
                [obsids, meta, metafits, uvfits, params.dical_args]
            }
            | hypCalSol

        // - hypCalSol gives multiple solutions, flatMap gives 1 tuple per solution.
        hypCalSol.out.flatMap { obsid, meta, solns, _ ->
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
        //                     startTime = logDateFmt.parse(startMatch[0][1])
        //                 }
        //                 convergedMatch = (diCalLog.getText() =~ /([\d: -]+) INFO  All timesteps: (\d+)\/(\d+)/)
        //                 if (startTime && convergedMatch) {
        //                     convergedTime = logDateFmt.parse(convergedMatch[0][1])
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

        // phase fit on calibration solutions
        obsMetafits.cross(
                eachCal.map { obsid, meta, soln -> [[obsid, meta.name], meta, soln] }
                    .groupTuple(by: 0)
                    .map{ obsName, metas, solns ->
                        def (obsid, _) = obsName
                        [obsid, coerceList(metas)[0], solns]
                    }
            )
            .map { obsMetafits_, hypCalSol_ ->
                def (obsid, metafits) = obsMetafits_
                def (_, meta, solns) = hypCalSol_
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
                def results = stats.RESULTS?[0]?:[]
                def (nans, convergences) = results.split { it == "NaN" }
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
                storeDir: "${results_dir}${params.cal_suffix}"
            )
            // display output path and number of lines
            | view { [it, it.readLines().size()] }

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
        calQA.out | plotCalQA
        calQA.out.map { _, meta, json ->
                name = meta.name
                [name, json] }
            .groupTuple(by: 0)
            .map { name, jsons ->
                [name, jsons.sort(false)]
            }
            | plotCalJsons

        // collect calQA results as .tsv
        calQA.out
            // form row of tsv from json fields we care about
            .map { obsid, meta, json ->
                stats = parseJson(json)
                bad_ants = stats.BAD_ANTS?:[]
                // convg_var = stats.CONVERGENCE_VAR
                [
                    obsid,
                    meta.subobs?:'',
                    meta.name,
                    stats.STATUS?:'',
                    bad_ants.size,
                    displayInts(bad_ants, delim=' '),
                    (stats.PERCENT_UNUSED_BLS?:0) / 100,
                    (stats.PERCENT_NONCONVERGED_CHS?:0) / 100,
                    stats.RMS_CONVERGENCE?:'',
                    stats.SKEWNESS?:'',
                    stats.RECEIVER_VAR?:'',
                    stats.DFFT_POWER?:'',
                    stats.FAILURE_REASON?:'',
                ].join("\t")
            }
            .collectFile(
                name: "cal_metrics${params.cal_suffix}.tsv", newLine: true, sort: true,
                seed: [
                    "OBS", "SUBOBS", "CAL NAME", "STATUS",
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
                .map { obsid, meta, json -> [obsid, meta, parseJson(json)] }
                .filter { _, __, stats -> stats != null }
                .map { obsid, meta, stats ->
                    ([ obsid, meta.subobs?:'', meta.name ] + getMetric(stats)).join("\t")
                }
                .collectFile(
                    name: "calqa_${metric}.tsv", newLine: true, sort: true,
                    seed: [
                        "OBS", "SUBOBS", "CAL NAME", metric,
                    ].join("\t"),
                    storeDir: "${results_dir}${params.cal_suffix}"
                )
                // display output path and number of lines
                | view { [it, it.readLines().size()] }
        }

        allTiles = obsMetaMetafitsVis.flatMap { obs, meta, metafits, vis ->
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
                    phaseFits.out.combine(allTiles).map { obsid, meta, tsv, _, tiles ->
                            phaseFits = parseCsv(coerceList(tsv)[0], true, 0, '\t')
                            ([obsid, meta.lst] + tiles.collect { tile ->
                                tilePhaseFits = phaseFits.find { "${it['tile_id']}" == "${tile}" }?:[:]
                                tilePhaseFits[key]?:''
                            }).join("\t")
                        }
                        .toSortedList()
                        .flatten()
                ).collectFile(
                    name: "phase_fits_${key}.tsv", newLine: true, sort: false,
                    storeDir: "${results_dir}${params.cal_suffix}"
                )
                .view { [it, it.readLines().size()] }
        }

        // channel of obsids and names that pass qa. tuple(obsid, name)
        // - take tuple(obsid, cal_name, json) from calQA.out
        // - filter on json.STATUS == "PASS"
        // - take obsid and name
        obsMetaPass = calQA.out
            .map { obsid, meta, json -> [obsid, meta, parseJson(json)] }
            .filter { _, __, stats -> stats != null }
            .map { obsid, meta, stats ->
                def (fail_code, reason) = calqa_pass(stats)
                def newMeta = [fail_code: fail_code, reason: reason];
                if (stats.BAD_ANTS) {
                    def prepFlags = (meta.prepFlags?:[]) as Set
                    def calFlags = ([]) as Set
                    if (!params.noCalFlags) {
                        calFlags.addAll(stats.BAD_ANTS?:[])
                    }
                    def newflags = (calFlags - prepFlags) as ArrayList
                    if (newflags) {
                        newMeta.calFlags = deepcopy(newflags.sort(false))
                    }
                }
                [obsid, mapMerge(meta, newMeta), stats]
            }
            .tap { obsMetaStats }
            // TODO: reintroduce status filter
            // .filter { _, meta, stats ->
            //     stats.STATUS == null || stats.STATUS == "PASS"
            // }
            .filter { _, meta, stats ->
                meta.fail_code == 0x00
            }
            .map { obsid, meta, stats -> [obsid, meta] }

        fail_codes = obsMetaStats.map{ obs, meta, stats ->
                [obs, failCodes[meta.fail_code?:0x00], meta.reason?:""]
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
                    name: "reasons_calqa.tsv", newLine: true, sort: true,
                    seed: ([ "OBSID", "FAIL CODE", "REASON" ]).join("\t"),
                    storeDir: "${results_dir}"
                )
                | view { it.readLines().size() }

    emit:
        fail_codes = fail_codes.filter { obsid, fail_code, _ -> fail_code != failCodes[0x00] }
        // channel of calibration solutions that pass qa. tuple(obsid, name, cal)
        // - take tuple(obsid, meta, soln) from allCal
        // - match with obsMetaPass on (obsid, cal_name)
        obsMetaCalPass = allCal
            .cross(obsMetaPass) {def (obsid, meta) = it; [obsid, meta.name]}
            .map{ allCal_, obsMetaPass_ ->
                def (obsid, _, soln) = allCal_
                def (__, meta) = obsMetaPass_
                [obsid, meta, soln]
            }
        // channel of files to archive, and their buckets
        archive =
            hypCalSol.out
                .flatMap { obsid, meta, solns, _ ->
                    coerceList(solns).collect { soln ->
                        [ "soln", soln ]
                    }
                }
                .mix(
                    calQA.out.map { _, __, json -> ["calqa", json]}
                )

        // channel of video name and frames to convert
        frame = plotCalQA.out.mix(plotSols.out)
            .flatMap { _, meta, pngs ->
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
            .mix(phaseFits.out.flatMap { _, meta, __, pngs ->
                pngs.collect { png ->
                    suffix = png.baseName.split(' ')[-1]
                    ["phasefits_${suffix}", png]
                }
            })
            .groupTuple()
        // channel of files to zip
        zip = calQA.out.map { _, meta, json -> ["calqa_${meta.name}", json] }
            .mix(solJson.out.map { obsid, meta, json -> ["soljson_${meta.name}", json] })
            .mix(phaseFits.out.map { obsid, meta, tsv, _ -> ["phasefits_${meta.name}", tsv] })
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
                    displayInts(db_ants, delim=' '),
                    displayInts(manual_ants, delim=' '),
                    displayInts(prep_ants, delim=' '),
                    displayInts(cal_ants, delim=' '),
                    stats.NPOOR_ANTS?:'',
                    displayInts(vis_ants, delim=' '),
                    stats.NPOOR_BLS?:'',
                    displayInts(new_ants as ArrayList, delim=' '),
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

        visQA.out | plotVisQA

        if (params.noeor) {
            obsMetaUVEoR = obsMetaUV
        } else {
            obsMetaUVEoR = obsMetaUV
                .filter { _, meta, __ -> (meta.eorband != null && meta.eorfield != null) }
        }

        // ps_metrics
        if (params.nopsmetrics) {
            obsMetaUVPass_ = obsMetaUVEoR
            fail_codes = channel.empty()
        } else {
            obsMetaUVEoR | psMetrics

            // collect psMetrics as a .dat
            psMetrics.out
                // read the content of each ps_metrics file including the trailing newline
                .map { obsid, vi_name, dat -> dat.getText() }
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
                    [obsid, mapMerge(meta, newMeta)]
                }

            passReasons = psMeta
                // get cmt reduced metrics failures
                .groupTuple(by: 0)
                .flatMap { obsid, metas ->
                    def nosubMeta = metas.find { it.sub == null} ?: [:]
                    def (nosub_fail_code, nosubReason) = cmt_ps_metrics_pass(nosubMeta)
                    def subMetas = metas.findAll { it.sub != null } ?: []
                    def subReasons = subMetas.collect { subMeta -> cmt_ps_metrics_pass_sub(nosubMeta, subMeta) }
                    def (asub_fail_code, asubReason) = subReasons.find { fail_code, reason -> fail_code != 0x00 }?:[null, null]
                    if (nosub_fail_code == 0x00 && asubReason != null) {
                        nosub_fail_code = asub_fail_code
                        nosubReason = asubReason
                    }
                    return [
                        [obsid, nosubMeta, failCodes[nosub_fail_code], nosubReason]
                    ] + [subMetas, subReasons].transpose().collect { subMeta, subReason ->
                        def (sub_fail_code, sub_reason) = subReason
                        [obsid, subMeta, failCodes[sub_fail_code], sub_reason?:nosubReason]
                    }
                }

            passReasons
                // .filter { _, __, fail_code, ___ -> fail_code != failCodes[0x00] }
                .map { obsid, meta, fail_code, reason ->
                    [obsid, meta.name, meta, fail_code, reason]
                }
                // .groupTuple( by: 0 )
                .map { obsid, name, meta, fail_code, reason ->
                    [
                        obsid,
                        name,
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
                        "P_WINDOW", "P_WEDGE", "NUM_CELLS",
                        "FAIL CODE", "REASON"
                    ].join("\t"),
                    storeDir: "${results_dir}${params.cal_suffix}"
                )
                // display output path and number of lines
                | view { [it, it.readLines().size()] }


            fail_codes = passReasons.map { obsid, meta, fail_code, reason ->
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
                    storeDir: "${results_dir}${params.cal_suffix}"
                )
                | view { it.readLines().size() }

            obsMetaPass = passReasons
                .filter { obsid, meta, fail_code, reason -> reason == null }
                .map { obsid, meta, _, __ -> [obsid, meta] }

            obsMetaUVPass_ = obsMetaUV.cross(obsMetaPass) { def (obsid, meta) = it; [obsid, meta.name] }
                .map { obsMetaUV_ , obsMetaPass_ ->
                    def (obsid, _, vis) = obsMetaUV_
                    def (__, meta) = obsMetaPass_
                    [obsid, meta, vis]
                }
        }

        // delay spectrum
        obsMetaUVEoR | delaySpec

    emit:
        // channel of files to archive, and their buckets
        archive = visQA.out.map { _, __, json -> ["visqa", json]}
        // channel of video name and frames to convert
        frame = plotVisQA.out
            .flatMap { _, meta, pngs ->
                coerceList(pngs).collect { png ->
                    def suffix = png.baseName.split('_')[-1]
                    [''+ "visqa_${meta.name}_${suffix}", png]
                }
            }
            .mix(uvPlot.out.flatMap {_, name, pngs -> coerceList(pngs).collect { png ->
                def pol = png.baseName.split('_')[-2..-1].join('_')
                [''+"visqa_${name}_${pol}", png]
            }})
            .mix(delaySpec.out.map { _, meta, png -> ["visqa_dlyspec_${meta.name}", png] })
            .groupTuple()
        // channel of files to zip
        zip = visQA.out.map { _, meta, json -> [''+"visqa_${meta.name}", json] }
            .groupTuple()
        obsMetaUVPass = obsMetaUVPass_
        fail_codes = fail_codes.filter { obsid, fail_code -> fail_code != failCodes[0x00] }
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
def wscleanDConvParams = mapMerge(wscleanParams, [
    args: "${params.wsclean_args} ${params.wsclean_dconv_args}",
    // args: params.wsclean_dconv_args
    niter: params.img_niter,
    minor_clean_gain: params.img_minor_clean_gain,
    major_clean_gain: params.img_major_clean_gain,
    auto_threshold: params.img_auto_threshold,
    auto_mask: params.img_auto_mask,
    mwa_path: params.img_mwa_path,
])


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

def polModes = [
    "I": [glob: "-I", pol: "I", idg: true], // Total intensity science, i.e., not interested in QUV: only image I with -pol I.
    "Q": [glob: "-Q", pol: "Q", idg: true],
    "U": [glob: "-U", pol: "U", idg: true],
    "V": [glob: "-V", pol: "V", idg: true],
    "IQUV": [glob: "-{I,Q,U,V}", pol: "IQUV", idg: true], // Total intensity science, but “nice to have QUV” for e.g. sensivity analysis
    "IQUV_join": [glob: "-{I,Q,U,V}", pol: "IQUV", join_pols: true, idg: true], // Interested in all stokes parameter, cleaning each polarization in a joined way
    "IV_join": [glob: "-{I,V}", pol: "i,v", join_pols: true, ],
    "IQUV_lqu": [glob: "-{Q,U}", pol: "IQUV", link_pols: "q,u", idg: true ], // Interested in rotation measure synthesis
    "QU_sq": [glob: "-{Q,U}", pol: "q,u", join_pols: true, sq_chan_join: true ],
    "XXYY_join": [glob: "-{XX,YY}", pol: "xx,yy", join_pols: true, ],
    "XXXYXX_lxxyy": [glob: "-{XX,YY,XY,XYi}", pol: "xx,xy,yx,yy", link_pols: "xx,yy"],
]

// image visibilities and QA images
workflow img {
    take:
        // tuple of (obsid, meta, vis)
        obsMetaVis
    main:

        // wsclean: make deconvolved images
        if (params.img_split_intervals) {
            // add -t???? suffix to name, to split wsclean over multiple jobs
            splitObsMetaVis = obsMetaVis.flatMap { obsid, meta, vis ->
                (0..(meta.ntimes-1)).collect { i ->
                    [obsid, mapMerge(meta, [interval: [i, i+1], inter_tok: sprintf("-t%04d", i)]), vis]
                }
            }
        } else {
            splitObsMetaVis = obsMetaVis
        }

        // print("params.nodeconv: ${params.nodeconv}")
        if (params.nodeconv) {
            splitObsMetaVis.map {obsid, meta, vis ->
                    def imgParams = [ img_name: ''+"${meta.cal_prog}_${obsid}${meta.subobs?:''}_${meta.name}${meta.inter_tok?:''}" ]
                    [obsid, meta, vis, mapMerge(wscleanParams, imgParams)]
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
            //         newImgParams = mapMerge(imgParams, wscleanDConvParams)
            //         [obsid, meta, vis, newImgParams, dirtyImgs]
            //     }

            splitObsMetaVis.map { obsid, meta, vis ->
                    def imgParams = [ img_name: ''+"${meta.cal_prog}_${obsid}${meta.subobs?:''}_${meta.name}${meta.inter_tok?:''}" ]
                    [obsid, meta, vis, mapMerge(wscleanDConvParams, imgParams), []]
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
                        if (multiinterval && !meta.inter_tok) {
                            imgParams.glob += "-t????"
                        }
                        if (multichannel) {
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
            .branch { obsid, meta, img ->
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
                    def high = parseCsv(hist, true, 2)
                        .find { row -> Float.compare(parseFloatOrNaN(row.quantile), params.thumbnail_quantile as Float) == 0 }
                    high = parseFloatOrNaN(high == null ? null : high.value)
                    def low = parseCsv(hist, true, 2)
                        .find { row -> Float.compare(parseFloatOrNaN(row.quantile), (1-params.thumbnail_quantile) as Float) == 0 }
                    low = parseFloatOrNaN(low == null ? null : low.value)
                    // subobs = meta.subobs?:''
                    [obsid, meta.inter_tok, meta.name, meta.inter_suffix, high, low]
                }
                .filter { _obsid, _interval, _name, _suff, high, low -> !(isNaN(high) || isNaN(low)) }
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
                    def _suffLimits = [
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
                    storeDir: "${results_dir}${params.img_suffix}${params.cal_suffix}"
                )
                // display output path and number of lines
                | view { [it, it.readLines().size()] }

            // value channel containing a map from img suffix to max limit excluding outliers
            // suffLimits = channel.from([:])
            suffLimits = imgLimits.combine(suffs_)
                .flatMap{ obsid, interval, name, limits, suffs ->
                    suffs.collect { suff ->
                        def hilo = limits[suff]?: [Float.NaN, Float.NaN]
                        [new Tuple(name, suff)] + hilo
                    }
                }
                .groupTuple(by: 0)
                .map { group, highs_, lows_ ->
                    thumbnail_vmax = params.thumbnail_vmax?: {
                        def highs = coerceList(highs_)
                            .findAll { !isNaN(it) }
                            .sort(false)
                        def high_index = ((highs.size() - 1) * params.thumbnail_quantile).round() as Integer
                        highs[high_index]
                    }()

                    thumbnail_vmin = params.thumbnail_vmin ?: {
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
                            (imgs, highs, _)  = order.collect { pol -> polImgLimits[pol] }.transpose()
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
            thumbnail.out.flatMap { obs, meta, png ->
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

        if (params.krvis) {
            obsMetaImgMfsPass
                .map { obsid, meta, img -> [obsid, meta.name, meta.inter_tok?:'', meta, img] }
                .groupTuple(by: 0..2)
                .filter { obsid, name, _, metas, __ -> metas.find { it.pol == 'I' } != null }
                .map { obsid, name, interval, metas, imgs ->
                    def iMeta = metas.find { it.pol == 'I' }
                    def (_, iquvImgs) = [metas, imgs].transpose().findAll { meta, img ->
                        ['I', 'Q', 'U', 'V'].contains(meta.pol)
                    }.transpose()
                    [obsid, iMeta, iquvImgs, file('/pawsey/mwa/mwaeor/dev/telescope_data_visualisation/romaO.npy')]
                }
                | krVis
        } else {
            channel.empty() | krVis
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
            fail_codes = channel.empty()
        } else {
            obsMetaImgGroup
                .map { obsid, meta, imgMetas, imgs -> [obsid, meta, imgs] }
                .filter { _, meta, imgs ->
                    // imgQA looks at pks flux, which needs to be deconvolved, only works with eor0
                    // meta.prod == 'image' &&
                    meta.eorfield == 0 &&
                    imgs.size() == 3
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
                        xx.RMS_ALL?:"", xx.RMS_BOX?:"", xx_pks.PEAK_FLUX?:"", xx_pks.INT_FLUX?:"",
                        yy.RMS_ALL?:"", yy.RMS_BOX?:"", yy_pks.PEAK_FLUX?:"", yy_pks.INT_FLUX?:"",
                        v.RMS_ALL?:"", v.RMS_BOX?:"", v_pks.PEAK_FLUX?:"", v_pks.INT_FLUX?:"",
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
                    def (asub_fail_code, asubReason) = subReasons.find { fail_code, reason -> fail_code != 0x00 }?:[null, null]
                    if (nosub_fail_code == 0x00 && asubReason != null) {
                        nosub_fail_code = asub_fail_code
                        nosubReason = asubReason
                    }
                    return [
                        [obsid, nosubMeta, failCodes[nosub_fail_code], nosubReason]
                    ] + [subMetas, subReasons].transpose().collect { subMeta, subReason ->
                        def (sub_fail_code, sub_reason) = subReason
                        [obsid, subMeta, failCodes[sub_fail_code], sub_reason?:nosubReason]
                    }
                }
                .tap { obsMetaReasons }

            obsMetaReasons
                .map { obsid, meta, fail_code, reason ->
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
                    storeDir: "${results_dir}${params.img_suffix}${params.cal_suffix}"
                )
                | view { it.readLines().size() }

            obsMetaReasons
                // .filter { _, __, fail_code, ___ -> fail_code != failCodes[0x00] }
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
                    storeDir: "${results_dir}${params.img_suffix}${params.cal_suffix}"
                )
                // display output path and number of lines
                | view { [it, it.readLines().size()] }

            obsMetaPass = obsMetaReasons
                .filter { obsid, meta, fail_code, reason -> reason == null }
                .map { obsid, meta, _, __ -> [obsid, meta] }

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
                [''+"imgqa_${meta.name}_polcomp_${meta.prod}_${meta.orderName}", png]
            }
            .mix(polMontage.out.map { _, meta, png ->
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
            .mix(krVis.out.flatMap {_, meta, pngs ->
                pngs.collect { png ->
                    def suffix = png.baseName.split('-')[-2..-1].join('-');
                    ["krvis_${meta.name}_${suffix}", png]
                }
            })
            .groupTuple()
        // channel of files to zip
        zip = imgQA.out.map { _, meta, json -> ["imgqa_${meta.name}", json] }
            .mix(obsMetaImgMfsPass.map { _, meta, img -> ["img_${meta.name}", img] })
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
        //     //             newMeta.chan_tok = sprintf("-ch%03d", ch)
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
                def (obsid, obsMeta, vis) = obsMetaUV_
                def (__, chunkMeta, chunk) = chunkMetaObs_
                [chunk, chunkMeta.name, chunkMeta, vis]
            }
            .groupTuple(by: 0..1)
            // .view { it -> "\n -> imgCombine obsMetaUV x chunkMetaObs ${it}\n" }
            .map { chunk, _, chunkMetas, viss ->
                [chunk, coerceList(chunkMetas)[0], viss.flatten(), deepcopy(wscleanParams)]
            }
            // .view { it -> "\n -> imgCombine before filter ${it}\n" }
            .filter { chunk, chunkMeta, viss, imgParams -> chunkMeta.name }

        if (params.img_split_coarse_chans && !params.noimgcombine) {
            if (params.sky_chans == null || params.sky_chans.size == 0) {
                throw new Exception("img_split_coarse_chans is enabled but params.sky_chans=${params.sky_chans}")
            }
            nCoarseChans = params.sky_chans.size
            print(" -> deleteme nCoarseChans=${nCoarseChans}")
            fineChansPerCoarse = obsMetaUV.map { obsid, meta, uvfits ->
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
                        def chan_tok = sprintf("-ch%03d", ch)
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
                    [chunk, meta, vis, mapMerge(wscleanDConvParams, imgParams), []]
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
            .filter { _, meta, img ->
                meta.prod !=~ /uv-.*/ \
                && (meta.chan?:-1 == -1 || params.thumbnail_all_chans) \
                && (meta.prod == (params.nodeconv ? "dirty" : "image"))
            }
            | thumbnail

        if (params.img_split_coarse_chans) {
            frame = thumbnail.out
                .filter { group, meta, png ->
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

    channel.of(obsids_file)
        .splitCsv()
        .filter { line -> !line[0].startsWith('#') }
        .map { line ->
            def (obsid, _) = line
            def vis = file("${params.outdir}/${obsid}/cal${params.cal_suffix}/${baseMeta.cal_prog}_${obsid}_${baseMeta.name}.uvfits")
            [ obsid, deepcopy(baseMeta), vis ]
        }
        .filter { obs, __, vis -> obs.size() > 0 && vis.exists() }
        .tap { obsMetaVis }
        // .view { "before uvMeta ${it}"}
        | uvMeta

    obsVis = obsMetaVis.map { obs, _, vis -> [obs, vis] }
    obsVis.join(uvMeta.out)
        .map { obs, vis, meta_, metaJson ->
            def uvmeta = parseJson(metaJson)
            def newMeta = [
                lowfreq: uvmeta.freqs[0],
                freq_res: (uvmeta.freqs[1] - uvmeta.freqs[0]),
                nchans: (uvmeta.freqs?:[]).size(),
                ntimes: (uvmeta.times?:[]).size(),
            ]
            ['config', 'eorband', 'num_ants', 'total_weight'].each { key ->
                newMeta[key] = uvmeta[key]
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
                "nchans", "eorband", "eorfield", "lowfreq", "freq_res",
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
        .flatMap { group, sorts, obss, gmetas, viss ->
            [sorts, obss, gmetas].transpose()
                .sort { it -> it[0] }
                .collate(params.chunkSize, params.chunkRemainder)
                .take(params.chunkCount)
                .collect { chunk ->
                    def (chunkSorts_, chunkObss, chunkMetas) = chunk.transpose()
                    def meta = coerceList(chunkMetas)[0]
                    def chunkSorts = coerceList(chunkSorts_)
                    // meta.obsids = chunkObss.sort(false)
                    obsids = coerceList(chunkObss).sort(false)
                    def newMeta = [
                        nobs: obsids.size(),
                        hash: obsids.join(' ').md5()[0..7],
                        sort_bounds: [chunkSorts[0], chunkSorts[-1]],
                    ]
                    ["${group}_${meta.hash}", mapMerge(meta, newMeta), obsids]
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
            .map { obsid, meta, vis_tot, vis_diff, noise_tot, noise_diff, weights ->
                [obsid, meta.name?:'', meta, [vis_tot, vis_diff, noise_tot, noise_diff, weights]]
            }
            .groupTuple(by: 0..1)
            .map { obsid, name, metas, grids ->
                [obsid, coerceList(metas)[0], grids.flatten()]
            }
            .cross(
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
                def newMeta = [ ext: ''+"${chunk}_${chunkMeta.name}"]
                // print("\n -> chunkMetas 0: ${chunkMeta}\n")
                [chunk, mapMerge(chunkMeta, newMeta), exts, grids.flatten()]
            }
            // .view { it -> "\n -> before chipsCombine ${it}\n" }
            .filter { chunk, chunkMeta, exts, grid -> chunkMeta.name }
            .combine(pols)
            .map { chunk, chunkMeta, exts, grid, pol ->
                [chunk, mapMerge(chunkMeta, [pol:pol]), exts, grid.findAll { it.baseName.contains(pol) }]
            }
            .filter { chunk, chunkMeta, exts, grids -> grids.size() > 0 }
            | chipsCombine

        chipsCombine.out.map { group, meta, _, dats -> [ group, meta, dats ] }
            | chipsLssa
            | chips1d_tsv

        chips1d_tsv.out
            .map { println("chips1d_tsv out ${it}") }

        all_k_modes = chips1d_tsv.out.flatMap { chunk, meta, tsv ->
                tsv.readLines().collect {
                    // try { it.split("\t")[0].toFloat() } catch(NumberFormatException e) { null }
                    it.split("\t")[0]
                }
                .findAll {
                    try {
                        it.split("\t")[0].toFloat()
                    } catch(NumberFormatException) {
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
                        table = parseCsv(coerceList(tsv)[0], true, 0, '\t')
                        def (deltas, noises) = k_modes.collect { k_mode ->
                            row = table.find { it.k_modes == k_mode }
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
                name: ''+"chips1d_delta_${params.lssa_bin}.tsv", newLine: true, sort: false,
                storeDir: "${results_dir}${params.img_suffix}${params.cal_suffix}"
            )
            .tap { chips1d_delta_tsv }
            | view { [it, it.readLines().size()] }

        // chips1d_delta_tsv.out.map { }

        singles1D = chipsLssa.out.map { chunk, meta, grid ->
                def newMeta = [
                    ptype: '1D',
                    pol: 'both',
                    title: ''+"crosspower\\n${chunk}\\n${meta.name}",
                    plot_name: 'chips1d',
                    max_power: 1e15,
                    min_power: 1e3,
                    tags: [],
                ]
                [chunk, mapMerge(meta, newMeta), grid]
            }

        singles1DDelta = chipsLssa.out.map { chunk, meta, grid ->
                def newMeta = [
                    ptype: '1D',
                    pol: 'both',
                    plot_delta: true,
                    title: ''+"crosspower\\n${chunk}\\n${meta.name}",
                    plot_name: 'chips1d',
                    max_power: 1e15,
                    min_power: 1e3,
                    tags: [],
                ]
                [chunk, mapMerge(meta, newMeta), grid]
            }

        singles2D = chipsLssa.out.map { chunk, meta, grid ->
                def newMeta = [
                    ptype: '2D',
                    pol: 'both',
                    title: ''+"crosspower\\n${chunk}\\n${meta.name}",
                    plot_name: 'chips2d',
                    max_power: 1e15,
                    min_power: 1e3,
                    tags: [],
                ]
                [chunk, mapMerge(meta, newMeta), grid]
            }

        comps1D = chipsLssa.out
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
                    max_power: 1e15,
                    min_power: 1e3,
                ]
                [chunk, mapMerge(nosubMeta, newMeta), grids.flatten()]
            }
        diffs2D_sub = chipsLssa.out
            // group grids from vis name together
            .groupTuple(by: 0)
            .filter { _, metas, __ -> metas.size() >= 2 }
            .map { chunk, metas, grids ->
                def nosubMeta = metas.find { it.sub == null} ?: [:]
                def subMeta = metas.find { it.sub == 'sub'} ?: [:]
                def newMeta = [
                    ptype: '2D_diff',
                    pol: 'both',
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
                [chunk, mapMerge(nosubMeta, newMeta), grids.flatten()]
            }
        diffs2D_iono = chipsLssa.out
            // group grids from vis name together
            .groupTuple(by: 0)
            .filter { _, metas, __ -> metas.size() >= 2 & metas.find { it.sub == null} != null & metas.find { it.sub == 'ionosub'} != null }
            .map { chunk, metas, grids ->
                def nosubMeta = metas.find { it.sub == null} ?: [:]
                def ionosubMeta = metas.find { it.sub == 'ionosub'} ?: [:]
                def newMeta = [
                    ptype: '2D_diff',
                    pol: 'both',
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
                [chunk, mapMerge(nosubMeta, newMeta), grids.flatten()]
            }

        diffs2D_sub_ionosub = chipsLssa.out
            // group grids from vis name together
            .groupTuple(by: 0)
            .filter { _, metas, __ -> metas.size() >= 2 & metas.find { it.sub == 'sub'} != null & metas.find { it.sub == 'ionosub'} != null }
            .map { chunk, metas, grids ->
                def ionosubMeta = metas.find { it.sub == 'ionosub'} ?: [:]
                def subMeta = metas.find { it.sub == 'sub'} ?: [:]
                def newMeta = [
                    ptype: '2D_diff',
                    pol: 'both',
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
                [chunk, mapMerge(ionosubMeta, newMeta), grids.flatten()]
            }

        ratios2D_sub = chipsLssa.out
            // group grids from vis name together
            .groupTuple(by: 0)
            .filter { _, metas, __ -> metas.size() >= 2 & metas.find { it.sub == 'sub'} != null & metas.find { it.sub == null} != null }
            .map { chunk, metas, grids ->
                def nosubMeta = metas.find { it.sub == null} ?: [:]
                def subMeta = metas.find { it.sub == 'sub'} ?: [:]
                def newMeta = [
                    ptype: '2D_ratio',
                    pol: 'both',
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
                [chunk, mapMerge(nosubMeta, newMeta), grids.flatten()]
            }

        ratios2D_ionosub = chipsLssa.out
            // group grids from vis name together
            .groupTuple(by: 0)
            .filter { _, metas, __ -> metas.size() >= 2 & metas.find { it.sub == null} != null & metas.find { it.sub == 'ionosub'} != null }
            .map { chunk, metas, grids ->
                def nosubMeta = metas.find { it.sub == null} ?: [:]
                def ionosubMeta = metas.find { it.sub == 'ionosub'} ?: [:]
                def newMeta = [
                    ptype: '2D_ratio',
                    pol: 'both',
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
                [chunk, mapMerge(nosubMeta, newMeta), grids.flatten()]
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
                    [''+"${meta.plot_name}_${meta.name}${suffix}", png]
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
        asvoPrep.out.map { _, __, vis -> ["prep", vis] }
            | archive
    }
}

// default entrypoint: get preprocessed vis from asvo and run qa
workflow {
    // get obsids from csv
    obsCSV = channel.of(obsids_file)
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
        prep(channel.empty())
    } else {
        prep(ws.out.obsMeta.join(ws.out.obsMetafits))
    }

    // channel of obsids that pass the flag gate
    prep.out.subobsMetaVisPass
        .map { obsid, meta, vis -> [obsid, meta.subobs?:'', meta.name?:'', vis].join('\t') }
        .collectFile(
            name: "prep_subobs_name_pass.csv", newLine: true, sort: true,
            storeDir: "${results_dir}"
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
        .filter { _, fail_code -> fail_code == failCodes[0x00] }
        .join(ws.out.obsMeta)
        .map { obsid, fail_code, meta ->
            [obsid, meta.lst.round().intValue()]
        }
        .groupTuple(by: 1)
        .map { obsids, lst ->
            [
                sprintf("%+03d", lst),
                obsids.size(),
            ].join("\t")
        }
        .collectFile(
            name: "lst_counts_all.tsv", newLine: true, sort: true,
            seed: ([ "LST", "COUNT" ]).join("\t"),
            storeDir: "${results_dir}${params.img_suffix}${params.cal_suffix}"
        )
        .view { it.readLines().size() }
        .map { tsv ->
            meta = [
                name: "lst_counts_all", title: "lst counts after qa",
                x: "LST", y: "COUNT"
            ]
            [meta, tsv]
        }
        | tsvScatterPlot

    // do ewp counts
    all_fail_codes
        .filter { _, fail_code -> fail_code == failCodes[0x00] }
        .join(ws.out.obsMeta)
        .map { obsid, fail_code, meta ->
            [obsid, meta.ew_pointing]
        }
        .groupTuple(by: 1)
        .map { obsids, ewp ->
            [
                ewp == null ? "" : sprintf("%+1d", ewp),
                obsids.size(),
            ].join("\t")
        }
        .collectFile(
            name: "ewp_counts_all.tsv", newLine: true, sort: true,
            seed: ([ "LST", "COUNT" ]).join("\t"),
            storeDir: "${results_dir}${params.img_suffix}${params.cal_suffix}"
        )

    all_fail_codes.groupTuple(by: 1)
        .map { obsids, fail_code ->
            [
                fail_code,
                obsids.size(),
            ].join("\t")
        }
        .collectFile(
            name: "fail_counts_all.tsv", newLine: true, sort: true,
            seed: ([ "FAIL CODE", "COUNT" ]).join("\t"),
            storeDir: "${results_dir}${params.img_suffix}${params.cal_suffix}"
        )
}

// given a channel of tuple (obs, metafits, vis), calibrate and analyse
workflow qaPrep {
    take:
        subobsMetaVis
        obsMetafits

    main:
        obsMetaMetafitsVis = obsMetafits.cross(subobsMetaVis).map{ obsMetafits_, subobsMetaVis_ ->
                def (obsid, metafits) = obsMetafits_
                def (_, meta, vis) = subobsMetaVis_
                [obsid, meta, metafits, vis]
            }
            .filter { _, meta, metafits, vis -> metafits != null && vis != null }
        obsMetafitsVis = obsMetaMetafitsVis.map { obsid, meta, metafits, vis -> [obsid, metafits, vis] }
        obsMeta = obsMetaMetafitsVis.map { obsid, meta, metafits, vis -> [obsid, meta] }

        // get sourcelists for each obs (only currently used in subtraction, not calibration)
        obsMetaMetafitsVis.map { obsid, meta, metafits, vis -> [obsid, metafits ] }
            .unique()
            | (hypSrclistAO & hypSrclistYaml)

        // obsMetafitsSrclist: channel of tuple(obsid, metafits, srclist)
        // cluster unless --nocluster
        if (params.nocluster) {
            channel.empty() | rexCluster
            obsMetafitsSrclist = obsMetafits.cross(hypSrclistAO.out)
                .map { obsMetafits_, hypSrclistAO_ ->
                    def (obsid, metafits) = obsMetafits_
                    def (_, srclist) = hypSrclistAO_
                    [obsid, metafits, srclist] }
        } else {
            hypSrclistAO.out | rexCluster
            obsMetafitsSrclist = obsMetafits.cross(rexCluster.out)
                .map { obsMetafits_, rexCluster_ ->
                    def (obsid, metafits) = obsMetafits_
                    def (_, srclist) = rexCluster_
                    [obsid, metafits, srclist] }
        }

        // calibrate each obs that passes flag gate unless --nocal:
        if (params.nocal) {
            // empty channel disables a process
            channel.empty() | cal
        } else {
            obsMetaMetafitsVis | cal
        }

        cal.out.obsMetaCalPass
            .map { def (obsid, meta, soln) = it;
                def calFlags = deepcopy(meta.calFlags?:[])
                ([obsid, meta.name, meta.subobs?:'', displayInts(meta.prepFlags?:[]), displayInts(calFlags)]).join("\t")
            }
            .collectFile(
                name: "pass_cal.tsv", newLine: true, sort: true,
                seed: ([ "OBS", "NAME", "SUBOBS", "PREPFLAGS", "CALFLAGS" ]).join("\t"),
                storeDir: "${results_dir}${params.cal_suffix}"
            )
            | view { [it, it.readLines().size()] }

        // channel of arguments for hypApply{UV,MS}
        // - take tuple(obsid, meta, soln) from cal.out.obsMetaCalPass
        // - match with tuple(obsid, metafits, prepUVFits) by obsid

        // apply calibration solutions to uvfits and ms unless --nouv or --noms
        obsMetafitsVis.cross(cal.out.obsMetaCalPass).map { obsMetafitsVis_, obsMetaCalPass_ ->
                def (obsid, metafits, prepUVFits) = obsMetafitsVis_
                def (_, meta, soln) = obsMetaCalPass_
                def newMeta = [
                    time_res: params.apply_time_res,
                    freq_res: params.apply_freq_res,
                    nodut1: params.nodut1,
                    apply_args: params.apply_args?:'',
                ]
                def calFlags = deepcopy(meta.calFlags?:[])
                if (calFlags) {
                    newMeta.apply_args = ''+"${params.apply_args?:''} --tile-flags ${calFlags.join(' ')}"
                }
                // calFlags = [8,9,10,11,12,13,15,27,30,39,40,42,44,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,79,88,89,90,91,92,93,94,95,117]
                // newMeta.apply_args = "${newMeta.apply_args} --tile-flags ${calFlags.join(' ')}"
                newMeta.apply_name = '' + "${newMeta.time_res}s_${newMeta.freq_res}kHz"
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
                def (obsid, meta, _, offsets, __) = hypIonoSubUV_;
                def (___, srclist) = hypSrclistYaml_;
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
        //             storeDir: "${results_dir}${params.cal_suffix}"
        //         )
        //         | view { [it, it.readLines().size()] }
        // }

        // collect ionoqa results as .tsv
        cthulhuPlot.out
            // form row of tsv from json fields we care about
            .map { obsid, _, __, ___, json, ____ ->
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
                storeDir: "${results_dir}${params.cal_suffix}"
            )
            // display output path and number of lines
            | view { [it, it.readLines().size()] }

        // obsIonoPass = cthulhuPlot.out.map { obsid, _, __, ___, json, ____ ->
        //         stats = parseJson(json);
        //         [obsid, stats.QA?:null]
        //     }
        //     .filter { _, qa ->
        //         qa != null && (params.filter_max_hyp_ionoqa == null || qa < params.filter_max_hyp_ionoqa)
        //     }
        //     .map { obsid, qa -> [obsid] }

        // channel of calibrated, subtracted and ionosubtracted uvfits: tuple(obsid, name, uvfits)
        obsMetaUV = hypApplyUV.out.map { obsid, meta, vis, _ -> [obsid, meta, vis] }
            .mix(hypIonoSubUV.out.map { obsid, meta, vis, _, __ -> [obsid, meta, vis] })
            .mix(hypSubUV.out.map { obsid, meta, vis, _ -> [obsid, meta, vis] })

        // improve uvfits meta
        obsMetaUV | uvMeta
        obsMetaUVSmart = obsMetaUV.cross(uvMeta.out) {[it[0], it[1].name]}
            .map { obsMetaUV_, uvMeta_ ->
                def (obsid, meta, vis) = obsMetaUV_;
                def (_, __, json) = uvMeta_;
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
                .filter { _, meta, __ ->
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
                storeDir: "${results_dir}${params.img_suffix}${params.cal_suffix}"
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
                        ["${group}_${hash}", mapMerge(groupMeta, chunkMeta), obs_list, metas]
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
            .map { chunk, name, obsids, metas ->
                [chunk, metas[0], obsids]
            }
            // finally, we have a channel of visibilities from the same obs group and vis type
            // .view { it -> "\n -> chunkMetaPass: ${it}"}
            .tap { chunkMetaPass }

        groupChunkMetaPass
            .map { chunk, chunkMeta, obsids, all_metas ->
                sort_bounds = chunkMeta.sort_bounds?:[Float.NaN, Float.NaN]
                (
                    sort_bounds.collect { sprintf("%9.9f", it) }
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
            .unique()
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

        uvfits.out.obsMetaUVPass.cross(obsMetaPass) { def (obsid, meta) = it; [obsid, meta.name] }
            .map { obsMetaUVPass_, obsMetaPass_ ->
                def (obsid, _, uvfits) = obsMetaUVPass_;
                def (__, meta) = obsMetaPass_;
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
        obsMetaUVPass
        frame = cal.out.frame
            .mix(uvfits.out.frame)
            .mix(img.out.frame)
            .mix(chips.out.frame)
            .mix(cthulhuPlot.out.flatMap {
                def (_, meta, pngs, __, ___, tecs) = it;
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
                hypIonoSubUV.out.map { it ->
                    def (obsid, meta, _, offsets) = it;
                    [''+"offsets_${meta.name}", offsets]
                }
                .groupTuple()
            )

        fail_codes = uvfits.out.fail_codes.join(img.out.fail_codes, remainder: true)
            .map { obsid, uv_code, img_code ->
                codes = [uv_code, img_code].findAll { it != null }
                aFailCode = coerceList(firstFail([codes]))[0]
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
                def cachebust = ''+"${latest}_x" + sprintf("%04d", files.size())
                // sorted = files.collect { path -> file(deepcopy(path.toString())) }.sort(false)
                [name, files, cachebust]
            }
            | tarchive
            | view { zip_, cachebust -> [zip_, zip_.size()] }
}

// entrypoint: get externally preprocessed uvfits files from csv file and run qa
workflow extPrep {

    def name = params.visName ?: "ssins"
    def cal_prog = params.cal_prog ?: "hyp"

    channel.of(obsids_file)
        .splitCsv()
        .filter { line -> !line[0].startsWith('#') }
        .map { line ->
            def (obsid, _) = line
            obsid
        }
        .unique()
        .map { obsid_ ->
            def obsid = coerceList(obsid_)[0]
            def meta = [name:name, cal_prog:cal_prog, obsid: obsid]
            def vis = file("${params.outdir}/${obsid}/prep/birli_${obsid}_2s_40kHz.${name}.uvfits")
            [ obsid, deepcopy(meta), vis ]
        }
        .filter { _, __, vis -> vis.exists() }
        .tap { obsMetaVis }
        .map { obsid, meta, vis -> obsid }
        .tap { obsids }

    obsids | ws

    obsWsmetaVis = obsMetaVis.join(ws.out.obsMeta)
        .map { obsid, meta, vis, wsMeta ->
            [ deepcopy(obsid), mapMerge(meta, wsMeta), vis ]
        }

    flag(obsWsmetaVis, ws.out.obsMetafits)

    obsFlagmetaVis = flag.out.subobsMetaPass.map { obsid, meta -> [[obsid, meta.subobs?:''], meta] }
            .join(obsWsmetaVis.map { obsid, meta, uvfits -> [[obsid, meta.subobs?:''], uvfits] })
            // .join(flag.out.subobsMetaRxAnts.map { obsid, meta, rxAnts -> [[obsid, meta.subobs?:''], rxAnts] })
            .map { obsSubobs, meta, uvfits ->
                def (obsid, _) = obsSubobs
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
        .map { obsids, fail_code ->
            [ fail_code, obsids.size() ].join("\t")
        }
        .collectFile(
            name: "fail_counts_all.tsv", newLine: true, sort: true,
            seed: ([ "FAIL CODE", "COUNT" ]).join("\t"),
            storeDir: "${results_dir}${params.img_suffix}${params.cal_suffix}"
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

workflow tsvTest {
    chan = channel.of(
        [0, [a:2, b:3]],
        [1, [b:4, c:5]]
    )
    // N	a	b	c
    // 0	2	3
    // 1		4	5
    keys = chan.flatMap { number, hash -> hash.keySet() }
        .unique()

    channel.of("N").concat(keys)
        .toList()
        .map { it.join("\t") }
        .concat(
            chan.combine(keys.toSortedList().map{[it]}).map { number, hash, keys ->
                    ([number] + keys.collect { key -> hash[key]?:'' }).join("\t")
                }
                .toSortedList()
                .flatten()
        ).collectFile(
            name: "test.tsv", newLine: true, sort: false,
            storeDir: "${results_dir}"
        )
        .view { [it, it.readLines().size()] }
}
