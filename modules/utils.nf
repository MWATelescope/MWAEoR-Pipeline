def obsids_file() {
    obsids_file = file(params.obsids_path)
    if (params.obsids_suffix) {
        obsids_file = file("${obsids_file.parent}/${obsids_file.baseName}${params.obsids_suffix}.${obsids_file.extension}")
    }
    obsids_file
}

def results_dir() {
    "${params.resultsdir}/results${params.obsids_suffix}" + (params.visName ? ".${params.visName}" : "") + "${params.result_suffix}"
}

// whether imaging is configured for multiple channels
def is_multichannel() {
    def img_channels_out = (params.img_channels_out instanceof String
        ? params.img_channels_out.split(' ')[0]
        : params.img_channels_out)
    return (img_channels_out as int > 1)
}
// whether imaging is configured for multiple intervals
def is_multiinterval() {
    return (params.img_intervals_out as int > 1) || params.img_split_intervals
}

def coerceList(x) {
    if (x instanceof nextflow.util.ArrayBag) {
        x.toList()
    }
    else if (x instanceof List) {
        x
    }
    else {
        [x]
    }
}

def deepcopy(orig) {
    def bos = new ByteArrayOutputStream()
    def oos = new ObjectOutputStream(bos)
    try {
        oos.writeObject(orig)
        oos.flush()
    }
    catch (Exception e) {
        println("error deepcopying ${orig} ${e}")
        org.codehaus.groovy.runtime.StackTraceUtils.sanitize(e).printStackTrace()
        error(e)
    }
    def byteArray = bos.toByteArray()
    def ident = orig.getClass().getName()
    if (ident =~ /.*HashMap.*/) {
        ident += ' ' + orig.keySet().toString()
    }
    println("deepcopy ${byteArray.size()} bytes ${ident}")
    def bin = new ByteArrayInputStream(byteArray)
    def ois = new ObjectInputStream(bin)
    return ois.readObject()
}

def mapMerge(a, b) {
    // > a = [a:1]
    // > b = a + [a:2]
    // > println a
    // [a:1]
    return deepcopy(a + b)
}

def exitCodes() {
    // sources:
    // - man 7 signal
    // - https://www.gnu.org/software/wget/manual/html_node/Exit-Status.html
    return [
        1: "general or permission",
        2: "incorrect usage",
        3: "File I/O error",
        4: "Network failure",
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
        140: "SIGUSR2 - out of time"
    ]
}


// get birli argstr (asvo), argstr (cli), and filename suffix from params
def birli_argstr_suffix() {
    // asvo argstr reference: https://github.com/MWATelescope/manta-ray-client#conversion-job-options
    // birli cli reference: https://github.com/MWATelescope/Birli#usage
    def suffix = params.prep_suffix ?: ''
    def args_asvo = [output: "uvfits"]
    def args_cli = [:]
    if (params.van_vleck) {
        args_cli['van-vleck'] = null
        suffix += "_vv"
    }
    if (params.prep_time_res_s != null) {
        args_asvo.avg_time_res = params.prep_time_res_s
        args_cli['avg-time-res'] = params.prep_time_res_s
        suffix += "_${params.prep_time_res_s}s"
    }
    if (params.prep_freq_res_khz != null) {
        // validate this resolution is supported
        args_asvo.avg_freq_res = params.prep_freq_res_khz
        args_cli['avg-freq-res'] = params.prep_freq_res_khz
        suffix += "_${params.prep_freq_res_khz}kHz"
    }
    if (params.prep_rfi != null && !params.prep_rfi) {
        args_asvo.no_rfi = "true"
        args_cli['no-rfi'] = null
        suffix += '_norfi'
    }
    if (params.prep_pc != null && params.prep_pc.size() == 2) {
        // --phase-centre <RA> <DEC>
        args_cli['phase-centre'] = params.prep_pc

        args_asvo.centre = "custom"
        // <ra formatted as: 0.0 deg> ICRS (J2000.0).
        args_asvo.phase_centre_ra = String.format("%f", params.prep_edge_width[0])
        // <dec formatted as: +00.0 deg> ICRS (J2000.0).
        args_asvo.phase_centre_dec = String.format("%+f", params.prep_edge_width[1])
    }
    if (params.prep_edge_width != null) {
        args_cli['flag-edge-width'] = params.prep_edge_width
        args_asvo.flag_edge_width = params.prep_edge_width
        suffix += "_edg${params.prep_edge_width}"
    }
    def argstr_asvo = args_asvo.collect { k, v -> '' + "${k}=${v}" }.join(',')
    def argstr_cli = args_cli
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
    [argstr_asvo, argstr_cli, suffix]
}

def hyp_apply_name(time_res, freq_res) {
    def suffix = ''
    if (time_res != null) {
        suffix += "${time_res}s"
    }
    if (freq_res != null) {
        if (suffix) {
            suffix += "_"
        }
        suffix += "${freq_res}kHz"
    }
    suffix
}

//


def groovy2bashAssocArray(map, name) {
    // size = map.collect { _, v -> v.size() }.max()
    '' + "declare -A ${name}=(" + map.collect { k, v -> "[${k}]=" + '"' + "${v}" + '"' }.join(" ") + ")".toString()
}


// TODO: this is nasty and is probably a source of many memory
def parseJson(path) {
    // static final jslurp = new groovy.json.JsonSlurperClassic()
    def getSlurper = {
        new groovy.json.JsonSlurperClassic()
    }.memoize()

    // TODO: fix nasty hack to deal with NaNs
    try {
        // try reading the file first, it may not be ready yet.
        def size = file(path).size()
        print("reading ${size} bytes from ${path}")
    }
    catch (Exception _e) {
        Thread.sleep(5)
    }
    try {
        def text = path.getText().replaceAll(/(NaN|-?Infinity)/, '"$1"')
        def parsed = getSlurper().parseText(text)
        return deepcopy(parsed)
    }
    catch (Exception e) {
        println("error parsing ${path} ${e}")
    }
}

def parseCsv(path, header=true, skipBytes=0, delimeter=',') {
    def allLines = path.readLines()
    if (!header) {
        return allLines.collect { it.split(delimeter) }
    }
    else {
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
    }
    else if (s instanceof Float) {
        return s
    }
    else if (s instanceof java.math.BigDecimal) {
        return s.floatValue()
    }
    try {
        return Float.parseFloat(s)
    }
    catch (NumberFormatException _e) {
        return Float.NaN
    }
}

// def logDateFmt() {
//     new java.text.SimpleDateFormat("yyyy-MM-dd HH:mm:ss")
// }

def isNaN(n) {
    if (n == null || n == "") {
        return true
    }
    // def cls = null
    try {
        // cls = n.getClass();
        return (n as float).isNaN()
    }
    catch (MissingMethodException e) {
        org.codehaus.groovy.runtime.StackTraceUtils.sanitize(e).printStackTrace()
        error("${e}")
        return true
    }
}

// display a long list of ints, replace bursts of consecutive numbers with ranges
// start, end, delim
def displayRange(s, e, d=',') {
    return s == e ? "${s}${d}" : s == e - 1 ? "${s}${d}${e}${d}" : "${s}-${e}${d}"
}

def displayInts(l_, delim=',') {
    def max_ints_display = 50
    def mid_ints_display = (max_ints_display - 1).intdiv(2)
    def l = (l_ as ArrayList).sort(false).unique()
    if (l.size() == 0) {
        return ""
    }
    if (l.size() == 1) {
        return "${l[0]}"
    }
    def (sb, start, end) = ['' << '', l[0], l[0]]
    l[1..-1].each { i ->
        if (i == end + 1) {
            // Continue the current range
            end = i
        }
        else {
            // End current range and start new one
            sb << displayRange(start, end, delim)
            start = i
            end = i
        }
    }
    def result = (sb << displayRange(start, end, delim))[0..-2].toString()
    def size = result.size()
    if (size > max_ints_display) {
        try {
            return result.substring(0, mid_ints_display) + "..." + result.substring(size - mid_ints_display)
        }
        catch (StringIndexOutOfBoundsException e) {
            println("error: ${e} ${result} ${size} ${mid_ints_display}")
            error(e)
        }
    }
    else {
        return result
    }
}

// collapse a list into its contiguous ranges
def contigRanges(l) {
    if (l.size() == 0) {
        return []
    }
    if (l.size() == 1) {
        return [(l[0]..l[0])]
    }
    def (sb, start, end) = [[], l[0], l[0]]
    l
        .tail()
        .each { i ->
            if (i == end + 1) {
                // Continue current range
                end = i
            }
            else {
                // End current range, start new one
                sb << displayRange(start, end)
                start = i
                end = i
            }
        }
    return (sb << (start..end))[0..-1]
}

def get_seconds(time, unit) {
    def unit_conversions = [
        s: 1,
        ms: 0.001,
        µs: 0.000001,
        ns: 1E-9
    ]
    if (unit_conversions[unit] == null) {
        println("unknown duration unit ${unit} for time ${time}")
        time
    }
    else {
        time * unit_conversions[unit]
    }
}

def prepqa_pass(flagMeta) {
    def reasons = []
    def fail_code = 0
    // no error
    if (params.flag_occupancy_threshold != null && flagMeta.total_occ > params.flag_occupancy_threshold) {
        reasons += "total_occ(${String.format('%.2f', flagMeta.total_occ)})>${params.flag_occupancy_threshold}"
        fail_code = fail_code == 0 ? 49 : fail_code
    }
    if (params.rfi_occupancy_threshold != null && flagMeta.total_non_preflagged_bl_occ > params.rfi_occupancy_threshold) {
        reasons += "rfi_occ(${String.format('%.2f', flagMeta.total_non_preflagged_bl_occ)})>${params.rfi_occupancy_threshold}"
        fail_code = fail_code == 0 ? 50 : fail_code
    }
    if (params.ssins_occupancy_threshold != null && flagMeta.ssins_total > params.ssins_occupancy_threshold) {
        reasons += "ssins_total(${String.format('%.2f', flagMeta.ssins_total)})>${params.ssins_occupancy_threshold}"
        fail_code = fail_code == 0 ? 51 : fail_code
    }
    if (params.ssins_streak_threshold != null && flagMeta.ssins_streak > params.ssins_streak_threshold) {
        reasons += "ssins_streak(${String.format('%.2f', flagMeta.ssins_streak)})>${params.ssins_streak_threshold}"
        fail_code = fail_code == 0 ? 52 : fail_code
    }
    if (params.ssins_narrow_threshold != null && flagMeta.ssins_narrow_total > params.ssins_narrow_threshold) {
        reasons += "ssins_narrow_total(${String.format('%.2f', flagMeta.ssins_narrow_total)})>${params.ssins_narrow_threshold}"
        fail_code = fail_code == 0 ? 53 : fail_code
    }
    if (params.ssins_dab_threshold != null && flagMeta.ssins_dab_total > params.ssins_dab_threshold) {
        reasons += "ssins_dab_total(${String.format('%.2f', flagMeta.ssins_dab_total)})>${params.ssins_dab_threshold}"
        fail_code = fail_code == 0 ? 54 : fail_code
    }
    if (!params.noprepqafilter && (flagMeta.prep_status ?: "GOOD") != "GOOD") {
        reasons += "prep_status=${flagMeta.prep_status}"
        fail_code = fail_code == 0 ? 55 : fail_code
    }

    if (params.filter_bad_tile_frac != null && flagMeta.n_tiles != null && (flagMeta.prepFlags ?: []).size() > 0) {
        // print("prepqa_pass filter_bad_tile_frac=${params.filter_bad_tile_frac} n_tiles=${flagMeta.n_tiles} flagAnts=${flagMeta.flagAnts?:[]} prepFlags=${flagMeta.prepFlags?:[]}")
        def n_tiles = flagMeta.n_tiles
        def n_bad_tiles = ((((flagMeta.flagAnts ?: []) as Set) + ((flagMeta.prepFlags ?: []) as Set)) as ArrayList).size()
        if (n_bad_tiles > params.filter_bad_tile_frac * n_tiles) {
            reasons += "n_bad_tiles(${n_bad_tiles}) > ${params.filter_bad_tile_frac} * n_tiles(${n_tiles})"
            fail_code = fail_code == 0 ? 63 : fail_code
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
    def fail_code = 0
    if (params.filter_max_rms_convg != null && meta.rms_convg != null) {
        if (meta.rms_convg > params.filter_max_rms_convg) {
            reasons += "rms_convg(${meta.rms_convg}) > ${params.filter_max_rms_convg}"
            fail_code = fail_code == 0 ? 64 : fail_code
        }
    }
    if (params.filter_max_unused_bl_frac != null && meta.unused_bl_frac != null) {
        if (meta.unused_bl_frac > params.filter_max_unused_bl_frac) {
            reasons += "unused_bl_frac(${meta.unused_bl_frac}) > ${params.filter_max_unused_bl_frac}"
            fail_code = fail_code == 0 ? 65 : fail_code
        }
    }
    if (params.filter_max_unconvg_chs != null && meta.unconvg_chs != null) {
        if (meta.unconvg_chs > params.filter_max_unconvg_chs) {
            reasons += "unconvg_chs(${meta.unconvg_chs}) > ${params.filter_max_unconvg_chs}"
            fail_code = fail_code == 0 ? 66 : fail_code
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
    def fail_code = 0
    def is_unsub = ((meta.sub ?: "") == "")
    if (params.filter_max_ps_window_unsub != null && meta.p_window != null) {
        if (meta.p_window > params.filter_max_ps_window_unsub) {
            return [96, "${meta.sub ?: ""}_p_win(${meta.p_window}) > max_ps_window_unsub(${params.filter_max_ps_window_unsub})"]
        }
    }
    if (is_unsub && params.filter_min_ps_window_unsub != null && meta.p_window != null) {
        if (meta.p_window < params.filter_min_ps_window_unsub) {
            return [97, "p_win(${meta.p_window}) < min_ps_window_unsub(${params.filter_min_ps_window_unsub})"]
        }
    }
    return [fail_code, null]
}

// POWER	Normal P_win/P_wg	< 0.1	Window power a small fraction of wedge power
// POWER	P_wg_sub/P_wg	< 0.3	More than 70% wedge power subtracted
// POWER	P_win_sub/P_win	< 1.0, >0.1	Window power not crazy after subtraction
def cmt_ps_metrics_pass_sub(nosubMeta, subMeta) {
    def (fail_code, subReason) = cmt_ps_metrics_pass(subMeta)
    if (fail_code != 0) {
        return [fail_code, subReason]
    }

    // filter ps_ratio unsub
    if (nosubMeta.p_window != null && nosubMeta.p_wedge != null) {
        def p_win_p_wg = nosubMeta.p_window / nosubMeta.p_wedge
        if (params.filter_max_ps_ratio_unsub != null && p_win_p_wg > params.filter_max_ps_ratio_unsub) {
            return [98, "p_win:p_wg(${p_win_p_wg}) > max_ps_ratio_unsub(${params.filter_max_ps_ratio_unsub})"]
        }
        if (params.filter_min_ps_ratio_unsub != null && p_win_p_wg < params.filter_min_ps_ratio_unsub) {
            return [99, "p_win:p_wg(${p_win_p_wg}) < min_ps_ratio_unsub(${params.filter_min_ps_ratio_unsub})"]
        }
    }

    // filter ps_window sub ratio
    if (subMeta.p_window != null && nosubMeta.p_window != null) {
        def sub_p_win = subMeta.p_window / nosubMeta.p_window
        if (params.filter_max_ps_win_ratio_sub != null && sub_p_win > params.filter_max_ps_win_ratio_sub) {
            return [100, "sub_p_win:p_win(${sub_p_win}) > max_ps_win_ratio_sub(${params.filter_max_ps_win_ratio_sub})"]
        }
        if (params.filter_min_ps_win_ratio_sub != null && sub_p_win < params.filter_min_ps_win_ratio_sub) {
            return [101, "sub_p_win:p_win(${sub_p_win}) < min_ps_win_ratio_sub(${params.filter_min_ps_win_ratio_sub})"]
        }
    }

    // filter ps_wedge sub ratio
    if (subMeta.p_wedge != null && nosubMeta.p_wedge != null) {
        def sub_p_wg = subMeta.p_wedge / nosubMeta.p_wedge
        if (params.filter_max_ps_wedge_ratio_sub != null && sub_p_wg > params.filter_max_ps_wedge_ratio_sub) {
            return [102, "sub_p_wg:p_wg(${sub_p_wg}) > max_ps_wedge_ratio_sub{${params.filter_max_ps_wedge_ratio_sub}}"]
        }
        if (params.filter_min_ps_wedge_ratio_sub != null && sub_p_wg < params.filter_min_ps_wedge_ratio_sub) {
            return [103, "sub_p_wg:p_wg(${sub_p_wg}) < min_ps_wedge_ratio_sub{${params.filter_min_ps_wedge_ratio_sub}}"]
        }
    }

    return [fail_code, subReason]
}

// rules that apply to both subMeta and nosubMeta
def cmt_imgqa_pass(_meta) {
    def fail_code = 0

    // if (params.filter_max_pks_int_diff != null && meta.xx_pks_int != null && meta.yy_pks_int != null) {
    //     pks_int_diff = (meta.xx_pks_int - meta.yy_pks_int).abs()
    //     if (pks_int_diff > params.filter_max_pks_int_diff) {
    //         fail_code = ((meta.sub?:"")=="") ? 0x74 : 0x75
    //         return [fail_code, "${meta.sub?:""}_pks_int |xx-yy| (${pks_int_diff}) > max_pks_int_diff(${params.filter_max_pks_int_diff})"]
    //     }
    // }
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
    def fail_code = 0

    // filter vrms_box_nosub
    if (nosubMeta != null && nosubMeta.v_rms_box != null) {
        if (params.filter_max_vrms_box_nosub && nosubMeta.v_rms_box > params.filter_max_vrms_box_nosub) {
            return [112, "unsub_v_rms_box(${nosubMeta.v_rms_box}) > max_vrms_box(${params.filter_max_vrms_box_nosub})"]
        }
        if (params.filter_min_vrms_box_nosub && nosubMeta.v_rms_box < params.filter_min_vrms_box_nosub) {
            return [113, "unsub_v_rms_box(${nosubMeta.v_rms_box}) < max_vrms_box(${params.filter_min_vrms_box_nosub})"]
        }
    }

    // filter pks_int_v_ratio_nosub
    if (nosubMeta != null && nosubMeta.xx_pks_int != null && nosubMeta.yy_pks_int != null && nosubMeta.v_pks_int != null) {
        def pks_int_v_ratio_nosub = nosubMeta.v_pks_int / (nosubMeta.xx_pks_int + nosubMeta.yy_pks_int)
        if (params.filter_max_pks_int_v_ratio_nosub != null && pks_int_v_ratio_nosub > params.filter_max_pks_int_v_ratio_nosub) {
            return [114, "pks_int ${nosubMeta.sub ?: ''} v/(xx+yy) (${pks_int_v_ratio_nosub}) > max_pks_int_v_ratio_nosub(${params.filter_max_pks_int_v_ratio_nosub})"]
        }
        if (params.filter_min_pks_int_v_ratio_nosub != null && pks_int_v_ratio_nosub < params.filter_min_pks_int_v_ratio_nosub) {
            return [115, "pks_int ${nosubMeta.sub ?: ''} v/(xx+yy) (${pks_int_v_ratio_nosub}) < min_pks_int_v_ratio_nosub(${params.filter_min_pks_int_v_ratio_nosub})"]
        }
    }

    // filter pks_int_diff_sub
    if (subMeta != null && subMeta.xx_pks_int != null && subMeta.yy_pks_int != null) {
        def pks_int_diff_sub = (subMeta.xx_pks_int - subMeta.yy_pks_int).abs()
        if (params.filter_max_pks_int_diff_sub != null && pks_int_diff_sub > params.filter_max_pks_int_diff_sub) {
            return [116, "${subMeta.sub ?: ""}_pks_int |xx-yy| (${pks_int_diff_sub}) > max_pks_int_diff_sub(${params.filter_max_pks_int_diff_sub})"]
        }
        if (params.filter_min_pks_int_diff_sub != null && pks_int_diff_sub < params.filter_min_pks_int_diff_sub) {
            return [117, "${subMeta.sub ?: ""}_pks_int |xx-yy| (${pks_int_diff_sub}) < min_pks_int_diff_sub(${params.filter_min_pks_int_diff_sub})"]
        }
    }

    // filter pks_int_sub_ratio xx
    if (nosubMeta != null && subMeta != null && nosubMeta.xx_pks_int != null && subMeta.xx_pks_int != null) {
        def pks_int_sub_ratio_xx = subMeta.xx_pks_int / nosubMeta.xx_pks_int
        if (params.filter_max_pks_int_sub_ratio_xx != null && pks_int_sub_ratio_xx > params.filter_max_pks_int_sub_ratio_xx) {
            return [118, "${subMeta.sub ?: ""}_xx_pks_int:xx_pks_int (${pks_int_sub_ratio_xx}) > max_pks_int_sub_ratio_xx(${params.filter_max_pks_int_sub_ratio_xx})"]
        }
        if (params.filter_min_pks_int_sub_ratio_xx != null && pks_int_sub_ratio_xx < params.filter_min_pks_int_sub_ratio_xx) {
            return [119, "${subMeta.sub ?: ""}_xx_pks_int:xx_pks_int (${pks_int_sub_ratio_xx}) < min_pks_int_sub_ratio_xx(${params.filter_min_pks_int_sub_ratio_xx})"]
        }
    }
    // filter pks_int_sub_ratio yy
    if (nosubMeta != null && subMeta != null && nosubMeta.yy_pks_int != null && subMeta.yy_pks_int != null) {
        def pks_int_sub_ratio_yy = subMeta.yy_pks_int / nosubMeta.yy_pks_int
        if (params.filter_max_pks_int_sub_ratio_yy != null && pks_int_sub_ratio_yy > params.filter_max_pks_int_sub_ratio_yy) {
            return [120, "${subMeta.sub ?: ""}_yy_pks_int:yy_pks_int (${pks_int_sub_ratio_yy}) > max_pks_int_sub_ratio_yy(${params.filter_max_pks_int_sub_ratio_yy})"]
        }
        if (params.filter_min_pks_int_sub_ratio_yy != null && pks_int_sub_ratio_yy < params.filter_min_pks_int_sub_ratio_yy) {
            return [121, "${subMeta.sub ?: ""}_yy_pks_int:yy_pks_int (${pks_int_sub_ratio_yy}) < min_pks_int_sub_ratio_yy(${params.filter_min_pks_int_sub_ratio_yy})"]
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
    if (is_multichannel()) {
        def chan_tok = tokens.removeLast()
        if (chan_tok != "MFS") {
            meta.chan_tok = chan_tok
            try {
                meta.chan = (chan_tok as int)
            }
            catch (NumberFormatException e) {
                print("error parsing channel ${chan_tok} from ${img}. ${meta}")
                error(e)
            }
        }
    }
    // suffix without interval
    meta.inter_suffix = [meta.chan_tok, meta.pol, meta.prod].join('-')
    meta.suffix = meta.inter_suffix
    def inter_tok = tokens.removeLast()
    if (is_multiinterval() && inter_tok =~ 't[0-9]{4}') {
        meta.inter_tok = inter_tok
        meta.inter = inter_tok[1..-1] as int
        meta.suffix = "${inter_tok}-${meta.inter_suffix}"
    }
    deepcopy(meta)
}

def groupMeta(meta) {
    def newMeta = [:]
    newMeta.sort = 1
    if (meta.p_window ?: Float.NaN != Float.NaN && meta.p_wedge ?: Float.NaN != Float.NaN) {
        newMeta.sort *= parseFloatOrNaN(meta.p_window) / parseFloatOrNaN(meta.p_wedge)
    }
    if (meta.total_weight ?: Float.NaN != Float.NaN) {
        def weight = parseFloatOrNaN(meta.total_weight) / 1E+9
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
        group_tokens << String.format("ewp%+1d", meta.ew_pointing)
    }
    if (params.groupByLst && meta.lst != null) {
        def nearest_lst = ((meta.lst.round().intValue() + 180) % 360 - 180)
        group_tokens << String.format("lst%+03d", nearest_lst)
    }

    newMeta.group = group_tokens.join('_')

    [
        "cal_prog",
        "time_res",
        "freq_res",
        "lowfreq",
        "nchans",
        "eorband",
        "eorfield",
        "config",
        "name",
        "sub"
    ].each { key ->
        if (meta[key] != null) {
            newMeta[key] = meta[key]
        }
    }
    deepcopy(newMeta)
}

def getFailReason(code) {
    def failCodes = [
        0: "000 - clean",
        1: "001 - phase centre",
        2: "002 - channel selection",
        3: "003 - pointing",
        4: "004 - sun elevation",
        5: "005 - capture mode",
        6: "006 - array configuration",
        16: "016 - bad tiles",
        17: "017 - dead dipoles",
        32: "032 - data quality",
        33: "033 - large iono qa",
        34: "034 - null iono qa",
        35: "035 - no files",
        49: "049 - total occupancy",
        50: "050 - rfi occupancy",
        51: "051 - ssins occupancy",
        52: "052 - ssins streak",
        53: "053 - ssins narrow",
        54: "054 - ssins dab",
        55: "055 - prep status",
        63: "063 - prepqa bad tiles",
        64: "064 - high calqa rms convg",
        65: "065 - high calqa unused bls",
        66: "066 - high calqa unconvg chs",
        68: "068 - high calqa skewness",
        69: "069 - low calqa skewness",
        70: "070 - high calqa rx_var",
        71: "071 - low calqa rx_var",
        72: "072 - high calqa dfft_pow",
        73: "073 - low calqa dfft_pow",
        79: "079 - calqa bad tiles",
        96: "096 - large unsub p_win",
        97: "097 - small unsub p_win",
        98: "098 - large unsub p_win:p_wg",
        99: "099 - small unsub p_win:p_wg",
        100: "100 - large sub_p_win:p_win",
        101: "101 - small sub_p_win:p_win",
        102: "102 - large sub_p_wg:p_wg",
        103: "103 - small sub_p_wg:p_wg",
        112: "112 - large unsub v_rms_box",
        113: "113 - small unsub v_rms_box",
        114: "114 - large unsub pks_int v/(xx+yy)",
        115: "115 - small unsub pks_int v/(xx+yy)",
        116: "116 - large unsub pks_int |xx-yy|",
        117: "117 - small unsub pks_int |xx-yy|",
        118: "118 - large sub_xx_pks_int:xx_pks_int",
        119: "119 - small sub_xx_pks_int:xx_pks_int",
        120: "120 - large sub_yy_pks_int:yy_pks_int",
        121: "121 - small sub_yy_pks_int:yy_pks_int"
    ]
    failCodes[code] ?: "unknown fail code ${code}"
}

// pass in a list of codes and (optianally) other associated items,
// return the first non-clean code with its associated items
//
// example:
// codes = [0x00, 0x31]
// reasons = ["clean", "total_occ"]
// firstFail([codes, reasons]) => [0x31, "total_occ"]
def firstFail(pt) {
    if (pt == null || pt.size() == 0) {
        return [0, null]
    }
    def tp = pt.transpose()
    if (tp == null || tp.size() == 0) {
        return [0, null]
    }
    try {
        def failures = tp.findAll { it -> it[0] != getFailReason(0) }
        if (failures == null || failures.size() == 0) {
            return tp[0]
        }
        else {
            return failures[0]
        }
    }
    catch (Exception e) {
        print("pt=${pt}")
        org.codehaus.groovy.runtime.StackTraceUtils.sanitize(e).printStackTrace()
        e
    }
}

def wrap_angle(a) {
    return (Float.valueOf(a) + 180) % 360 - 180
}

def wsSummarize(_obsid, wsJson, filesJson, tapJson, quality_update, manualAnts) {
    def wsStats = parseJson(wsJson)
    def fileStats = parseJson(filesJson)
    def tapStats = parseJson(tapJson)

    def obs_name = wsStats.obsname
    def groupid = wsStats.groupid
    // def _projectid = wsStats.projectid

    def ra_phase_center = wsStats.ra_phase_center
    if (ra_phase_center == null) {
        ra_phase_center = (wsStats.metadata ?: [:]).ra_pointing ?: Float.NaN
    }
    ra_phase_center = wrap_angle(ra_phase_center)

    def dec_phase_center = wsStats.dec_phase_center
    if (dec_phase_center == null) {
        dec_phase_center = (wsStats.metadata ?: [:]).dec_pointing ?: Float.NaN
    }
    dec_phase_center = wrap_angle(dec_phase_center)

    def az_pointing = (wsStats.metadata ?: [:]).azimuth_pointing ?: Float.NaN
    az_pointing = wrap_angle(az_pointing)
    def el_pointing = (wsStats.metadata ?: [:]).elevation_pointing
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
    if ([-90, 0, 90].contains(az_pointing.round() as int)) {
        ew_pointing = ((90 - el_pointing) / 7).round() as int
        if (az_pointing > 0) {
            ew_pointing *= -1
        }
    }
    else {
        println("unknown azel pointing az=${az_pointing} el=${el_pointing}")
    }

    // mapping from eor1 pointing to sweet pointing
    def eor12sweet = [
        0: 0,
        1: 2,
        2: 4,
        3: 10,
        4: 12,
        5: 26,
        6: 28,
        7: 51,
        8: 54,
        9: 83,
        10: 87
    ]

    def gridpoint_number = (wsStats.metadata ?: [:]).gridpoint_number
    def gridpoint_name = (wsStats.metadata ?: [:]).gridpoint_name
    def sweet_pointing = null
    if (gridpoint_name == "sweet") {
        sweet_pointing = gridpoint_number
    }
    else if (gridpoint_name == "EOR1") {
        if (eor12sweet[gridpoint_number] != null) {
            sweet_pointing = eor12sweet[gridpoint_number]
        }
        else {
            println("unknown EOR1 gridpoint_number ${gridpoint_number}")
        }
    }
    else {
        println("unknown gridpoint_name ${gridpoint_name}")
    }
    def lst = wrap_angle((wsStats.metadata ?: [:]).local_sidereal_time_deg.floatValue())

    def eorfield = null
    def eorband = null
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
    }
    else if (nearest_ra == 60 && nearest_dec == -30) {
        eorfield = 1
    }
    else if (nearest_ra == 155 && nearest_dec == -10) {
        eorfield = 2
    }
    else if (nearest_ra == 15 && nearest_dec == -27) {
        eorfield = 3
    }
    else {
        println("unknown eor field for ${nearest_ra} ${nearest_dec}")
    }

    // eor bands
    // |   band | cent          | range                |
    // | ------ | ------------- | -------------------- |
    // | 0 low  |  120 (154MHz) | 109-132 (139-170MHz) |
    // | 1 high |  142 (182MHz) | 131-154 (167-198MHz) |
    // | 2 ulow |   70  (90MHz) |  59-88   (75-113MHz) |

    def coarse_chans = ((wsStats.rfstreams ?: [:])["0"] ?: [:]).frequencies ?: []
    def center_chan = null
    if (tapStats.center_channel_number != null) {
        center_chan = tapStats.center_channel_number as int
    }
    else if (coarse_chans != null && coarse_chans.size() > 0) {
        center_chan = coarse_chans[Math.max(coarse_chans.size() / 2 as int - 1, 0)]
    }

    if (center_chan == 142) {
        eorband = 1
    }
    else if (center_chan == 120) {
        eorband = 0
    }
    else if (center_chan == 70) {
        eorband = 2
        print("unknown eor band for ${center_chan}")
    }

    def nscans = ((wsStats.stoptime ?: 0) - (wsStats.starttime ?: 0)) / (wsStats.int_time ?: 1)
    def delays = (wsStats.alldelays ?: [:]).values().flatten()
    def quality = wsStats.quality ?: [:]
    def tiles = wsStats.tdict ?: [:]
    def tile_nums = tiles.collect { k, _v -> k as int }
    def tile_names = tiles.collect { _k, v -> v[0] }
    def tile_rxs = tiles.collect { _k, v -> v[1] }
    def n_tiles = tile_names.size()
    // def n_lb = tile_names.count { it =~ /(?i)lb/ }
    // def n_hex = tile_names.count { it =~ /(?i)hex/ }

    def bad_tiles = wsStats.bad_tiles ?: [:]
    def n_bad_tiles = bad_tiles.size()
    def n_good_tiles = n_tiles - n_bad_tiles
    def bad_tile_frac = n_bad_tiles / n_tiles
    def n_dead_dipoles = delays.count { it == 32 }
    def dead_dipole_frac = null
    if (delays.size() != null && delays.size() > 0) {
        dead_dipole_frac = n_dead_dipoles / delays.size()
    }
    def dataquality = Float.valueOf(wsStats.dataquality ?: 0)
    def dataqualitycomment = wsStats.dataqualitycomment ?: ''
    def manual_dataquality = Float.valueOf(quality_update.dataquality ?: 0)
    if (quality_update.dataquality != null && dataquality != manual_dataquality) {
        dataquality = manual_dataquality
        dataqualitycomment = "manual: ${quality_update.dataqualitycomment ?: ''}"
    }
    def faults = wsStats.faults ?: [:]
    def badstates = (faults.badstates ?: [:]).values().flatten()
    def badpointings = (faults.badpointings ?: [:]).values().flatten()
    def badfreqs = (faults.badfreqs ?: [:]).values().flatten()
    def badgains = (faults.badgains ?: [:]).values().flatten()
    def badbeamshape = (faults.badbeamshape ?: [:]).values().flatten()
    def fail_reasons = []
    def capture_mode = wsStats.mode

    def config = tapStats.mwa_array_configuration
    def sun_elevation = Float.valueOf(tapStats.sun_elevation ?: 'NaN')
    def sun_pointing_distance = Float.valueOf(tapStats.sun_pointing_distance ?: 'NaN')

    def bad_ants = bad_tiles.collect { tile_nums.indexOf(it) + 1 }

    // fail codes
    def fail_code = 0
    // no error

    // 0x0X - observational
    if (params.filter_eorfield != null && eorfield != params.filter_eorfield) {
        fail_reasons += [String.format("phase_radec(%+.1f,%+.1f)!=eor%d", ra_phase_center, dec_phase_center, params.filter_eorfield)]
        fail_code = fail_code == 0 ? 1 : fail_code
    }
    if (params.filter_ra != null && (ra_phase_center - params.filter_ra).abs() > 0.1) {
        fail_reasons += [String.format("phase_ra(%+.1f)!=ra(%+.1f)", ra_phase_center, params.filter_ra)]
        fail_code = fail_code == 0 ? 1 : fail_code
    }
    if (params.filter_dec != null && (dec_phase_center - params.filter_dec).abs() > 0.1) {
        fail_reasons += [String.format("phase_dec(%+.1f)!=dec(%+.1f)", dec_phase_center, params.filter_dec)]
        fail_code = fail_code == 0 ? 1 : fail_code
    }
    if (params.filter_eorband != null && eorband != params.filter_eorband) {
        fail_reasons += ["center_chan=${center_chan}"]
        fail_code = fail_code == 0 ? 2 : fail_code
    }
    if (params.filter_sweet_pointings != null && !params.filter_sweet_pointings.contains(sweet_pointing)) {
        fail_reasons += ["sweet_pointing=${sweet_pointing}"]
        fail_code = fail_code == 0 ? 3 : fail_code
    }
    if (params.filter_ew_pointings != null && !params.filter_ew_pointings.contains(ew_pointing)) {
        fail_reasons += ["ew_pointing=${ew_pointing}"]
        fail_code = fail_code == 0 ? 3 : fail_code
    }
    if (params.filter_sun_elevation != null && sun_elevation != null && sun_elevation > params.filter_sun_elevation) {
        fail_reasons += [String.format("sun_elevation(%+.1f)>%+.1f", Float.valueOf(sun_elevation), Float.valueOf(params.filter_sun_elevation))]
        fail_code = fail_code == 0 ? 4 : fail_code
    }
    if (params.filter_min_sun_pointing_distance != null && sun_pointing_distance != null && sun_pointing_distance > params.filter_min_sun_pointing_distance) {
        fail_reasons += [String.format("sun_pointing_distance(%+.1f)>%+.1f", Float.valueOf(sun_pointing_distance), Float.valueOf(params.filter_min_sun_pointing_distance))]
        fail_code = fail_code == 0 ? 4 : fail_code
    }
    if (capture_mode == "NO_CAPTURE") {
        fail_reasons += ["no cap fr fr"]
        fail_code = fail_code == 0 ? 5 : fail_code
    }
    if (params.filter_config != null && config != null) {
        // e.g. config = "Phase II Compact", params.filter_config = ["Phase I"]
        if (!params.filter_config.contains(config)) {
            fail_reasons += [String.format("config(%s) not in %s", config, params.filter_config)]
            fail_code = fail_code == 0 ? 6 : fail_code
        }
    }

    // 0x1X - runtime
    if (params.filter_bad_tile_frac != null && bad_tile_frac > params.filter_bad_tile_frac) {
        fail_reasons += ["bad_tiles(${bad_tiles.size()})=${displayInts(bad_tiles)}"]
        fail_code = fail_code == 0 ? 16 : fail_code
    }
    if (params.filter_dead_dipole_frac != null && dead_dipole_frac > params.filter_dead_dipole_frac) {
        fail_reasons += ["dead_dipole_frac=${dead_dipole_frac}"]
        fail_code = fail_code == 0 ? 17 : fail_code
    }

    // 0x2X - quality
    if (params.filter_quality != null && dataquality > params.filter_quality) {
        fail_reasons += ["dataquality=${dataquality} (${dataqualitycomment})"]
        fail_code = fail_code == 0 ? 32 : fail_code
    }
    if (params.filter_ionoqa && quality.iono_qa != null && quality.iono_qa > params.filter_ionoqa) {
        fail_reasons += [String.format("rts_iono_qa(%.1f)>%.1f", Float.valueOf(quality.iono_qa), Float.valueOf(params.filter_ionoqa))]
        fail_code = fail_code == 0 ? 33 : fail_code
    }
    if (params.filter_ionoqa && quality.iono_qa == null) {
        fail_reasons += ["rts_iono_qa is null"]
        fail_code = fail_code == 0 ? 34 : fail_code
    }
    if (fileStats.num_data_files < 2) {
        fail_reasons += ["no data files"]
        fail_code = fail_code == 0 ? 35 : fail_code
    }

    fail_code = getFailReason(fail_code)

    [
        fail_code: fail_code,
        fail_reasons: fail_reasons,
        obs_name: obs_name,
        groupid: groupid,
        corrmode: wsStats.mode,
        delaymode: wsStats.delaymode_name,
        starttime_mjd: parseFloatOrNaN(tapStats.starttime_mjd),
        starttime_utc: tapStats.starttime_utc,
        sun_elevation: sun_elevation,
        sun_pointing_distance: sun_pointing_distance,
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
        freq_res: wsStats.freq_res,
        coarse_chans: coarse_chans,
        eorband: eorband,
        centre_freq: ((coarse_chans[0] + coarse_chans[-1]) * 1.28E+6 / 2),
        int_time: wsStats.int_time,
        nscans: nscans,
        config: config,
        n_tiles: n_tiles,
        n_good_tiles: n_good_tiles,
        tile_nums: tile_nums,
        tile_rxs: tile_rxs,
        bad_tiles: bad_tiles,
        bad_ants: bad_ants,
        manual_ants: (manualAnts as ArrayList),
        iono_magnitude: quality.iono_magnitude,
        iono_pca: quality.iono_pca,
        iono_qa: quality.iono_qa,
        dataquality: dataquality,
        dataqualitycomment: dataqualitycomment,
        bad_tile_frac: bad_tile_frac,
        n_bad_tiles: n_bad_tiles,
        dead_dipole_frac: dead_dipole_frac,
        n_dead_dipoles: n_dead_dipoles,
        badstates: badstates.size(),
        badpointings: badpointings.size(),
        badfreqs: badfreqs.size(),
        badgains: badgains.size(),
        badbeamshape: badbeamshape.size(),
        fault_str: faults.shortstring.replaceAll(java.util.regex.Pattern.compile('\\s+'), '|'),
        num_data_files: fileStats.num_data_files,
        num_data_files_archived: fileStats.num_data_files_archived
    ]
}

def wscleanParams() {
    [
        suffix: params.img_suffix,
        weight: params.img_weight,
        size: params.img_size,
        scale: params.img_scale,
        channels_out: params.img_channels_out,
        intervals_out: params.img_intervals_out,
        split_intervals: params.img_split_intervals,
        pol: params.img_pol,
        args: params.wsclean_args
    ]
}

def wscleanDConvParams() {
    mapMerge(
        wscleanParams(),
        [
            args: "${params.wsclean_args} ${params.wsclean_dconv_args}",
            niter: params.img_niter,
            minor_clean_gain: params.img_minor_clean_gain,
            major_clean_gain: params.img_major_clean_gain,
            auto_threshold: params.img_auto_threshold,
            auto_mask: params.img_auto_mask,
            mwa_path: params.img_mwa_path
        ]
    )
}

def openWithDelay(f) {
    // open f, read a byte, then wait for a random time between 50 and 500 ms
    f.withReader {
        it.read(0)
        Thread.sleep(50 + new Random().nextInt(450))
    }
}
