#!/usr/bin/env python3
"""Per-file ionospheric statistics for hyperdrive peel JSONs (local paths).

Adapted from the iono_stats.py embedded in https://projects.pawsey.org.au/low0.peel/index.html
(section "HOW THE IONOSPHERIC STATISTICS WORK"); same definitions, but reads local files and
additionally reports per-timestep unfit fractions and the unfit fraction by source rank.

usage: iono_stats_local.py OUT.jsonl FILE.json [FILE.json ...]
"""
import json, sys
import numpy as np
EXCLUDE_FIRST = 0
EXCLUDE_LAST = 1
REF_LAMBDA_M = 2.0          # 150 MHz
RAD_M2_TO_ARCSEC = REF_LAMBDA_M**2 * 206264.8

def one(path):
    with open(path) as fh:
        d = json.load(fh)
    names = list(d.keys())
    nts = max((len(s.get("alphas", [])) for s in d.values()), default=0)
    a = np.full((len(names), nts), np.nan); b = np.full((len(names), nts), np.nan)
    g = np.full((len(names), nts), np.nan)
    for i, s in enumerate(d.values()):
        al = np.array([x if x is not None else np.nan for x in s.get("alphas", [])], dtype=float)
        be = np.array([x if x is not None else np.nan for x in s.get("betas", [])], dtype=float)
        ga = np.array([x if x is not None else np.nan for x in s.get("gains", [])], dtype=float)
        a[i, :len(al)] = al; b[i, :len(be)] = be; g[i, :len(ga)] = ga
    mag = np.hypot(a, b)
    nan_mask = ~np.isfinite(mag)
    default_mask = (a == 0) & (b == 0) & (g == 1)
    bad = nan_mask | default_mask; mag[bad] = np.nan
    with np.errstate(all="ignore"):
        ts_med = np.nanmedian(mag, axis=0) if mag.size else np.array([])
        ts_p90 = np.nanpercentile(mag, 90, axis=0) if mag.size else np.array([])
    ts_default_frac = default_mask.mean(axis=0) if mag.size else np.array([])
    lo = min(EXCLUDE_FIRST, nts); hi = max(lo, nts - EXCLUDE_LAST)
    magt = mag[:, lo:hi]; at = a[:, lo:hi]; bt = b[:, lo:hi]
    nan_t = nan_mask[:, lo:hi]; def_t = default_mask[:, lo:hi]; bad_t = bad[:, lo:hi]
    ts_med_used = ts_med[lo:hi]
    with np.errstate(all="ignore"):
        a_ok = np.where(bad_t, np.nan, at); b_ok = np.where(bad_t, np.nan, bt)
        nfin = np.sum(~bad_t, axis=1)
        src_rms = np.sqrt(np.nanvar(a_ok, axis=1) + np.nanvar(b_ok, axis=1)) if nts else np.array([])
        src_rms = np.where(nfin >= 2, src_rms, np.nan)
        fin_rms = src_rms[np.isfinite(src_rms)]
        tvar_med = np.median(fin_rms) if fin_rms.size else np.nan
        tvar_p90 = np.percentile(fin_rms, 90) if fin_rms.size else np.nan
        src_mean_mag = np.nanmean(magt, axis=1) if nts else np.array([])
        # per-source: fraction of (trimmed) timesteps unfit; sources are in sky-model (brightness) order
        src_def_frac = def_t.mean(axis=1) if nts else np.array([])
        # unfit fraction by brightness-rank decile of the source list
        n = len(names); rank_bins = {}
        for lo_r, hi_r in [(0,10),(10,100),(100,300),(300,600),(600,1000)]:
            hi_r = min(hi_r, n)
            if lo_r < hi_r:
                rank_bins[f"{lo_r}-{hi_r}"] = float(f"{def_t[lo_r:hi_r].mean():.4f}")
        # median offset of fitted sources by rank bin (arcsec at 150 MHz)
        rank_med = {}
        for lo_r, hi_r in [(0,10),(10,100),(100,300),(300,600),(600,1000)]:
            hi_r = min(hi_r, n)
            if lo_r < hi_r:
                v = magt[lo_r:hi_r]; v = v[np.isfinite(v)]
                rank_med[f"{lo_r}-{hi_r}"] = float(f"{np.median(v)*RAD_M2_TO_ARCSEC:.3g}") if v.size else None
    f = lambda x: None if (x is None or not np.isfinite(x)) else float(f"{x:.4g}")
    finite_ts = ts_med_used[np.isfinite(ts_med_used)]
    nfit_src = int(np.sum(src_def_frac < 1.0)) if nts else 0   # sources fitted at >=1 timestep
    return {"key": path.split("/")[-1], "path": path, "nsrc": n, "nts": int(nts),
            "ts_excl_first": lo, "ts_excl_last": int(nts - hi),
            "nan_frac": float(f"{nan_t.mean():.4f}") if mag.size else 1.0,
            "default_frac": float(f"{def_t.mean():.4f}") if mag.size else 0.0,
            "nfit_src": nfit_src,
            "ts_default_frac": [float(f"{x:.4f}") for x in ts_default_frac],
            "default_frac_by_rank": rank_bins,
            "med_arcsec_by_rank": rank_med,
            "ts_med": [f(x) for x in ts_med], "ts_p90": [f(x) for x in ts_p90],
            "med": f(np.median(finite_ts)) if finite_ts.size else None,
            "med_arcsec": f(np.median(finite_ts) * RAD_M2_TO_ARCSEC) if finite_ts.size else None,
            "max": f(finite_ts.max()) if finite_ts.size else None,
            "tvar_med": f(tvar_med), "tvar_p90": f(tvar_p90)}

if __name__ == "__main__":
    out = sys.argv[1]
    with open(out, "w") as fo:
        for p in sys.argv[2:]:
            try:
                r = one(p)
            except Exception as e:
                r = {"key": p.split("/")[-1], "path": p, "error": str(e)}
            fo.write(json.dumps(r) + "\n"); fo.flush()
            print(r["key"], r.get("error") or f"nsrc={r['nsrc']} nts={r['nts']} unfit={r['default_frac']:.3f} med={r['med_arcsec']}\"")
