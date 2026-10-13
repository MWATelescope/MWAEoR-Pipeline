#!/usr/bin/env python3
"""Compare CHIPS 1D power spectra (chips1D_tsv.py output) between peel variants.

usage: compare_1d.py OUT_PREFIX baseline=PATH.tsv variant1=PATH.tsv [variant2=PATH.tsv ...]
The first pair is the baseline. Writes OUT_PREFIX.tsv (per-k table) and prints the
EoR-window (0.1 < k < 1 h/Mpc) summary: ratio of window-summed P(k), and median per-bin ratio.
"""
import sys, math
import numpy as np

def read_tsv(path):
    k, p, n = [], [], []
    with open(path) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        ik, ip, inoise = hdr.index("k_modes"), hdr.index("power"), hdr.index("noise")
        for line in fh:
            cols = line.rstrip("\n").split("\t")
            if len(cols) <= max(ik, ip): continue
            try:
                kk = float(cols[ik]); pp = float(cols[ip]) if cols[ip].strip() else math.nan
                nn = float(cols[inoise]) if len(cols) > inoise and cols[inoise].strip() else math.nan
            except ValueError:
                continue
            k.append(kk); p.append(pp); n.append(nn)
    return np.array(k), np.array(p), np.array(n)

out = sys.argv[1]
pairs = [a.split("=", 1) for a in sys.argv[2:]]
names = [p[0] for p in pairs]
data = {nm: read_tsv(path) for nm, path in pairs}
k0, p0, n0 = data[names[0]]
for nm in names[1:]:
    assert np.allclose(data[nm][0], k0), f"k bins differ for {nm}"
win = (k0 > 0.1) & (k0 < 1.0) & np.isfinite(p0) & (p0 > 0)
with open(out + ".tsv", "w") as fo:
    fo.write("k_modes\t" + "\t".join(f"P_{nm}" for nm in names) + "\t" + "\t".join(f"ratio_{nm}" for nm in names[1:]) + "\n")
    for i in range(len(k0)):
        row = [f"{k0[i]:.4f}"] + [f"{data[nm][1][i]:.6g}" for nm in names]
        row += [f"{data[nm][1][i]/p0[i]:.4f}" if (np.isfinite(p0[i]) and p0[i] != 0) else "nan" for nm in names[1:]]
        fo.write("\t".join(row) + "\n")
print(f"# {out}: {win.sum()} k bins in EoR window 0.1<k<1 (of {len(k0)})")
print("variant\tsum_ratio_window\tmedian_bin_ratio_window\tmin_bin_ratio\tmax_bin_ratio")
for nm in names[1:]:
    p = data[nm][1]
    r = p[win] / p0[win]
    print(f"{nm}\t{p[win].sum()/p0[win].sum():.4f}\t{np.nanmedian(r):.4f}\t{np.nanmin(r):.4f}\t{np.nanmax(r):.4f}")
