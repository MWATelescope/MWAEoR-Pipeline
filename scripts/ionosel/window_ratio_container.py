#!/usr/bin/env python3
"""Run INSIDE the chips_wrappers image. EoR-window power ratios variant/baseline computed from the
same 2D arrays plotchips_all.py plots (ChipsDataProducts.read_data_and_create_2Darray), with the
same horizon wedge line it draws (0.5*pi*wedge_factor*k_perp, wedge_factor from cosmology).

usage: window_ratio_container.py BASEDIR GROUP OUT.tsv
Windows (all: weights>0 in both, finite, 0.1 < |k| < 1 h/Mpc):
  W1 1D-cuts        : 0.025 < k_perp <= 0.045, k_par >= 0.11            (the pipeline's chips1d_tsv window)
  W2 window_kperp01 : k_par > horizon + 0.05 h/Mpc, k_perp <= 0.10
  W3 window_all     : k_par > horizon + 0.05 h/Mpc, any k_perp
Also reports the whole wedge (k_par <= horizon) for completeness.
Statistics per window: weighted-mean power ratio sum(P1*w)/sum(P0*w) (w = baseline weights),
median per-cell ratio, fraction of cells where the variant is lower than baseline.
"""
import sys, numpy as np
from chips_wrappers.setup import make_args
from chips_wrappers.ps_methods.chips_data import ChipsDataProducts

basedir, group, out = sys.argv[1:4]
NAME = "ionosub_ssins_30l_src8k_300it_8s_80kHz_sel-{v}_i1000"
common = (f"plotchips_all.py --basedir {basedir}/ --bias_mode 11 --lowerfreq_orig 138900000.0 --chan_width 80000 "
          f"--umax 300 --N_chan_orig 384 --kperp_max 0.045 --kperp_min 0.025 --kparra_min 0.11 --num_k_edges 80 "
          f"--plot_type 2D --polarisation xx --chips_tag {group}_{NAME.format(v='none')}").split()
sys.argv = common
args = make_args.check_args(make_args.get_args())
cd = ChipsDataProducts(args)

rows = []
for pol in ("xx", "yy"):
    P = {}; W = {}
    for v in ("none", "pred4", "meas4"):
        arr, _ = cd.read_data_and_create_2Darray(pol, f"{group}_{NAME.format(v=v)}")
        P[v] = np.ma.filled(np.ma.masked_invalid(arr), np.nan).astype(float)
        W[v] = np.array(cd.weights, dtype=float)
    cd.get_horizon_and_beam_lines()
    kper, kpa = np.array(cd.kper), np.array(cd.kpa)
    KP, KL = np.meshgrid(kper, kpa)
    K = np.hypot(KP, KL)
    horiz = np.broadcast_to(cd.line_horiz, KP.shape)
    ok = (W["none"] > 0) & np.isfinite(P["none"]) & (P["none"] > 0)
    shell = ok & (K > 0.1) & (K < 1.0)
    wins = {
        "W1_1Dcuts":        shell & (KP > 0.025) & (KP <= 0.045) & (KL >= 0.11),
        "W2_window_kperp<=0.1": shell & (KL > horiz + 0.05) & (KP <= 0.10),
        "W3_window_all":    shell & (KL > horiz + 0.05),
        "wedge":            ok & (KL <= horiz) & (K < 1.0),
    }
    for v in ("pred4", "meas4"):
        for name, m in wins.items():
            mm = m & (W[v] > 0) & np.isfinite(P[v])
            w = W["none"][mm]
            r = np.sum(P[v][mm] * w) / np.sum(P["none"][mm] * w)
            cell = P[v][mm] / P["none"][mm]
            rows.append((pol, v, name, int(mm.sum()), r, float(np.nanmedian(cell)), float(np.mean(cell < 1))))
print(f"wedge_factor={cd.wedge_factor:.4f}  central_freq={cd.central_freq/1e6:.2f} MHz  z={cd.z:.3f}")
with open(out, "w") as fo:
    fo.write("pol\tvariant\twindow\tncell\tratio_wmean\tmedian_cell_ratio\tfrac_cells_lower\n")
    for r in rows:
        line = f"{r[0]}\t{r[1]}\t{r[2]}\t{r[3]}\t{r[4]:.4f}\t{r[5]:.4f}\t{r[6]:.3f}"
        fo.write(line + "\n"); print(line)
