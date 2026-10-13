# Iono-sub source selection (hyperdrive PR 93): power-spectrum comparison, 2026-09-03 EoR0 low band

Run 2026-10-13 on DUG (not Pawsey: this checkout lives on DUG, which has A100 GPUs, so the CUDA
image was used directly and no HIP build was needed). Branch `iono-sub-selection-2026-09-03`
(pushed), commits 6d14a6c and b6f5f32 on top of e8f3998.

- container: `docker://mwatelescope/hyperdrive:0.8.1-bmetrics-cuda12.5.1-ubuntu24.04` pulled to
  `/data/curtin_mwaeor/singularity/hyperdrive_0.8.1-bmetrics_cuda12.5.1-ubuntu24.04.sif`;
  verified `mwa_hyperdrive 0.8.1`, 4 `iono-sub-max` matches in `peel --help`, `bmetrics` present,
  and the peel logs show the A100 (capability 8.0, CUDA 12.5) and the PR 92 convention lines.
  The dug profile's `withLabel: hyperdrive` now uses `params.hyperdrive_sif`.
- variants: `eor0low-ionosel-{none,pred4,meas4}.conf`, copies of `eor0low-ph3.conf` differing only in
  `peel_suffix` and `hyp_peel_args`; `ionosub_nsrcs = 1000`, `sub_nsrcs = 8000`, same sky model
  (`srclist_pumav3_EoR0LoBES_EoR1pietro_CenA-GP_2023-11-07.yaml`), same `hyp_srclist_args`, same
  beam file. Archiving disabled, `bucket_prefix = mwaeor:low0-sel` (mwaeor:low0.peel untouched).
- obsids: `obsids-eor0low20260903.csv` (66), subsets `-quiet` (25, ≤16:59 UTC) and `-event` (23, 17:01–17:49 UTC).
- commands: COMMANDS.md (this directory). Preprocessing and calibration were fully cached from the
  earlier eor0low20261002_1 run; only peel, ionosub grids, combine, lssa, 1D and plots ran.
- CHIPS: lssa_fg_general, cal_suffix _norfi, bias mode 11, kperp 0.025–0.045, kparra_min 0.11.
  Group dirs: full night eor0low_phase3_be56c57b, quiet eor0low_phase3_25c2db54, event eor0low_phase3_ac453249
  under /data/curtin_mwaeor/nfdata/<group>/lssa_fg_general_norfi_b11/ionosub_ssins_30l_src8k_300it_8s_80kHz_sel-<v>_i1000/.

## Did the selection do something? (peel JSONs, 66 per variant, stats script = port of the low0.peel index.html one)

| variant | unfit entries (α=0, β=0, g=1) | sources fitted at ≥1 timestep (median/obs) | median offset of fitted (″ @150 MHz) | per-source temporal RMS (″) |
|---|---|---|---|---|
| A baseline | 1.8% (1.1–2.5%) | 995 | 43.1 | 22.8 |
| B predicted, ratio 4 | 59.4% (56.2–62.4%) | 410 (379–444) | 29.0 | 13.2 |
| C measured, ratio 4 after 3 timeblocks | 40.6% (18.2–56.4%) | 994 at start, ≈370–460 by the end | 33.0 | 21.7 |

B keeps apparent flux ≥ (median of brightest 10%)/4 ≈ 1.2 Jy: ranks 1–300 almost entirely, 62% of
301–600, none of 601–1000. That is ~410 sources, more than the 100–300 the brief expected, so the
ratio-8 fallback was not needed. C drops ~490 sources at timeblock 3 and keeps trimming (unfit by
timestep: 0.02 ×3, then 0.40 0.47 0.50 … 0.54). The baseline reproduces the brief's offset-vs-rank
pattern (13.6″ ranks 1–10 → 64″ ranks 601–1000) and the event (brightest-10 offset ≈8″ at 16:05,
peaking 47″ at 17:47, back to ≈10″ by 18:15). Details and per-obsid table: peel_stats.md.

## Main table: variant x polarisation (full night, 66 obs)

| variant | pol | unfit % (peel JSON) | median offset of fitted sources (″ @150 MHz) | 1D EoR-window P ratio to baseline (sum over 0.1<k<1, k⊥ 0.025–0.045, k∥≥0.11) | 2D EoR-window ratio (weighted mean, above horizon+0.05, k⊥≤0.1) | 2D wedge ratio (weighted mean / median cell) |
|---|---|---|---|---|---|---|
| A baseline (sel-none) | XX | 1.8% | 43.1 | 1 (reference) | 1 (reference) | 1 (reference) |
| A baseline (sel-none) | YY | 1.8% | 43.1 | 1 (reference) | 1 (reference) | 1 (reference) |
| B predicted (sel-pred4) | XX | 59.4% | 29.0 | 0.9992 | 0.9997 | 1.0070 / 1.011 |
| B predicted (sel-pred4) | YY | 59.4% | 29.0 | 0.9990 | 1.0000 | 1.0035 / 1.029 |
| C measured (sel-meas4) | XX | 40.6% | 33.0 | 1.0017 | 0.9996 | 0.9982 / 1.049 |
| C measured (sel-meas4) | YY | 40.6% | 33.0 | 1.0013 | 0.9998 | 1.0019 / 1.024 |

The peel-JSON columns are per variant (identical for XX and YY; the fit is not per-polarisation).

## Quiet (25 obs, ≤16:59 UTC) and event (23 obs, 17:01–17:49 UTC) subsets

| subset | variant | pol | 1D window ratio | 2D window ratio (W2) | 2D wedge ratio (wmean / median cell) |
|---|---|---|---|---|---|
| quiet | B predicted (sel-pred4) | XX | 1.0002 | 0.9999 | 1.0101 / 1.014 |
| quiet | B predicted (sel-pred4) | YY | 0.9988 | 1.0000 | 1.0055 / 1.024 |
| quiet | C measured (sel-meas4) | XX | 1.0046 | 1.0001 | 1.0187 / 1.059 |
| quiet | C measured (sel-meas4) | YY | 1.0029 | 0.9998 | 1.0051 / 1.028 |
| event | B predicted (sel-pred4) | XX | 0.9974 | 0.9996 | 1.0017 / 1.014 |
| event | B predicted (sel-pred4) | YY | 0.9990 | 0.9997 | 1.0055 / 1.036 |
| event | C measured (sel-meas4) | XX | 1.0081 | 0.9997 | 0.9721 / 1.060 |
| event | C measured (sel-meas4) | YY | 0.9992 | 0.9995 | 0.9853 / 1.032 |

Window definitions. 1D: the pipeline's chips1d_tsv bins, ratio of summed P(k) over bins with
0.1 < k < 1 h/Mpc and positive baseline power (16 XX / 11 YY bins; many higher-k bins are negative,
i.e. noise-dominated, see 1D_table.md). 2D: from the same crosspower/outputweights arrays
plotchips_all.py plots, weighted-mean power ratio Σ P_var w / Σ P_base w over cells with k∥ above the
CHIPS horizon line (wedge_factor 3.557 at z = 8.69) + 0.05 h/Mpc, k⊥ ≤ 0.1, 0.1 < |k| < 1 (≈270 cells);
"wedge" is k∥ ≤ horizon, |k| < 1 (≈1250 cells). Per-k tables: */1D_ratios_{xx,yy}.tsv; per-window
rows incl. median cell ratio and fraction of cells lower: */window_ratios_*.tsv.

## 2D plots

full-night/ (and quiet/, event/): chips2D_{xx,yy}_sel-{none,pred4,meas4}.png (crosspower) and
chips2D_ratio_{xx,yy}_sel-{pred4,meas4}_over_sel-none.png (variant / baseline, colour range 0.5–2).
What they show: the EoR window is white (ratio 1 to within the per-cell noise) for both variants in
both polarisations. Both variants are blue (10–30% more power) in the wedge at k⊥ ≈ 0.07–0.15 h/Mpc,
k∥ < 0.1. C additionally shows a blue band of excess power hugging the horizon line from
k⊥ ≈ 0.05 up to k⊥ ≈ 0.3 (k∥ ≈ 0.2–1), i.e. power leaking along the wedge edge; B does not.

## Conclusion

Restricting the ionospheric fit did not reduce foreground residual power. In the EoR window the
variant/baseline ratios are 0.9995–1.0001 (2D weighted mean) and within ±0.3% (1D sums) for both
variants, both polarisations, the full night and the quiet and event subsets: no detectable change
at this depth. In the wedge the selection makes things slightly worse: 0.3–1.9% more
weighted-mean power and 1–6% higher median cell ratio, concentrated at k⊥ ≈ 0.1, and variant C
adds excess power along the horizon line. So, at least for this night and ratio 4, the fitted
(α, β) of the 300–1000th sources are doing net useful work in the subtraction even though their
scatter is large; not shifting them leaves more foreground in the wedge and gains nothing above it.
Of the two, B (predicted, up-front) is the better behaved: no horizon-line excess, smaller wedge
penalty, half the peel runtime (15 vs 31 min per obsid on an A100). Neither is an improvement
over the baseline for the power spectrum. If this is pursued, the next things to try are a looser
ratio (8–16, which keeps more of ranks 300–600 where the wedge penalty originates) or the opposite
test (fit only, do not subtract unfitted sources) to separate "noisy offset" from "faint model".

## Flags

- Not on Pawsey: DUG, `-profile dug`. The brief's Garrawarla/Setonix container instructions were
  replaced by the dug-profile change above. No HIP image was built or needed.
- The local reference peel for these obsids in /data/curtin_mwaeor/nfdata was made with hyperdrive
  0.7.0 (the old dug sif), whereas the Acacia low0.peel set is 0.8.1. The baseline here (A) was
  re-run with the new 0.8.1 image, so A/B/C differ only in the selection flags.
- Each pipeline run reports 2 failed tasks: `chipsPlot … 2D_ratio.xx+yy_11`, the pipeline's
  pre-existing nosub-vs-ionosub combined-polarisation ratio plot. Error is ignored by the pipeline,
  unrelated to this test, and present identically in all 9 runs (3 variants × full/quiet/event).
- The three variants (and later the six subset runs) were launched concurrently, so the shared
  unsubtracted (`ssins_30l_src8k_300it_8s_80kHz`) combine/lssa was computed by each run rather than
  cached across them; harmless, just wasted CPU.
- `--nopsmetrics --nodelayspec --noplotvisqa --nothumbnail --nopolcomp --novideo` were set to run
  only what was needed; `nopowerspec`/`nochips` were not set.
- Kept as-is but worth questioning separately: the EoR0 low conf uses the 167–197 MHz interpolated
  beam file (`MWA_embedded_element_pattern_rev2_interp_167_197MHz.h5`) for 138.9–169.5 MHz data.
- Peel wall time: 31 min/obsid baseline, ~15 min/obsid for B and C (A100, 16 cpus, 374 GB).
