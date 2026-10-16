# AGENTS.md

Notes for an agent (or a new human) operating this checkout. The README describes what the
pipeline does; this file is the non-obvious operational knowledge that is not derivable from the
code. Keep it current: when you learn something the hard way, add it here, and delete anything
that stops being true.

## 1. Site notes

Site-specific notes (paths, network, queues, containers) live in one file per cluster. For DUG
(`-profile dug`, this checkout) read [AGENTS-DUG.md](AGENTS-DUG.md); it is local and not committed.
Pawsey (Setonix, Garrawarla) profiles remain in `nextflow.config` but are not maintained.

## 2. Running the pipeline

Standard launch, from the checkout, in `tmux`. The exact command for your site is in the site
file; the shape is:

```bash
nextflow run main.nf -profile <profile> -c <field>.conf --obsids_suffix=-<name> -name <run name>
```

- `-with-tower` needs `TOWER_ACCESS_TOKEN` and `TOWER_WORKSPACE_ID` in the environment; drop it if unsure.
- `--obsids_suffix=-foo` reads `obsids-foo.csv` in the checkout: the suffix is inserted before the
  extension of `params.obsids_path`. One obsid per line, optional second column, `#` comments.
- The results directory is `results<obsids_suffix>[.<visName>]<result_suffix>`, and most per-stage
  tsvs land in a sibling with `<img_suffix><cal_suffix>` appended (so `_norfi` for the ph3 confs).
  The param is `result_suffix`, not `results_suffix`.
- Live confs: `eor0high-ph3.conf`, `eor0low-ph3.conf`, `eor1high-ph3.conf`, `cenA-ph3.conf`. The
  other `*.conf` files are earlier experiments kept for reference. Each conf's first comment line is
  its launch command.
- Several runs may share the checkout, `nfwork` and `nfdata` at once; give each its own `-name` and
  log file. The script is parsed at launch, so editing `main.nf` after a run has started does not
  affect it.

### Entry points (`-entry`)

| entry | starts from | notes |
|---|---|---|
| (default) | obsid list | web-service metadata, ASVO preprocessed download, flag/SSINS QA, then `qaPrep` (cal, vis, img, chips) |
| `asvoRawFlow` | obsid list | downloads raw gpubox files and runs Birli locally (`birli_args`); used for Phase III nights |
| `extRaw` | raw files already in `nfdata/<obsid>/raw/` | |
| `extPrep` | prep uvfits already in `nfdata/<obsid>/prep/` | reads `birli_<obsid><prep_suffix>.<visName>.uvfits` with `visName` defaulting to `ssins`; use after pulling prep from Acacia; flagQA then measures occupancy on the SSINS-flagged file |
| `extChips` | existing calibrated/peeled uvfits | `visName` names the vis |
| `extPSMetrics` | `external.csv` of obsid,vis | |
| `archivePrep` | obsid list | preprocess and push prep to Acacia only |

### `storeDir` rules (the most common source of confusion)

- A task is skipped only when **every** declared output exists in `storeDir`, including
  `optional: true` logs. Products copied in by hand (e.g. pulled from Acacia) need placeholder
  files beside them: an empty `hyp_di-cal_<obsid>_<name>.log`, `hyp_apply_<name>.log`,
  `hyp_vis-<name>_uv.log`, and `{}` in `hyp_peel_<obsid>_<name>_uv.json`. Otherwise the step reruns
  and silently overwrites the pulled files. The `.nextflow.log` says "is unable to find [UnixPath]".
- Outputs are keyed by **name**, not by parameters. If you change a parameter without changing
  the name that reaches the filename (`prep_suffix`, `cal_suffix`, `dical_name`, `apply_*`,
  `peel_suffix`, `img_suffix`, `ionosub_nsrcs`), the old product is reused. Example: every EoR0
  low product made before 2026-10-14 was calibrated with the 167-197 MHz interpolated beam (wrong
  for 139-170 MHz); fixing `beam_path` alone would have reused them.
- Product names, for grepping `nfdata`:
  - prep: `prep/birli_<obsid><prep_suffix>.uvfits` and `.ssins.uvfits` (SSINS-flagged copy)
  - solution: `cal<cal_suffix>/hyp_soln_<obsid>_ssins_<dical_name>.fits` (+ `hyp_di-cal_*.log`)
  - calibrated vis: `cal<cal_suffix>/hyp_<obsid>_ssins_<dical_name>_<apply_time>s_<apply_freq>kHz.uvfits`
  - iono-subtracted vis: `..._ionosub_ssins_<dical_name>_8s_80kHz<peel_suffix>_i<ionosub_nsrcs>.uvfits` + `hyp_peel_..._uv.json`
  - QA: `cal_qa<cal_suffix>/`, `vis_qa<cal_suffix>/`, `img_qa<img_suffix><cal_suffix>/`, `ps_metrics<cal_suffix>/`
  - CHIPS: `nfdata/<eorNband>_<config>[_<date>][_ewp][_lst][_<groupSuffix>]_<chunkhash>/{ps_metrics,lssa_fg_general}<cal_suffix>[_b<bias>]/<vis name>/`
    with a `manifest_<group>.csv` listing the obsids in the chunk.

### Running a comparison / variant

Copy the field conf, keep everything identical except the names that separate the products, and
turn archiving off so the reference buckets stay clean:

```groovy
dical_name = "30l_src8k_300it_40s"   // or peel_suffix = '_sel-meas4', or result_suffix = '-fnxa'
archive = false
bucket_prefix = "mwaeor:low0-<variant>"
```

Put the launch command and the list of intended differences in the conf header. Past examples:
`eor0low-solint40s.conf` (di-cal solution interval), the removed `eor0low-ionosel-*.conf` at
commit `e6d8af0` (peel args), the Fornax A model test (`dical_name` + `result_suffix`).

### Command line gotchas

- A number in scientific notation on the command line stays a String
  (`--filter_max_rms_convg=2e-5` fails with "Cannot compare java.lang.Float ... String"). Put such
  values in a conf file.
- `--noimg` sets the flag true; use `--noimg=false` to turn an option back on.
- `-resume` is rarely needed: `storeDir` already makes completed steps skip.

### Gates and fail codes

Each stage can reject an obs; the first failing code wins and is written to
`fail_counts_all.tsv`, `reasons_*.tsv` and the per-stage metrics tsvs in the results dir. Codes are
grouped by stage (`modules/utils.nf: getFailReason`):

| range | stage | params |
|---|---|---|
| 001-006 | metadata (`wsGate`) | `filter_eorfield/eorband/ra/dec`, `filter_sweet_pointings`, `filter_ew_pointings`, `filter_sun_elevation`, `filter_min_sun_pointing_distance`, `filter_config` |
| 016-017 | tiles / dipoles | `filter_bad_tile_frac`, `filter_min_good_tiles`, `filter_dead_dipole_frac` |
| 032-035 | data quality, iono QA, no files | `filter_quality`, `quality-updates.csv`, `filter_ionoqa` |
| 049-055 | flag occupancy, RFI, SSINS | `flag_occupancy_threshold`, `rfi_occupancy_threshold`, `ssins_{occupancy,streak,narrow,dab}_threshold` |
| 063 | prepVisQA bad tiles | `filter_bad_tile_frac`, `noprepqafilter` |
| 064-079 | calQA | `filter_max_rms_convg`, `filter_max_unused_bl_frac`, `filter_max_unconvg_chs`, `filter_max_cal_amp_rms` |
| 096-103 | ps_metrics window/wedge | `filter_{max,min}_ps_window_unsub`, `filter_{max,min}_ps_ratio_unsub`, `filter_{max,min}_ps_{win,wedge}_ratio_sub` |
| 112-121 | imgQA | `filter_{max,min}_vrms_box_nosub`, `filter_*_pks_int_*` |

A gate tuned for one configuration can block a variant: 40 s calibration timeblocks hold a third of
the data, so their RMS convergence sits above the 1e-5 `filter_max_rms_convg` of the ph3 confs
(1.06e-5 vs 9.4e-6 for 2 min). Loosen it in the variant conf rather than silently dropping obs.

### Debugging a run

- Progress lines show the task hash: `[b4/e70585] qaPrep:cal:hypCalSol (1472486760)` means
  `<workDir>/b4/e70585*/` holds `.command.sh` (the rendered script), `.command.log`
  and `.exitcode`.
- `nextflow log <run name> -f workdir,exit,status,process,tag,duration,realtime,peak_vmem`
  lists tasks; `nextflow log -q` lists run names; `.nextflow/history` has every launch command.
- Interesting lines in the log file given by `-log`:
  `grep -E '(COMPLETED; exit: [^0]|WARN|ERROR|unable to find)' <log>`.
- Retries: GPU and prep tasks retry on exit 137-140 (OOM/timeout) with time scaled by
  `task.attempt`; most other failures are `ignore`d so the obs just drops out. An obs missing from a
  later stage usually failed a gate, not a task: check `reasons_*.tsv` first.
- `storeDir` outputs are copied after the task; a task can be `COMPLETED` while the copy is still
  in flight on a slow filesystem. Wait for the file size to settle before reading it.
- `maxForks` (and the `rate_limit_*` labels) cap concurrency **per process**, not per run: a run
  with peel, di-cal and gridding all capped at 2 can still hold 6 GPU/CPU jobs. Nor do they add up
  across runs or users (curtin_ericj writes to the same `nfdata`; `sacct -a` only shows your own jobs,
  so count concurrency from the hyperdrive logs in `nfdata` instead).
- Hyperdrive cost is set by the beam file, not by I/O: di-cal with the 2 GB
  `MWA_embedded_element_pattern_rev2_interp_167_197MHz.h5` spends 35-55 min in "reading input data
  and sky modelling" against ~9 min with `mwa_full_embedded_element_pattern.h5` (same 8000 sources,
  60 timesteps, 768 channels); the peel loop is ~85-97 vs ~27 min. Reading the input uvfits from /data
  takes under a minute. Only the interpolated beam slows with concurrency (above ~24 hyperdrive jobs),
  hence `maxForks = 12` for `hypCalSol` and `hypIonoSubUV` in the dug profile (2026-10-16).
- `prepVisQA` crashes (plotting) on data whose autocorrelations are zero, e.g. simulation-only
  visibilities; use `--noprepqa` for those.
- `-entry extChips` passes the wrong field to `gridvisdiff` for EoR1 data.
- Group imaging (`imgCombine`, label `nvme_full`) converts every obs in the group to MS in scratch;
  on DUG the default RAM-disk scratch fills ("Error writing to temporary data file"). Run it with a
  conf that sets `process { withLabel: 'nvme_full' { scratch = false } }`.
- With `--noimg=false`, CHIPS and image combining only receive obs that pass imgQA, and the imgQA
  filter only passes EoR0 (`meta.eorfield == 0`). On EoR1 add `--noimgqa` or both stages get nothing.
- hyperdrive applies a solution file as V_cal = J·V·Jᴴ; DP3 DDECal gains G calibrate as G⁻¹VG⁻ᴴ, so
  converting DDECal solutions to hyperdrive format means writing J = G⁻¹.

## 3. Data flow in and out

- **ASVO**: `asvoPrep` uses `params.giant_squid` with `MWA_ASVO_API_KEY`
  from the environment, submits a conversion job and exits 75 (ignored) until the job is ready.
  Rerun the pipeline later to pick the files up. `prep_time_res_s`/`prep_freq_res_khz` null means
  native resolution (ph3).
- **Acacia** (Pawsey object store) via the rclone remote `mwaeor:`; buckets are
  `${bucket_prefix}.{prep,soln,uvfits,peel,img,imgqa,visqa,grids,...}` with `bucket_prefix`
  `mwaeor:high0`, `mwaeor:low0`, `mwaeor:high1`. `archive = false` by default; `pullPrep` and
  `pullCalSol` fetch from the buckets inside tasks. Not everything is archived (e.g. `high1.prep`
  has only the `.ssins.uvfits`, EoR1-high has no peel bucket). From the login node rclone does
  slow; bulk pulls belong in a batch job (array over obsids).
- **MWA TAP / web services**: `scripts/query_eor{0high,0low,1high}.py` list recent obsids
  (need `pyvo` and network access); `wsMeta`, `tapMeta`
  and `wsMetafits` fetch metadata per obs into `nfdata/<obsid>/meta/` and `raw/`.
- **Manual flags**: `tile-updates.csv` (`startgps,endgps,tiles,comment`, tile indices as in the
  metafits TILE_DATA HDU, pipe separated) and `quality-updates.csv`
  (`obsid,dataquality,comment`) are read every run.

## 4. Field specifics

- EoR0 low band is 139-170 MHz (coarse channels 109-130); EoR0/EoR1 high is 167-197 MHz
  (131-154). `MWA_embedded_element_pattern_rev2_interp_167_197MHz.h5` only covers the high band:
  use `mwa_full_embedded_element_pattern.h5` for low.
- Phase III confs (`*-ph3.conf`) use `cal_suffix = "_norfi"` (no AOFlagger RFI flagging, SSINS
  instead), `dical_name = "30l_src8k_300it"`, apply to 8 s / 80 kHz, iono-sub of 1000 sources,
  `lssa_bias_mode = 11`, `nosub = true` (plain model subtraction off).
- `scripts/combine_dug.sh` is the manual CHIPS combine (sbatch, large memory) for when
  `chipsCombine` must be driven by hand; `scripts/rfistats.sh` summarises nightly RFI metrics.

## 5. Working style expected here

- Report outcomes with evidence paths (results dir, job id, product filename) and an ETA for
  anything still queued.
- Do not cancel, delete or overwrite other runs' products; `squeue -u $USER` shows all of them.
- Save non-obvious discoveries here (and in memory), not in the README.
