# AGENTS-DUG.md

Operational notes for running this pipeline on DUG (`-profile dug`, Slurm partition
`curtin_mwaeor`). Local to this checkout (listed in `.gitignore`); the general notes are in
[AGENTS.md](AGENTS.md). Site details that are not public (proxy address, hostnames, node inventory)
are deliberately left out: read them from `nextflow.config` or the environment.

## Paths

| What | Where |
|---|---|
| per-obsid products (`params.outdir`, Nextflow `storeDir`) | `/data/curtin_mwaeor/nfdata/<obsid>/` |
| CHIPS group products | `/data/curtin_mwaeor/nfdata/<group>_<chunkhash>/` |
| Nextflow work dir | `/data/curtin_mwaeor/nfwork/<ab>/<hash>/` |
| aggregate results (`params.resultsdir`) | `/data/curtin_mwaeor/nfresults/results<obsids_suffix><result_suffix>[<img_suffix>]<cal_suffix>/` |
| singularity images | `/data/curtin_mwaeor/singularity/` |
| other software (giant-squid, jq, conda env, old sifs) | `/data/curtin_mwaeor/sw/` |
| beam and source lists | `/data/curtin_mwaeor/mwa_full_embedded_element_pattern.h5`, `/data/curtin_mwaeor/srclist_*.yaml` |
| session scratch (login node only, compute nodes cannot see it) | `/localData/...` |

Compute nodes do not mount `/localData`, so anything a Slurm job must read or write (logs, scripts,
sbatch `-o`) goes under `/data`.

## Launching

```bash
nextflow -log /data/curtin_mwaeor/nfresults/<results-dir>/logs/nextflow/<name>.nextflow.log \
  run main.nf -profile dug -with-tower -c eor0low-ph3.conf \
  --obsids_suffix=-eor0low20260903 --nothumbnail --nopolcomp --novideo -name <name>
```

Run inside `tmux` (`nohup ... &` also works). Java, Nextflow, Singularity and rclone come from
environment modules and are normally already loaded in a login shell. `-with-tower` needs
`TOWER_ACCESS_TOKEN` and `TOWER_WORKSPACE_ID` in the environment; drop it if they are not set.
Every non-interactive bash prints a harmless "No such file or directory" for the module init
script; filter it out rather than chasing it.

## Network

Outbound traffic needs the site proxy. Processes that download or upload prepend
`params.proxy_prelude` (see the `dug` profile). On the login node export the same
`http_proxy`/`https_proxy` values yourself before `curl`, `git push`, `rclone` or
`singularity pull`. Do not push unless asked.

## Compute

- GPU work (`hypCalSol`, `hypSubUV`, `hypIonoSubUV` and the MS variants) goes to the A100 nodes via
  `params.gpu_constraint`. There is no HIP build; the CUDA images are the right ones.
- CPU work uses `params.cpu_constraint`. Leave it alone: zen (rome) nodes give illegal-instruction
  errors and KNL nodes are far too slow.
- Everything you launch shares one account queue. The usual pending reason, "Nodes required for job
  are DOWN, DRAINED or reserved", just means waiting; `QOSGrpJobsLimit` means the account cap.
  Typical step times for ETAs: a 2 minute obs di-cal about 11 min on an A100, an iono-sub peel
  about 17 min. Estimate from queue position plus those, and give ETAs in status reports.
- Jobs are named `nf-<process>.<tag>`, so `squeue -u $USER | grep <obsid>` finds an obs. Another
  run of yours may be working on the same obsid (e.g. a different peel variant); check the job's
  `WorkDir` with `scontrol show job <id>` before touching it.
- Processes run in `scratch = $TMPDIR_SHM` (RAM disk), or `$TMPDIR` for the `nvme` label, and
  outputs are **copied** to `storeDir`. The work dir keeps only `.command.*` and `.exitcode`, so
  look for products in `nfdata`, not `nfwork`. The filesystem is slow under load, so job time
  limits are scaled by `scratchFactor` in the confs.
- Task hash `[b4/e70585]` in the progress output means `/data/curtin_mwaeor/nfwork/b4/e70585*/`.
- rclone from the login node is slow (~20 MB/s); bulk Acacia pulls belong in a Slurm job (array
  over obsids, logs under `/data`).

## Containers

- Singularity only, no Docker. The profile binds `/data` into containers; when you run an image by
  hand on the login node add `--bind /data` yourself, or files will "not exist".
- Default process image is the mwa-demo CUDA image; `mwa_qa` tasks use the mwa_qa image; GPU
  hyperdrive tasks use `params.hyperdrive_sif`, set inside the `dug` profile; CPU hyperdrive tasks
  (`solutions-apply`, `solutions-plot`) use `docker://mwatelescope/hyperdrive:autos-dev` from the
  cache dir.
- Container selection is resolved when the config is parsed. A `-c` conf overrides ordinary params
  at run time (its `beam_path` wins over the profile's), but changing which image a process uses
  has to be done in the profile block, or in a separate git worktree with its own profile edit.
- New images: `singularity pull` (with the proxy) into `/data/curtin_mwaeor/singularity/`, keep the
  old file, and point the profile at the new name. Check which hyperdrive wrote a solution with
  the `SOFTWARE` keyword in the soln FITS primary header.

## Misc

- `scripts/combine_dug.sh` is the manual CHIPS combine (sbatch, large memory) for when
  `chipsCombine` must be driven by hand.
- ASVO conversion uses `/data/curtin_mwaeor/sw/bin/giant-squid`; the TAP query scripts run in the
  conda env `/data/curtin_mwaeor/sw/conda/dev` and need the proxy.
