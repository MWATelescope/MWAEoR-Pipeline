# MWA EoR Nextflow Pipeline

## Flow

<!-- docker run --rm -v ${PWD}:/home/node/data minlag/mermaid-cli:latest -i /home/node/data/flow.mmd -o /home/node/data/flow.png --scale 3 --backgroundColor transparent --width 1000 --height 1800 -->

```mermaid
%%{init: { 'theme': 'base', 'themeVariables': {
  'darkMode': false,
  'background': '#00ffffff',
  'primaryColor': '#ff0000',
  'primaryBorderColor': '#002b36',
  'primaryTextColor': '#002b36',
  'secondaryColor': '#eee8d5',
  'secondaryBorderColor': '#002b36',
  'secondaryTextColor': '#002b36',
  'tertiaryColor': '#fdf6e3',
  'tertiaryBorderColor': '#002b36',
  'tertiaryTextColor': '#002b36',
  'lineColor': '#002b36',
  'textColor': '#002b36'
}  } }%%

%% ff0000
graph TD;
  classDef in fill:#2aa198;
  classDef out fill:#d33682;
  classDef file fill:#268bd2;
  classDef proc fill:#b58900;
  classDef decision fill:#cb4b16;
  classDef incompleteDecision fill:#EEC3B1;

    subgraph meta
      obsids>fa:fa-file obsids.csv ]; class obsids in
      wsMeta[[fa:fa-download wsMeta]]; class wsMeta proc
      ws([fa:fa-server MWA web services ]); class ws in
      wsmetaJson[/fa:fa-file wsmeta.json /]; class wsmetaJson file
      %% wsfilesJson[/fa:fa-file files.json /]; class wsfilesJson file
      metafits[/fa:fa-file metafits /]; class metafits file

      obsids --> wsMeta --> ws --> wsmetaJson & metafits
      wsGate{"fa:fa-check wsGate &nbsp;&nbsp;"}; class wsGate decision;
      nfConfigMeta>fa:fa-file-code nextflow.config ]; class nfConfigMeta in
      qualityUpdates>fa:fa-file-code quality-updates.csv ]; class qualityUpdates in
      qualityUpdates --"quality updates"-----> wsGate
      wsmetaJson --"pointing/tile/dipole/quality info"--> wsGate
      nfConfigMeta --"pointing/tile/dipole/quality filters"-----> wsGate
    end

    subgraph prep
      asvoPrep[[fa:fa-bolt asvoPrep]]; class asvoPrep proc
      asvo([fa:fa-server MWA ASVO ]); class asvo in
      prepUVFits[/fa:fa-file prep uvfits /]; class prepUVFits file

      wsGate --"passing obs, metafits"--> asvoPrep --> asvo --> prepUVFits

      flagQA[["fa:fa-gem flagQA &nbsp;"]]; class flagQA proc;
      occupancyJson[/fa:fa-file-code occupancy.json /]; class occupancyJson file
      prepUVFits --> flagQA --> occupancyJson;
      %% metafits -.-> flagQA;

      prepVisQA[["fa:fa-gem prepVisQA &nbsp;"]]; class prepVisQA proc;
      prepVisJson[/fa:fa-file-code prepVis.json /]; class prepVisJson file
      prepUVFits --> prepVisQA --> prepVisJson;

      flagGate{"fa:fa-flag flagGate &nbsp;&nbsp;"}; class flagGate decision;
      occupancyJson --> flagGate

      nfConfigPrep>fa:fa-file-code nextflow.config ]; class nfConfigPrep in
      nfConfigPrep --"occupancy_threshold"----> flagGate;

      tileUpdates>fa:fa-file-code tile-updates.csv ]; class tileUpdates in
      tile_flags[/fa:fa-file-code tile_flags.csv /]; class tile_flags file
      tileUpdates --"manual flags"----> tile_flags
      prepVisJson --"prep flags"--> tile_flags

    end


    subgraph cal
      %% goodUVFits[/fa:fa-file good uvfits /]; class goodUVFits file
      calSol[/fa:fa-file-excel cal solutions /]; class calSol file
      hypCalSol[["fa:fa-wrench hypCalSol &nbsp;"]]; class hypCalSol proc
      %% dicalLog[/fa:fa-file dical log /]; class dicalLog file

      flagGate --"passing obs, metafits, uvfits, flags"--> hypCalSol --> calSol
      %% hypCalSol --> dicalLog
      %% dicalArgs --> hypCalSol
      %% metafits -...-> hypCalSol

      polyCal[/fa:fa-file-excel poly solutions /]; class polyCal file
      polyFit[["fa:fa-chart-line polyFit &nbsp;"]]; class polyFit proc
      calSol --> polyFit --> polyCal

      nfConfigCal>fa:fa-file-code nextflow.config ]; class nfConfigCal in
      hypCalSol --"dical args"--- nfConfigCal

      %% end
      %% subgraph cal_qa

      calQA[["fa:fa-gem calQA &nbsp;"]]; class calQA proc;
      calMetricsJson[/"fa:fa-file cal_metrics.json "/]; class calMetricsJson file
      calSol & polyCal --> calQA --> calMetricsJson
      %% metafits -...-> calQA

      calGate{"fa:fa-check calGate &nbsp;&nbsp;"}; class calGate decision;
      calMetricsJson --> calGate

      acaciaCal([fa:fa-box-open acacia]); class acaciaCal out
      calSol & polyCal & calMetricsJson --> acaciaCal

      plotSolutions[["fa:fa-gem plotSolutions "]]; class plotSolutions proc
      plotSol[/"fa:fa-file-image plot_{phases,amps}.png "/]; class plotSol file
      polyCal & calSol --> plotSolutions --> plotSol
      %% metafits -...-> plotSolutions
    end

    subgraph uvfits
      %% nfConfigApply>fa:fa-file-code nextflow.config ]; class nfConfigApply in
      calUVFits[/fa:fa-file cal uvfits /]; class calUVFits file

      hypApplyUV[["fa:fa-times-circle hypApplyUV &nbsp;"]]; class hypApplyUV proc
      calGate --"uvfits"--> hypApplyUV --> calUVFits
      %% metafits -...-> hypApplyUV
      nfConfigCal --"apply args"--> hypApplyUV

      %% end
      %% subgraph sub
      subUVFits[/fa:fa-file sub uvfits /]; class subUVFits file
      hypSubUV[["fa:fa-minus-circle hypSubUV &nbsp;"]]; class hypSubUV proc
      calUVFits --> hypSubUV --> subUVFits

      %% end
      %% subgraph ps_metrics
      psMetricsDat[/fa:fa-file ps_metrics.dat /]; class psMetricsDat file
      psMetrics[["fa:fa-gem ps_metrics &nbsp;"]]; class psMetrics proc;
      subUVFits & calUVFits --> psMetrics --> psMetricsDat

      %% end
      %% subgraph vis_qa
      visMetricsJson[/fa:fa-file vis_metrics.json /]; class visMetricsJson file
      visQA[["fa:fa-gem visQA &nbsp;"]]; class visQA proc;
      calUVFits --> visQA --> visMetricsJson

      visGate{"fa:fa-check visGate &nbsp;&nbsp;"}; class visGate incompleteDecision;
      visMetricsJson & psMetricsDat -.-> visGate

      acaciaUVFits([fa:fa-box-open acacia]); class acaciaUVFits out
      visMetricsJson & calUVFits & subUVFits --> acaciaUVFits
    end

    subgraph img
      calMS[/fa:fa-file cal ms /]; class calMS file
      hypApplyMS[["fa:fa-times-circle hypApplyMS &nbsp;"]]; class hypApplyMS proc
      visGate --"uvfits"--> hypApplyMS --> calMS
      %% calGate --"uvfits"--> hypApplyMS --> calMS
      %% applyArgs -...-> hypApplyMS
      %% metafits -...-> hypApplyMS
      subMS[/fa:fa-file sub ms /]; class subMS file
      hypSubMS[["fa:fa-minus-circle hypSubMS &nbsp;"]]; class hypSubMS proc
      calMS --> hypSubMS --> subMS

      imgDirty[/"fa:fa-file-image image-{XX,YY,V}-dirty.fits "/]; class imgDirty file
      wscleanDirty[["fa:fa-image wscleanDirty &nbsp;"]]; class wscleanDirty proc
      calMS & subMS --> wscleanDirty --> imgDirty

      %% end
      %% subgraph img_qa
      imgMetricsJson[/"fa:fa-file img_metrics.json &nbsp;"/]; class imgMetricsJson file
      imgQA[[fa:fa-gem imgQA]]; class imgQA proc;
      imgDirty --> imgQA --> imgMetricsJson

      acaciaImg([fa:fa-box-open acacia]); class acaciaImg out
      imgMetricsJson & imgDirty --> acaciaImg

      imgGate{"fa:fa-check imgGate &nbsp;&nbsp;"}; class imgGate incompleteDecision;
      imgMetricsJson -.-> imgGate

      nfConfigCal --"apply args"--> hypApplyMS
    end

    subgraph "Power Spectrum"
      chips[[fa:fa-bolt chips]]; class chips proc
      powerSpectrum[/"fa:fa-file-image power_spectrum "/]; class powerSpectrum file
      imgGate -.-> chips -.-> powerSpectrum
    end

    %% subgraph archive
      %% acacia([fa:fa-box-open acacia]); class acacia out
      %% calUVFits --> acacia
      %% calMetricsJson --> acacia
    %% end
  %% end
```

## Components

main tasks:

- obsid → **`asvoRaw`** → obsid, metafits, \*gpufits
  - if obsid raw not stored, schedule ASVO download job to Accacia with
    [Giant Squid](github.com/mwaTelescope/giant-squid), wait,
    download with `wget`, untar with `tar`
  - ASVO wait times can be between a few minutes and a few days, so job will
    wait with ASVO socket open for an hour, then exponential backoff for $2^a$
    hours for attempt number $a$, up to 5 attempts
  - store: `${obsid}/raw`
  - resources: mem
- metafits, \*gpufits → **`birliPrep`** → prepUVFits, \*mwaf, birliLog
  - if prepUVFits for obsid not stored, preprocess and flag with `birli`
  - store: `${obsid}/prep`
  - resources: mem, cpu
- metafits, mwaf → **`flagQA`** → occupancy
  - get total flag occupancy, and occupancy for each coarse channel
  - reject obs flag occupancy is above threshold
- metafits, prepUVFits, dicalArgs → **`hypCalSol`** → \*calSol, \*dicalLog
  - if calSols not stored, `hyperdrive di-calibrate` with dicalArgs
  - store: `${obsid}/cal${params.cal_suffix}`
  - resources: mem, gpu
- metafits, prepUVFits, calSol, visName, applyArg → **`hypApplyUV`** → calUVFits, applyLog
  - if calUVFits for (obsid × visName) not stored, `hyperdrive solutions-apply` with applyArg
  - store: `${obsid}/cal${params.cal_suffix}`
  - resources: mem, cpu
- name, metafits, calUVFits → **`wscleanDirty`** → img{XX,YY,V}Dirty
  - dirty images of stokes XX,YY,V via [wsclean](https://gitlab.com/aroffringa/wsclean)

qa tasks:

- obsid, metafits → **`metaJson`** → obsid, metafitsJson
  - get flagged inputs from metafits
  - size of raw vis for missing hdu fraction
- obsid, prepUVfits, birliLog → **`prepStats`** → obsid, prepStatsJson
  - from prep uvfits: count timesteps, channels, baselines
  - from birli log: collect
  - store: `${obsid}/prep`
- obsid, name, metafits, calSol → **`calQA`** (#3) → obsid, calMetricsXJson, calMetricsYJson
  - [mwa_qa](https://github.com/Chuneeta/mwa_qa) `scripts/run_valqa.py`
  - needs a unique name for each calibration
- obsid, name, metafits, calSol → **`plotSolutions`** → obsid, phasesPng, ampsPng
  - plot calibration solution gains,phases with `hyperdrive solutions-plot`
- obsid, name, calUVFits → **`visQA`** → obsid, visMetrics
  - [mwa_qa](https://github.com/Chuneeta/mwa_qa) `scripts/run_visqa.py`
  - needs a unique name for each vis
- obsid, name, calUVFits → **`psMetrics`** (#1/#2) → obsid, psMetricsDat
  - window/wedge/total power/subtraced power iono proxy via [chips](https://github.com/cathryntrott/chips) `src/ps_power_metrics.c`
- obsid, name, imgXXDirty, imgYYDirty, imgVDirty → **`imgQA`** (#4) → imgMetricsJson
  - [mwa_qa](https://github.com/Chuneeta/mwa_qa) `scripts/run_valqa.py`
  - flux density ratio of stokes, ratio of corner thermal noise RMS, variance, ratio stokes V/XX V/YY
  - needs a unique name for each image

## Output Directory Structure:

- `${params.outdir}/` (e.g. `/data/curtin_mwaeor/data/`)
  - `${obsid}/`
    - `raw/` - raw gpubox fits, metafits, metafits json
    - `prep/` - preprocessed birli vis uvfits (weights encode flags), flag mwaf, birli log,
      preprocessing stats json
    - `cal${params.cal_suffix}/` - calibrated vis uvfits, di-cal solution fits, hyperdrive di-cal,
      apply log
    - `img/` - wscleaned dirty XX,YY,V image fits
    - `flag_qa/` - flag occupancy json
    - `cal_qa/` - calibration metrics json, hyperdrive solution amp,phase plot png
    - `vis_qa/` - visibility metrics json
    - `img_qa/` - image metrics json
    - `ps_metrics/` - power spectrum metrics dat, log
  - `results/` - aggregate metrics for spreadsheet import
    - `raw_stats.tsv` - info from raw stage: raw/metafits size on disk, raw file count
    - `metafits_stats.tsv` - info from metafits: raw dimensions, flagged inputs
    - `prep_stats.tsv` - info from prep stage: prep uvfits dimensions, missing HDUs, uvfits / mwaf size on disk
    - `occupancy.tsv` - mwaf occupancy: total, by coarse channel
    - `cal_metrics.tsv` - calQA metrics for each di calibration
    - `vis_metrics.tsv` - visQA metrics for each vis
    - `img_metrics.tsv` - imgQA metrics for each image
    - `ps_metrics.tsv` - power spectrum metrics for each vis

## Configuration

populate `obsids.csv`:

```
obsid1
obsid2
```

populate `dical_args.csv` as `dical_name,dical_args`, e.g. :

```
30l_src4k,--uvw-min 30l -n 4000
50l_src4k,--uvw-min 50l -n 4000
```

populate `apply_args.csv` as `dical_name,apply_name,apply_args`:

```
30l_src4k,8s_80kHz,--time-average 8s --freq-average 80kHz
50l_src4k,8s_80kHz,--time-average 8s --freq-average 80kHz
30l_src4k,4s_80kHz,--time-average 4s --freq-average 80kHz
```

all other config is in `nextflow.config`. See: <https://www.nextflow.io/docs/latest/config.html>

parameters can also be specified in the newflow command line for each run too. e.g. set `params.obsids_path` with `--obsids_path=...`

## usage

load modules - conflicts with singularity for dumb Java reasons.

```
module unload singularity
module load nextflow
```

run everything (use `obsids.csv` by default)

```
nextflow run main.nf -profile dug -with-timeline -with-report -with-dag --asvo_api_key=$MWA_ASVO_API_KEY
```

run only a subset of obsids:

```
nextflow run main.nf -profile dug -with-timeline -with-report -with-dag --asvo_api_key=$MWA_ASVO_API_KEY \
  --obsids_path=obsids-1090X.csv
```

## Handy DuG commands

### get an interactive session with gpus

```bash
salloc --partition=curtin_mwaeor --constraint=v100 --time=1:00:00
srun --pty bash
```

### run hyperdrive di-cal

needs GPU, see above

```bash
module load cuda/11.3.1
module load gcc-rt/9.2.0
export HYPERDRIVE_CUDA_COMPUTE=... # 80 for a100s, 70 for v100s
export MWA_BEAM_FILE="/data/curtin_mwaeor/data/mwa_full_embedded_element_pattern.h5"
export CUDA_VISIBLE_DEVICES=... # optional, restrict CUDA devices.
export obsid=...
export cal_name=... # optional, e.g. "30l_src4k", "50l_src4k"
export apply_args=... # optional, e.g. "--uvw-min 30l -n 4000", "-v"
export sourcelist="/data/curtin_mwaeor/data/srclist_pumav3_EoR0LoBESv2_fixedEoR1pietro+ForA_phase1+2.yaml"
cd /data/curtin_mwaeor/data/${obsid}/
/usr/bin/time /data/curtin_mwaeor/sw/bin/hyperdrive di-calibrate ${cal_args} \
  --data "raw/${obsid}.metafits" "prep/birli_${obsid}.uvfits" \
  --beam "${MWA_BEAM_FILE}" \
  --source-list "${sourcelist}" \
  --outputs "hyp_soln_${obsid}_${cal_name}.fits" \
  | tee hyp_di-cal_${obsid}_${cal_name}_test.log
```

### run hyperdrive apply solutions

```bash
module load gcc-rt/9.2.0
export obsid=...
export cal_name=... # either "30l_src4k" or "50l_src4k"
export apply_args=... # optional, e.g. "--time-average 16s --freq-average 160kHz", "-v"
cd /data/curtin_mwaeor/data/${obsid}/
/usr/bin/time /data/curtin_mwaeor/sw/bin/hyperdrive solutions-apply ${apply_args} \
  --data "raw/${obsid}.metafits" "prep/birli_${obsid}.uvfits" \
  --solutions "cal/hyp_soln_${obsid}_${cal_name}.fits" \
  --outputs "hyp_${obsid}_${cal_name}_test.uvfits" | tee "hyp_${obsid}_${cal_name}_test.log"
```

### run wsclean

```bash
# salloc --partition=curtin_mwaeor --constraint='knl&nogpu' --time=1:00:00
# srun --pty bash
module use /data/curtin_mwaeor/sw/modules
module load wsclean/3.1-knl-gcc
export obsid=...
cd /data/curtin_mwaeor/data/${obsid}
OMP_NUM_THREADS=60 /data/curtin_mwaeor/sw/wsclean/3.1-knl-gcc/bin/wsclean -j 60 \
    -weight briggs -1.0 \
    -name wsclean_hyp_${obsid}_sub_30l_src4k_8s_80kHz \
    -size 2048 2048 \
    -scale 40asec \
    -pol xx,yy,v \
    -channels-out 24 \
    -niter 0 \
    cal/hyp_${obsid}_30l_src4k_4s_80kHz.ms
```

### get times of jobs for a run

```bash
export run="elegant_euclid"
nextflow log $run -F 'process=="calQA"' -f workdir,exit,status,process,duration
```

### get all lines matching pattern from all scripts executed recently

```bash
export process="metaJson"
export first_run="$(nextflow log -q | head -n 1)"
export first_run="scruffy_bhaskara"
echo -n $'${start} ${workdir} ${script.split(\'\\n\').findAll {it =~ /.*out.*/}[0]}' | tee template.md
nextflow log -after $first_run -F 'process=="'${process}'"&&exit==0' -t template.md
```

### Updating singularity images

e.g. Birli or mwa_qa

```bash
module load singularity
# export container="birli"
export container="mwa_qa"
# export url="docker://mwatelescope/${container}:latest"
export url="docker://d3vnull0/${container}:latest"
cd /data/curtin_mwaeor/sw/singularity/
if [ ! -d $container ]; then mkdir $container; fi
singularity pull --force --dir $container "$url"
```

optional: update `profiles.dug.params.birli` in `nextflow.config`

### Cancel all failed jobs

```bash
squeue --states SE --format %A -h | sort | xargs scancel
```

### get disk usage

```bash
for stage in cal cal_qa flag_qa img img_qa prep ps_metrics raw vis_qa; do
  find /data/curtin_mwaeor/data -type d -name "${stage}" | xargs du -csh | tee "du_${stage}.tsv"
done
```

### dump gpufits

for legacy: `corr_type="MWA_ORD"`
for mwax: `corr_type="MWAX"`

```bash
module load python/3.9.7
export obsid=...;
export corr_type=...;
cd /data/curtin_mwaeor/data/${obsid}/raw/;
for gpufits in $(ls *.fits); do \
  python /data/curtin_mwaeor/src/Birli/tests/data/dump_gpufits.py \
    ${gpufits} \
    --corr-type $corr_type --timestep-limit 7 \
    | tee ${gpufits}.txt
done
```

### Carta

no config

```bash
singularity exec --bind /data/curtin_mwaeor/data:/images --contain /data/curtin_mwaeor/sw/singularity/carta/carta_latest.sif carta --no_user_config --no_system_config --no_browser /images
```

user config

```bash
singularity exec --bind /data/curtin_mwaeor/data:/images /data/curtin_mwaeor/sw/singularity/carta/carta_latest.sif carta --no_system_config --debug_no_auth /images
```

open `all_imgs`, sort by name, unix search `wsclean_hyp_*_sub_poly_30l_src4k_8s_80kHz-MFS-XX-dirty.fits`

## make csv file of all obsids at various stages

```bash
# downloaded, preprocessed.
find /data/curtin_mwaeor/data/ -type f -path '*/prep/birli_*.uvfits' -print | sed -e 's!.*/\([0-9]\{10\}\)/.*!\1!g' | sort | uniq | tee obsids-downloaded.csv
# calibrated
find /data/curtin_mwaeor/data/ -type f -path '*/cal/hyp_*.uvfits' -print | sed -e 's!.*/\([0-9]\{10\}\)/.*!\1!g' | sort | uniq | tee obsids-calibrated.csv
```
