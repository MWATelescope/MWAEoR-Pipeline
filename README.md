# MWA EoR Nextflow Pipeline

## Flow

The MWAEoR pipeline workflow takes a list of obsids and produces detailed metrics about the quality of the data by analysing the data in various stages. At the end of each stage is a gate, which prevents observations from passing to the next stage if an issue is found.

### Metadata

```mermaid
%%{init:{'theme':'base','themeVariables':{'primaryBorderColor':'#002b36','secondaryColor':'#eee8d5','tertiaryColor':'#fdf6e3','tertiaryBorderColor':'#002b36','lineColor':'#002b36','textColor':'#002b36'}}}%%
graph TD;
  classDef in fill:#2aa198;
  classDef out fill:#d33682;
  classDef file fill:#268bd2;
  classDef proc fill:#b58900;
  classDef decision fill:#cb4b16;

  subgraph
    obsids --> wsMeta --> ws --> wsmetaJson;
    ws ---> metafits
    qualityUpdates --"quality updates"-----> wsGate;
    wsmetaJson --"pointing, tile, dipole, quality"--> wsGate;
    nfConfigMeta --"filters"-----> wsGate;

    obsids>fa:fa-file obsids.csv ]; class obsids in;
    wsMeta[[fa:fa-download wsMeta ]]; class wsMeta proc;
    ws([fa:fa-server MWA web services ]); class ws in;
    wsmetaJson[/fa:fa-file wsmeta.json /]; class wsmetaJson file;
    %% wsfilesJson[/fa:fa-file files.json /]; class wsfilesJson file;
    metafits[/fa:fa-file metafits /]; class metafits file;
    wsGate{fa:fa-check wsGate }; class wsGate decision;
    nfConfigMeta>fa:fa-file-code nextflow.config ]; class nfConfigMeta in
    qualityUpdates>fa:fa-file-code quality-updates.csv ]; class qualityUpdates in
  end
```

The metadata stage gathers information about each observation from MWA Web Services. This ranges from scheduling information to information about faults that occurred during the observation, as well as information about the files archived for the observation. This information is used to filter observations based on a set of criteria that prevents observations with too many faults from passing through to further stages.

### Preprocessing

```mermaid
%%{init:{'theme':'base','themeVariables':{'primaryBorderColor':'#002b36','secondaryColor':'#eee8d5','tertiaryColor':'#fdf6e3','tertiaryBorderColor':'#002b36','lineColor':'#002b36','textColor':'#002b36'}}}%%
graph TD;
  classDef in fill:#2aa198;
  classDef out fill:#d33682;
  classDef file fill:#268bd2;
  classDef proc fill:#b58900;
  classDef decision fill:#cb4b16;

  subgraph
    nfConfig --"spectral, temporal resolution"--> asvoPrep
    asvoPrep --> asvo --> prepUVFits;
    prepUVFits & metafits --> flagQA --> occupancyJson --> flagGate;
    prepUVFits --> prepVisQA --> prepVisJson;
    nfConfig --"occupancy threshold"-----> flagGate;
    tileUpdates --"manual flags"-----> tile_flags;
    prepVisJson --"prep flags"--> tile_flags;

    asvoPrep[[fa:fa-bolt asvoPrep ]]; class asvoPrep proc;
    asvo([fa:fa-server MWA ASVO ]); class asvo in;
    prepUVFits[/fa:fa-file prep uvfits /]; class prepUVFits file;
    metafits[/fa:fa-file metafits /]; class metafits file;
    flagQA[[fa:fa-gem flagQA]]; class flagQA proc;
    occupancyJson[/fa:fa-file-code occupancy.json /]; class occupancyJson file;
    prepVisQA[[fa:fa-gem prepVisQA ]]; class prepVisQA proc;
    prepVisJson[/fa:fa-file-code prepVis.json /]; class prepVisJson file;
    flagGate{fa:fa-flag flagGate }; class flagGate decision;
    nfConfig>fa:fa-file-code nextflow.config ]; class nfConfig in;
    tileUpdates>fa:fa-file-code tile-updates.csv ]; class tileUpdates in;
    tile_flags[/fa:fa-file-code tile_flags.csv /]; class tile_flags file;
  end
```

The Preprocessing stage produces and analyses preprocessed uvfits visibilities. If the preprocessed visibilities for an observation are not found in storage, the preprocessing stage schedules a [Birli](https://github.com/MWATelescope/Birli) conversion job on [MWA ASVO](https://asvo.mwatelescope.org/) via [Acacia](https://pawsey.org.au/systems/acacia/) with [Giant Squid](github.com/mwaTelescope/giant-squid). The Birli conversion job includes flagging with [AOFlagger](https://gitlab.com/aroffringa/aoflagger). Further preprocessing parameters can be configured in nextflow,

ASVO wait times can be between a few minutes and a few days, so an exponential backoff is used to wait $2^a$ hours after attempt number $a$ for up to 5 attempts until the conversion job is ready. Finally the archive is downloaded with `wget`, hash-validated, and inflated with `tar`.

The [prep vis qa script](https://github.com/Chuneeta/mwa_qa/blob/main/scripts/run_prepvisqa.py) is used to determine bad tiles using outlier analysis.

Further bad tiles can be specified manually using a tile updates file

Total flag occupancy of the preprocessed visibilities, as well as the occupancy for each coarse channel is determined by the flagQA step.

Finally, observations whose flag occupancy is above the defined occupancy threshold are rejected.

### Direction Independent Calibration

```mermaid
%%{init:{'theme':'base','themeVariables':{'primaryBorderColor':'#002b36','secondaryColor':'#eee8d5','tertiaryColor':'#fdf6e3','tertiaryBorderColor':'#002b36','lineColor':'#002b36','textColor':'#002b36'}}}%%
graph TD;
  classDef in fill:#2aa198;
  classDef out fill:#d33682;
  classDef file fill:#268bd2;
  classDef proc fill:#b58900;
  classDef decision fill:#cb4b16;

  subgraph
    nfConfig --"dical args"--> hypCalSol --> calSol
    prepUVFits & metafits --> hypCalSol
    %% hypCalSol --> dicalLog
    calSol --> polyFit --> polyCal
    polyCal & calSol & metafits --> calQA --> calMetricsJson
    calMetricsJson --> calGate
    calSol & polyCal & calMetricsJson --> acaciaCal
    polyCal & calSol & metafits --> plotSolutions --> plotSol

    calSol[/fa:fa-file-excel cal solutions/]; class calSol file
    hypCalSol[[fa:fa-wrench hypCalSol]]; class hypCalSol proc
    %% dicalLog[/fa:fa-file dical log /]; class dicalLog file
    nfConfig>fa:fa-file-code nextflow.config]; class nfConfig in
    metafits[/fa:fa-file metafits/]; class metafits file
    prepUVFits[/fa:fa-file prep uvfits/]; class prepUVFits file
    polyCal[/fa:fa-file-excel poly solutions /]; class polyCal file
    polyFit[[fa:fa-chart-line polyFit]]; class polyFit proc
    calQA[[fa:fa-gem calQA]]; class calQA proc;
    calMetricsJson[/"fa:fa-file cal_metrics.json "/]; class calMetricsJson file
    calGate{fa:fa-check calGate}; class calGate decision;
    acaciaCal([fa:fa-box-open acacia]); class acaciaCal out
    plotSolutions[[fa:fa-gem plotSolutions]]; class plotSolutions proc
    plotSol[/"fa:fa-file-image plot_{phases,amps}.png "/]; class plotSol file
  end
```

[Hyperdrive](github.com/mwaTelescope/mwa_hyperdrive) is used to generate one or more direction independent calibration solutions using the model file and several calibration parameter sets defined in the nextflow config. The statistical properties of the calibration solutions are analysed with the [cal qa](https://github.com/Chuneeta/mwa_qa/blob/main/scripts/run_calqa.py). If the variance of an observation's calibration solutions is too high, then the observation is rejected.

### Calibrated Visibility Analysis

```mermaid
%%{init:{'theme':'base','themeVariables':{'primaryBorderColor':'#002b36','secondaryColor':'#eee8d5','tertiaryColor':'#fdf6e3','tertiaryBorderColor':'#002b36','lineColor':'#002b36','textColor':'#002b36'}}}%%
graph TD;
  classDef in fill:#2aa198;
  classDef out fill:#d33682;
  classDef file fill:#268bd2;
  classDef proc fill:#b58900;
  classDef decision fill:#cb4b16;

  subgraph
    metafits & prepUVFits --> hypApplyUV --> calUVFits
    nfConfig --"apply args"--> hypApplyUV
    calUVFits --> hypSubUV --> subUVFits
    subUVFits & calUVFits --> psMetrics --> psMetricsDat
    calUVFits --> visQA --> visMetricsJson
    visMetricsJson & psMetricsDat --> visGate
    visMetricsJson & calUVFits & subUVFits --> acaciaUVFits

    nfConfig>fa:fa-file-code nextflow.config ]; class nfConfig in
    calUVFits[/fa:fa-file cal uvfits /]; class calUVFits file
    hypApplyUV[[fa:fa-times-circle hypApplyUV ]]; class hypApplyUV proc
    metafits[/fa:fa-file metafits /]; class metafits file
    prepUVFits[/fa:fa-file prep uvfits /]; class prepUVFits file
    subUVFits[/fa:fa-file sub uvfits /]; class subUVFits file
    hypSubUV[[fa:fa-minus-circle hypSubUV ]]; class hypSubUV proc
    psMetricsDat[/fa:fa-file ps_metrics.dat /]; class psMetricsDat file
    psMetrics[[fa:fa-gem ps_metrics ]]; class psMetrics proc;
    visMetricsJson[/fa:fa-file vis_metrics.json /]; class visMetricsJson file
    visQA[[fa:fa-gem visQA ]]; class visQA proc;
    visGate{fa:fa-check visGate }; class visGate decision;
    acaciaUVFits([fa:fa-box-open acacia ]); class acaciaUVFits out
  end
```

- obsid, name, calUVFits → **`visQA`** → obsid, visMetrics
  - [mwa_qa](https://github.com/Chuneeta/mwa_qa) `scripts/run_visqa.py`
  - needs a unique name for each vis
- obsid, name, calUVFits → **`psMetrics`** (#1/#2) → obsid, psMetricsDat
  - window/wedge/total power/subtraced power iono proxy via [chips](https://github.com/cathryntrott/chips) `src/ps_power_metrics.c`
- metafits, prepUVFits, calSol, visName, applyArg → **`hypApplyUV`** → calUVFits, applyLog
  - if calUVFits for (obsid × visName) not stored, `hyperdrive solutions-apply` with applyArg
  - store: `${obsid}/cal${params.cal_suffix}`
  - resources: mem, cpu

### Image Analysis

```mermaid
%%{init:{'theme':'base','themeVariables':{'primaryBorderColor':'#002b36','secondaryColor':'#eee8d5','tertiaryColor':'#fdf6e3','tertiaryBorderColor':'#002b36','lineColor':'#002b36','textColor':'#002b36'}}}%%
graph TD;
  classDef in fill:#2aa198;
  classDef out fill:#d33682;
  classDef file fill:#268bd2;
  classDef proc fill:#b58900;
  classDef decision fill:#cb4b16;

  subgraph
    metafits & prepUVFits --> hypApplyMS --> calMS --> hypSubMS --> subMS
    calMS & subMS --> wscleanDirty --> imgDirty --> imgQA --> imgMetricsJson
    imgMetricsJson & imgDirty --> acaciaImg
    imgMetricsJson -.-> imgGate
    nfConfig --"apply args"--> hypApplyMS

    nfConfig>fa:fa-file-code nextflow.config ]; class nfConfig in
    metafits[/fa:fa-file metafits /]; class metafits file
    prepUVFits[/fa:fa-file prep uvfits /]; class prepUVFits file
    calMS[/fa:fa-file cal ms /]; class calMS file
    hypApplyMS[[fa:fa-times-circle hypApplyMS ]]; class hypApplyMS proc
    subMS[/fa:fa-file sub ms /]; class subMS file
    hypSubMS[[fa:fa-minus-circle hypSubMS ]]; class hypSubMS proc
    imgDirty[/"fa:fa-file-image image-{XX,YY,V}-dirty.fits "/]; class imgDirty file
    wscleanDirty[[fa:fa-image wscleanDirty ]]; class wscleanDirty proc
    imgMetricsJson[/fa:fa-file img_metrics.json /]; class imgMetricsJson file
    imgQA[[fa:fa-gem imgQA]]; class imgQA proc;
    acaciaImg([fa:fa-box-open acacia]); class acaciaImg out
    imgGate{fa:fa-check imgGate }; class imgGate decision;

  end
```

- name, metafits, calUVFits → **`wscleanDirty`** → img{XX,YY,V}Dirty
  - dirty images of stokes XX,YY,V via [wsclean](https://gitlab.com/aroffringa/wsclean)
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
# imaged
```

# WIP

### hard link images (easier for carta)

```bash
find /data/curtin_mwaeor/data -path '*/img/wsclean*.fits' -exec sh -c 'ln {} all_imgs/$(basename {})' \;
```

## copy all metafits to hopper

```bash
cd /data/curtin_mwaeor/data
find . -regextype posix-extended  -regex './[0-9]{10}/raw/[0-9]{10}.metafits.json' | while read -r p; do cp $p /data/curtin_mwaeor/FRB_hopper/; done
```

```bash
for obsid in 1060539904 1060540152 1060540272 1061319352 1061320080 1061320328 1061320448 1061320568 1065196248 1065196368 1065196856 1065196976 1065197344 1185484472 1259411080 1259411320 1259411560 1259412160 1259412280 1322135592 1322135712 1322135832 1322307760 1322307880 1322308000 1322308240 1322480408 1322480528 1322654736 1322826184 1322826424 1322826664 1322826784; do
  echo $'\n'$obsid
  ls /data/curtin_mwaeor/data/${obsid}/cal_qa/*src4k_{X,Y}.json
done
```

## show the script of last few runs

`head -n -1` only if a pipeline is currently running

```bash
while read run; do
  echo $run
  nextflow log $run \
    -f workdir,exit,status,process,duration,script \
    -filter 'process == "wscleanDirty" && exit == 0'
done < <(nextflow log -q | tail -n 50 | head -n -1) | tee runs.txt
```

```bash
while read run; do
  echo $run
  nextflow log $run \
    -f tag,exit,status,process,realtime \
    -filter 'process == "cal:hypCalSol"'
  # nextflow log $run \
  #   -f workdir,exit,status,process,duration,script \
  #   -filter 'process == "wscleanDirty" && exit == 0'
done < <(nextflow log -q | tail -n 50 | head -n 1) | tee runs.txt
cat runs.txt | grep -E $'(-name|work)' | tee runs_clean.txt
```

## dump img_qa

```bash
find /data/curtin_mwaeor/data -type f -path '*/img_qa/*MFS.json' \
  | gawk 'match()
```

## dump flags

```bash
module load python/3.9.7
export obsid=1059505936
for i in {01..12}; do \
  python /data/curtin_mwaeor/src/Birli/tests/data/dump_mwaf.py \
    /data/curtin_mwaeor/data/${obsid}/prep/${obsid}_${i}.mwaf --summary
done
```

## clear blank logs

```bash
find /data/curtin_mwaeor/data/ -name '*.log' -size 0 -delete
```

## cleanup folders

```
cd /data/curtin_mwaeor/data
for obsid in ...; do
  [ -d /data/curtin_mwaeor/data/$obsid ] && rm -rf /data/curtin_mwaeor/data/$obsid
done
```

## cleanup files

```bash
cd /data/curtin_mwaeor/data
for obsid in $(ls); do
  # export path=${obsid}/raw/${obsid}.metafits.json
  # export path=${obsid}/prep/${obsid}_prep_stats.json
  # export path=${obsid}/{vis_qa,img_qa,cal_qa}/*.json
  # export path=${obsid}/vis_qa/*.json
  export path=${obsid}/img_qa/*.json
  # export path=${obsid}/cal_qa/*.json
  for path in $path; do
    [ -f "$path" ] && echo "$path" && rm -rf "$path"
  done
done
find /data/curtin_mwaeor/data -path '*/*_qa/*.json' # -delete
find /data/curtin_mwaeor/data -path '*/cal/*.ms' -exec sh -c 'echo rm -rf {}' \;
```

```
ls /data/curtin_mwaeor/data/*/img/*.fits | cut -d / -f 5 | sort | uniq | tee obsids-with-img.csv | tee /dev/stderr | wc -l
```

```bash
for obsid in \
  1259412400\
  1259412280\
  1259412160\
  1259412040\
  1259411920\
  1259411800\
  1259411680\
  1259411560\
  1259411440\
  1259411320\
  1259411200\
  1259411080\
  1065197464\
  1065197344\
  1065197224\
  1065197096\
  1065196976\
  1065196856\
  1065196736\
  1065196608\
  1065196488\
  1065196368\
  1065196248\
  1065196120\
  1065196000\
; do
  echo $obsid
  # find /data/curtin_mwaeor/data/$obsid -type f -path './img_qa/*.png' #-delete
  rm -f /data/curtin_mwaeor/data/${obsid}/img_qa/*.png
done
```

## don't delete:

1090008640
1094490328

### dump metafits

not needed any more, just look in `${obsid}.metafits.json`

```bash
module load python/3.9.7
for obsid in ... ; do \
  python /data/curtin_mwaeor/src/Birli/tests/data/dump_metafits.py \
    /data/curtin_mwaeor/data/${obsid}/raw/${obsid}.metafits \
    | tee /data/curtin_mwaeor/data/${obsid}/raw/${obsid}.metafits.txt
done
```

## add prefix to existing file stage

```bash
bash -c 'unset PATH; /bin/find /data/curtin_mwaeor/data -type d -name img_qa -execdir sh -c "pwd && mv {} img_qa-dci1000" \;'
```

## ipython via conda

create env

```bash
conda create --prefix /data/curtin_mwaeor/sw/conda/dev astropy ipython scipy 'numpy<1.23.0' scipy matplotlib pandas seaborn six ipykernel
```

add new modules

```bash
conda install --prefix /data/curtin_mwaeor/sw/conda/dev astropy ipython scipy 'numpy<1.23.0' scipy matplotlib pandas seaborn six ipykernel
```

```bash
module load miniconda/4.10.3
conda activate /data/curtin_mwaeor/sw/conda/dev
ipython
```

```python
from astropy.io import fits
hdus=fits.open('/data/curtin_mwaeor/data/1060539904/raw/1060539904.metafits')
hdus=fits.open('/data/curtin_mwaeor/data/test/1343457784/raw/1343457784.metafits')
hdus=fits.open('/data/curtin_mwaeor/data/1322827024/raw/1322827024.metafits')
print('\n'.join(map(lambda row: ','.join(map(str, row)), hdus['TILEDATA'].data['Delays'].tolist())))
hdus=fits.open('/data/curtin_mwaeor/data/1060539904/cal/hyp_soln_1060539904_30l_src4k.fits')
hdus=fits.open('/data/curtin_mwaeor/data/1322827024/cal/hyp_soln_1322827024_30l_src4k.fits')
hdus=fits.open('/data/curtin_mwaeor/data/1090094680/cal/hyp_soln_1090094680_30l_src4k.fits')
print('\n'.join(map(lambda row: ','.join(map(str, row)), hdus['TILES'].data['DipoleGains'].tolist())))

hdus=fits.open('/data/curtin_mwaeor/data/1060540392/prep/birli_1060540392.uvfits')
```

### mwa_qa debugging

```python
import sys
from astropy.io import fits
sys.path.insert(0, '/data/curtin_mwaeor/src/mwa_qa')
from mwa_qa.read_uvfits import UVfits
# hdus=fits.open('/data/curtin_mwaeor/data/1061319840/cal/hyp_1061319840_sub_30l_src4k_8s_80kHz.uvfits')
uvf = UVfits('/data/curtin_mwaeor/data/1061319840/cal/hyp_1061319840_sub_30l_src4k_8s_80kHz.uvfits')
```

### run json query on metrics

```bash
find /data/curtin_mwaeor/data/ -type f -path '*/cal_qa/*X.json' -print -exec jq -r '(.CONVERGENCE//[])[0]|length' {} \; | tee has_convergence.txt
```

### rendering image

```python
import matplotlib.pyplot as plt
from astropy.visualization import astropy_mpl_style, make_lupton_rgb, simple_norm, SqrtStretch, ImageNormalize, MinMaxInterval
from mpl_toolkits.axes_grid1.axes_rgb import RGBAxes
from astropy.wcs import WCS
from astropy.io import fits
import numpy as np

prefix = '/data/curtin_mwaeor/data/1259412400/img-dci1000/wsclean_hyp_1259412400_poly_30l_src4k_8s_80kHz-MFS-'
plt.style.use(astropy_mpl_style)

with fits.open(f'{prefix}XX-dirty.fits') as hdus:
  header = hdus[0].header
wcs = WCS(header)

img_data_xx = fits.getdata(f'{prefix}XX-dirty.fits', ext=0)[0,0,:,:]
img_data_yy = fits.getdata(f'{prefix}YY-dirty.fits', ext=0)[0,0,:,:]
img_data_v = fits.getdata(f'{prefix}V-dirty.fits', ext=0)[0,0,:,:]
# normalize to 0-1 using 95th percentile
i_percentile = np.percentile(np.stack((img_data_xx, img_data_yy)), 95)
img_data_xx /= i_percentile
img_data_yy /= i_percentile
v_percentile = np.percentile(img_data_v, 95)
img_data_v /= v_percentile
fig = plt.gcf()
fig.set_size_inches(20, 20)
plt.subplot(projection=wcs, slices=('x', 'y', 0, 0))
plt.imshow(make_lupton_rgb(img_data_xx, img_data_v, img_data_yy, Q=1))
plt.savefig('wsclean_hyp_1259412400_poly_30l_src4k_8s_80kHz-MFS-polarimetry-dirty.png')
# plt.imshow(hdu.data, vmin=-2.e-5, vmax=2.e-4, origin='lower')
# ny, nx = img_data_xx.shape
# ax=RGBAxes(plt.figure(), [0.1, 0.1, 0.8, 0.8], pad=0.0)
# ax.imshow(img_data_xx, img_data_v, img_data_yy)
# rgb = np.zeros((ny, nx, 3))
# rgb[:,:,0] = img_data_xx
# rgb[:,:,1] = img_data_v
# rgb[:,:,2] = img_data_yy
# plt.imshow(img_data_xx, cmap='Reds', extent=extent)
# plt.imshow(img_data_yy, cmap='Blues', extent=extent)
```

### tap stuff

```python
import pyvo as vo
tap = vo.dal.TAPService("http://vo.mwatelescope.org/mwa_asvo/tap")
results = tap.search("select obs_id, gpubox_files_archived, gpubox_files_total, total_archived_data_bytes \
  from mwa.observation \
  where obs_id IN (1087250776)")
```

### render movies

#### prepvisqa

```bash
module load ffmpeg
# ,crop=in_w-3600:in_h:800:in_h
ffmpeg -y -framerate 5 -pattern_type glob -i '/data/curtin_mwaeor/data/??????????/prep/prepvis_metrics_??????????_rms.png' -c:v libx264 -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2" -pix_fmt yuv420p "/data/curtin_mwaeor/data/results/prepvis_metrics_rms.mp4"
```

#### polcomps

```bash
module load ffmpeg
# ,crop=in_w-3600:in_h:800:in_h
for name in 30l sub_30l poly_30l sub_poly_30l; do
  ffmpeg -y -framerate 5 -pattern_type glob -i '/data/curtin_mwaeor/data/??????????/img_qa/??????????_'$name'_*.png' -c:v libx264 -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2" -pix_fmt yuv420p "/data/curtin_mwaeor/data/results/polcomp_${name}.mp4"
done
for name in 30l sub_30l; do
  ffmpeg -y -framerate 5 -pattern_type glob -i '/data/curtin_mwaeor/data/12????????/img_qa/??????????_'$name'_*.png' -c:v libx264 -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2" -pix_fmt yuv420p "/data/curtin_mwaeor/data/results/polcomp_12_${name}.mp4"
done
```

#### calQA

```bash
module load ffmpeg
# ,crop=in_w-3600:in_h:800:in_h
for name in 30l_src4k_dlyspectrum poly_30l_src4k_dlyspectrum 30l_src4k_fft poly_30l_src4k_fft 30l_src4k_variance poly_30l_src4k_variance; do
  ffmpeg -y -framerate 5 -pattern_type glob -i '/data/curtin_mwaeor/data/??????????/cal_qa/calmetrics_??????????_'$name'.png' -c:v libx264 -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2" -pix_fmt yuv420p "/data/curtin_mwaeor/data/results/calmetrics_${name}.mp4"
done
```

#### plot solutions

```bash
module load ffmpeg
# ,crop=in_w-3600:in_h:800:in_h
for name in 30l_src4k_phases poly_30l_src4k_phases 30l_src4k_amps poly_30l_src4k_amps; do
  ffmpeg -y -framerate 5 -pattern_type glob -i '/data/curtin_mwaeor/data/??????????/cal_qa/hyp_soln_??????????_'$name'.png' -c:v libx264 -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2" -pix_fmt yuv420p "/data/curtin_mwaeor/data/results/hyp_soln_${name}.mp4"
done
for name in 30l_src4k_phases 30l_src4k_amps; do
  ffmpeg -y -framerate 5 -pattern_type glob -i '/data/curtin_mwaeor/data/12????????/cal_qa/hyp_soln_??????????_'$name'.png' -c:v libx264 -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2" -pix_fmt yuv420p "/data/curtin_mwaeor/data/results/hyp_soln_12_${name}.mp4"
done
```

#### redo imagqa from acacia

```bash
export name="sub_30l_src4k_8s_80kHz"
while read -r obsid; do rclone copy mwaeor:high0.imgqa/wsclean_hyp_${obsid}_${name}.json .; jq -r '['${obsid}',"'${name}'",.XX.RMS_ALL,.XX.RMS_BOX,.XX.PKS0023_026?.PEAK_FLUX,.XX.PKS0023_026?.INT_FLUX,.YY.RMS_ALL,.YY.RMS_BOX,.YY.PKS0023_026.PEAK_FLUX,.YY.PKS0023_026.INT_FLUX,.V.RMS_ALL,.V.RMS_BOX,.V.PKS0023_026?.PEAK_FLUX,.V.PKS0023_026?.INT_FLUX,.V_XX.RMS_RATIO,.V_XX.RMS_RATIO_BOX]|@csv' wsclean_hyp_${obsid}_${name}.json; done < obsids-stage2imgqa.csv | tee img_metrics.csv
```

#### Fast flagged vs unflagged

```bash
for name in 30l_src4_fast_amps 30l_src4_fast_phases; do
  ffmpeg -y -framerate 5 -pattern_type glob -i '/data/curtin_mwaeor/data/??????????/cal_qa-flagged-fast/hyp_soln_??????????_'$name'.png' -c:v libx264 -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2" -pix_fmt yuv420p "/data/curtin_mwaeor/data/results-test/hyp_soln_${name}-flagged.mp4"
  ffmpeg -y -framerate 5 -pattern_type glob -i '/data/curtin_mwaeor/data/??????????/cal_qa-unflagged-fast/hyp_soln_??????????_'$name'.png' -c:v libx264 -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2" -pix_fmt yuv420p "/data/curtin_mwaeor/data/results-test/hyp_soln_${name}-unflagged.mp4"
done
for name in 30l_src4_fast_variance 30l_src4_fast_fft 30l_src4_fast_dlyspectrum; do
  ffmpeg -y -framerate 5 -pattern_type glob -i '/data/curtin_mwaeor/data/??????????/cal_qa-flagged-fast/calmetrics_??????????_'$name'.png' -c:v libx264 -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2" -pix_fmt yuv420p "/data/curtin_mwaeor/data/results-test/calmetrics_${name}-flagged.mp4"
  ffmpeg -y -framerate 5 -pattern_type glob -i '/data/curtin_mwaeor/data/??????????/cal_qa-unflagged-fast/calmetrics_??????????_'$name'.png' -c:v libx264 -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2" -pix_fmt yuv420p "/data/curtin_mwaeor/data/results-test/calmetrics_${name}-unflagged.mp4"
done
```

### Hollow out measurement sets

```bash
find /data/curtin_mwaeor/data/ -type d -path '*/cal/*.ms' -exec rm -rf {}/* \;
```

### (Dev) Upload flags

```bash
find . -name '*.mwaf' -exec sh -c 'rclone copyto {} "dev:eor0high.mwaf/$(basename -- {})"' \;
```

### (Dev) exfil ms

```bash
export obsid=1092338912
tar -czvg hyp_${obsid}_30l_src4k_2s_80kHz.ms.tar.gz -C /data/curtin_mwaeor/FRB_hopper hyp_${obsid}_30l_src4k_2s_80kHz.ms
rclone copy hyp_${obsid}_30l_src4k_2s_80kHz.ms.tar.gz "dev:eor0high.ms"
```
