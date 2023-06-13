# MWA EoR Nextflow Pipeline

## Flow

The MWAEoR pipeline workflow takes a list of obsids and produces detailed metrics about the quality of the data by analysing the data in various stages. At the end of each stage is a gate, which prevents observations from passing to the next stage if an issue is found. Finally, power spectra are produced from the best performing observations.

### Metadata

```mermaid
%%{init:{'theme':'base','themeVariables':{'primaryBorderColor':'#002b36','secondaryColor':'#eee8d5','tertiaryColor':'#fdf6e3','tertiaryBorderColor':'#002b36','lineColor':'#002b36','textColor':'#002b36'}}}%%
graph TD;
  classDef in fill:#2aa198;
  classDef out fill:#d33682;
  classDef file fill:#268bd2;
  classDef proc fill:#b58900;
  classDef decision fill:#cb4b16;

  subgraph " "
    obsids --> wsMeta --> ws --> wsmetaJson;
    ws --> metafits
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

The metadata stage gathers information about each observation from MWA Web Services. This ranges from scheduling information to information about faults that occurred during the observation, as well as the status the files archived for the observation. This is used to filter observations based on a set of criteria that prevents observations with too many faults from being downloaded in the first place.

### Preprocessing

```mermaid
%%{init:{'theme':'base','themeVariables':{'primaryBorderColor':'#002b36','secondaryColor':'#eee8d5','tertiaryColor':'#fdf6e3','tertiaryBorderColor':'#002b36','lineColor':'#002b36','textColor':'#002b36'}}}%%
graph TD;
  classDef in fill:#2aa198;
  classDef out fill:#d33682;
  classDef file fill:#268bd2;
  classDef proc fill:#b58900;
  classDef decision fill:#cb4b16;

  subgraph " "
    nfConfig --"spectral, temporal resolution"--> asvoPrep
    asvoPrep --> asvo --> prepUVFits;
    prepUVFits & metafits --> flagQA --> occupancyJson --> flagGate;
    prepUVFits --> prepVisQA --> prepVisJson;
    nfConfig --"occupancy threshold"-----> flagGate;
    tileUpdates --"manual flags"-------> tile_flags;
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

The Preprocessing stage produces and analyses uvfits visibilities that have been preprocessed, flagged and averaged by [Birli](https://github.com/MWATelescope/Birli) via [MWA ASVO](https://asvo.mwatelescope.org/). The flag occupancy of the uvfits files are analysed by `flagQA` for the presence of RFI over each coarse channel, and the RMS amplitude of the autocorrelations are measured across frequency and antenna by [`prepVisQA`](https://github.com/Chuneeta/mwa_qa/blob/main/scripts/run_prepvisqa.py). These measures are used to produce a list of outlier antennas, which are flagged in later stages.

If the preprocessed files are not present, the `asvoPrep` process schedules a conversion job on ASVO using [Giant Squid](github.com/mwaTelescope/giant-squid). When the files are ready, they are downloaded from [Acacia](https://pawsey.org.au/systems/acacia/).

Finally, observations whose `flagQA` flag occupancy is above the configured threshold threshold are rejected.

### Direction Independent Calibration

```mermaid
%%{init:{'theme':'base','themeVariables':{'primaryBorderColor':'#002b36','secondaryColor':'#eee8d5','tertiaryColor':'#fdf6e3','tertiaryBorderColor':'#002b36','lineColor':'#002b36','textColor':'#002b36'}}}%%
graph TD;
  classDef in fill:#2aa198;
  classDef out fill:#d33682;
  classDef file fill:#268bd2;
  classDef proc fill:#b58900;
  classDef decision fill:#cb4b16;

  subgraph " "
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

[Hyperdrive](github.com/mwaTelescope/mwa_hyperdrive) is used to generate one or more direction independent calibration solutions using the [MWA Long Baseline Epoch of Reionisation Survey (LoBES)](https://arxiv.org/abs/2110.08400) catalogue.

The pipeline has the facilities to compare multiple independent calibration solutions for each obsid. A "smoothed" version of the calibration solutions is generated by `polyFit`. Each calibration solution, as well as the smoothed solution are analysed independently.

The statistical properties of the calibration solutions are analysed with [`calQA`](https://github.com/Chuneeta/mwa_qa/blob/main/scripts/run_calqa.py). If the variance of an observation's calibration solutions is too high, then the observation is rejected.

### Calibrated Visibility Analysis

```mermaid
%%{init:{'theme':'base','themeVariables':{'primaryBorderColor':'#002b36','secondaryColor':'#eee8d5','tertiaryColor':'#fdf6e3','tertiaryBorderColor':'#002b36','lineColor':'#002b36','textColor':'#002b36'}}}%%
graph TD;
  classDef in fill:#2aa198;
  classDef out fill:#d33682;
  classDef file fill:#268bd2;
  classDef proc fill:#b58900;
  classDef decision fill:#cb4b16;

  subgraph " "
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

Calibrated UVFits visibilities are obtained by applying calibration solutions to the preprocessed visibilities. An analysis of redundant baseline groups is performed by [`visQA`](https://github.com/Chuneeta/mwa_qa/blob/main/scripts/run_visqa.py) to detect anomalous baselines that do not conform to the behaviour of their group.

An analysis of power spectrum metrics is performed by [`ps_metrics`](https://github.com/cathryntrott/chips) on these visibilities, as well as their residual, which is obtained by subtracting the simulated sky model from the calibrated visibilities. In addition to the un-subtracted visibilities (`nosub` for short), and residual (`sub` for short) visibilities, an ionospheric subtraction also produces `ionosub` visibilities. All three types of visibilities are analysed separately. EoR Power spectrum window and wedge measurements are used to gate observations with too much contamination.

### Image Analysis

```mermaid
%%{init:{'theme':'base','themeVariables':{'primaryBorderColor':'#002b36','secondaryColor':'#eee8d5','tertiaryColor':'#fdf6e3','tertiaryBorderColor':'#002b36','lineColor':'#002b36','textColor':'#002b36'}}}%%
graph TD;
  classDef in fill:#2aa198;
  classDef out fill:#d33682;
  classDef file fill:#268bd2;
  classDef proc fill:#b58900;
  classDef decision fill:#cb4b16;

  subgraph " "
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

Calibrated and residual visibilities from the calibrated visibility analysis stage are analysed in measurement set format. [wsclean](https://gitlab.com/aroffringa/wsclean) is used to make dirty deconvolved images of Stokes XX,YY and V polarisations, which are analysed by [`imgQA`](https://github.com/Chuneeta/mwa_qa/blob/main/scripts/run_imgqa.py).

### Power Spectra

```mermaid
%%{init:{'theme':'base','themeVariables':{'primaryBorderColor':'#002b36','secondaryColor':'#eee8d5','tertiaryColor':'#fdf6e3','tertiaryBorderColor':'#002b36','lineColor':'#002b36','textColor':'#002b36'}}}%%
graph TD;
  classDef in fill:#2aa198;
  classDef out fill:#d33682;
  classDef file fill:#268bd2;
  classDef proc fill:#b58900;
  classDef decision fill:#cb4b16;

  subgraph "chipsSpec"
    chipsGrid[[fa:fa-times-circle CHIPS: gridvisdiff & prepare_diff ]]; class chipsGrid proc;
    chipsGridOut[/"fa:fa-file {vis_tot,vis_diff,noise_tot,noise_diff,weights}_${pol}_${ext}.dat"/]; class chipsGridOut file;
    chipsCombine[[fa:fa-times-circle CHIPS: combine_data]]; class chipsCombine proc;
    %% chipsCombineOut[/"fa:fa-file {vis_tot,vis_diff,noise_tot,noise_diff,weights}_${pol}_${group}.dat"/]; class chipsCombineOut file;
    chipsLssa[[fa:fa-times-circle CHIPS: lssa_fg_simple]]; class chipsLssa proc;
    chipsLssaOut[/"fa:fa-file {crosspower,residpower,residpowerimag,totpower,flagpower,fg_num,outputweights}*.dat"/]; class chipsLssaOut file;
    %% chipsLssaOut[/"fa:fa-file {crosspower,residpower,residpowerimag,totpower,flagpower,fg_num,outputweights}_${pol}_${group}.dat"/]; class chipsLssaOut file;

    uvfits --> chipsGrid --> chipsGridOut --> chipsCombine --> chipsGridOut
    chipsGridOut --> chipsLssa --> chipsLssaOut

    class uvfits file;
  end
```

Finally, observations that have passed all QA steps are gridded with CHIPS, and combined into integrated power spectra.

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

### obsids

populate `obsids.csv` with each obsid to be processed, one per line

```
obsid1
obsid2
```

### filters

these can be filtered by `params.filter_pointings`, `params.filter_bad_tile_frac`, `params.filter_dead_dipole_frac` and `params.filter_quality` in `nextflow.config`

### quality updates

if you discover several obsids are of poor quality, but you haven't yet updated the data quality number in the mwa metadata database, you can populate `quality-updates.csv` with `obsid,dataquality,dataqualitycomment` to filter these obsids

### tile updates path

if you discover tiles that need to be flagged manually over a range of obsids, you can populate `tile-updates.csv` with `startgps,endgps,tiles,comment` to manually flag these tiles during calibration. `tiles` are pipe separated tile indices as used in the TILE_DATA metafits HDU, and Hyperdrive calibraiton solution plots (tile61 == tile id #40 != tile id #61), e.g.

```txt
1125766160,1125767872,40|41|42|43|44|45|46|47,chaotic fringes in high channels on recv06
1126921360,1126923520,55,Tile98 fringes
1321443064,1321448344,4,Tile15 fringes
```

### di-cal args

multiple hyperdrive calibrations per obsid can be obtained by populating `params.dical_args` as a map of `[dical_name:dical_args]` in `nextflow.config`, where `dical_args` are passed to `hyperdrive di-cal`, e.g.

```
dical_args = [
  "30l_src4k": "--uvw-min 30l -n 4000",
  "50l_src1k": "--uvw-min 50l -n 1000",
]
```

all other config is in `nextflow.config`. See: <https://www.nextflow.io/docs/latest/config.html>

parameters can also be specified in the newflow command line for each run too. e.g. set `params.obsids_path` with `--obsids_path=...`

## usage

- Ensure a profile has been set up in `nextflow.config` for your cluster or local machine (see: `profiles.dug`, `profiles.garrawarla`).
- Ensure `nextflow` version 21 is available (more on this later).
- Ensure the `export MWA_ASVO_API_KEY=...` environment is set.
- Tip: You'll want to run this inside of `tmux` or `screen`, otherwise all your jobs will be cancelled when your ssh connection is interrupted. `nohup` doesn't seem to work.
- Tip: I like to work out of a vscode remote ssh session on hpc-data, since the head node is much slower. just bear in mind you need to `ssh garrawarla` to run the pipeline or see `squeue`.

```bash
nextflow run main.nf -profile <profile> -with-timeline -with-report -with-dag --asvo_api_key=$MWA_ASVO_API_KEY
```

if you would like to run run only a subset of obsids, you can specify them in a file, e.g. `obsids-test.csv`

```bash
nextflow run main.nf --obsids_path=obsids-test.csv
```

you can also modify any of the params in `nextflow.config` from the command line, e.g.

```bash
nextflow run main.nf --img_suffix='-briggs+0.5' --img_weight='briggs +0.5'
```

### pawsey quirks

Running this pipeline on Pawsey requires the `cluster_edit` nextflow module. Get the available versions of nextflow with `module avail nextflow` and load the module with (e.g.)

```bash
module load nextflow_cluster_edit/22.04.3_cluster_edit
```

then run nextflow with

```bash
nextflow run main.nf -profile garrawarla
```

architecture reference:
- garra workq & gpuq: `Intel(R) Xeon(R) Gold 6230` (cascade lake)
- garra head: `Intel(R) Xeon(R) Silver 4215` - (cascade lake)
- hpc-data{1..6} (zeus): `AMD EPYC 7351`

### tracing

```bash
tail -n 999999 -F .nextflow.log | grep -E 'COMPLETED; exit: [^0]'
```


### dug quirks

nextflow module conflicts with singularity for dumb Java reasons.

```
module unload singularity
module load nextflow
nextflow run main.nf -profile dug
```

## handy Nextflow commands

### get times of jobs of the last run

```bash
export run_name="$(nextflow log -q | tail -n 1)"
export filter='process=~/.*hypCalSol/'
nextflow log "${run_name}" -F "${filter}" -f "workdir,exit,status,process,tag,duration,realtime,memory,peak_vmem"
```

### get all lines matching pattern from all scripts executed recently

```bash
export process="hypCalSol"
export first_run="$(nextflow log -q | head -n 1)"
echo -n $'${start} ${workdir} ${script.split(\'\\n\').findAll {it =~ /.*wsclean.*/}[0]}' | tee template.md

# or for csv
echo -n $'${process},${tag},${realtime},${memory},${read_bytes},${write_bytes},${rss},${vmem},${cpus},${pcpu},${pmem}' | tee template.csv
nextflow log -after $first_run -F 'exit==0' -t template.csv | tee stats.csv
```

### show the script of last few runs

`head -n -1` only if a pipeline is currently running

```bash
while read run; do
  echo $run
  nextflow log $run \
    -f workdir,exit,status,process,duration,script \
    -filter 'process == "wscleanDirty" && exit == 0'
done < <(nextflow log -q | tail -n 50 | head -n -1) | tee runs.txt
```

### suffixes

```
# this will put results in $outdir/test/results-a-b-c-d
nextflow run main.nf -profile garrawarla --outdir=$outdir/test/ --img_suffix=-a --cal_suffix=-b --obsids_suffix=-c --results_suffix=-d
```

### Run for all clusters

```bash
for c in obsids-rn-*.csv; do \
  c=${c#obsids}
  nextflow run main.nf -profile garrawarla -w /astro/mwaeor/dev/nfwork --obsids_suffix=${c%.csv} --noprep
done
```

## Handy Singularity commands

### Updating singularity images

```bash
module load singularity
# - specify the container to update:
export container="mwa_qa" # or "birli"
# - specify the url of the docker image
export url="docker://mwatelescope/${container}:latest" # or export url="docker://d3vnull0/${container}:latest"
# - cd into singularity directory
cd /pawsey/mwa/singularity/ # or /data/curtin_mwaeor/sw/singularity/
# - make a directory for the container if it doesn't exist
[ -d $container ] || mkdir $container
# - create the singularity image
singularity pull --force --dir $container "$url"
```

optional: update `profiles.<profile>.params.<container>_sif` in `nextflow.config` if things have changed

## Handy slurm commands

### Cancel all failed jobs

```bash
squeue -u $USER --format %A -h --states SE | sort | xargs scancel
```

### get info for jobid: standard out

```bash
export jobid=3405303
squeue --json | /pawsey/mwa/dev/bin/jq -r '.jobs[]|select(.job_id==$jobid)|.standard_output' --argjson jobid 3459463
```

### live updating queue for current user

```bash
while true; do echo ""; date -Is; squeue -u $USER -t R --format "%.7i %.9P %60j %.2t %.8M %.8L %.3c %.4m %.4d %.34Z"; squeue -u $USER -t PD --format "%.7i %.9P %60j %.2t %.8M %.8L %.3c %.4m %.4d %.34Z" | head -n 20; sleep 15; done
```

### get job directory

```bash
squeue -u $USER -t R --format "%.7i %.9P %30j %.2t %.8M %.8L %.36Z"
```

## Handy Giant Squid Commands

### schedule conversion jobs

run this before any pipeline runs which download new obsids

[conversion options](https://wiki.mwatelescope.org/pages/viewpage.action?pageId=5963859#MWAASVO:CommandLineClient(mantarayclient)-ConversionJobOptions)

```bash
giant-squid submit-conv -p output=uvfits,no_rfi=true, obsids-eclipse-2023-04-20.csv
```

### get obsids that are ready to download

```bash
giant-squid list --states=Ready -j | jq -r '.[]|[.obsid]|@csv' | tee obsids-ready.csv | tee /dev/stderr | wc -l
```

### reschedule conversion jobs where some failed.

```bash
giant-squid list --states=Error -j | jq -r '.[]|[.obsid]|@csv' \
  | while read -r obsid; do giant-squid submit-conv \
    -p timeres=2,freqres=40,conversion=uvfits,preprocessor=birli $obsid; done
```

## handy storage commands

### get disk usage of each stage

```bash
export outdir="/data/curtin_mwaeor/data"
# export outdir="/astro/mwaeor/dev/nfdata"
for stage in cal cal_qa flag_qa img img_qa prep ps_metrics raw vis_qa; do
  find $outdir -type d -name "${stage}" | xargs du -csh | tee "du_${stage}.tsv"
done
```

### make csv file of all obsids at various stages

```bash
export outdir="/data/curtin_mwaeor/data"
# export outdir="/astro/mwaeor/dev/nfdata"
# downloaded, preprocessed.
find $outdir -type f -path '*/prep/birli_*.uvfits' -print | sed -e 's!.*/\([0-9]\{10\}\)/.*!\1!g' | sort | uniq | tee obsids-downloaded.csv | tee /dev/stderr | wc -l
# calibrated
find $outdir -type f -path '*/cal/hyp_*.uvfits' -print | sed -e 's!.*/\([0-9]\{10\}\)/.*!\1!g' | sort | uniq | tee obsids-calibrated.csv | tee /dev/stderr | wc -l
# iono peeled
find $outdir -type f -path '*/cal/hyp_peel*.json' -print | sed -e 's!.*/\([0-9]\{10\}\)/.*!\1!g' | sort | uniq | tee obsids-peeled.csv | tee /dev/stderr | wc -l
# images
ls $outdir/*/img/*.fits | cut -d / -f 5 | sort | uniq | tee obsids-imgd.csv | tee /dev/stderr | wc -l
```

### get files not downloaded

```bash
comm -13 <(sort obsids-downloaded.csv) <(sort obsids-stage1.csv) > obsids-stage1-undownloaded.csv
```

### clear blank logs

```bash
find $outdir -name '*.log' -size 0 -delete
```

### add suffix to existing file stage

```bash
bash -c 'unset PATH; /bin/find ${outdir} -type d -name img_qa -execdir sh -c "pwd && mv {} {}-dci1000" \;'
```

add suffix but keep original cal

```bash
export suffix="-oldhyp"
for dir in ??????????/{cal_qa,vis_qa,img_qa,img-briggs+0.5_4096px,img_qa-briggs+0.5_4096px}; do
  [ -d "${dir}" ] && mv "${dir}" "${dir}${suffix}";
done
for dir in ??????????/cal; do
  [ -d "${dir}" ] && mv "${dir}" "${dir}-kill";
done
for dir in ??????????/cal-kill; do
  mkdir -p "${dir/-kill/}";
  for file in ${dir}/hyp_*.fits ${dir}/hyp_di-cal*.log; do
    [ -f "${file}" ] && mv "${file}" "${file/-kill/}";
  done
  rm -rf "${dir}";
done
for dir in ??????????/*-oldhyp; do
  rm -rf "${dir}";
done
```

### cleanup folders

```
cd ${outdir}
for obsid in ...; do
  [ -d ${outdir}/$obsid ] && rm -rf ${outdir}/$obsid
done
```

### cleanup files

```bash
cd ${outdir}
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
find ${outdir} -path '*/*_qa/*.json' # -delete
find ${outdir} -path '*/cal/*.ms' -exec sh -c 'echo rm -rf {}' \;
```

```bash
for obsid in ...; do
  echo $obsid
  # find ${outdir}/$obsid -type f -path './img_qa/*.png' #-delete
  rm -f ${outdir}/${obsid}/img_qa/*.png
done
```

### download calibration solutions

```bash
module load rclone
cd $outdir
while read -r obsid; do
  mkdir -p ${obsid}/cal-noflag
  rclone copy mwaeor:high0.soln/hyp_soln_${obsid}_30l_src4k.fits ${obsid}/cal-noflag/
  touch ${obsid}/cal-noflag/hyp_di-cal_${obsid}_30l_src4k_8s_80kHz.log
done < obsids.csv
```

### Bootstrap visibility files

```bash
salloc --nodes=1 --mem=350G --gres=gpu:1 --time=8:00:00 --clusters=garrawarla --partition=gpuq -c 18 --tmp=500G
cd /nvmetmp/
module use pawsey/mwa/software/python3/modulefiles
module load rclone hyperdrive/peel

export MWA_BEAM_FILE="/pawsey/mwa/mwa_full_embedded_element_pattern.h5"
export metafits="${obsid}.metafits"
export prep_uvfits="birli_${obsid}_2s_40kHz.uvfits"
export hyp_srclist="/pawsey/mwa/software/python3/srclists/master/srclist_pumav3_EoR0LoBESv2_fixedEoR1pietro+ForA_phase1+2_edit.txt"
export hyp_soln="hyp_soln_${obsid}_30l_src4k.fits"
export cal_uvfits="hyp_${obsid}_30l_src4k_8s_80kHz.uvfits"
export sub_uvfits="hyp_${obsid}_sub_30l_src4k_8s_80kHz.uvfits"
export ionosub_uvfits="hyp_${obsid}_ionosub_30l_src4k_8s_80kHz.uvfits"
export ionosub_json="hyp_${obsid}_ionosub_30l_src4k_8s_80kHz.json"
export tile_flags="1,2,3,4..."

wget -O "${metafits}" "http://ws.mwatelescope.org/metadata/fits?obs_id=${obsid}"
rclone copy mwaeor:high0.prep/${prep_uvfits} .
rclone copy mwaeor:high0.soln/${hyp_soln} .

hyperdrive solutions-apply \
        --data "${metafits}" "${prep_uvfits}" \
        --solutions "${hyp_soln}" \
        --tile-flags ${tile_flags} \
        --outputs "${cal_uvfits}"

hyperdrive vis-subtract \
        --data "${cal_uvfits}" \
        --source-list "${hyp_srclist}" \
        --outputs "${sub_uvfits}" \
        --invert \
        --num-sources 8000

hyperdrive peel \
        --data "${cal_uvfits}" \
        --source-list "${hyp_srclist}" \
        --outputs "${ionosub_uvfits}" "${ionosub_json}" \
        --sub 8000 \
        --iono-sub 1000

# copy files off /nvmetmp before the interactive job ends
```

## Misc handy commands

### rclone

list buckets

```bash
rclone lsd mwaeor:
```

### Carta

```bash
export carta_sif=/data/curtin_mwaeor/sw/singularity/carta/carta_latest.sif
export carta_sif=/pawsey/mwa/singularity/carta/carta_latest.sif
```

no config

```bash
singularity exec --bind $PWD:/images --contain $carta_sif carta --no_user_config --no_system_config --no_browser /images
```

user config

```bash
singularity exec --bind $outdir:/images /data/curtin_mwaeor/sw/singularity/carta/carta_latest.sif carta --no_system_config --debug_no_auth /images
```

open `all_imgs`, sort by name, unix search `wsclean_hyp_*_sub_poly_30l_src4k_8s_80kHz-MFS-XX-dirty.fits`

### export image snippet

first enable code snippets and custom title text

```js
for (const frame of app.frames) {
  frame.setTitleCustomText(frame.filename);
  await app.setRasterScalingMatchingEnabled(frame, true);
  await app.setActiveFrame(frame);
  await new Promise(r => setTimeout(r, 1000));
  console.log(app.activeFrame.filename);
  await app.exportImage(1);
  await new Promise(r => setTimeout(r, 1000));
};
```

### hard link images (easier for carta)

```bash
find ${outdir} -path '*/img/wsclean*.fits' -exec sh -c 'ln {} all_imgs/$(basename {})' \;
```

### copy all metafits to hopper

```bash
find ${outdir} -regextype posix-extended  -regex './[0-9]{10}/raw/[0-9]{10}.metafits.json' | while read -r p; do cp $p /data/curtin_mwaeor/FRB_hopper/; done
```

### ipython via conda

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
from mwa_qa.read_uvfits import UVfits
from mwa_qa.read_calfits import CalFits
# hdus=fits.open('1061319840/cal/hyp_1061319840_sub_30l_src4k_8s_80kHz.uvfits')
uvf = UVfits('1061319840/cal/hyp_1061319840_sub_30l_src4k_8s_80kHz.uvfits')
cf = CalFits('1251908712/cal/hyp_soln_1251908712_30l_src4k.fits')
```

### run json query on metrics

```bash
find /data/curtin_mwaeor/data/ -type f -path '*/cal_qa/*X.json' -print -exec jq -r '(.CONVERGENCE//[])[0]|length' {} \; | tee has_convergence.txt
```


### tap stuff

```bash
module load singularity
singularity exec --bind $PWD --home $PWD --cleanenv         /pawsey/mwa/singularity/tap/tap_latest.sif python
```

```python
import pyvo as vo
tap = vo.dal.TAPService("http://vo.mwatelescope.org/mwa_asvo/tap")
results = tap.search("select obs_id, gpubox_files_archived, gpubox_files_total, total_archived_data_bytes \
  from mwa.observation \
  where obs_id IN (1087250776)")
```

```python
__import__('pyvo').dal.TAPService("http://vo.mwatelescope.org/mwa_asvo/tap").search("""
SELECT TOP 100000 obs_id, projectid, projectshortname, starttime_mjd, starttime_utc, ra_pointing, dec_pointing, ra_phase_center, dec_phase_center, gridpoint_number, gridpoint_name, local_sidereal_time_deg, obsname, freq_res, int_time, total_tiles, bad_tiles, tiles_flagged, dataquality, total_archived_data_bytes/1024/1024/1024 as size_gb
FROM mwa.observation
WHERE total_archived_data_bytes > 0
  AND dataquality <= 1
  AND ra_phase_center = 0
  AND dec_phase_center = -27
  AND first_channel_number = 131
  AND last_channel_number = 154
ORDER BY obs_id ASC
""").to_table().to_pandas().to_csv('eor0high_obs.tsv', sep='\t', index=False)
```

### combine images

```bash
module load singularity
mkdir -p "${outdir}/combined"
for name in sub_30l ionosub_30l; do
  for pol in XX YY V; do
    for prefix in 109; do
      export result="${outdir}/combined/wsclean_hyp_${prefix}_${name}_src4k_8s_80kHz-MFS-${pol}-image.fits"
      singularity exec /pawsey/mwa/singularity/imagemagick/imagemagick_latest.sif convert \
        ${outdir}/${prefix}???????/img/wsclean_hyp_??????????_${name}_src4k_8s_80kHz-MFS-${pol}-image.fits \
        -average -auto-gamma -auto-level $result
      echo $result
    done
  done
done

for name in sub_30l ionosub_30l; do
  for pol in XX YY V; do
    for prefix in 10; do
    export result="${outdir}/combined/avg_hyp_${prefix}_${name}_src4k_8s_80kHz_MFS-${pol}-image.png";
    singularity exec /pawsey/mwa/singularity/imagemagick/imagemagick_latest.sif convert \
      ${outdir}/${prefix}????????/img_qa/??????????_${name}_src4k_8s_80kHz_MFS-${pol}-image.png \
      -average -auto-gamma -auto-level $result;
    echo $result;
    done;
  done;
done
```

### manually render movies

#### render prepvisqa

```bash
# singularity exec /pawsey/mwa/singularity/ffmpeg/ffmpeg_latest.sif ffmpeg
module load ffmpeg
# ,crop=in_w-3600:in_h:800:in_h
ffmpeg -y -framerate 5 -pattern_type glob \
 -i "${outdir}/??????????/prep/prepvis_metrics_??????????_rms.png" \
  -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2" \
  -pix_fmt yuv420p \
  "${outdir}/results/prepvis_metrics_rms.mp4"
```

#### render polcomps

```bash
module load ffmpeg
# ,crop=in_w-3600:in_h:800:in_h
for name in 30l sub_30l; do
  ffmpeg -y -framerate 5 -pattern_type glob\
  -i "${outdir}/??????????/img_qa/??????????_${name}_*.png" \
  -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2" -pix_fmt yuv420p \
  "${outdir}/results/polcomp_${name}.mp4"
done
for range in 122 125 132; do
  for name in 30l sub_30l; do
    ffmpeg -y -framerate 5 -pattern_type glob\
    -i "${outdir}/${range}???????/img_qa-briggs+0.5_4096px/??????????_${name}_*.png" \
    -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2" \
    -pix_fmt yuv420p \
    "${outdir}/results/polcomp_${range}_${name}_briggs+0.5_4096px.mp4"
  done
done
```

#### render visQA

```bash
module load ffmpeg
# ,crop=in_w-3600:in_h:800:in_h
ffmpeg -y -framerate 5 -pattern_type glob\
  -i "${outdir}/??????????/vis_qa/hyp_??????????_30l_src4k_8s_80kHz_vis_metrics_rms.png" \
  -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2" \
  -pix_fmt yuv420p \
  "${outdir}/results/vis_metrics_rms.mp4"
```

#### render calQA

```bash
module load ffmpeg
# ,crop=in_w-3600:in_h:800:in_h
for name in 30l_src4k_dlyspectrum poly_30l_src4k_dlyspectrum 30l_src4k_fft poly_30l_src4k_fft 30l_src4k_variance poly_30l_src4k_variance; do
  ffmpeg -y -framerate 5 -pattern_type glob\
  -i "${outdir}/??????????/cal_qa/calmetrics_??????????_${name}.png" \
  -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2" \
  -pix_fmt yuv420p \
  "${outdir}/results/calmetrics_${name}.mp4"
done
```

#### plot solutions

```bash
module load ffmpeg
# ,crop=in_w-3600:in_h:800:in_h
for name in 30l_src4k_phases poly_30l_src4k_phases 30l_src4k_amps poly_30l_src4k_amps; do
  ffmpeg -y -framerate 5 -pattern_type glob\
  -i "${outdir}/??????????/cal_qa/hyp_soln_??????????_${name}.png" \
  -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2" \
  -pix_fmt yuv420p \
  "${outdir}/results/hyp_soln_${name}.mp4"
done
for name in 30l_src4k_phases 30l_src4k_amps; do
  ffmpeg -y -framerate 5 -pattern_type glob\
  -i "${outdir}/12????????/cal_qa/hyp_soln_??????????_${name}.png" \
  -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2" \
  -pix_fmt yuv420p \
  "${outdir}/results/hyp_soln_12_${name}.mp4"
done
```

### kickstart obsid from acacia

```bash
export obsid=...
cd $outdir
mkdir $obsid
cd $obsid
mkdir -p raw
wget -O raw/${obsid}.metafits "http://ws.mwatelescope.org/metadata/fits?obs_id=${obsid}&include_ppds=1"
mkdir -p prep
rclone copy mwaeor:high0.prep/birli_${obsid}_2s_40kHz.uvfits prep/
mkdir -p cal
touch cal/hyp_di-cal_${obsid}_30l_src4k.log
rclone copy mwaeor:high0.soln/hyp_soln_${obsid}_30l_src4k.fits cal/
touch cal/hyp_apply_30l_src4k_8s_80kHz.log
rclone copy mwaeor:high0.uvfits/hyp_${obsid}_30l_src4k_8s_80kHz.uvfits cal/
touch cal/hyp_vis-sub_30l_src4k_8s_80kHz_uv.log
rclone copy mwaeor:high0.uvfits/hyp_${obsid}_sub_30l_src4k_8s_80kHz.uvfits cal/
mkdir -p img
touch img/wsclean_30l_src4k_8s_80kHz.log
rclone copy mwaeor:high0.img/wsclean_hyp_${obsid}_30l_src4k_8s_80kHz-MFS-$'{XX,YY,V}'-dirty.fits img/
touch img/wsclean_sub_30l_src4k_8s_80kHz.log
```

#### redo imagqa from acacia

```bash
export name="sub_30l_src4k_8s_80kHz"
while read -r obsid; do rclone copy mwaeor:high0.imgqa/wsclean_hyp_${obsid}_${name}.json .; jq -r '['${obsid}',"'${name}'",.XX.RMS_ALL,.XX.RMS_BOX,.XX.PKS0023_026?.PEAK_FLUX,.XX.PKS0023_026?.INT_FLUX,.YY.RMS_ALL,.YY.RMS_BOX,.YY.PKS0023_026.PEAK_FLUX,.YY.PKS0023_026.INT_FLUX,.V.RMS_ALL,.V.RMS_BOX,.V.PKS0023_026?.PEAK_FLUX,.V.PKS0023_026?.INT_FLUX,.V_XX.RMS_RATIO,.V_XX.RMS_RATIO_BOX]|@csv' wsclean_hyp_${obsid}_${name}.json; done < obsids-stage2imgqa.csv | tee img_metrics.csv
```

#### copy from backup

```bash
for o in ??????????; do
  f=cal/hyp_soln_${o}_30l_src4k.fits;
  # f=raw/${o}.metafits;
  [ -f $o/$f ] && echo $o | tee -a /pawsey/mwa/mwaeor/dev/MWAEoR-Pipeline/obsids-calib.csv

  # echo -n $f - ;
  # && [ -f $o/$f ] && echo -n e \
  # && [ ! -f ../nfdata/$o/$f ] && echo -n r \
  # && mkdir -p ../nfdata/$o/${f%/*} && mv $o/$f ../nfdata/$o/$f && echo -n m; \
  # echo '';
done
```

#### Fast flagged vs unflagged

```bash
for name in 30l_src4_fast_amps 30l_src4_fast_phases; do
  ffmpeg -y -framerate 5 -pattern_type glob -i '/data/curtin_mwaeor/data/??????????/cal_qa-flagged-fast/hyp_soln_??????????_'$name'.png' -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2" -pix_fmt yuv420p "/data/curtin_mwaeor/data/results-test/hyp_soln_${name}-flagged.mp4"
  ffmpeg -y -framerate 5 -pattern_type glob -i '/data/curtin_mwaeor/data/??????????/cal_qa-unflagged-fast/hyp_soln_??????????_'$name'.png' -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2" -pix_fmt yuv420p "/data/curtin_mwaeor/data/results-test/hyp_soln_${name}-unflagged.mp4"
done
for name in 30l_src4_fast_variance 30l_src4_fast_fft 30l_src4_fast_dlyspectrum; do
  ffmpeg -y -framerate 5 -pattern_type glob -i '/data/curtin_mwaeor/data/??????????/cal_qa-flagged-fast/calmetrics_??????????_'$name'.png' -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2" -pix_fmt yuv420p "/data/curtin_mwaeor/data/results-test/calmetrics_${name}-flagged.mp4"
  ffmpeg -y -framerate 5 -pattern_type glob -i '/data/curtin_mwaeor/data/??????????/cal_qa-unflagged-fast/calmetrics_??????????_'$name'.png' -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2" -pix_fmt yuv420p "/data/curtin_mwaeor/data/results-test/calmetrics_${name}-unflagged.mp4"
done
```

### sinfo get available resources

```bash
# name, partition, state, cpus (avail/idle/other/total), tmp(mb), free mem
sinfo --Node --format="%9N %10P %6t %14C %9d %e" | tee sinfo_idle.txt
```

### ionpeel?

```bash
module use /pawsey/mwa/software/python3/modulefiles
module load hyperdrive
module load singularity
export obsid=1254670880
cd ${outdir}
mkdir -p "${obsid}/iono"
hyperdrive srclist-by-beam --metafits "${obsid}/raw/${obsid}.metafits" \
  -n 400 -o ao \
  /astro/mwaeor/dev/calibration/srclist_pumav3_EoR0LoBESv2_fixedEoR1pietro+ForA_phase1+2_edit.txt \
  "${obsid}/iono/ionpeel_model_400.txt"
singularity exec \
  -B /pawsey/mwa:/usr/lib/python3/dist-packages/mwapy/data \
  /pawsey/mwa/singularity/mwa-reduce/mwa-reduce.img \
  cluster \
    "${obsid}/iono/ionpeel_model_400.txt" \
    "${obsid}/iono/ionpeel_clustered_model_400_50.txt" \
    50
singularity exec \
  -B /pawsey/mwa:/usr/lib/python3/dist-packages/mwapy/data \
  /pawsey/mwa/singularity/mwa-reduce/mwa-reduce.img \
  ionpeel \
  "${obsid}/cal/hyp_${obsid}_30l_src4k_8s_80kHz.ms" \
  "${obsid}/iono/ionpeel_clustered_model_400_50.txt" \
  "${obsid}/iono/ionpeel_clustered_soln_400_50.bin"

# applyion wsclean-image.fits ioncorrected.fits clustered-model.txt ionsolutions.bin
singularity exec \
  -B /pawsey/mwa:/usr/lib/python3/dist-packages/mwapy/data \
  /pawsey/mwa/singularity/mwa-reduce/mwa-reduce.img \
  applyion
```

```bash
module load wsclean
export range=125
mkdir -p "${outdir}/${range}XXXXXXX/img_qa-briggs+0.5_4096px"
cd "${outdir}/${range}XXXXXXX/img_qa-briggs+0.5_4096px"
wsclean \
  -weight briggs +0.5 \
  -name wsclean_hyp_${range}XXXXXXX_sub_30l_src4k_8s_80kHz \
  -size 4096 4096 \
  -scale 40asec \
  -pol xx,yy,v \
  -channels-out 68 \
  -niter 0 \
  -parallel-reordering 8 \
  ../../${range}???????/cal/hyp_${range}???????_sub_30l_src4k_8s_80kHz.ms
```

### quick apply MS with different res

```bash
salloc --nodes=1 --mem=350G --gres=gpu:1 --time=8:00:00 --clusters=garrawarla --partition=gpuq --account=mwaeor -c 18
export obsid=1255528600
export time_res=2
export freq_res=80
module load hyperdrive/peel
mkdir /dev/shm/deleteme
cd /dev/shm/deleteme
cp /astro/mwaeor/dev/nfdata/${obsid}/prep/birli_${obsid}_2s_40kHz.uvfits .
cp /astro/mwaeor/dev/nfdata/${obsid}/cal/hyp_soln_${obsid}_30l_src4k.fits .
cp /astro/mwaeor/dev/nfdata/$obsid/raw/${obsid}.metafits .
hyperdrive solutions-apply \
  --time-average "${time_res}s" --freq-average "${freq_res}kHz" \
  --data "${obsid}.metafits" "birli_${obsid}_2s_40kHz.uvfits" \
  --solutions "hyp_soln_${obsid}_30l_src4k.fits" \
  --outputs "hyp_${obsid}_30l_src4k_${time_res}s_${freq_res}kHz.ms"
tar -zcvf "hyp_${obsid}_30l_src4k_${time_res}s_${freq_res}kHz.ms.tar.gz" \
  "hyp_${obsid}_30l_src4k_${time_res}s_${freq_res}kHz.ms"
cp "hyp_${obsid}_30l_src4k_${time_res}s_${freq_res}kHz.ms.tar.gz" /astro/mwaeor/dev/frbhopper/

hyperdrive vis-sub \
  --data "${obsid}.metafits" "hyp_${obsid}_30l_src4k_${time_res}s_${freq_res}kHz.ms" \
  --beam "/pawsey/mwa/mwa_full_embedded_element_pattern.h5" \
  --source-list "/astro/mwaeor/dev/calibration/srclist_pumav3_EoR0LoBESv2_fixedEoR1pietro+ForA_phase1+2_edit.txt" \
  --invert --num-sources 4000 \
  --outputs "hyp_${obsid}_sub_30l_src4k_${time_res}s_${freq_res}kHz.ms"

hyperdrive peel \
  --data "${obsid}.metafits" "hyp_${obsid}_30l_src4k_${time_res}s_${freq_res}kHz.ms" \
  --beam "/pawsey/mwa/mwa_full_embedded_element_pattern.h5" \
  --source-list "/astro/mwaeor/dev/calibration/srclist_pumav3_EoR0LoBESv2_fixedEoR1pietro+ForA_phase1+2_edit.txt" \
  --iono-sub 1000 \
  --sub 4000 \
  --outputs "hyp_${obsid}_ionosub_30l_src4k_${time_res}s_${freq_res}kHz.ms"

```

### wsclean

```bash
salloc --nodes=1 --mem=350G --time=8:00:00 --clusters=garrawarla --partition=workq --account=mwaeor --ntasks=20 --tmp=500G
module load singularity
mkdir /dev/shm/deleteme
cd /dev/shm/deleteme
export obsid=1090012424
# export suffix=""
# export suffix="_sub"
export suffix="_ionosub"
singularity exec --bind /astro --bind /pawsey --bind $PWD:/tmp --writable-tmpfs --pwd /tmp --home $PWD --cleanenv \
    /pawsey/mwa/singularity/casa/casa.img \
    casa -c "importuvfits('/astro/mwaeor/dev/nfdata/${obsid}/cal/hyp_${obsid}${suffix}_30l_src4k_8s_80kHz.uvfits', 'hyp_${obsid}${suffix}_30l_src4k_8s_80kHz.ms')"
export img_weight='briggs -1.0'
export img_size=2048
export img_scale='40asec'
export img_pol='xx,yy,v'
export img_channels_out=24
export img_intervals_out=15
wsclean \
    -name hyp_${obsid}${suffix}_30l_src4k_8s_80kHz \
    -weight ${img_weight} \
    -size ${img_size} ${img_size} \
    -scale ${img_scale} \
    -pol ${img_pol} \
    -channels-out ${img_channels_out} \
    -niter 0 \
    -intervals-out ${img_intervals_out} \
    hyp_${obsid}${suffix}_30l_src4k_8s_80kHz.ms
```

### fits2png

```bash
salloc --nodes=1 --mem=350G --time=8:00:00 --clusters=garrawarla --partition=workq --account=mwaeor --ntasks=20 --tmp=500G
module load singularity
for fits in *109*.fits; do
  singularity exec --cleanenv --home /astro/mwaeor/dev/mplhome \
    /pawsey/mwa/singularity/mwa_qa/mwa_qa_latest.sif python \
    /pawsey/mwa/mwaeor/dev/MWAEoR-Pipeline/templates/thumbnail.py \
      --fits=${fits} \
      --output_name=${fits%.fits}.png \
      --plot_title=${fits%.fits} \
      --vmin_quantile=0.3 --vmax_quantile=0.99
  echo $fits
done
```

### chips combine

```bash
salloc --nodes=1 --mem=300G --time=3:00:00 --clusters=garrawarla --partition=workq --account=mwaeor --ntasks=1 --cpus-per-task=30 --tmp=500G
mkdir /nvmetmp/deleteme
cd /nvmetmp/deleteme


export DATADIR="$PWD"
export INPUTDIR="$PWD/"
export OUTPUTDIR="$PWD/"
export OBSDIR="$PWD/"
export OMP_NUM_THREADS="$(nproc)"

module load chips/cmt singularity
export name="ionosub_30l_src4k_8s_80kHz"

for pol in xx yy; do
  echo -n "" > combine_${pol}_${name}.txt
  for chunk in \
    "eor0high_phase1-128T_p0_0341889a" \
    "eor0high_phase1-128T_p0_118103a9" \
    "eor0high_phase1-128T_p0_36a89d4c" \
    "eor0high_phase1-128T_p0_44de4e7a" \
    "eor0high_phase1-128T_p0_92e7fdb3" \
    "eor0high_phase1-128T_p0_a3dbf06d" \
    "eor0high_phase1-128T_p0_b0d4ed7b" \
    "eor0high_phase1-128T_p0_bac147ae" \
    "eor0high_phase1-128T_p0_ca42f782" \
    "eor0high_phase1-128T_p0_dc1f6b77" \
  ; do
    echo ${pol}.${chunk}_${name} >> combine_${pol}.${name}.txt
    cp /astro/mwaeor/dev/nfdata/${chunk}/ps_metrics/${name}/*_${pol}.${chunk}_${name}.dat .
  done

  combine_data combine_${pol}.${name}.txt 384 "${pol}.${name}" 1
  cp *_${pol}.${name}.dat /astro/mwaeor/dev/nfdata/combined/ps_metrics/

  lssa_fg_simple ${name} 384 80 ${pol} 300 ${name} 0 1
  cp crosspower_${pol}_*.iter.${name}.dat /astro/mwaeor/dev/nfdata/combined/ps_metrics/

  singularity exec --cleanenv -B /nvmetmp --home /astro/mwaeor/dev/mplhome /pawsey/mwa/singularity/ssins/ssins_latest.sif python \
      /pawsey/mwa/mwaeor/dev/MWAEoR-Pipeline/templates/jline_plotchips.py \
      --basedir ./ \
      --polarisation $pol \
      --plot_type 2d \
      --chips_tag $name \
      --title $'crosspower\n'$name \
      --min_power=1e3 \
      --max_power=1e15

  singularity exec --cleanenv -B /nvmetmp --home /astro/mwaeor/dev/mplhome /pawsey/mwa/singularity/ssins/ssins_latest.sif python \
    /pawsey/mwa/mwaeor/dev/MWAEoR-Pipeline/templates/jline_plotchips.py \
    --title $name \
    --basedir ./  \
    --chips_tag $name  \
    --plot_type 1D  \
    --polarisation $pol  \
    --min_power 1E+3 \
    --max_power 1E+15 \
    --lowerfreq 166995000 \
    --umax 300 \
    --N_chan 384 \
    --num_k_edges 80

  cp *.png /astro/mwaeor/dev/nfdata/combined/ps_metrics/
done

```

### manual ionoqa

```bash
salloc --nodes=1 --mem=350G --time=8:00:00 --clusters=garrawarla --partition=workq --account=mwaeor --ntasks=20 --tmp=500G
module load singularity
mkdir /dev/shm/deleteme
cd /dev/shm/deleteme
export obsid=1061320200
cp /astro/mwaeor/dev/nfdata/${obsid}/cal/hyp_peel_${obsid}_ionosub_30l_src4k_8s_80kHz_uv.json offsets.json
singularity exec --bind \$PWD --cleanenv --home /astro/mwaeor/dev/mplhome /pawsey/mwa/singularity/mwa_qa/mwa_qa_latest.sif python \
    /pawsey/mwa/mwaeor/dev/MWAEoR-Pipeline/templates/ionoqa.py \
    --offsets=offsets.json \
    --obsid=${obsid} \
    --output_name="/pawsey/mwa/mwaeor/dev/MWAEoR-Pipeline/test/ionoqa_${obsid}.json"
```