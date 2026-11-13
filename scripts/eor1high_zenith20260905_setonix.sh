#!/bin/bash
#SBATCH --account=mwaeor
#SBATCH --partition=work
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=24:00:00
#SBATCH --job-name=nf-eor1zen
# Nextflow head job for the 16 zenith obs of 2026-09-05 (EoR1 high) on Setonix, with the recipe of
# eor1high-ph3-setonix.conf, one run per CHIPS band. Submit from the checkout:
#   sbatch -o /scratch/mwaeor/dev/claude-eor/logs/nf_head.%j.log scripts/eor1high_zenith20260905_setonix.sh
# Finished steps are skipped (storeDir), so it is safe to submit again. One line per stage is
# appended to /scratch/mwaeor/dev/claude-eor/STATUS.md (scripts/eor1high_status_setonix.sh).
set -eu
BASE=/scratch/mwaeor/dev/claude-eor
# Birli measurement sets (ms0: no flags, 8 s / 80 kHz) and metafits of the night, read only
SRC=/scratch/mwaeor/dev/birli-dbuild
HYP=/scratch/mwaeor/dev/hyp-bfi
SRCLIST=srclist_pumav3_EoR0LoBES_EoR1pietro_CenA-GP_2023-11-07_no-shapelets_FnxA-v2-20260905.yaml
SUFFIX=-eor1high-zenith20260905
cd "${SLURM_SUBMIT_DIR:-$PWD}"
module load nextflow/25.10.4 singularity/4.3.5-nohost

# containers named in the conf, pulled once
mkdir -p $BASE/logs $BASE/containers $BASE/sing_cache $BASE/sing_tmp
export SINGULARITY_CACHEDIR=$BASE/sing_cache SINGULARITY_TMPDIR=$BASE/sing_tmp
TAG=bmetrics-92d3a5c-rocm6.3.3-ubuntu24.04
[ -e $BASE/containers/hyperdrive_$TAG.sif ] \
    || singularity pull $BASE/containers/hyperdrive_$TAG.sif docker://mwatelescope/hyperdrive:$TAG
[ -e $BASE/containers/mwa-demo_f55fd52-hyperdrive-$TAG.sif ] \
    || singularity pull $BASE/containers/mwa-demo_f55fd52-hyperdrive-$TAG.sif docker://mwatelescope/mwa-demo:f55fd52-hyperdrive-$TAG

# inputs: link, never copy over anything that is already there
[ -e $BASE/$SRCLIST ] || cp $HYP/$SRCLIST $BASE/
for o in $(grep -v '^#' obsids$SUFFIX.csv); do
    mkdir -p $BASE/nfdata/$o/prep $BASE/nfdata/$o/raw
    ms=$BASE/nfdata/$o/prep/birli_${o}_8s_80kHz_norfi.ms
    [ -e $ms ] || ln -s $SRC/ms0/$o.ms $ms
    [ -e $BASE/nfdata/$o/raw/$o.metafits ] || ln -s $SRC/meta/$o.metafits $BASE/nfdata/$o/raw/$o.metafits
done

id=${SLURM_JOB_ID:-$$}
status() { bash scripts/eor1high_status_setonix.sh "head $id: $1" "$2" || true; }
trap 'status "head job exited with status $?" "read logs/nf_head.$id.log; submit again to continue"' EXIT
# lssa_bias_mode 11, 12, 10 = channels 0-191, 96-287, 192-383 of 384
for b in 11 12 10; do
    status "starting band b$b" "nextflow run zen20260905_b${b}_$id, log logs/nf_b$b.$id.log"
    nextflow -log $BASE/logs/nf_b$b.$id.log run main.nf -profile setonix -c eor1high-ph3-setonix.conf \
        -entry extRaw --obsids_suffix=$SUFFIX --lssa_bias_mode=$b -name zen20260905_b${b}_$id "$@"
done
