#!/bin/bash
#SBATCH --account=mwaeor
#SBATCH --partition=work
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:10:00
#SBATCH --job-name=nf-eor1zen-status
# Append one line of product counts to claude-eor/STATUS.md: tail -3 STATUS.md shows the state.
#   scripts/eor1high_status_setonix.sh "<stage>" ["<what comes next>"]
BASE=/scratch/mwaeor/dev/claude-eor
D=$BASE/nfdata
n() { ls $1 2>/dev/null | wc -l; }
echo "$(date -u +%FT%TZ) | ${1:-check} | of 16 obs: prep $(n "$D/*/prep/*.uvfits") soln $(n "$D/*/cal_norfi/hyp_soln_*.fits") cal $(n "$D/*/cal_norfi/hyp_1*_8s_80kHz.uvfits") ionosub $(n "$D/*/cal_norfi/hyp_1*_ionosub_*.uvfits") | of 64 obs grids $(n "$D/1*/ps_metrics_norfi/*/grid_*.tar.gz") | of 4 combined grids $(n "$D/eor*/ps_metrics_norfi/*/grid_*.tar.gz") | lssa crosspower b11 $(n "$D/eor*/lssa_fg_general_norfi_b11/*/crosspower_*") b12 $(n "$D/eor*/lssa_fg_general_norfi_b12/*/crosspower_*") b10 $(n "$D/eor*/lssa_fg_general_norfi_b10/*/crosspower_*") of 4 each | plots $(n "$D/eor*/lssa_fg_general_norfi_b*/*/*.png") | slurm jobs queued or running: $(squeue -u $USER -h | wc -l) | next: ${2:-}" >> $BASE/STATUS.md
