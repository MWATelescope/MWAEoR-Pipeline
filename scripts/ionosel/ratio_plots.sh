#!/bin/bash
# Make 2D ratio (variant / baseline) and 2D single plots from the pipeline's CHIPS lssa outputs,
# using the same plotchips_all.py (chips_wrappers image) as the pipeline's chipsPlot process.
# usage: ratio_plots.sh GROUPDIR OUTDIR     e.g. /data/curtin_mwaeor/nfdata/eor0low_phase3_xxxx/lssa_fg_general_norfi_b11
set -eu
GROUPDIR=$1; OUTDIR=$2; mkdir -p "$OUTDIR/dats"; rm -f "$OUTDIR"/*.png
IMG=/data/curtin_mwaeor/singularity/d3vnull0-chips_wrappers-latest.img
BASE=ionosub_ssins_30l_src8k_300it_8s_80kHz_sel-none_i1000
for v in none pred4 meas4; do
  d=$GROUPDIR/ionosub_ssins_30l_src8k_300it_8s_80kHz_sel-${v}_i1000
  ln -sf $d/*.dat "$OUTDIR/dats/"
done
GROUP=$(basename $(dirname $GROUPDIR))
# parameters exactly as the pipeline's chipsPlot/chips1d_tsv used them for this group (from the .command.sh of those tasks)
NCHAN=384; MAXU=300; LOWFREQ=138900000.0; CW=80000; BIAS=11
echo "group=$GROUP nchan=$NCHAN maxu=$MAXU lowfreq=$LOWFREQ chanwidth=$CW bias=$BIAS"
COMMON="--basedir $OUTDIR/dats/ --bias_mode $BIAS --lowerfreq_orig $LOWFREQ --chan_width $CW --umax $MAXU --N_chan_orig $NCHAN --kperp_max 0.045 --kperp_min 0.025 --kparra_min 0.11 --num_k_edges 80 --plot_mode png"
cd "$OUTDIR"
for pol in xx yy; do
  for v in none pred4 meas4; do
    tag=${GROUP}_ionosub_ssins_30l_src8k_300it_8s_80kHz_sel-${v}_i1000
    singularity exec -B /data,/localData $IMG plotchips_all.py $COMMON --plot_type 2D --polarisation $pol --chips_tag $tag --chips_tag_label "sel-$v" --outputdir ./ >/dev/null 2>&1 || echo "2D $v $pol failed"
    f=$(ls -t *.png | head -1); mv "$f" "chips2D_${pol}_sel-${v}.png"
    mkdir -p npz; cp 2D_coords_and_power.npz npz/${pol}_sel-${v}.npz
  done
  for v in pred4 meas4; do
    tag=${GROUP}_ionosub_ssins_30l_src8k_300it_8s_80kHz_sel-${v}_i1000
    singularity exec -B /data,/localData $IMG plotchips_all.py $COMMON --plot_type 2D_Ratio --polarisation $pol --chips_tag_one $tag --chips_tag_two ${GROUP}_${BASE} --chips_tag_one_label "sel-$v" --chips_tag_two_label "sel-none" --max_power 2 --min_power 0.5 --outputdir ./ >/dev/null 2>&1 || echo "ratio $v $pol failed"
    f=$(ls -t *.png | head -1); mv "$f" "chips2D_ratio_${pol}_sel-${v}_over_sel-none.png"
  done
done
ls -la *.png
