#!/usr/bin/env python

from astropy.io import fits
from astropy.visualization import astropy_mpl_style, simple_norm
import astropy.units as u
import numpy as np
import matplotlib.pyplot as plt
from mwa_qa.read_uvfits import UVfits
from argparse import ArgumentParser
import sys
import pandas as pd
import copy
from math import ceil
import shlex

"""
example:

```bash
salloc --nodes=1 --mem=350G --time=8:00:00 --clusters=garrawarla --partition=workq --account=mwaeor --ntasks=20 --tmp=500G
module load singularity
cd /nvmetmp/
export obsid=1365977896
export prep=birli_${obsid}_2s_40kHz.uvfits
export prep=birli_${obsid}_0.5s_10kHz.uvfits
cp /astro/mwaeor/dev/nfdata/${obsid}/prep/\${prep} .
cp /astro/mwaeor/dev/nfdata/${obsid}/raw/${obsid}.metafits .
singularity exec --cleanenv --home /astro/mwaeor/dev/mplhome -B /nvmetmp /pawsey/mwa/singularity/mwa_qa/mwa_qa_latest.sif python \
    /pawsey/mwa/mwaeor/dev/MWAEoR-Pipeline/templates/autoplot.py \
    --uvfits=\${prep} \
    --metafits=${obsid}.metafits \
    --output_name="/pawsey/mwa/mwaeor/dev/MWAEoR-Pipeline/${obsid}_autoplot.png" \
    --plot_title="${obsid}" \
    --transparent \
    --highlight_ants 14 24 25
singularity exec --cleanenv --home /astro/mwaeor/dev/mplhome -B /nvmetmp /pawsey/mwa/singularity/mwa_qa/mwa_qa_latest.sif python \
    /pawsey/mwa/mwaeor/dev/MWAEoR-Pipeline/templates/autoplot.py \
    --uvfits=\${prep} \
    --metafits=${obsid}.metafits \
    --dpi=300 \
    --output_name="/pawsey/mwa/mwaeor/dev/MWAEoR-Pipeline/${obsid}_autoplot_rx18.png" \
    --plot_title="${obsid} rx18" \
    --sel_ants {128..135}
singularity exec --cleanenv --home /astro/mwaeor/dev/mplhome -B /nvmetmp /pawsey/mwa/singularity/mwa_qa/mwa_qa_latest.sif python \
    /pawsey/mwa/mwaeor/dev/MWAEoR-Pipeline/templates/autoplot.py \
    --uvfits=\${prep} \
    --metafits=${obsid}.metafits \
    --dpi=300 \
    --output_name="/pawsey/mwa/mwaeor/dev/MWAEoR-Pipeline/${obsid}_autoplot_rx06.png" \
    --plot_title="${obsid}" \
    --sel_ants {104..111}
```
"""


def get_parser():
    parser = ArgumentParser(
        description="Plot autocorrelations in 3D for each pol.")

    parser.add_argument('--uvfits', help='uvfits file to plot')
    parser.add_argument(
        '--metafits', help='metafits file for antenna metadata')
    plot_group = parser.add_argument_group('PLOTTING OPTIONS')
    plot_group.add_argument('--output_name', default=None,
                            help='Name for output plot file, defaults to using the input uvfits name')
    plot_group.add_argument('--plot_title', default=None,
                            help="Optional title for the plot")
    plot_group.add_argument('--dpi', default=300, type=int,
                            help="dots per inch for the plot")
    plot_group.add_argument('--plot_style', default="lines",
                            choices=["lines", "imshow", "scatter"])
    plot_group.add_argument('--log_scale', default=False, action="store_true")
    plot_group.add_argument(
        '--transparent', default=False, action="store_true")

    flag_group = parser.add_argument_group('FLAGGING OPTIONS')
    flag_group.add_argument('--flag_style', default="metafits",
                            choices=["metafits", "weights"])
    flag_group.add_argument('--edge_width', default=80, type=float,
                            help="width of the edge channels to be flagged if metafits style [kHz]")
    flag_group.add_argument('--no_flag_centre', default=False, action="store_true",
                            help="flag the center channel")
    flag_group.add_argument('--sel_ants', default=None, nargs='+',
                            help="only plot tiles at these indices")
    flag_group.add_argument('--highlight_ants', default=[], nargs='+',
                            help="only plot tiles at these indices")
    flag_group.add_argument('--sel_times', default=None, nargs='+',
                            help="only plot times at these indices")

    return parser


def split_strip_filter(str):
    return list(filter(None, map(lambda tok: tok.strip(), str.split(','))))


def autoplot(args):
    with fits.open(args.metafits) as metafits:
        metafits.info()

        mpr = metafits['PRIMARY']
        sky_chans = [*map(int, split_strip_filter(mpr.header['CHANNELS']))]
        num_cchans = len(sky_chans)
        print(f"sky_chans={sky_chans} ({num_cchans})")
        int_time = mpr.header['INTTIME'] * u.s
        quack_time = mpr.header['QUACKTIM'] * u.s
        # freq_res = mpr.header['FINECHAN'] * u.kHz
        bandwidth = mpr.header['BANDWDTH'] * u.MHz

        mtd = metafits['TILEDATA']
        metafits_ants = pd.DataFrame(dict([
            (col, list(mtd.data.field(col)[:]))
            for col in ['Antenna', 'Rx', 'Slot', 'Pol', 'TileName', 'Flag']
        ])).drop_duplicates('Antenna').sort_values('Antenna')

    print(metafits_ants.to_string(max_rows=len(metafits_ants)), file=sys.stderr)

    uv = UVfits(args.uvfits)

    freq_res = (uv.channel_width * u.Hz).to(u.kHz)
    first_times = uv.unique_times[0:2]
    if len(first_times) > 1:
        int_time = ((first_times[1] - first_times[0]) * u.day).to(u.s)
        int_time = int_time.round(2)
    print(f"freq_res={freq_res}, int_time={int_time}")

    plt.style.use([astropy_mpl_style, 'dark_background'])

    blt_idxs = np.where(
        uv.ant_1_array - uv.ant_2_array == 0,
    )[0]

    sel_ant_names = metafits_ants['TileName'].values
    sel_ants = metafits_ants['Antenna'].values
    if args.sel_ants is not None:
        sel_ants = list(map(int, args.sel_ants))
        sel_ant_names = sel_ant_names[sel_ants]
        blt_auto_ants = uv.ant_1_array[blt_idxs]
        blt_idxs = blt_idxs[np.where(np.isin(blt_auto_ants, sel_ants))[0]]

    auto_ants = np.unique(uv.ant_1_array[blt_idxs])
    print(f"auto_ants={auto_ants}", file=sys.stderr)

    freqs = (uv.freq_array * u.Hz).to(u.MHz)
    nchans = len(freqs)

    with fits.open(uv.uvfits_path) as hdus:
        vis_hdu = hdus['PRIMARY']
        ncplx = vis_hdu.data.data.shape[-1]
        data = vis_hdu.data.data[blt_idxs, 0, 0, :, :, :].reshape(
            (uv.Ntimes, -1, nchans, uv.Npols, ncplx))

        print(data.shape, file=sys.stderr)

    # reals, imaginaries squard
    reals2 = data[:, :, :, :, 0] ** 2
    imags2 = data[:, :, :, :, 1] ** 2
    # auto amplitudes squared
    autos2 = reals2 + imags2

    if args.flag_style == "metafits":
        # flag quack time, edge channels, center channels
        quack_scans = round((quack_time / int_time).decompose().value)
        edge_chans = round(
            ((float(args.edge_width) * u.kHz) / freq_res).decompose().value)
        print(f"quack_scans={quack_scans}, edge_chans={edge_chans}")
        autos2[:quack_scans, :, :, :] = np.nan

        flagged_tile_idxs = np.array([], dtype=bool)
        if args.sel_ants is None:
            flagged_tile_idxs = np.where(metafits_ants['Flag'].values == 1)[0]
        autos2[:, flagged_tile_idxs, :, :] = np.nan

        cchan_bandwidth = bandwidth / num_cchans
        # number of fine chans per coarse
        num_fchans = round((cchan_bandwidth / freq_res).decompose().value)
        num_cchans = nchans // num_fchans
        center_fine_chan = round(num_fchans/2)
        if edge_chans > 0:
            chan_flags = np.full((num_cchans, num_fchans), True)
            chan_flags[:, edge_chans:-edge_chans] = False
        else:
            chan_flags = np.full((num_cchans, num_fchans), False)
        chan_flags[:, center_fine_chan] = not args.no_flag_centre
        # print("chan_flags:", chan_flags)
        chan_idxs = np.where(chan_flags.flatten())[0]
        # print("chan_idxs: ", chan_idxs)
        autos2[:, :, chan_idxs, :] = np.nan
    elif args.flag_style == "weights":
        # flag based on uvfits weights
        wghts = data[:, :, :, :, 2]
        autos2[np.where(wghts < 0)] = np.nan

    # pols = ["XX", "YY"]
    pols = ["XX", "YY"]
    uv_pol_order = ["XX", "YY", "XY", "YX"]
    fig, axs = plt.subplots(len(pols), sharex=True, sharey=True)
    legend_done = False
    for pol, ax in zip(pols, axs.flatten()):
        pol_idx = uv_pol_order.index(pol)

        # rms across time, for each antenna, freq
        rms_ant_freq = np.sqrt(np.nanmean(autos2[:, :, :, pol_idx], axis=0))
        # rms across time, freq for each antenna
        rms_ant = np.sqrt(np.nanmean(rms_ant_freq ** 2, axis=1))
        # median, std for this pol
        pol_median = np.nanmedian(rms_ant)
        pol_std = np.nanstd(rms_ant)
        pol_low_cutoff = max(0, pol_median - 3 * pol_std)
        pol_high_cutoff = pol_median + 3 * pol_std
        print(
            f"{pol} median={pol_median}, std={pol_std}, high={pol_high_cutoff}, low={pol_low_cutoff}")

        norm = simple_norm(
            rms_ant_freq,
            # 'log' if args.log_scale else 'linear',
            min_cut=pol_low_cutoff,
            max_cut=pol_high_cutoff,
            clip=False
        )

        ax.set_title(f"{pol}")
        if args.plot_style == "imshow":
            cmap = copy.copy(plt.get_cmap('viridis'))
            # cmap.set_over('fuchsia')
            # cmap.set_under('red')
            ax.imshow(rms_ant_freq, interpolation='none',
                      norm=norm, cmap=cmap, aspect=6)
            ax.set_ylabel("Antenna")
        elif args.plot_style == "lines":
            for ant_idx, ant_name, line in zip(sel_ants, sel_ant_names, rms_ant_freq):
                ax.plot(freqs, np.log(line) if args.log_scale else line,
                        alpha=(0.2 if args.highlight_ants else 0.5),
                        label=f"{ant_idx}|{ant_name}")
            if args.highlight_ants:
                highlight_idxs = [*map(int, args.highlight_ants)]
                for line in rms_ant_freq[highlight_idxs, :]:
                    ax.plot(freqs, np.log(line) if args.log_scale else line,
                            alpha=1, linewidth=4, color='yellow')
            ax.set_ylabel("log(RMS)" if args.log_scale else "RMS")
        elif args.plot_style == "scatter":
            for ant_idx, ant_name, line in zip(sel_ants, sel_ant_names, rms_ant_freq):
                ax.scatter(
                    freqs, np.log(line) if args.log_scale else line,
                    label=f"{ant_idx}|{ant_name}",
                    s=1,
                    edgecolor='none',
                    marker='.',
                #    alpha=(0.2 if args.highlight_ants else 0.5),
                )
            if args.highlight_ants:
                highlight_idxs = [*map(int, args.highlight_ants)]
                for line in rms_ant_freq[highlight_idxs, :]:
                    ax.scatter(freqs, np.log(line) if args.log_scale else line,
                               alpha=1, s=4, color='yellow')
            ax.set_ylabel("log(RMS)" if args.log_scale else "RMS")
        ax.grid(None)
        # ax.set_xticks([])
        ax.set_yticks([])
        if args.plot_style == "lines" and not legend_done:
            ax.legend(fontsize=6)
            legend_done = True
    if args.plot_style == "imshow":
        axs[-1].set_ylabel("Antenna")
    elif args.plot_style == "lines":
        axs[-1].set_xlabel("Freq [MHz]")

    # fig.add_subplot(111, frameon=False)
    # plt.tick_params(labelcolor='none', which='both', top=False, bottom=False, left=False, right=False)
    # plt.grid(None)
    # plt.tight_layout()
    plot_title = args.plot_title
    if args.plot_title is None:
        plot_title = f"{args.uvfits}"
    plt.suptitle(plot_title)
    print(f"title={plot_title}")
    output_name = args.output_name
    if args.output_name is None:
        output_name = args.uvfits.replace(".uvfits", "_autoplot.png")
    print(f"output_name={output_name}")
    plt.savefig(output_name, bbox_inches='tight', dpi=args.dpi, transparent=args.transparent)


def main():

    parser = get_parser()

    if len(sys.argv) > 1:
        args = parser.parse_args()
    else:
        # is being called directly from nextflow
        args = parser.parse_args([
            "--uvfits=${uvfits}",
            "--metafits=${metafits}",
            "--output_name=${autoplot}",
        ] + shlex.split("${args}"))

    print("args:", args)
    autoplot(args)


if __name__ == '__main__':
    main()
