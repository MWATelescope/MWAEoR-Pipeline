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
import os

"""
example:

```bash
salloc --nodes=1 --mem=350G --time=8:00:00 --clusters=garrawarla --partition=workq --account=mwaeor --ntasks=20 --tmp=500G
module load singularity
mkdir /dev/shm/deleteme
cd /dev/shm/deleteme
export obsid=1320408968
cp /astro/mwaeor/dev/nfdata/${obsid}/prep/birli_${obsid}_2s_40kHz.uvfits birli_${obsid}_2s_40kHz.uvfits
cp /astro/mwaeor/dev/nfdata/${obsid}/raw/${obsid}.metafits ${obsid}.metafits
singularity exec --cleanenv --home /astro/mwaeor/dev/mplhome /pawsey/mwa/singularity/mwa_qa/mwa_qa_latest.sif python \
    /astro/mwaeor/dev/MWAEoR-Pipeline/templates/autoplot.py \
    --uvfits=birli_${obsid}_2s_40kHz.uvfits \
    --metafits=${obsid}.metafits \
    --output_name="/astro/mwaeor/dev/MWAEoR-Pipeline/${obsid}_autoplot.png" \
    --plot_title="${obsid}" \
    --transparent \
    --highlight_ants 14 24 25
```
"""


def get_parser():
    parser = ArgumentParser(
        description="Plot autocorrelations in 3D for each pol.")

    parser.add_argument('--uvfits', help='uvfits file to plot')
    parser.add_argument(
        '--metafits', help='metafits file for antenna metadata')
    plot_group = parser.add_argument_group('PLOTTING OPTIONS')
    plot_group.add_argument('--output_name', default=False,
                            help='Name for output plot file, defaults to using the input uvfits name')
    plot_group.add_argument('--plot_title', default=None,
                            help="Optional title for the plot")
    plot_group.add_argument('--plot_dpi', default=100,
                            help="dots per inch for the plot")
    plot_group.add_argument('--plot_style', default="lines",
                            choices=["lines", "imshow"])
    plot_group.add_argument('--log_scale', default=False, action="store_true")
    plot_group.add_argument(
        '--transparent', default=False, action="store_true")

    flag_group = parser.add_argument_group('FLAGGING OPTIONS')
    flag_group.add_argument('--flag_style', default="metafits",
                            choices=["metafits", "weights"])
    flag_group.add_argument('--edge_width', default=80,
                            help="width of the edge channels to be flagged if metafits style [kHz]")
    flag_group.add_argument('--highlight_ants', default=[], nargs='+',
                            help="only plot tiles at these indices")

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
        freq_res = mpr.header['FINECHAN'] * u.kHz
        bandwidth = mpr.header['BANDWDTH'] * u.MHz

        mtd = metafits['TILEDATA']
        metafits_ants = pd.DataFrame(dict([
            (col, list(mtd.data.field(col)[:]))
            for col in ['Antenna', 'Rx', 'Slot', 'Pol', 'TileName', 'Flag']
        ])).drop_duplicates('Antenna').sort_values('Antenna')

    print(metafits_ants, file=sys.stderr)

    uv = UVfits(args.uvfits)

    plt.style.use([astropy_mpl_style, 'dark_background'])

    blt_idxs = np.where(
        uv.ant_1_array - uv.ant_2_array == 0,
    )[0]

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
        quack_scans = ceil((quack_time / int_time).decompose().value)
        edge_chans = ceil(
            ((args.edge_width * u.kHz) / freq_res).decompose().value)
        print(f"quack_scans={quack_scans}, edge_chans={edge_chans}")
        autos2[:quack_scans, :, :, :] = np.nan

        flagged_tile_idxs = np.where(metafits_ants['Flag'].values == 1)[0]
        autos2[:, flagged_tile_idxs, :, :] = np.nan

        cchan_bandwidth = bandwidth / num_cchans
        # number of fine chans per coarse
        num_fchans = ceil((cchan_bandwidth / freq_res).decompose().value)
        num_cchans = nchans // num_fchans
        center_fine_chan = ceil(num_fchans/2)
        chan_flags = np.full((num_cchans, num_fchans), True)
        chan_flags[:, edge_chans:-edge_chans] = False
        chan_flags[:, center_fine_chan] = True
        chan_idxs = np.where(chan_flags.flatten())[0]
        autos2[:, :, chan_idxs, :] = np.nan
    elif args.flag_style == "weights":
        # flag based on uvfits weights
        wghts = data[:, :, :, :, 2]
        autos2[np.where(wghts < 0)] = np.nan

    plt.style.use([astropy_mpl_style, 'dark_background'])

    pols = ["XX", "YY"]
    uv_pol_order = ["XX", "YY", "XY", "YX"]
    fig, axs = plt.subplots(len(pols), sharex=True, sharey=True)
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
            for line in rms_ant_freq:
                ax.plot(freqs, np.log(line) if args.log_scale else line,
                        alpha=(0.2 if args.highlight_ants else 0.5))
            if args.highlight_ants:
                highlight_idxs = [*map(int, args.highlight_ants)]
                for line in rms_ant_freq[highlight_idxs, :]:
                    ax.plot(freqs, np.log(line) if args.log_scale else line,
                            alpha=1, linewidth=4, color='yellow')

            ax.set_ylabel("log(RMS)" if args.log_scale else "RMS")
        ax.grid(None)
        # ax.set_xticks([])
        ax.set_yticks([])
    if args.plot_style == "imshow":
        axs[-1].set_ylabel("Antenna")
    elif args.plot_style == "lines":
        axs[-1].set_xlabel("Freq [MHz]")

    # fig.add_subplot(111, frameon=False)
    # plt.tick_params(labelcolor='none', which='both', top=False, bottom=False, left=False, right=False)
    # plt.grid(None)
    # plt.tight_layout()
    plt.suptitle(args.plot_title)
    print(f"title={args.plot_title}")
    plt.savefig(args.output_name, bbox_inches='tight',
                dpi=args.plot_dpi, transparent=args.transparent)


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
            "--plot_title=${title}",
        ] + shlex.split("${args}"))

    autoplot(args)


if __name__ == '__main__':
    main()
