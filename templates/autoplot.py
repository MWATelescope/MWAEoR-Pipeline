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
from os.path import realpath
# Install finufft if not available: pip install finufft
try:
    import finufft
except ImportError:
    print("finufft is required for NUFFT. Please install with 'pip install finufft'")
    sys.exit(1)

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
    plot_group.add_argument(
        '--reduction', default="rms", choices=["rms", "mean"],
        help="Reduction method to use for the plot")
    plot_group.add_argument(
        '--ft', default=False, action="store_true",
        help="Fourier transform along frequency axis")
    plot_group.add_argument(
        '--min-delay', default=10, type=int,
        help="max x axis range for Fourier transform plot [ns]")
    plot_group.add_argument(
        '--max-delay', default=1000, type=int,
        help="max x axis range for Fourier transform plot [ns]")
    plot_group.add_argument('--sel_ants', default=None, nargs='+',
                            help="only plot tiles at these indices")
    plot_group.add_argument('--highlight_ants', default=[], nargs='+',
                            help="highlight tiles at these indices")
    plot_group.add_argument('--sel_times', default=None, nargs='+',
                            help="only plot times at these indices")

    flag_group = parser.add_argument_group('FLAGGING OPTIONS')
    flag_group.add_argument('--flag_style', default="metafits",
                            choices=["metafits", "weights"])
    flag_group.add_argument('--edge_width', default=80, type=float,
                            help="width of the edge channels to be flagged if metafits style [kHz]")
    flag_group.add_argument('--no_flag_centre', default=False, action="store_true",
                            help="flag the center channel")

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
    # visibility amplitudes
    amps = np.sqrt(reals2 + imags2)

    if args.flag_style == "metafits":
        # flag quack time, edge channels, center channels
        quack_scans = round((quack_time / int_time).decompose().value)
        edge_chans = round(
            ((float(args.edge_width) * u.kHz) / freq_res).decompose().value)
        print(f"quack_scans={quack_scans}, edge_chans={edge_chans}, no_flag_centre={args.no_flag_centre}")
        amps[:quack_scans, :, :, :] = np.nan

        flagged_tile_idxs = np.array([], dtype=bool)
        if args.sel_ants is None:
            flagged_tile_idxs = np.where(metafits_ants['Flag'].values == 1)[0]
        amps[:, flagged_tile_idxs, :, :] = np.nan

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
        amps[:, :, chan_idxs, :] = np.nan
    elif args.flag_style == "weights":
        # flag based on uvfits weights
        wghts = data[:, :, :, :, 2]
        amps[np.where(wghts <= 0)] = np.nan

    # pols = ["XX", "YY"]
    pols = ["XX", "YY"]
    uv_pol_order = ["XX", "YY", "XY", "YX"]
    fig, axs = plt.subplots(len(pols), sharex=True, sharey=True, figsize=(36, 36))
    legend_done = False
    for pol, ax in zip(pols, axs.flatten()):
        pol_idx = uv_pol_order.index(pol)

        # reduce across time, for each antenna, freq
        if args.reduction == "rms":
            rdx_ant_freq = np.sqrt(np.nanmean(amps[:, :, :, pol_idx] ** 2, axis=0))
        elif args.reduction == "mean":
            rdx_ant_freq = np.nanmean(amps[:, :, :, pol_idx], axis=0)

        # Apply Fourier transform along frequency axis if requested
        x_label = "Freq [MHz]"
        x_values = freqs
        if args.ft:
            # Use original complex visibilities for proper delay transform
            freq_hz = freqs.to(u.Hz).value.astype(np.float64)
            min_delay_ns, max_delay_ns = float(args.min_delay), float(args.max_delay)

            nants = rdx_ant_freq.shape[0]
            nfreqs = rdx_ant_freq.shape[1]

            # Use amplitude data with zero-padding for better delay resolution
            print(f"Computing delay transform: {min_delay_ns}-{max_delay_ns} ns", file=sys.stderr)

            # Zero-pad to improve delay resolution - target ~1 ns resolution
            freq_res_hz = freq_hz[1] - freq_hz[0]
            target_delay_res_ns = 1.0
            pad_factor = max(1, int(1e9 / (target_delay_res_ns * freq_res_hz * nfreqs)))
            nfreqs_padded = nfreqs * pad_factor

            print(f"Zero-padding: {nfreqs} -> {nfreqs_padded} for {1e9/(nfreqs_padded*freq_res_hz):.1f} ns resolution", file=sys.stderr)

            rdx_ant_delay = np.full((nants, nfreqs_padded), np.nan, dtype=np.float64)

            for ant in range(nants):
                # Use the original amplitude data that's already been reduced
                data_freq = rdx_ant_freq[ant].copy()
                valid_mask = ~np.isnan(data_freq)
                if not np.any(valid_mask):
                    continue

                # Replace NaN with 0 for FFT
                data_freq[~valid_mask] = 0.0

                # Apply window to reduce spectral leakage
                window = np.blackman(len(data_freq))
                data_windowed = data_freq * window

                # Zero-pad the data
                data_padded = np.zeros(nfreqs_padded)
                data_padded[:nfreqs] = data_windowed

                # FFT to delay domain
                delay_spectrum = np.fft.fft(data_padded)
                rdx_ant_delay[ant] = np.abs(delay_spectrum)

            # Calculate delay axis with improved resolution
            delay_res_s = 1.0 / (nfreqs_padded * freq_res_hz)
            delays_s = np.arange(nfreqs_padded) * delay_res_s
            delays_ns = delays_s * 1e9

            print(f"Actual delay range: 0 to {delays_ns.max():.1f} ns, res={delay_res_s*1e9:.2f} ns", file=sys.stderr)
            print(f"Channel BW = {freq_res_hz/1e6:.3f} MHz -> delay period = {1e9/freq_res_hz:.1f} ns", file=sys.stderr)
            print(f"Delay data shape: {rdx_ant_delay.shape}, has finite data: {np.any(np.isfinite(rdx_ant_delay))}", file=sys.stderr)
            print(f"Data range: {np.nanmin(rdx_ant_delay):.2e} to {np.nanmax(rdx_ant_delay):.2e}", file=sys.stderr)

            # Filter to requested delay range
            delay_mask = (delays_ns >= min_delay_ns) & (delays_ns <= max_delay_ns)
            print(f"Requested range {min_delay_ns}-{max_delay_ns} ns has {np.sum(delay_mask)} points", file=sys.stderr)

            if np.any(delay_mask) and np.sum(delay_mask) > 5:
                rdx_ant_delay = rdx_ant_delay[:, delay_mask]
                delays_ns = delays_ns[delay_mask]
                print(f"Using filtered range: {delays_ns.min():.1f} to {delays_ns.max():.1f} ns", file=sys.stderr)
            else:
                # Show a broader range if requested range is too narrow
                if delays_ns.max() < max_delay_ns:
                    # Show all computed delays
                    print(f"Showing full computed range: 0 to {delays_ns.max():.1f} ns", file=sys.stderr)
                else:
                    # Show first part of delays up to max_delay_ns
                    subset_end = np.where(delays_ns <= max_delay_ns)[0]
                    if len(subset_end) > 0:
                        subset_end = subset_end[-1] + 1
                        rdx_ant_delay = rdx_ant_delay[:, :subset_end]
                        delays_ns = delays_ns[:subset_end]
                        print(f"Showing subset: 0 to {delays_ns.max():.1f} ns", file=sys.stderr)
                    else:
                        print(f"No delays <= {max_delay_ns} ns found!", file=sys.stderr)

            # For NUFFT, we already have the desired delay range
            rdx_ant_freq = rdx_ant_delay
            x_values = delays_ns
            x_label = "Delay [ns]"

            print(f"NUFFT result: {x_values.min():.1f} to {x_values.max():.1f} ns, {len(x_values)} points", file=sys.stderr)

            ax.set_xlim(min_delay_ns, max_delay_ns)




        ax.set_title(f"{pol}")
        quantity = "amps"
        if args.reduction != "none":
            quantity = f"{args.reduction} {quantity}"
        if args.log_scale:
            quantity = f"log({quantity})"
        if args.plot_style == "imshow" or args.ft:
            # for normalization, reduce across time, freq for each antenna
            rms_ant = np.sqrt(np.nanmean(rdx_ant_freq ** 2, axis=1))
            # median, std for this pol
            pol_median = np.nanmedian(rms_ant)
            pol_std = np.nanstd(rms_ant)
            pol_low_cutoff = max(0, pol_median - 3 * pol_std)
            pol_high_cutoff = pol_median + 3 * pol_std
            print(
                f"{pol} median={pol_median}, std={pol_std}, high={pol_high_cutoff}, low={pol_low_cutoff}")

            norm = simple_norm(
                rdx_ant_freq,
                # 'log' if args.log_scale else 'linear',
                min_cut=pol_low_cutoff,
                max_cut=pol_high_cutoff,
                clip=False
            )

            cmap = copy.copy(plt.get_cmap('viridis'))
            # cmap.set_over('fuchsia')
            # cmap.set_under('red')

            # Adjust aspect ratio for delay vs frequency plots
            if args.ft:
                aspect_ratio = len(x_values) / len(sel_ants) * 0.5  # Better for delay plots
                extent = [x_values.min(), x_values.max(), 0, len(sel_ants)]
            else:
                aspect_ratio = 6  # Original aspect for frequency plots
                extent = None

            im = ax.imshow(rdx_ant_freq, interpolation='none',
                          norm=norm, cmap=cmap, aspect=aspect_ratio,
                          extent=extent, origin='lower')
            ax.set_ylabel("Antenna")
        elif args.plot_style == "lines":
            for ant_idx, ant_name, line in zip(sel_ants, sel_ant_names, rdx_ant_freq):
                ax.plot(x_values, np.log(line) if args.log_scale else line,
                        alpha=(0.2 if args.highlight_ants else 0.5),
                        label=f"{ant_idx}|{ant_name}")
            if args.highlight_ants:
                highlight_idxs = [*map(int, args.highlight_ants)]
                for line in rdx_ant_freq[highlight_idxs, :]:
                    ax.plot(x_values, np.log(line) if args.log_scale else line,
                            alpha=1, linewidth=4, color='yellow')
            ax.set_ylabel(quantity)
        elif args.plot_style == "scatter":
            for ant_idx, ant_name, line in zip(sel_ants, sel_ant_names, rdx_ant_freq):
                ax.scatter(
                    x_values, np.log(line) if args.log_scale else line,
                    label=f"{ant_idx}|{ant_name}",
                    s=1,
                    edgecolor='none',
                    marker='.',
                #    alpha=(0.2 if args.highlight_ants else 0.5),
                )
            if args.highlight_ants:
                highlight_idxs = [*map(int, args.highlight_ants)]
                for line in rdx_ant_freq[highlight_idxs, :]:
                    ax.scatter(x_values, np.log(line) if args.log_scale else line,
                               alpha=1, s=4, color='yellow')
            ax.set_ylabel(quantity)
        ax.grid(None)
        # ax.set_xticks([])
        # ax.set_yticks([])
        if args.plot_style != "imshow" and not legend_done:
            ax.legend(fontsize=6)
            legend_done = True
    if args.plot_style == "imshow":
        axs[-1].set_ylabel("Antenna")
    elif args.plot_style in ["lines", "scatter"]:
        axs[-1].set_xlabel(x_label)

    # fig.add_subplot(111, frameon=False)
    # plt.tick_params(labelcolor='none', which='both', top=False, bottom=False, left=False, right=False)
    # plt.grid(None)
    # plt.tight_layout()
    output_name = args.output_name
    if args.output_name is None:
        suffix = f"{args.plot_style}"
        if args.ft:
            suffix += "_ft"
        output_name = args.uvfits.replace(".uvfits", f"_autoplot_{suffix}.png")

    plot_title = args.plot_title
    if args.plot_title is None:
        plot_title = output_name.split('/')[-1].replace('.png', '')
    plt.suptitle(plot_title)
    plt.savefig(output_name, bbox_inches='tight', dpi=args.dpi, transparent=args.transparent)
    print(realpath(output_name))


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
