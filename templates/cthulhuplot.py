#!/usr/bin/env python

import sys

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from scipy.interpolate import griddata
from collections import OrderedDict
import astropy.units as u
from astropy.time import Time, TimeDelta
from astropy.coordinates import EarthLocation, SkyCoord, AltAz, ICRS
from astropy.wcs import WCS
from scipy import signal

import json
import yaml
from yaml import CSafeLoader as SafeLoader
from argparse import ArgumentParser
import shlex
from math import pi, radians
import os
import pandas as pd

from cthulhu.reconstruct import Obsid
from cthulhu.plot_tools import setup_subplot, plot_tec, plot_vector_arrows, generate_diagnostic_figure
from cthulhu.rts_log_tools import lm2radec, MWA_LOCATION
from pyGrad2Surf.g2s import g2s, g2s_dirichlet


"""
example:

```bash
salloc --nodes=1 --mem=350G --time=8:00:00 --clusters=garrawarla --partition=workq --account=mwaeor --ntasks=20 --tmp=500G
module load singularity
mkdir /dev/shm/deleteme
cd /dev/shm/deleteme
export obsid=1065539256
export obsid=1065542792
export obsid=1095694800
cp /astro/mwaeor/dev/nfdata/${obsid}/cal/${obsid}_reduced_n8000.yaml srclist.yaml
# cp /astro/mwaeor/dev/nfdata/${obsid}/cal/hyp_peel_${obsid}_ionosub_30l_src4k_8s_80kHz_uv.json offsets.json
cp /astro/mwaeor/dev/nfdata/${obsid}/cal/hyp_peel_${obsid}_ionosub_30l_src4k_8s_80kHz_i4000_uv.json offsets.json
eval singularity exec --bind \$PWD --cleanenv --home /astro/mwaeor/dev/mplhome /pawsey/mwa/singularity/mwa_qa/mwa_qa_latest.sif python \
    /pawsey/mwa/mwaeor/dev/MWAEoR-Pipeline/templates/cthulhuplot.py \
    --srclist=srclist.yaml \
    --offsets=offsets.json \
    --obsid=${obsid} \
    --plot="cthulhuplot_${obsid}.png" \
    --tec="tec_${obsid}.png" \
    --plot-altaz 90 0 \
    --offset-type arrow \
    --obs-frequency 154000000 \
    --obs-radec 0 -27 \
    --show-axes
cp *.png /pawsey/mwa/mwaeor/dev/MWAEoR-Pipeline/test/
```
"""


def get_parser():
    parser = ArgumentParser(
        description="render TEC plot for offsets.")

    parser.add_argument('--srclist', help='source list file (yaml)')
    parser.add_argument('--offsets', help='offset file (json|yaml)')
    parser.add_argument('--obsid', default=0, type=int)
    parser.add_argument('--time-res', default=8, type=int,
                        help='time resolution (seconds)')
    parser.add_argument('--time_offset', default=0, type=int,
                        help='offset from obsid to first timestep (seconds)')
    plot_group = parser.add_argument_group('OUTPUT OPTIONS')
    plot_group.add_argument('--plot',
                            help='Name of output diagnostic plot file', default=None)
    plot_group.add_argument('--tec',
                            help='Name of output tec file', default=None)
    plot_group.add_argument('--json',
                            help='Name of output json file', default=None)
    plot_group.add_argument('--csv',
                            help='Name of output csv file', default=None)
    plot_group.add_argument('--average', action='store_true', default=False,
                            help="Average offsets together")

    tec_group = parser.add_argument_group('TEC OPTIONS')
    tec_group.add_argument('--resolution', default=1024, type=int,
                           help="Number of pixels in the TEC map")
    tec_group.add_argument('--dpi', default=50, type=float,
                           help="dots per inch for the TEC map")
    tec_group.add_argument('--cell-size', default=240/60/60, type=float,
                           help="pixel size in degrees for the TEC map")
    # tec_group.add_argument('--max-offset', default=1.2e-4, type=float,
    #                        help="max alpha/beta for TEC plot")
    tec_group.add_argument('--interpolation', default='cubic',
                           help="interpolation for alpha/beta in TEC plot")
    tec_group.add_argument('--obs-radec', default=[], nargs=2, type=float,
                           help="ra/dec for centre in the observation")
    tec_group.add_argument('--plot-radec', default=[], nargs=2, type=float,
                           help="ra/dec for centre of TEC plot")
    tec_group.add_argument('--plot-altaz', default=[], nargs=2, type=float,
                           help="altitude / azimuth angle for centre of TEC plot")
    tec_group.add_argument('--blur', default=None,
                           help="size of gaussian [degrees] to blur alpha/beta TEC plot")
    tec_group.add_argument('--obs-frequency', default=None, type=int,
                           help="frequency [Hz] of the observation for Beam FWHM calculation")
    tec_group.add_argument('--plot-frequency', default=200, type=int,
                           help="frequency [Hz] for TEC plot offsets")
    tec_group.add_argument('--offset-type', default="none", choices=["scatter", "arrow", "none"],
                           help="offset plot type")
    tec_group.add_argument('--show-axes', default=False, action="store_true",
                           help="show axes on plot")

    return parser


# MWA_LOCATION = EarthLocation(lat=-26.756528*u.deg, lon=116.670810*u.deg)

def main():

    parser = get_parser()

    if len(sys.argv) > 1:
        args = parser.parse_args()
    else:
        # is being called directly from nextflow
        args = parser.parse_args([
            # "--srclist=${srclist}",
            "--offsets=${offsets}",
            "--obsid=${obsid}",
            "--plot=${plot}",
            "--tec=${tec}",
            "--csv=${csv}",
            "--json=${json}",
        ] + shlex.split(extra) if (extra := "${extra}") else [])

    print(f"{vars(args)=}")

    # read source list
    # with open(args.srclist, "r") as h:
    #     srclist = yaml.load(h, Loader=SafeLoader)
    #     # srclist = json.load(h)

    # calculate offsets at a given freq
    obs_wavelength = 299792458 / args.obs_frequency
    fwhm = np.degrees(obs_wavelength / 4.5)
    print(f"{fwhm=}")
    obs_wavelength_2 = obs_wavelength**2

    obs_ra, obs_dec = args.obs_radec

    # plot_wavelength = 299792458 / args.plot_frequency
    # plot_wavelength_2 = plot_wavelength**2

    # read iono constants
    ras = []
    decs = []
    src_names = []
    _, ext = os.path.splitext(args.offsets)
    if ext == ".json":
        with open(args.offsets, "r") as h:
            iono_consts = json.load(h)
        alphas = []
        betas = []
        for src_name, consts in iono_consts.items():
            pos = consts["weighted_catalogue_pos_j2000"]
            ras.append(pos["ra"])
            decs.append(pos["dec"])
            alphas.append(consts["alphas"])
            betas.append(consts["betas"])
            src_names.append(src_name)

        alphas = np.array(alphas)
        betas = np.array(betas)
        n_times = min(alphas.shape[1], betas.shape[1])
        l_shifts = alphas * obs_wavelength_2
        m_shifts = betas * obs_wavelength_2
        n_times = min(l_shifts.shape[1], m_shifts.shape[1])
        ras = np.array(ras)
        decs = np.array(decs)
        ras[np.where(ras > 180)] -= 360
        ra_shifts = None
        dec_shifts = None

    elif ext == ".yaml":
        with open(args.offsets, "r") as h:
            iono_consts = yaml.load(h, Loader=SafeLoader)
        l_shifts = []
        m_shifts = []
        ra_shifts = []
        dec_shifts = []
        for src in iono_consts['sources'].values():
            src_names.append(src['name'])
            ras.append(src['ra'])
            decs.append(src['dec'])
            # l_shifts.append(src['l_shifts'])
            # m_shifts.append(src['m_shifts'])
            ra_shifts.append(src['ra_shifts'])
            dec_shifts.append(src['dec_shifts'])
        # l_shifts = np.array(l_shifts)
        # m_shifts = np.array(m_shifts)
        ras = np.array(ras)
        decs = np.array(decs)
        ra_shifts = np.array(ra_shifts)
        dec_shifts = np.array(dec_shifts)
        n_times = min(ra_shifts.shape[1], dec_shifts.shape[1])
        ra_shifts, dec_shifts = -1 * ra_shifts, -1 * dec_shifts
        ras[np.where(ras > 180)] -= 360

    start_time = Time(args.obsid, format="gps", scale="utc", location=MWA_LOCATION)
    centroid_time = start_time + TimeDelta(args.time_offset + (n_times * args.time_res / 2), format="sec")

    if args.plot_altaz:
        centroid_radec = SkyCoord(az=args.plot_altaz[1], alt=args.plot_altaz[0], unit=(u.deg, u.deg), frame="altaz",
                    obstime=centroid_time, location=MWA_LOCATION).transform_to(ICRS)
        # aaframe = AltAz(location=MWA_LOCATION, obstime=centroid_time)
        # centroid_radec = SkyCoord(radians(args.plot_altaz[1]), radians(args.plot_altaz[0]), unit="rad", frame=aaframe).transform_to(ICRS)
        centroid_ra, centroid_dec = centroid_radec.ra.deg, centroid_radec.dec.deg
    elif args.plot_radec:
        centroid_ra, centroid_dec = args.plot_radec
    else:
        centroid_ra, centroid_dec = (np.mean(ras), np.mean(decs))

    if ra_shifts is None or dec_shifts is None:
        ra_shifts = np.zeros_like(l_shifts)
        dec_shifts = np.zeros_like(m_shifts)
        for i in range(n_times):
            ra_shifts_ts, dec_shifts_ts = lm2radec(ras, decs, np.rad2deg(l_shifts[:, i]), np.rad2deg(m_shifts[:, i]), {"primary_beam_pointing_centre": [obs_ra, obs_dec]})
            ra_shifts[:, i] = np.array(ra_shifts_ts)
            dec_shifts[:, i] = np.array(dec_shifts_ts)
        # ra_shifts = np.array(ra_shifts)
        # dec_shifts = np.array(dec_shifts)

    # print(l_shifts[0:5])
    av_ra_shifts = np.mean(ra_shifts, axis=1)
    av_dec_shifts = np.mean(dec_shifts, axis=1)
    av_total_shifts = np.sqrt(av_ra_shifts**2 + av_dec_shifts**2)
    max_total_shift_idx = np.argmax(av_total_shifts)
    print(av_ra_shifts[0:5])
    o = Obsid((ras, decs, av_ra_shifts, av_dec_shifts), args.obsid, radius=fwhm/2)
    o.ra_centre = obs_ra
    o.dec_centre = obs_dec
    o.pca()
    o.obsid_metric()

    print(len(src_names), ras.shape, ra_shifts.shape)

    if args.json:
        print(f"saving json to {args.json}")
        data = OrderedDict()
        data['PCA_EIGENVALUE'] = o.pca_variance[0]
        p = data['PCA_EIGENVALUE']
        data['MED_ABS_OFFSET'] = o.s3[0]
        m = data['MED_ABS_OFFSET']
        data['QA'] = 25 * m
        if p > 0.6:
            data['QA'] += 64 * p * ( p - 0.6 )

        data['MAX_SHIFT_SRC'] = str(src_names[max_total_shift_idx])
        data['MAX_SHIFT_RA'] = ras[max_total_shift_idx]
        data['MAX_SHIFT_DEC'] = decs[max_total_shift_idx]
        data['N_TIMES'] = n_times
        data['N_SRCS'] = len(src_names)
        with open(args.json, 'w') as f:
            json.dump(data, f)

    if args.csv:
        print(f"saving csv to {args.csv}")
        df = pd.DataFrame({
            'name': src_names,
            'ra': ras,
            'dec': decs,
            'avg_ra_shift': av_ra_shifts,
            'avg_dec_shift': av_dec_shifts
        })
        df.set_index('name', inplace=True)
        df.to_csv(args.csv)

    if args.plot or args.tec:
        print(f"saving plots to {args.plot=} {args.tec=}")

        if args.blur:
            blur_pixels = args.blur / args.cell_size
            kernel_size = int(min(args.resolution / 2, blur_pixels * 8))
            kernel = signal.windows.gaussian(kernel_size, blur_pixels)
            kernel = np.outer(kernel, kernel)

        x_axis = np.arange(args.resolution)
        y_axis = np.arange(args.resolution)
        grid_x, grid_y = np.meshgrid(x_axis, y_axis)

        if args.average:
            # plot the average of shifts over time

            wcs = WCS({
                'NAXIS': 2, 'NAXIS1': args.resolution, 'NAXIS2': args.resolution,
                'CTYPE1': 'RA---SIN', 'CRPIX1': float(args.resolution/2 + 1), 'CRVAL1': centroid_ra, 'CDELT1': -args.cell_size, 'CUNIT1': 'deg',
                'CTYPE2': 'DEC--SIN', 'CRPIX2': float(args.resolution/2 + 1), 'CRVAL2': centroid_dec, 'CDELT2': args.cell_size, 'CUNIT2': 'deg',
            })
            pixels_x, pixels_y = wcs.all_world2pix(ras, decs, 0)
            visible_idxs = np.where((pixels_x > 1) & (pixels_x < args.resolution - 1) & (pixels_y > 1) & (pixels_y < args.resolution - 1))[0]
            pixels = np.vstack((pixels_x, pixels_y)).T
            o.reconstruct_tec(resolution=args.resolution, frequency=args.plot_frequency/1e6, interp_method=args.interpolation)
            o.tec_power_spectrum()
            o.metrics.append(
                [av_total_shifts[max_total_shift_idx], f"max total shift"])
            o.metric_weights.append(0)
            o.metrics.append([0, f"{src_names[max_total_shift_idx]}"])
            o.metric_weights.append(0)
            o.metrics.append([ras[max_total_shift_idx], "max RA"])
            o.metric_weights.append(0)
            o.metrics.append([decs[max_total_shift_idx], "max Dec"])
            o.metric_weights.append(0)
            path = args.plot
            generate_diagnostic_figure(
                o, filename=os.path.basename(path),
                directory=os.path.dirname(path) or '.',
                overwrite=True
            )
        else:
            # plot each time separately
            for t in range(n_times):
                local_time = start_time + TimeDelta(args.time_offset + (t * args.time_res), format="sec")
                obsid = int(local_time.gps)
                local_utc = local_time.to_value(format='iso')

                if args.plot_altaz:
                    aaframe = AltAz(location=MWA_LOCATION, obstime=local_time)
                    radec = SkyCoord(radians(args.plot_altaz[1]), radians(args.plot_altaz[0]), unit="rad", frame=aaframe).transform_to(ICRS)
                    ra, dec = radec.ra.deg, radec.dec.deg
                elif args.plot_radec:
                    ra, dec = args.plot_radec
                else:
                    ra, dec = (np.mean(ras), np.mean(decs))

                print(f"{obsid=} {local_utc=} {ra=} {dec=}")

                wcs = WCS({
                    'NAXIS': 2, 'NAXIS1': args.resolution, 'NAXIS2': args.resolution,
                    'CTYPE1': 'RA---SIN', 'CRPIX1': float(args.resolution/2 + 1), 'CRVAL1': ra, 'CDELT1': -args.cell_size, 'CUNIT1': 'deg',
                    'CTYPE2': 'DEC--SIN', 'CRPIX2': float(args.resolution/2 + 1), 'CRVAL2': dec, 'CDELT2': args.cell_size, 'CUNIT2': 'deg',
                })
                pixels_x, pixels_y = wcs.all_world2pix(ras, decs, 0)
                visible_idxs = np.where((pixels_x > 1) & (pixels_x < args.resolution - 1) & (pixels_y > 1) & (pixels_y < args.resolution - 1))[0]
                pixels = np.vstack((pixels_x, pixels_y)).T
                assert not np.any(np.isnan(pixels))

                t_ra_shifts = ra_shifts[:, t]
                t_dec_shifts = dec_shifts[:, t]
                t_l_shifts = l_shifts[:, t]
                t_m_shifts = m_shifts[:, t]
                assert len(t_ra_shifts) == len(t_dec_shifts) == len(t_l_shifts) == len(t_m_shifts) == len(ras) == len(decs)
                t_total_shifts = np.sqrt(t_ra_shifts**2 + t_dec_shifts**2)
                # t_lm_shifts = np.sqrt(t_l_shifts**2 + t_m_shifts**2)

                max_total_shift_idx = np.argmax(t_total_shifts)

                t_shift_median = np.median(t_total_shifts)
                t_shift_std = np.std(t_total_shifts)
                t_shift_cutoff = t_shift_median + 3 * t_shift_std
                # plt.hist(t_total_shifts, bins=100)
                # plt.axvline(t_shift_cutoff, color='k', linestyle='dashed', linewidth=1)
                # plt.savefig(f'hist_shifts_{obsid}.png')
                # plt.clf()

                pixel_offsets_x, pixel_offsets_y = wcs.all_world2pix(ras + t_ra_shifts, decs + t_dec_shifts, 0)
                pixel_offsets_x -= pixels_x
                pixel_offsets_y -= pixels_y

                o = Obsid((ras, decs, t_ra_shifts, t_dec_shifts), obsid, radius=fwhm/2)
                o.ra_centre = obs_ra
                o.dec_centre = obs_dec

                sane_idxs = np.where(t_total_shifts < t_shift_cutoff)[0]
                print(f"{len(visible_idxs)=} {len(sane_idxs)=} {len(t_total_shifts)=}")
                o.fra = o.ra[sane_idxs]
                o.fdec = o.dec[sane_idxs]
                o.fra_shifts = o.ra_shifts[sane_idxs]
                o.fdec_shifts = o.dec_shifts[sane_idxs]

                # left_radec = wcs.all_pix2world(0, args.resolution/2, 0)
                # o.radius = np.sqrt((left_radec[0] - ra)**2 + (left_radec[1] - dec)**2)
                # print(f"{left_radec=} {radius=} {o.radius=}")

                o.pca()
                o.obsid_metric()
                # o.reconstruct_tec(resolution=args.resolution, frequency=args.plot_frequency/1e6, interp_method=args.interpolation)
                o.reconstruct_tec(frequency=args.plot_frequency/1e6, interp_method=args.interpolation)

                # o.grid_dra = griddata(pixels[sane_idxs], o.fra_shifts, (grid_x, grid_y), method=args.interpolation, fill_value=0)
                # o.grid_ddec = griddata(pixels[sane_idxs], o.fdec_shifts, (grid_x, grid_y), method=args.interpolation, fill_value=0)
                # if args.blur:
                #     o.grid_dra = signal.fftconvolve(o.grid_dra, kernel, mode='same')
                #     o.grid_ddec = signal.fftconvolve(o.grid_ddec, kernel, mode='same')
                # o.tec = g2s(x_axis, y_axis, o.grid_dra, o.grid_ddec)
                # o.tec -= np.min(o.tec)

                o.tec_power_spectrum()
                o.metrics.append(
                    [t_total_shifts[max_total_shift_idx], f"max total shift"])
                o.metric_weights.append(0)
                o.metrics.append([0, f"{src_names[max_total_shift_idx]}"])
                o.metric_weights.append(0)
                o.metrics.append([ras[max_total_shift_idx], "max RA"])
                o.metric_weights.append(0)
                o.metrics.append([decs[max_total_shift_idx], "max Dec"])
                o.metric_weights.append(0)
                output_name, ext = os.path.splitext(args.plot)
                path = f"{output_name}-t{t:04d}{ext}"

                plot_idxs = np.intersect1d(visible_idxs, sane_idxs)
                o.fra = o.ra[plot_idxs]
                o.fdec = o.dec[plot_idxs]
                o.fra_shifts = o.ra_shifts[plot_idxs]
                o.fdec_shifts = o.dec_shifts[plot_idxs]
                t_l_shifts = t_l_shifts[plot_idxs]
                t_m_shifts = t_m_shifts[plot_idxs]

                generate_diagnostic_figure(
                    o, filename=os.path.basename(path),
                    directory=os.path.dirname(path) or '.',
                    overwrite=True
                )

                if not args.tec:
                    continue

                assert not np.any(np.isnan(t_ra_shifts))
                assert not np.any(np.isnan(t_dec_shifts))

                with plt.style.context('dark_background'):
                    img_fig = plt.figure(dpi=args.dpi, figsize=(args.resolution/args.dpi, args.resolution/args.dpi))
                    ax = img_fig.add_subplot(1, 1, 1, projection=wcs)
                    ra_bounds, dec_bounds = wcs.all_pix2world(
                        np.array([0, 0, 0, 1, 1, 1, 2, 2, 2]) * args.resolution / 2,
                        np.array([0, 1, 2, 0, 1, 2, 0, 1, 2]) * args.resolution / 2,
                        0
                    )
                    ax.scatter(ra_bounds, dec_bounds, s=0.000001, c='k', marker='o', transform=ax.get_transform('world'))

                    ax.set_ylabel("", visible=False)
                    ax.set_xlabel("", visible=False)
                    for spine in ax.spines.values():
                        spine.set_visible(False)
                    if args.show_axes:
                        ax.tick_params(axis="y", direction="in", pad=-50, horizontalalighment="left", verticalalignment="bottom")
                        ax.tick_params(axis="x", direction="in", pad=-
                                    20, verticalalignment="bottom")
                        for coord in ax.coords:
                            coord.set_ticklabel_visible(True)
                        plt.setp(ax.get_xticklabels()[0], visible=False)
                        plt.setp(ax.get_yticklabels()[-1], visible=False)
                        ax.coords.grid(True, color='gray', ls='--')
                    else:
                        for coord in ax.coords:
                            coord.set_ticks_visible(False)
                            coord.set_ticklabel_visible(False)
                            coord.set_axislabel('')

                    title = f"GPS {obsid} UTC {local_utc}"
                    ax.set_title(title, y=0.95, fontdict={'verticalalignment': 'top', 'family': 'monospace'})

                    # tec_vmin, tec_vmax, tec_sigma = np.min(tec_grid), np.max(tec_grid), np.std(tec_grid)
                    # tec_v = max(abs(tec_vmin), abs(tec_vmax))
                    # tec_pix_extent = [
                    #     *wcs.all_pix2world(o.tec_extent[0], o.tec_extent[2], 0),
                    #     *wcs.all_pix2world(o.tec_extent[1], o.tec_extent[3], 0)
                    # ]
                    # tec_pix_extent = [
                    #     tec_pix_extent[0], tec_pix_extent[2],
                    #     tec_pix_extent[1], tec_pix_extent[3]
                    # ]
                    # print(f"{tec_pix_extent=}")
                    # im = ax.imshow(o.tec, cmap='plasma', extent=tec_pix_extent)
                    # angles = (np.arctan2(o.fra_shifts, o.fdec_shifts) + np.pi) / (2 * np.pi)
                    angles = (np.arctan2(pixel_offsets_y, pixel_offsets_x) + np.pi) / (2 * np.pi)
                    if args.offset_type == "scatter":
                        # ax.scatter(pixels_x[plot_idxs], pixels_y[plot_idxs], c=angles[plot_idxs], cmap="rainbow", s=t_total_shifts[plot_idxs] * 1e5)
                        ax.scatter(o.fra, o.fdec, c=cmap(angles), s=t_total_shifts[plot_idxs] * 1e5)
                    elif args.offset_type == "arrow":
                        cmap = plt.cm.hsv
                        scale = args.resolution / 100
                        ax.quiver(pixels_x[plot_idxs], pixels_y[plot_idxs], np.arctan(scale * pixel_offsets_x[plot_idxs]), np.arctan(scale * pixel_offsets_y[plot_idxs]),
                                  angles='xy', scale_units='xy', color=cmap(angles[plot_idxs]))
                        # scale = args.resolution/20000
                        # ax.quiver(o.fra, o.fdec, np.arctan(o.fra_shifts), np.arctan(o.fdec_shifts),
                        #           transform=ax.get_transform('world'),
                        #           color=cmap(angles),
                        #         #   width=0.002, headwidth=2, headlength=3, headaxislength=2, minlength=3
                        #           )

                    base, ext = os.path.splitext(args.tec)

                    plt.savefig(f"{base}-t{t:04d}{ext}", bbox_inches='tight', pad_inches=0, dpi=args.dpi)

if __name__ == '__main__':
    main()