#!/usr/bin/env python

from astropy.io import fits
from astropy.visualization import astropy_mpl_style, make_lupton_rgb, simple_norm, SqrtStretch, imshow_norm
from astropy.visualization import ImageNormalize, PowerDistStretch, SinhStretch, LogStretch, AsinhStretch
from astropy.visualization.lupton_rgb import AsinhMapping, LinearMapping
from astropy.visualization.wcsaxes import WCSAxesSubplot, WCSAxes
from astropy.wcs import WCS
from matplotlib.patches import Rectangle, Circle
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.axes_divider import AxesDivider, HBoxDivider, make_axes_locatable, Size
from mpl_toolkits.axes_grid1.axes_rgb import make_rgb_axes
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import numpy as np
from os.path import basename, exists
from glob import glob
import argparse
import sys


def get_parser():
    parser = argparse.ArgumentParser(description="Render polarimetric composite for fits images.")
    parser.add_argument('fits', nargs='+', help='fits files to render')
    
    plot_group = parser.add_argument_group('PLOTTING OPTIONS')
    plot_group.add_argument('--output_name',
        help='Name of output plot file', default=None
    )
    plot_group.add_argument('--plot_title', default=None,
        help="Optional title for the plot")
    plot_group.add_argument('--plot_dpi', default=100)
    plot_group.add_argument('--norm_percentile', default=99.5)
    plot_group.add_argument('--xy_percentile', default=None)
    plot_group.add_argument('--v_percentile', default=None)
    
    return parser

def main():
    
    parser = get_parser()
    
    if len(sys.argv) > 1:
        args = parser.parse_args()
    else:
        # is being called directly from nextflow
        args = parser.parse_args([
            "--output_name=${polcomp}",
            "--plot_title=${title}",
            "--xy_percentile=${xy_percentile}",
            "--v_percentile=${v_percentile}",
        ] + "${fits.join(' ')}".split(' '))
    
    header = None
    paths = {}
    for pol in ["XX", "YY", "V"]:
        for path in args.fits:
            if f"-{pol}-" in path:
                paths[pol] = path
                break
        if not pol in paths:
            raise ValueError(f"No {pol} fits file found")
    data = {}
    for pol, path in paths.items():
        with fits.open(path) as hdus:
            if not header:
                header = hdus[0].header
            data[pol] = hdus[0].data[0,0,:,:]

    if not data:
        raise Exception("No data found")

    wcs = WCS(header)
    img_size = header['NAXIS1'], header['NAXIS2']
    # TODO: beam stuff
    # cell_size = header['CDELT1'], header['CDELT2']
    # beam_shape = header['BMAJ'], header['BMIN'], header['BPA']

    # normalize 0-percentile to 0-1, use same percentile for XX and YY
    
    xy_percentile = float(args.xy_percentile) if args.xy_percentile else \
        np.percentile(np.stack((data["XX"], data["YY"])), args.norm_percentile)
    
    v_percentile = float(args.v_percentile) if args.v_percentile else \
        np.percentile(data["V"], args.norm_percentile)

    subplot_kw = {"projection": wcs, "slices": ('x', 'y', 0, 0)}

    if any((k not in data) for k in ["XX", "YY", "V"]):
        print("could not make polcomp, XX, YY, or V missing")
        exit(0)

    rgbMap = AsinhMapping(0., stretch=3., Q=5.)
    rgb = rgbMap.make_rgb_image(
        data["XX"] / xy_percentile, 
        data["V"] / v_percentile, 
        data["YY"] / xy_percentile
    )

    # sources = np.array([
    #     [0, -27], # eor0
    #     [6.4549166666666675, -26.04], # PKS0023_026
    #     # [45, -26.04]
    # ])

    plt.style.use([astropy_mpl_style, 'dark_background'])
    
    img_fig = plt.figure(
        dpi=args.plot_dpi, 
        figsize=(img_size[0]/args.plot_dpi, 4*img_size[1]/3/args.plot_dpi)
    )
    # axis setup

    axd = img_fig.subplot_mosaic(
        [
            ["Comp", "XX"],
            ["Comp", "YY"],
            ["Comp", "V"],
        ],
        subplot_kw=subplot_kw,
        gridspec_kw={ "width_ratios": [3, 1], },
    )

    for ax in axd.values():
        ax.set_ylabel("", visible=False)
        ax.set_xlabel("", visible=False)
        ax.set_label("")

    divider = AxesDivider(axd['Comp'])
    locator = divider.new_locator(nx=0, ny=0)
    axd['Comp'].set_axes_locator(locator)

    sub_xsize = (1/3) * Size.AxesX(axd['Comp'])
    cb_xsize = (1/10) * sub_xsize
    ysize = (1./3) * Size.AxesY(axd['Comp'])

    divider.set_horizontal([Size.AxesX(axd['Comp']), sub_xsize, cb_xsize])
    divider.set_vertical([ysize, ysize, ysize])

    axd['Comp'].set_axes_locator(divider.new_locator(0, 0, ny1=-1))
    axd["Comp"].tick_params(axis="y", direction="in", pad=-30, horizontalalighment="left")
    axd["Comp"].tick_params(axis="x", direction="in", pad=-20, verticalalignment="bottom")
    axd["Comp"].imshow(rgb)
    axd["Comp"].set_title(args.plot_title, y=0.95, fontdict={'verticalalignment':'top'})

    for key, ny, neg_col, pos_col, vmax in [
        ['XX', 2, 'cyan', 'red', xy_percentile],
        ['YY', 1, 'yellow', 'blue', xy_percentile],
        ['V', 0, 'magenta', 'green', v_percentile]
    ]:
        ax = axd[key]
        ax.set_title(key, y=0.95, fontdict={'verticalalignment':'top'})
        ax.set_axes_locator(divider.new_locator(nx=1, ny=ny))
        cmap = LinearSegmentedColormap.from_list(f'{neg_col}-bl-{pos_col}', [neg_col, 'black', pos_col], gamma=1)
        norm=ImageNormalize(data[key], vmin=-vmax, vmax=vmax, clip=False)
        im = ax.imshow(data[key], cmap=cmap, norm=norm)
        for coord in ax.coords:
            coord.set_ticklabel_visible(False)
        cax = inset_axes(
            ax,
            width="5%",
            height="100%",
            loc='lower left',
            bbox_to_anchor=(1.0, 0., 1, 1),
            bbox_transform=ax.transAxes,
            borderpad=0,
            axes_kwargs={ "axisbelow": False}
        )
        img_fig.colorbar(im, cax=cax)

    # for ax in axd.values():
    #     ax.scatter(sources[:,0], sources[:,1], s=100, edgecolor='white', facecolor='none', transform=ax.get_transform('world'))
    #     ax.add_patch(Rectangle((0, 0), 100, 100, edgecolor='white', fill=False, zorder=999))

    #  hist_ax = inset_axes(
    #      axd["Comp"],
    #      width="100%",
    #      height="20%",
    #      loc='lower left',
    #      bbox_to_anchor=(0., -0.25, 1, 1),
    #      bbox_transform=axd['Comp'].transAxes,
    #      borderpad=0,
    #  )

    #  hist_ax.hist(data["XX"].flatten(), bins=1000, color='red', alpha=0.5, histtype='step')
    #  hist_ax.hist(data["YY"].flatten(), bins=1000, color='blue', alpha=0.5, histtype='step')
    #  hist_ax.hist(data["V"].flatten(), bins=1000, color='green', alpha=0.5, histtype='step')

    plt.savefig(args.output_name, bbox_inches='tight', dpi=args.plot_dpi)
    
if __name__ == '__main__':
    main()