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
import shlex


def get_parser():
    parser = argparse.ArgumentParser(description="Render polarimetric composite for fits images.")
    parser.add_argument('fits', nargs='+', help='fits files to render')
    
    plot_group = parser.add_argument_group('PLOTTING OPTIONS')
    plot_group.add_argument('--output_name',
        help='Name of output plot file', default=None
    )
    plot_group.add_argument('--plot_title', default=None,
        help="Optional title for the plot")
    plot_group.add_argument('--plot_dpi', default=100,
        help="dots per inch for the plot")
    plot_group.add_argument('--pol_order', nargs=3, default=["XX", "V", "YY"],
        help="Which pols are red, green and blue")
    plot_group.add_argument('--limits', nargs=3, default=[None, None, None],
        help="normalise each pol to a given brightness limit")
    plot_group.add_argument('--limit_percentile', default=99.5,
        help="percentile used to normalize pols if --limits is not specified")
    plot_group.add_argument('--red_blue_lock', default=True, 
        help="set red and blue max values to the max of the two")
    
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
        ] + \
            shlex.split("${args}") + \
            "${fits.join(' ')}".split(' '))
    
    header = None
    paths = {}
    for pol in args.pol_order:
        found = False
        for path in args.fits:
            if f"-{pol}-" in path:
                paths[pol] = path
                found = True
                break
        if not found:
            raise ValueError(f"No {pol} fits file found in paths: {', '.join(args.fits)}")
    data = {}
    for pol, path in paths.items():
        with fits.open(path) as hdus:
            if not header:
                header = hdus[0].header
            data[pol] = hdus[0].data[0,0,:,:]

    wcs = WCS(header)
    img_size = header['NAXIS1'], header['NAXIS2']
    # TODO: beam stuff
    # cell_size = header['CDELT1'], header['CDELT2']
    # beam_shape = header['BMAJ'], header['BMIN'], header['BPA']

    limits = {}
    for pol, norm in zip(args.pol_order, args.limits):
        if norm is None:
            limits[pol] = np.percentile(data[pol], args.limit_percentile)
        else:
            limits[pol] = float(norm)

    red_pol, green_pol, blue_pol = args.pol_order
    if args.red_blue_lock:
        limits[red_pol] = max(limits[red_pol], limits[blue_pol])
        limits[blue_pol] = limits[red_pol]

    subplot_kw = {"projection": wcs, "slices": ('x', 'y', 0, 0)}

    rgbMap = AsinhMapping(0., stretch=3., Q=5.)
    # normalize 0-norm to 0-1
    rgb = rgbMap.make_rgb_image(
        data[red_pol] / limits[red_pol], 
        data[green_pol] / limits[green_pol], 
        data[blue_pol] / limits[blue_pol]
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
            ["Comp", red_pol],
            ["Comp", blue_pol],
            ["Comp", green_pol],
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

    for pol, ny, neg_col, pos_col in [
        [red_pol, 2, 'cyan', 'red'],
        [blue_pol, 1, 'yellow', 'blue'],
        [green_pol, 0, 'magenta', 'green']
    ]:
        ax = axd[pol]
        ax.set_title(pol, y=0.95, fontdict={'verticalalignment':'top'})
        ax.set_axes_locator(divider.new_locator(nx=1, ny=ny))
        cmap = LinearSegmentedColormap.from_list(
            f'{neg_col}-bl-{pos_col}', [neg_col, 'black', pos_col], gamma=1
        )
        norm=ImageNormalize(data[pol], vmin=-limits[pol], vmax=limits[pol], clip=False)
        im = ax.imshow(data[pol], cmap=cmap, norm=norm)
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