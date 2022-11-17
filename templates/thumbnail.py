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
    parser = argparse.ArgumentParser(description="Render thumbnail for fits image.")

    parser.add_argument('--fits', default=False, help='fits file to render')

    plot_group = parser.add_argument_group('PLOTTING OPTIONS')
    plot_group.add_argument('--output_name',
        help='Name of output plot file', default=None
    )
    plot_group.add_argument('--plot_title', default=None,
        help="Optional title for the plot")
    plot_group.add_argument('--plot_dpi', default=100)

    return parser

def main():

    parser = get_parser()

    if len(sys.argv) > 1:
        args = parser.parse_args()
    else:
        # is being called directly from nextflow
        args = parser.parse_args([
            "--fits=${img}",
            "--output_name=${thumb}",
            "--plot_title=${obsid} - ${name} - ${pol}",
        ])

    plt.style.use(astropy_mpl_style)

    with fits.open(args.fits) as hdus:
        header = hdus[0].header
        data = hdus[0].data[0,0,:,:]

    wcs = WCS(header)
    img_size = header['NAXIS1'], header['NAXIS2']
    # TODO: beam stuff
    # cell_size = header['CDELT1'], header['CDELT2']
    # beam_shape = header['BMAJ'], header['BMIN'], header['BPA']

    # normalize 0-percentile to 0-1
    percentile = np.percentile(data, 99.5)

    with plt.style.context('dark_background'):
        img_fig = plt.figure(dpi=args.plot_dpi, figsize=(img_size[0]/args.plot_dpi, 4*img_size[1]/3/args.plot_dpi))

        subplot_kw = {"projection": wcs, "slices": ('x', 'y', 0, 0)}
        ax = img_fig.add_subplot(1, 1, 1, **subplot_kw)
        ax.set_ylabel("", visible=False)
        ax.set_xlabel("", visible=False)
        ax.set_label("")
        ax.tick_params(axis="y", direction="in", pad=-30, horizontalalighment="left")
        ax.tick_params(axis="x", direction="in", pad=-20, verticalalignment="bottom")

        for coord in ax.coords:
            coord.set_ticklabel_visible(True)
        ax.set_title(args.plot_title, y=0.95, fontdict={'verticalalignment':'top'})
        cmap = LinearSegmentedColormap.from_list(f'blue-bl-red', ['blue', 'black', 'red'], gamma=1)
        norm=ImageNormalize(data, vmin=-percentile, vmax=percentile, clip=False)
        ax.imshow(data, cmap=cmap, norm=norm)
        plt.savefig(args.output_name, bbox_inches='tight', dpi=args.plot_dpi)

if __name__ == '__main__':
    main()