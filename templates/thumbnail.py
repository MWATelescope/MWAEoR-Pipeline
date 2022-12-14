#!/usr/bin/env python

from astropy.io import fits
from astropy.visualization import astropy_mpl_style, SqrtStretch
from astropy.visualization import ImageNormalize, AsinhStretch, LinearStretch
from astropy.wcs import WCS
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
from argparse import ArgumentParser
import sys
import shlex

"""
example:

```bash
salloc --nodes=1 --mem=350G --time=8:00:00 --clusters=garrawarla --partition=workq --account=mwaeor --ntasks=20 --tmp=500G
module load singularity
mkdir /dev/shm/deleteme
cd /dev/shm/deleteme
# export obsid=1251909072
# export obsid=1090870648
export obsid=1090012424
# cp /astro/mwaeor/dev/nfdata/${obsid}/img-manualFlags/wsclean_hyp_${obsid}_30l_src4k_8s_80kHz-0000-V-uv-real.fits .
# cp /astro/mwaeor/dev/nfdata/${obsid}/img-4ch-768px/wsclean_hyp_${obsid}_ionosub_30l_src4k_8s_80kHz-0003-XX-uv-imag.fits img.fits
# cp /astro/mwaeor/dev/nfdata/${obsid}/img-4ch-768px/wsclean_hyp_${obsid}_ionosub_30l_src4k_8s_80kHz-0002-XX-image.fits img.fits
cp /astro/mwaeor/dev/nfdata/${obsid}/img-dev/wsclean-t0012-MFS-XX-dirty.fits img.fits
singularity exec --cleanenv --home /astro/mwaeor/dev/mplhome /pawsey/mwa/singularity/mwa_qa/mwa_qa_latest.sif python \
    /astro/mwaeor/dev/MWAEoR-Pipeline/templates/thumbnail.py \
    --fits=img.fits \
    --output_name="/astro/mwaeor/dev/MWAEoR-Pipeline/${obsid}_thumbtest.png" \
    --plot_title="${obsid}" \
    --transparent
```
"""


def make_fits_axis_array(header, axis):
    count = header[f"NAXIS{axis}"]
    crval = header[f"CRVAL{axis}"]
    cdelt = header[f"CDELT{axis}"]
    crpix = header[f"CRPIX{axis}"]
    return cdelt * (np.arange(count) - crpix) + crval


def get_parser():
    parser = ArgumentParser(
        description="Render thumbnail for fits image.")

    parser.add_argument('--fits', default=False, help='fits file to render')

    plot_group = parser.add_argument_group('PLOTTING OPTIONS')
    plot_group.add_argument('--output_name',
                            help='Name of output plot file', default=None
                            )
    plot_group.add_argument('--plot_title', default=None,
                            help="Optional title for the plot")
    plot_group.add_argument('--plot_dpi', default=100,
                            help="dots per inch for the plot")
    plot_group.add_argument('--limit', default=None,
                            help="normalise to a given brightness limit")
    plot_group.add_argument('--limit_percentile', default=99.5,
                            help="percentile used to normalize brightness if --limit is not specified")
    plot_group.add_argument('--transparent', default=False, action="store_true",
                            help="adds transparency to colormap")
    plot_group.add_argument('--symmetric', default=False, action="store_true",
                            help="make colormap symmetric")
    plot_group.add_argument('--cmap', default='plasma',
                            help="matplotlib colormap to use")

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
            "--plot_title=${title}"
        ] + shlex.split("${args}"))

    plt.style.use(astropy_mpl_style)

    with fits.open(args.fits) as hdus:
        header = hdus[0].header
        # print([*header.items()])
        data = hdus[0].data[0, 0, :, :]

    img_size = header['NAXIS1'], header['NAXIS2']
    # TODO: beam stuff
    # cell_size = header['CDELT1'], header['CDELT2']
    # beam_shape = header['BMAJ'], header['BMIN'], header['BPA']

    # normalize 0-percentile to 0-1
    if args.limit is None:
        limit = np.percentile(data, args.limit_percentile)
    else:
        limit = float(args.limit)

    with plt.style.context('dark_background'):
        img_fig = plt.figure(dpi=args.plot_dpi, figsize=(
            img_size[0]/args.plot_dpi, 4*img_size[1]/3/args.plot_dpi))

        subplot_kw = {}
        sky_coords = False
        if header['CTYPE1'] == 'RA---SIN':
            sky_coords = True
            subplot_kw["projection"] = WCS(header)
            subplot_kw["slices"] = ('x', 'y', 0, 0)
            stretch_a = 0.1
        else:
            stretch_a = 0.05
        ax = img_fig.add_subplot(1, 1, 1, **subplot_kw)

        nColors = 256
        if args.cmap:
            color_array = mpl.cm.get_cmap(args.cmap)(range(nColors))
        else:
            base_color_name = 'blue-bl-red' if args.symmetric else 'bl-red'
            base_colors = ['blue', 'black', 'red'] if args.symmetric \
                else ['black', 'red']
            color_array = LinearSegmentedColormap.from_list(
                base_color_name, base_colors)(range(nColors))
        if args.transparent:
            if args.symmetric:
                alphas = np.linspace(0, 1, nColors//2)
                alphas = np.stack((np.flip(alphas), alphas)).flatten()
            else:
                alphas = np.linspace(0, 1, nColors)
            color_array[:, -1] = alphas
        cmap = LinearSegmentedColormap.from_list(
            name='custom', colors=color_array)
        # Transforms normalized values [0,1] to [-1,1] before stretch and then back
        if args.symmetric:
            stretch = LinearStretch(
                slope=0.5, intercept=0.5) + AsinhStretch(stretch_a) + LinearStretch(slope=2, intercept=-1)
        else:
            stretch = SqrtStretch()

        vmin = -limit if args.symmetric else 0
        imshow_kw = {
            'cmap': cmap,
            'norm': ImageNormalize(data, vmin=vmin, vmax=limit, clip=False, stretch=stretch)
        }
        ax.set_label("")
        ax.set_title(args.plot_title, y=0.95, fontdict={
                     'verticalalignment': 'top'})
        if sky_coords:
            ax.set_ylabel("", visible=False)
            ax.set_xlabel("", visible=False)
            ax.tick_params(axis="y", direction="in", pad=-
                           50, horizontalalighment="left")
            ax.tick_params(axis="x", direction="in", pad=-
                           20, verticalalignment="bottom")
            for coord in ax.coords:
                coord.set_ticklabel_visible(True)
        else:
            ax.set_ylabel("V (wavenumbers)", visible=True)
            ax.set_xlabel("U (wavenumbers)", visible=True)
            ax.tick_params(axis="y", direction="in", pad=-30)
            ax.tick_params(axis="x", direction="in", pad=-20)
            ax.set_xticks([0])
            ax.set_yticks([0])
            # ax.axis('off')
            u = make_fits_axis_array(header, 1)
            v = make_fits_axis_array(header, 2)
            imshow_kw['extent'] = np.min(u), np.max(u), np.min(v), np.max(v)
        # ax.grid(ls='dotted')
        im = ax.imshow(data, **imshow_kw)
        # img_fig.colorbar(im)
        plt.savefig(args.output_name, bbox_inches='tight',
                    dpi=args.plot_dpi, transparent=args.transparent)


if __name__ == '__main__':
    main()
