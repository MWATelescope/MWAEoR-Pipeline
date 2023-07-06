#!/usr/bin/env python

from astropy.io import fits
from astropy.visualization import astropy_mpl_style, SqrtStretch
from astropy.visualization import ImageNormalize, AsinhStretch, LinearStretch
from astropy.visualization.mpl_normalize import simple_norm
from astropy.wcs import WCS
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
from argparse import ArgumentParser
import sys
import shlex
import os
import json
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

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
# cp /astro/mwaeor/dev/nfdata/\${obsid}/img-manualFlags/wsclean_hyp_\${obsid}_30l_src4k_8s_80kHz-0000-V-uv-real.fits .
# cp /astro/mwaeor/dev/nfdata/\${obsid}/img-4ch-768px/wsclean_hyp_\${obsid}_ionosub_30l_src4k_8s_80kHz-0003-XX-uv-imag.fits img.fits
# cp /astro/mwaeor/dev/nfdata/\${obsid}/img-4ch-768px/wsclean_hyp_\${obsid}_ionosub_30l_src4k_8s_80kHz-0002-XX-image.fits img.fits
cp /astro/mwaeor/dev/nfdata/\${obsid}/img-dev/wsclean-t0012-MFS-XX-dirty.fits img.fits
singularity exec --cleanenv -B \$PWD --home /astro/mwaeor/dev/mplhome /pawsey/mwa/singularity/mwa_qa/mwa_qa_latest.sif python \
    /pawsey/mwa/mwaeor/dev/MWAEoR-Pipeline/templates/thumbnail.py \
    --fits=img.fits \
    --thumb="/astro/mwaeor/dev/MWAEoR-Pipeline/\${obsid}_thumbtest.thumb" \
    --title="\${obsid}" \
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
    plot_group.add_argument('--thumb',
                            help='Name of output plot file', default=None)
    plot_group.add_argument('--title', default=None,
                            help="Optional title for the plot")
    plot_group.add_argument('--subtitle', default=None,
                            help="Optional subtitle for the plot")
    plot_group.add_argument('--dpi', default=100, type=int,
                            help="dots per inch for the plot")
    plot_group.add_argument('--vmin', default=None, type=float,
                            help="quantile used to normalize brightness if --limit is not specified")
    plot_group.add_argument('--vmax', default=None, type=float,
                            help="quantile used to normalize brightness if --limit is not specified")
    plot_group.add_argument('--vmin_quantile', default=0.01, type=float,
                            help="quantile used to normalize brightness if --vmin is not specified")
    plot_group.add_argument('--vmax_quantile', default=0.99, type=float,
                            help="quantile used to normalize brightness if --vmax is not specified")
    plot_group.add_argument('--transparent', default=False, action="store_true",
                            help="adds transparency to colormap")
    plot_group.add_argument('--symmetric', default=False, action="store_true",
                            help="make colormap symmetric")
    plot_group.add_argument('--cmap', default='plasma',
                            help="matplotlib colormap to use")
    plot_group.add_argument('--norm_args', default=None, help="json parameters for astropy.visualization.mpl_normalize.simple_norm")
    # plot_group.add_argument('--scale', default=1, type=float,
    #                         help="scale pixels by factor")

    return parser


def main():

    parser = get_parser()

    if len(sys.argv) > 1:
        args = parser.parse_args()
    else:
        # is being called directly from nextflow with args ${args}
        args = parser.parse_args(shlex.split('${argstr}'))

    print(f"{args=}")

    plt.style.use(astropy_mpl_style)

    subplot_kw = {}
    with fits.open(args.fits) as hdus:
        header = hdus[0].header
        # print([*header.items()])
        shape = hdus[0].data.shape
        print(shape)
        if len(shape) == 4:
            data = hdus[0].data[0, 0, :, :]
            slices = ('x', 'y', 0, 0)
        elif len(shape) == 2:
            data = hdus[0].data[:, :]
            slices = ('x', 'y')

        # header['CDELT1'] /= args.scale
        # header['CDELT2'] /= args.scale

    img_size = header['NAXIS1'], header['NAXIS2']
    # TODO: beam stuff
    # cell_size = header['CDELT1'], header['CDELT2']
    # beam_shape = header['BMAJ'], header['BMIN'], header['BPA']

    vmin = args.vmin
    if not vmin:
        vmin = np.quantile(data, args.vmin_quantile)

    vmax = args.vmax
    if not vmax:
        vmax = np.quantile(data, args.vmax_quantile)

    if args.symmetric:
        vmax = max(abs(vmin), abs(vmax))
        vmin = -vmax

    with plt.style.context('dark_background'):
        img_fig = plt.figure(dpi=args.dpi, figsize=(
            img_size[0]/args.dpi, 4*img_size[1]/3/args.dpi))

        imtype = None
        if header.get('CTYPE1') == 'RA---SIN':
            imtype = 'sky'
            subplot_kw["projection"] = WCS(header)
            subplot_kw["slices"] = slices
            stretch_a = 0.5
        elif header.get('CTYPE1') == 'U---WAV':
            imtype = 'uv'
            stretch_a = 0.05
        else:
            stretch_a = 0.1
        print("subplot_kw", subplot_kw)
        ax = img_fig.add_subplot(1, 1, 1, **subplot_kw)

        nColors = 256
        if not args.symmetric and args.cmap:
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

        # stretch = {
        #     'sqrt': SqrtStretch(),
        #     'linear': LinearStretch(),
        #     'asinh': AsinhStretch(stretch_a),
        # }[args.stretch]
        # Transforms normalized values [0,1] to [-1,1] before stretch and then back
        # if args.symmetric:
        #     stretch = LinearStretch(
        #         slope=0.5, intercept=0.5) + stretch + LinearStretch(slope=2, intercept=-1)

        norm_args = {}
        if args.norm_args:
            print(f"{args.norm_args!r}")
            norm_args = json.loads(args.norm_args)

        imshow_kw = {
            'cmap': cmap,
            # 'norm': ImageNormalize(data, vmin=vmin, vmax=vmax, clip=False, stretch=stretch)
            'norm': simple_norm(data, min_cut=vmin, max_cut=vmax, **norm_args)
        }

        ax.set_label("")
        title = args.title or header.get('OBJECT')
        subtitle = args.subtitle or header.get('DATE-OBS')
        title = chr(10).join(filter(None,[title, subtitle]))
        ax.set_title(title, y=0.95, fontdict={
                    'verticalalignment': 'top'})

        if imtype == 'sky':
            ax.set_ylabel("", visible=False)
            ax.set_xlabel("", visible=False)
            ax.tick_params(axis="y", direction="in", pad=-
                           50, horizontalalighment="left", verticalalignment="bottom")
            ax.tick_params(axis="x", direction="in", pad=-
                           20, verticalalignment="bottom")
            for coord in ax.coords:
                coord.set_ticklabel_visible(True)
        elif imtype == 'uv':
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
        else:
            ax.set_ylabel("Y (pixels)", visible=True)
            ax.set_xlabel("X (pixels)", visible=True)
            ax.set_xticks([])
            ax.set_yticks([])
        im = ax.imshow(data, **imshow_kw)
        ax.grid()
        cbaxes = inset_axes(ax, width="30%", height="3%", loc='upper right')
        img_fig.colorbar(im, cax=cbaxes, orientation="horizontal")
        fits_base, _ = os.path.splitext(os.path.basename(args.fits))
        thumb = args.thumb or f"{fits_base}.thumb"
        plt.savefig(thumb, bbox_inches='tight',
                    dpi=args.dpi, transparent=args.transparent, figsize=(10, 10))


if __name__ == '__main__':
    main()
