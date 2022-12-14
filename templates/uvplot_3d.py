#!/usr/bin/env python
from astropy.io import fits
from astropy.visualization import astropy_mpl_style
from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt
from mwa_qa.read_uvfits import UVfits
import numpy as np
import sys
import argparse


def get_parser():
    parser = argparse.ArgumentParser(description="Plot Visibilities in 3D for each pol.")

    parser.add_argument('--uvfits', default=False, help='uvfits file to plot')

    plot_group = parser.add_argument_group('PLOTTING OPTIONS')
    plot_group.add_argument('--output_name', default=False,
        help='Name for output plot file, defaults to using the input uvfits name'
    )
    plot_group.add_argument('--plot_title', default=None,
        help="Optional title for the plot")
    plot_group.add_argument('--plot_dpi', default=100,
        help="dots per inch for the plot")
    plot_group.add_argument('--z_max', default=None,
        help="Optional maximum value for the Z axis")
    plot_group.add_argument('--dpi', default=100,
        help="Plot DPI")

    return parser

def uvplot_3d(args):
    uv = UVfits(args.uvfits)

    Npairs = len([(ap[0], ap[1]) for ap in uv.antpairs if ap[0] != ap[1]])
    blt_idxs = np.where(
        uv.ant_1_array - uv.ant_2_array != 0,
    )[0]

    pols = [ "XX", "YY", "XY", "YX" ]

    plt.style.use(astropy_mpl_style)

    fig = plt.figure(dpi=args.plot_dpi)

    with fits.open(uv.uvfits_path) as hdus:
        vis_hdu = hdus['PRIMARY']

        for pol_idx, pol in enumerate(pols):
            plt.clf()
            ax = plt.axes(projection='3d')
            reals = vis_hdu.data.data[blt_idxs, 0, 0, :, pol_idx, 0].reshape(
                (uv.Ntimes, Npairs, uv.Nfreqs))
            imags = vis_hdu.data.data[blt_idxs, 0, 0, :, pol_idx, 1].reshape(
                (uv.Ntimes, Npairs, uv.Nfreqs))
            wghts = vis_hdu.data.data[blt_idxs, 0, 0, :, pol_idx, 2].reshape(
                (uv.Ntimes, Npairs, uv.Nfreqs))
            flag_idxs = np.where(wghts < 0)[0]
            reals[flag_idxs] = np.nan
            imags[flag_idxs] = np.nan

            xs = vis_hdu.data['UU'][blt_idxs].reshape(uv.Ntimes, Npairs)
            xs = np.nanmean(xs, axis=(0,))
            ys = vis_hdu.data['VV'][blt_idxs].reshape(uv.Ntimes, Npairs)
            ys = np.nanmean(ys, axis=(0,))

            means = np.nanmean(reals + 1j * imags, axis=(0,2))
            print(f"{pol} shapes. xs={xs.shape} ys={ys.shape} reals={reals.shape} means={means.shape}")

            zs = np.abs(means)
            cs = np.angle(means)

            z_max = float(args.z_max) if args.z_max else np.quantile(zs, 0.90)

            ax.scatter3D(xs, ys, zs, c=cs, cmap='rainbow', s=1)
            ax.set_xlabel("u", visible=True)
            ax.set_ylabel("v", visible=True)
            ax.set_zlabel("Amplitude", visible=True)
            ax.set_zlim3d(0, z_max)
            ax.set_title(f"{args.plot_title} - {pol}")
            plt.savefig(f"{args.output_name}_{pol}.png", bbox_inches='tight', dpi=args.plot_dpi)

def main():
    """
    example:

    ```bash
    singularity exec /pawsey/mwa/singularity/mwa_qa/mwa_qa_latest.sif python \
        /astro/mwaeor/dev/MWAEoR-Pipeline/templates/uvplot_3d.py \
        --uvfits=/astro/mwaeor/dev/nfdata/1061315448/cal/hyp_1061315448_sub_30l_src4k_8s_80kHz.uvfits \
        --output_name=/astro/mwaeor/dev/nfdata/1061315448/vis_qa/uvplot_1061315448_sub_30l_src4k_8s_80kHz \
        --plot_title="1061315448_sub_30l_src4k_8s_80kHz"
    ```
    """

    parser = get_parser()

    if len(sys.argv) > 1:
        args = parser.parse_args()
    else:
        # is being called directly from nextflow
        args = parser.parse_args([
            "--uvfits=${uvfits}",
            "--output_name=${uvplot}",
            "--plot_title=${title}",
            "--z_max=${z_max}",
        ])

    uvplot_3d(args)

if __name__ == '__main__':
    main()