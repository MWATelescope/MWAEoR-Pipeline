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
    plot_group.add_argument('--flag_ants', default=None, nargs='+',
                            help="only plot tiles at these indices")

    return parser

def uvplot_2d(args):
    uv = UVfits(args.uvfits)

    Npairs = len([(ap[0], ap[1]) for ap in uv.antpairs if ap[0] != ap[1]])
    blt_idxs = np.where(
        uv.ant_1_array - uv.ant_2_array != 0,
    )[0]
    if args.flag_ants is not None:
        flag_ants = list(map(int, args.flag_ants))
        blt_idxs = blt_idxs[np.where(np.logical_not(np.logical_or(
            np.isin(uv.ant_1_array[blt_idxs], flag_ants),
            np.isin(uv.ant_2_array[blt_idxs], flag_ants),
        )))[0]]
        Npairs = len([(ap[0], ap[1]) for ap in uv.antpairs if ap[0] != ap[1] and ap[0] not in flag_ants and ap[1] not in flag_ants])

    pols = [ "XX", "YY", "XY", "YX" ]

    plt.style.use(astropy_mpl_style)

    fig = plt.figure(dpi=args.plot_dpi)

    print("opening hdus")

    with fits.open(uv.uvfits_path) as hdus:
        vis_hdu = hdus['PRIMARY']

        plt.clf()
        ax = plt.axes()

        xs = vis_hdu.data['UU'][blt_idxs].reshape(uv.Ntimes, Npairs)
        xs = np.nanmean(xs, axis=(0,))
        ys = vis_hdu.data['VV'][blt_idxs].reshape(uv.Ntimes, Npairs)
        ys = np.nanmean(ys, axis=(0,))

        print(f"shapes. xs={xs.shape} ys={ys.shape}")

        ax.scatter(xs, ys, cmap='rainbow', s=1)
        ax.set_xlabel("u", visible=True)
        ax.set_ylabel("v", visible=True)
        ax.set_title(f"{args.plot_title}")
        plt.savefig(f"{args.output_name}.png", bbox_inches='tight', dpi=args.plot_dpi)

        # for pol_idx, pol in enumerate(pols):
        #     plt.clf()
        #     # ax = plt.axes(projection='3d')
        #     # reals = vis_hdu.data.data[blt_idxs, 0, 0, :, pol_idx, 0].reshape(
        #     #     (uv.Ntimes, Npairs, uv.Nchan))
        #     # imags = vis_hdu.data.data[blt_idxs, 0, 0, :, pol_idx, 1].reshape(
        #     #     (uv.Ntimes, Npairs, uv.Nchan))
        #     # wghts = vis_hdu.data.data[blt_idxs, 0, 0, :, pol_idx, 2].reshape(
        #     #     (uv.Ntimes, Npairs, uv.Nchan))
        #     # flag_idxs = np.where(wghts < 0)[0]
        #     # reals[flag_idxs] = np.nan
        #     # imags[flag_idxs] = np.nan

        #     xs = vis_hdu.data['UU'][blt_idxs].reshape(uv.Ntimes, Npairs)
        #     xs = np.nanmean(xs, axis=(0,))
        #     ys = vis_hdu.data['VV'][blt_idxs].reshape(uv.Ntimes, Npairs)
        #     ys = np.nanmean(ys, axis=(0,))

        #     # means = np.nanmean(reals + 1j * imags, axis=(0,2))
        #     # print(f"{pol} shapes. xs={xs.shape} ys={ys.shape} reals={reals.shape} means={means.shape}")
        #     print(f"{pol} shapes. xs={xs.shape} ys={ys.shape}")

        #     # zs = np.abs(means)
        #     # cs = np.angle(means)

        #     # z_max = float(args.z_max) if args.z_max else np.quantile(zs, 0.90)

        #     ax.scatter(xs, ys, cmap='rainbow', s=1)
        #     ax.set_xlabel("u", visible=True)
        #     ax.set_ylabel("v", visible=True)
        #     # ax.set_zlabel("Amplitude", visible=True)
        #     # ax.set_zlim3d(0, z_max)
        #     ax.set_title(f"{args.plot_title} - {pol}")
        #     plt.savefig(f"{args.output_name}_{pol}.png", bbox_inches='tight', dpi=args.plot_dpi)

def main():
    """
    example:

    ```bash
    salloc --nodes=1 --mem=350G --time=8:00:00 --clusters=garrawarla --partition=workq --account=mwaeor --ntasks=20 --tmp=500G

    # !!! NVMETMP IS TOO SLOW !!!
    mkdir -p /dev/shm/deleteme
    cd /dev/shm/deleteme
    module load singularity
    export obsid=1366000096
    cp /astro/mwaeor/dev/nfdata/${obsid}/prep/birli_${obsid}*.uvfits .
    export uvfits=$(ls -1 birli_${obsid}*.uvfits | head -n 1)
    singularity exec --cleanenv --home /astro/mwaeor/dev/mplhome -B /nvmetmp /pawsey/mwa/singularity/mwa_qa/mwa_qa_latest.sif python \
        /pawsey/mwa/mwaeor/dev/MWAEoR-Pipeline/templates/uvplot_2d.py \
        --uvfits=${uvfits} \
        --output_name=uvplot \
        --plot_title="uvplot"

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
            # "--z_max=${z_max}",
        ])

    uvplot_2d(args)

if __name__ == '__main__':
    main()