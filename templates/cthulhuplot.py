#!/usr/bin/env python

import sys

import numpy as np
import matplotlib.pyplot as plt

import json
import yaml
from yaml import CSafeLoader as SafeLoader
from argparse import ArgumentParser
import shlex
from math import pi
import os

from cthulhu.reconstruct import Obsid
from cthulhu.plot_tools import setup_subplot, plot_tec, plot_vector_arrows, generate_diagnostic_figure


"""
example:

```bash
salloc --nodes=1 --mem=350G --time=8:00:00 --clusters=garrawarla --partition=workq --account=mwaeor --ntasks=20 --tmp=500G
module load singularity
mkdir /dev/shm/deleteme
cd /dev/shm/deleteme
export obsid=1322482208
export obsid=1094230128
export obsid=1094230616
export obsid=1090012424
cp /astro/mwaeor/dev/nfdata/${obsid}/cal/${obsid}_reduced_n4000.yaml srclist.yaml
cp /astro/mwaeor/dev/nfdata/${obsid}/cal/hyp_peel_${obsid}_ionosub_30l_src4k_8s_80kHz_uv.json offsets.json
singularity exec --cleanenv --home /astro/mwaeor/dev/mplhome /pawsey/mwa/singularity/mwa_qa/mwa_qa_latest.sif python \
    /astro/mwaeor/dev/MWAEoR-Pipeline/templates/cthulhuplot.py \
    --srclist=srclist.yaml \
    --offsets=offsets.json \
    --obsid=${obsid} \
    --output_name="${obsid}_chulhuplot.png"
cp ${obsid}_chulhuplot*.png /astro/mwaeor/dev/MWAEoR-Pipeline/
```
"""


def get_parser():
    parser = ArgumentParser(
        description="render TEC plot for offsets.")

    parser.add_argument('--srclist', help='source list file (yaml)')
    parser.add_argument('--offsets', help='offset file (json)')
    parser.add_argument('--obsid', default=0, type=int)
    parser.add_argument('--time_res', default=1, type=int,
                        help='time resolution in seconds')
    plot_group = parser.add_argument_group('PLOTTING OPTIONS')
    plot_group.add_argument('--output_name',
                            help='Name of output plot file', default=None)
    plot_group.add_argument('--average', action='store_true', default=False,
                            help="Average offsets together")

    return parser


def main():

    parser = get_parser()

    if len(sys.argv) > 1:
        args = parser.parse_args()
    else:
        # is being called directly from nextflow
        args = parser.parse_args([
            "--srclist=${srclist}",
            "--offsets=${offsets}",
            "--obsid=${obsid}",
            "--output_name=${cthulhuplot}",
        ]  # + shlex.split("${args}")
        )

    # read source list
    with open(args.srclist, "r") as h:
        srclist = yaml.load(h, Loader=SafeLoader)

    # read iono constants
    with open(args.offsets, "r") as h:
        iono_consts = json.load(h)

    # print(srclist.keys())
    # print(iono_consts.keys())

    ras = []
    decs = []
    ra_shifts = []
    dec_shifts = []

    # calculate offsets at a given freq
    # freq = 200e6
    freq = 150e6
    wavelength = 299792458 / freq
    wavelength_2 = wavelength**2

    for src_name, consts in iono_consts.items():
        alphas = consts[0]
        betas = consts[1]

        # Get the position from the srclist
        src = srclist[src_name][0]
        ra = src["ra"]
        dec = src["dec"]

        ras.append(ra)
        decs.append(dec)
        ra_shifts.append(alphas)
        dec_shifts.append(betas)

    ras = np.array(ras)
    decs = np.array(decs)
    ra_shifts = np.array(ra_shifts)
    dec_shifts = np.array(dec_shifts)
    n_times = min(ra_shifts.shape[1], dec_shifts.shape[1])

    # Fix the branch cut.
    ras = np.rad2deg(np.unwrap(np.deg2rad(ras), discont=pi))

    if args.average:
        # plot the average of shifts over time
        t_ra_shifts = np.mean(ra_shifts, axis=1) * wavelength_2
        t_dec_shifts = np.mean(dec_shifts, axis=1) * wavelength_2
        o = Obsid((ras, decs, t_ra_shifts, t_dec_shifts), args.obsid)
        o.pca()
        o.obsid_metric()
        o.reconstruct_tec()
        o.tec_power_spectrum()
        generate_diagnostic_figure(
            o, filename=args.output_name, directory='.', overwrite=True)
    else:
        # plot each time separately
        for t in range(n_times):
            obsid = args.obsid + (t * args.time_res)
            t_ra_shifts = ra_shifts[:, t] * wavelength_2
            t_dec_shifts = dec_shifts[:, t] * wavelength_2
            o = Obsid((ras, decs, t_ra_shifts, t_dec_shifts), obsid)
            o.pca()
            o.obsid_metric()
            o.reconstruct_tec()
            o.tec_power_spectrum()
            output_name, ext = os.path.splitext(args.output_name)
            filename = f"{output_name}-t{t:04d}{ext}"
            generate_diagnostic_figure(
                o, filename=filename, directory='.', overwrite=True)

    # print(ras, decs)
    # print(args.obsid)
    # o = Obsid((ras, decs, ra_shifts, dec_shifts), obsid=args.obsid)
    # o.pca()
    # o.obsid_metric()
    # o.reconstruct_tec()
    # o.tec_power_spectrum()
    # generate_diagnostic_figure(
    #     o, filename=args.output_name, directory='.', overwrite=True)

    # fig, ax = plt.subplots()
    # setup_subplot(axis=ax)
    # plot_vector_arrows(axis=ax, obsid=o, scale=1e4)
    # ax.imshow(o.tec, extent=o.tec_extent, origin="lower") #, vmin=0, vmax=1)
    # ax.set_xlim(o.tec_extent[0], o.tec_extent[1])
    # ax.set_ylim(o.tec_extent[2], o.tec_extent[3])
    # plt.savefig(args.output_name, dpi=200)
    # # plt.show()
    print(f"saved to {args.output_name}")


if __name__ == '__main__':
    main()
