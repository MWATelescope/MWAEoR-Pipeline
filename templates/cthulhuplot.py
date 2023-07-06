#!/usr/bin/env python

import sys

import numpy as np
import matplotlib.pyplot as plt
from collections import OrderedDict

import json
import yaml
from yaml import CSafeLoader as SafeLoader
from argparse import ArgumentParser
import shlex
from math import pi
import os
import pandas as pd

from cthulhu.reconstruct import Obsid
from cthulhu.plot_tools import setup_subplot, plot_tec, plot_vector_arrows, generate_diagnostic_figure
from cthulhu.rts_log_tools import lm2radec


"""
example:

```bash
salloc --nodes=1 --mem=350G --time=8:00:00 --clusters=garrawarla --partition=workq --account=mwaeor --ntasks=20 --tmp=500G
module load singularity
mkdir /dev/shm/deleteme
cd /dev/shm/deleteme
export obsid=1061320200
cp /astro/mwaeor/dev/nfdata/${obsid}/cal/${obsid}_reduced_n8000.yaml srclist.yaml
cp /astro/mwaeor/dev/nfdata/${obsid}/cal/hyp_peel_${obsid}_ionosub_30l_src4k_8s_80kHz_uv.json offsets.json
singularity exec --bind \$PWD --cleanenv --home /astro/mwaeor/dev/mplhome /pawsey/mwa/singularity/mwa_qa/mwa_qa_latest.sif python \
    /pawsey/mwa/mwaeor/dev/MWAEoR-Pipeline/templates/cthulhuplot.py \
    --srclist=srclist.yaml \
    --offsets=offsets.json \
    --obsid=${obsid} \
    --output_name="${obsid}_chulhuplot.png"
cp ${obsid}_chulhuplot*.png /pawsey/mwa/mwaeor/dev/MWAEoR-Pipeline/test/
```
"""


def get_parser():
    parser = ArgumentParser(
        description="render TEC plot for offsets.")

    parser.add_argument('--srclist', help='source list file (yaml)')
    parser.add_argument('--offsets', help='offset file (json|yaml)')
    parser.add_argument('--obsid', default=0, type=int)
    parser.add_argument('--time_res', default=1, type=int,
                        help='time resolution (seconds)')
    parser.add_argument('--time_offset', default=0, type=int,
                        help='offset from obsid to first timestep (seconds)')
    plot_group = parser.add_argument_group('PLOTTING OPTIONS')
    plot_group.add_argument('--plot',
                            help='Name of output plot file', default=None)
    plot_group.add_argument('--json',
                            help='Name of output json file', default=None)
    plot_group.add_argument('--csv',
                            help='Name of output csv file', default=None)
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
            # "--srclist=${srclist}",
            "--offsets=${offsets}",
            "--obsid=${obsid}",
            "--plot=${plot}",
            "--csv=${csv}",
            "--json=${json}",
        ] + shlex.split(extra) if (extra := "${extra}") else [])

    # read source list
    # with open(args.srclist, "r") as h:
    #     srclist = yaml.load(h, Loader=SafeLoader)
    #     # srclist = json.load(h)

    # calculate offsets at a given freq
    freq = 200e6
    # freq = 150e6
    wavelength = 299792458 / freq
    wavelength_2 = wavelength**2

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

        l_shifts = np.array(alphas) * wavelength_2
        m_shifts = np.array(betas) * wavelength_2
        n_times = min(l_shifts.shape[1], m_shifts.shape[1])
        ras = np.array(ras)
        decs = np.array(decs)
        ras[np.where(ras > 180)] -= 360

        ra_shifts = np.zeros_like(l_shifts)
        dec_shifts = np.zeros_like(m_shifts)
        for i in range(n_times):
            ra_shifts_ts, dec_shifts_ts = lm2radec(ras, decs, np.rad2deg(l_shifts[:, i]), np.rad2deg(m_shifts[:, i]), {"primary_beam_pointing_centre": [0.0, -27.0]})
            ra_shifts[:, i] = np.array(ra_shifts_ts)
            dec_shifts[:, i] = np.array(dec_shifts_ts)
        ra_shifts = np.array(ra_shifts)
        dec_shifts = np.array(dec_shifts)
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
        ra_shifts, dec_shifts = -1 * ra_shifts, -1 * dec_shifts
        ras[np.where(ras > 180)] -= 360

    n_times = min(ra_shifts.shape[1], dec_shifts.shape[1])
    print(len(src_names), ras.shape, ra_shifts.shape)

    # print(l_shifts[0:5])
    av_ra_shifts = np.mean(ra_shifts, axis=1)
    av_dec_shifts = np.mean(dec_shifts, axis=1)
    av_total_shifts = np.sqrt(av_ra_shifts**2 + av_dec_shifts**2)
    max_total_shift_idx = np.argmax(av_total_shifts)
    print(av_ra_shifts[0:5])
    o = Obsid((ras, decs, av_ra_shifts, av_dec_shifts), args.obsid)
    o.pca()
    o.obsid_metric()

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

    if args.plot:
        print(f"saving plot to {args.plot}")
        if args.average:
            # plot the average of shifts over time
            o.reconstruct_tec()
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
                obsid = args.obsid + args.time_offset + (t * args.time_res)
                t_ra_shifts = ra_shifts[:, t]
                t_dec_shifts = dec_shifts[:, t]
                t_total_shifts = np.sqrt(t_ra_shifts**2 + t_dec_shifts**2)
                max_total_shift_idx = np.argmax(t_total_shifts)

                o = Obsid((ras, decs, t_ra_shifts, t_dec_shifts), obsid)
                o.pca()
                o.obsid_metric()
                o.reconstruct_tec()
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
                generate_diagnostic_figure(
                    o, filename=os.path.basename(path),
                    directory=os.path.dirname(path) or '.',
                    overwrite=True
                )

if __name__ == '__main__':
    main()