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
from collections import OrderedDict

from cthulhu.reconstruct import Obsid
from cthulhu.plot_tools import setup_subplot, plot_tec, plot_vector_arrows, generate_diagnostic_figure
from cthulhu.rts_log_tools import lm2radec



def get_parser():
    parser = ArgumentParser(
        description="render TEC plot for offsets.")

    parser.add_argument('--srclist', help='source list file (yaml)')
    parser.add_argument('--offsets', help='offset file (json)')
    parser.add_argument('--obsid', default=0, type=int)
    parser.add_argument('--time_offset', default=0, type=int,
                        help='offset from obsid to first timestep (seconds)')
    plot_group = parser.add_argument_group('PLOTTING OPTIONS')
    plot_group.add_argument('--output_name',
                            help='Name of output plot file', default=None)

    return parser


def main():

    parser = get_parser()

    if len(sys.argv) > 1:
        args = parser.parse_args()
    else:
        print('is being called directly from nextflow')
        args = parser.parse_args([
            "--offsets=${offsets}",
            "--obsid=${obsid}",
            "--output_name=${ionoqa}",
        ])
        # "--srclist=${srclist}",

    print(args)

    # read iono constants
    with open(args.offsets, "r") as h:
        iono_consts = json.load(h)

    ras = []
    decs = []
    alphas = []
    betas = []
    src_names = []

    # calculate offsets at a given freq
    freq = 200e6
    # freq = 150e6
    wavelength = 299792458 / freq
    wavelength_2 = wavelength**2

    for src_name, consts in iono_consts.items():
        # src_alphas = consts["alphas"]
        # src_betas = consts["betas"]

        # # Get the positions of all components for this source from the srclist.
        # src = srclist[src_name]
        # comp_ras = [comp["ra"] for comp in src]
        # comp_decs = [comp["dec"] for comp in src]

        # # Get the average position of the components. This should be weighted by
        # # brightness; for now, make all components equal weight.
        # ra = np.mean(comp_ras)
        # dec = np.mean(comp_decs)

        pos = consts["weighted_catalogue_pos_j2000"]
        ras.append(pos["ra"])
        decs.append(pos["dec"])
        alphas.append(consts["alphas"])
        betas.append(consts["betas"])
        src_names.append(src_name)

    ras = np.array(ras)
    decs = np.array(decs)
    l_shifts = np.array(alphas) * wavelength_2
    m_shifts = np.array(betas) * wavelength_2
    print(ras.shape, l_shifts.shape)

    # Fix the branch cut.
    # ras = np.rad2deg(np.unwrap(np.deg2rad(ras), discont=pi))
    ras[np.where(ras > 180)] -= 360

    n_times = min(l_shifts.shape[1], m_shifts.shape[1])

    ra_shifts = np.zeros_like(l_shifts)
    dec_shifts = np.zeros_like(m_shifts)
    for i in range(n_times):
        ra_shifts_ts, dec_shifts_ts = lm2radec(ras, decs, np.rad2deg(l_shifts[:, i]), np.rad2deg(m_shifts[:, i]), {"primary_beam_pointing_centre": [0.0, -27.0]})
        ra_shifts[:, i] = np.array(ra_shifts_ts)
        dec_shifts[:, i] = np.array(dec_shifts_ts)
    # ra_shifts, dec_shifts = lm2radec2(ras, decs, np.rad2deg(l_shifts), np.rad2deg(m_shifts))
    # ra_shifts = l_shifts
    # dec_shifts = m_shifts

    ra_shifts = np.array(ra_shifts)
    dec_shifts = np.array(dec_shifts)
    print(ra_shifts[0:5])
    print(l_shifts[0:5])

    # ra_shifts -= np.mean(ra_shifts)
    # dec_shifts -= np.mean(dec_shifts)

    data = OrderedDict()

    av_ra_shifts = np.mean(ra_shifts, axis=1)
    av_dec_shifts = np.mean(dec_shifts, axis=1)
    av_total_shifts = np.sqrt(av_ra_shifts**2 + av_dec_shifts**2)
    max_total_shift_idx = np.argmax(av_total_shifts)
    o = Obsid((ras, decs, av_ra_shifts, av_dec_shifts), args.obsid)
    o.pca()
    o.obsid_metric()

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

    print('data:')
    for key, value in data.items():
        print(f'{key} => {value} ({type(value)})')
    print(f"saving to {args.output_name}")
    with open(args.output_name, 'w') as f:
        json.dump(data, f)


if __name__ == '__main__':
    main()