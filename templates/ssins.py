#!/usr/bin/env python

import os
import json

from SSINS import Catalog_Plot as cp, util as ssins_util
from SSINS import SS, INS, MF, plot_lib
from matplotlib import pyplot as plt, cm
import numpy as np
import pylab
from pyuvdata import UVData, UVFlag, utils as uvutils
import sys


def get_parser():
    import argparse

    parser = argparse.ArgumentParser(
        description="Run SSINS on uvfits visibilities.")

    parser.add_argument('--uvfits', default=False,
                        help='uvfits file to flag')

    plot_group = parser.add_argument_group('PLOTTING OPTIONS')
    plot_group.add_argument('--output_prefix',
                            help='Prefix of output plot files', default=''
                            )
    plot_group.add_argument('--plot_title', default=None,
                            help="Optional title for the plot")
    plot_group.add_argument('--guard_width', default=0, type=int,
                            help="Guard width of RFI bands in Hz. Half a fine channel width is recommended.")

    return parser


def main():
    """
    example:

    ```bash
    singularity exec /pawsey/mwa/singularity/ssins/ssins_latest.sif python \
        /astro/mwaeor/dev/MWAEoR-Pipeline/templates/ssins.py \
        --uvfits=/astro/mwaeor/dev/nfdata/1061315448/prep/birli_1061315448_2s_40kHz.uvfits \
        --output_prefix=/astro/mwaeor/dev/nfdata/1061315448/prep/test \
        --plot_title="1061315448" \
        --guard_width=20000
    """

    parser = get_parser()

    if len(sys.argv) > 1:
        args = parser.parse_args()
    else:
        # is being called directly from nextflow
        args = parser.parse_args([
            "--uvfits=${uvfits}",
            "--plot_title=${base}",
            "--guard_width=${guard_width}",
            "--output_prefix=${base}_",
        ])

    ins_plot_args = {
        "file_ext": "png",
        "title": args.plot_title,
        "extent_time_format": "lst"
    }
    shape_dict = {
        "DAB-5A":   [1.74160e8 - args.guard_width, 1.75696e8 + args.guard_width],
        "DAB-5B":   [1.75872e8 - args.guard_width, 1.77408e8 + args.guard_width],
        "DAB-5C":   [1.77584e8 - args.guard_width, 1.79120e8 + args.guard_width],
        "DAB-5D":   [1.79296e8 - args.guard_width, 1.80832e8 + args.guard_width],
        "DAB-6A":   [1.81168e8 - args.guard_width, 1.82704e8 + args.guard_width],
        "DAB-6B":   [1.82880e8 - args.guard_width, 1.84416e8 + args.guard_width],
        "DAB-6C":   [1.84592e8 - args.guard_width, 1.86128e8 + args.guard_width],
        "DAB-6D":   [1.86304e8 - args.guard_width, 1.87840e8 + args.guard_width],
        "DAB-7A":   [1.88160e8 - args.guard_width, 1.89696e8 + args.guard_width],
        "DAB-7B":   [1.89872e8 - args.guard_width, 1.91408e8 + args.guard_width],
        "DAB-7C":   [1.91584e8 - args.guard_width, 1.93120e8 + args.guard_width],
        "DAB-7D":   [1.93296e8 - args.guard_width, 1.94832e8 + args.guard_width],
        "DAB-8A":   [1.95168e8 - args.guard_width, 1.96704e8 + args.guard_width],
        "DAB-8B":   [1.96880e8 - args.guard_width, 1.98416e8 + args.guard_width],
        "DAB-8C":   [1.98592e8 - args.guard_width, 2.00128e8 + args.guard_width],
        "DAB-8D":   [2.00304e8 - args.guard_width, 2.01840e8 + args.guard_width],
        "DAB-9A":   [2.02160e8 - args.guard_width, 2.03696e8 + args.guard_width],
        "DAB-9B":   [2.03872e8 - args.guard_width, 2.05408e8 + args.guard_width],
        "DAB-9C":   [2.05584e8 - args.guard_width, 2.07120e8 + args.guard_width],
        "DAB-9D":   [2.07296e8 - args.guard_width, 2.08832e8 + args.guard_width],
        "DAB-10A":  [2.09168e8 - args.guard_width, 2.10704e8 + args.guard_width],
        "DAB-10B":  [2.10880e8 - args.guard_width, 2.12416e8 + args.guard_width],
        "DAB-10C":  [2.12592e8 - args.guard_width, 2.14128e8 + args.guard_width],
        "DAB-10D":  [2.14304e8 - args.guard_width, 2.15840e8 + args.guard_width],
        "DAB-11A":  [2.16160e8 - args.guard_width, 2.17696e8 + args.guard_width],
        "DAB-11B":  [2.17872e8 - args.guard_width, 2.19408e8 + args.guard_width],
        "DAB-11C":  [2.19584e8 - args.guard_width, 2.21120e8 + args.guard_width],
        "DAB-11D":  [2.21296e8 - args.guard_width, 2.22832e8 + args.guard_width],
        "DAB-12A":  [2.23168e8 - args.guard_width, 2.24704e8 + args.guard_width],
        "DAB-12B":  [2.24880e8 - args.guard_width, 2.26416e8 + args.guard_width],
        "DAB-12C":  [2.26592e8 - args.guard_width, 2.28128e8 + args.guard_width],
        "DAB-12D":  [2.28304e8 - args.guard_width, 2.29840e8 + args.guard_width],
    }
    sig_thresh = {shape: 5 for shape in shape_dict}
    sig_thresh['narrow'] = 7
    sig_thresh['streak'] = 8

    ss = SS()
    ss.read(args.uvfits, read_data=False)

    # discaring the first and last tiemstamp
    times = np.unique(ss.time_array)[1:-1]
    ss.read(args.uvfits, read_data=True, times=times, diff=True)
    ss.apply_flags(flag_choice='original')

    # cp.VDH_plot(ss, 'ssins_vdh', file_ext='png',
    #     pre_flag=True, post_flag=True, pre_model=False, post_model=False,
    #     post_label='Post-Flag Data', pre_label='Pre-Flag Data',
    #     legend=True)

    # RFI in autos
    ins_autos = INS(ss, spectrum_type="auto")
    cp.INS_plot(ins_autos, f'{args.output_prefix}autos', **ins_plot_args)

    # RFI in crosses
    ins_cross = INS(ss, spectrum_type='cross')
    cp.INS_plot(ins_cross, f'{args.output_prefix}cross', **ins_plot_args)

    mf = MF(ins_cross.freq_array, sig_thresh,
            shape_dict=shape_dict, streak=True, narrow=True)
    ins_cross.metric_array[ins_cross.metric_array == 0] = np.ma.masked
    ins_cross.metric_ms = ins_cross.mean_subtract()
    ins_cross.sig_array = np.ma.copy(ins_cross.metric_ms)
    mf.apply_match_test(ins_cross)

    cp.INS_plot(ins_cross, f'{args.output_prefix}flagged', **ins_plot_args)
    occ_dict = ssins_util.calc_occ(ins_cross, mf, 0)
    with open(f'{args.output_prefix}ssins_occ.json', "w") as json_file:
        json.dump(occ_dict, json_file, indent=4)

    # with open(f'{args.output_prefix}ssins_events.json', "w") as events_file:
    #     events_obj = [
    #         {
    #             'time_bounds': [int(event[0].start), int(event[0].stop)],
    #             'freq_bounds': [int(event[1].start), int(event[1].stop)],
    #             'shape': event[2],
    #             'sig': event[3] if event[3] is None else float(event[3])
    #         } for event in ins_cross.match_events
    #     ]
    #     json.dump(events_obj, events_file, indent=4)

    ins_cross.write(f'{args.output_prefix}', output_type='mask', clobber=True)

    # flags = ins_cross.mask_to_flags()
    # uvd = UVData()
    # uvd.read(args.uvfits, times=times)
    # uvf = UVFlag(uvd, waterfall=True, mode='flag')
    # uvf = ins_cross.flag_uvf(uvf)
    # uvutils.apply_uvflag(uvd, uvf)

    # params.sky_chans.collect { ch ->
    #          "\"\${ch}\": [\${(ch-0.51)*1.28e8}, \${(ch+0.51)*1.28e8}]".toString()
    #      }


if __name__ == '__main__':
    main()
