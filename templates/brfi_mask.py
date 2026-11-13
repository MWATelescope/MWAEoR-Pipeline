#!/usr/bin/env python

"""
Whole-band time mask for one observation from the night-level bayesian-rfi fits (brfi_fit.py).

Rule (the VS mask): in a fitted coarse channel, a time-differenced sample is RFI if P(RFI) >=
threshold in both XX and YY. It flags the later of the two timesteps it spans. Timesteps are
grouped into steps of `--step-s` seconds from the start of the observation (the resolution the
visibilities are averaged to later), and a step is RFI in that channel if at least `--step-frac`
of its timesteps are (2 of 4 for 2 s data and 8 s steps). The whole band is flagged at a step if
any fitted coarse channel says RFI.

The mask is written in the layout of an SSINS mask (`Header/freq_array`, `Header/time_array`,
`Data/flag_array` [time, freq, pol]) on the time grid of the visibilities, so that ssins_apply.py
can apply it.

example:
    python templates/brfi_mask.py --label 1472670968 --ins birli_1472670968_INS.h5 --fits brfi_2026-09-05_ch*.npz --out birli_1472670968_brfi_mask.h5 --json birli_1472670968_brfi_occ.json
"""

import json
import shlex
import sys

import h5py
import numpy as np


def parse_args(argv=None):
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--label', required=True, help="this observation's label in the fits")
    parser.add_argument('--ins', required=True, help='the INS h5 of this observation')
    parser.add_argument('--fits', nargs="+", required=True, help="the night's brfi_fit.py outputs, one per coarse channel")
    parser.add_argument('--out', required=True, help='output mask h5')
    parser.add_argument('--json', required=True, help='output occupancy json')
    parser.add_argument('--thresh', type=float, default=0.9, help='P(RFI) needed in both XX and YY')
    parser.add_argument('--step-s', type=float, default=8.0, help='length of the steps flagged as a whole (s)')
    parser.add_argument('--step-frac', type=float, default=0.5, help='fraction of timesteps that flags a step')
    return parser.parse_args(argv)


def step_mask(p_xx, p_yy, thresh, per_step, step_frac):
    """Flags of the timesteps (one more than the differenced samples), whole steps at a time."""
    rfi = (np.nan_to_num(p_xx) >= thresh) & (np.nan_to_num(p_yy) >= thresh)
    samp = np.zeros(len(rfi) + 1, bool)
    samp[1:] |= rfi
    need = max(1, int(np.ceil(step_frac * per_step - 1e-9)))
    padded = np.pad(samp, (0, (-len(samp)) % per_step)).reshape(-1, per_step)
    steps = padded.sum(axis=1) >= need
    return np.repeat(steps, per_step)[:len(samp)], steps


def main():
    if len(sys.argv) > 1:
        args = parse_args()
    else:
        # is being called directly from nextflow
        args = parse_args(shlex.split("""${args}"""))
    print(vars(args))

    with h5py.File(args.ins, "r") as h5:
        mid_jd = np.sort(np.asarray(h5["Header/time_array"], float))
        freqs = np.asarray(h5["Header/freq_array"], float).reshape(-1)
        npols = len(np.asarray(h5["Header/polarization_array"]).reshape(-1))
    if len(mid_jd) > 1:
        int_jd = float(np.median(np.diff(mid_jd)))
    else:
        raise UserWarning(f"{args.ins} has fewer than two differenced samples")
    # the timesteps the differenced samples span
    times_jd = np.concatenate([[mid_jd[0] - int_jd / 2], mid_jd + int_jd / 2])
    per_step = max(1, int(round(args.step_s / (int_jd * 86400))))

    by_chan, flagged = {}, np.zeros(len(times_jd), bool)
    for path in sorted(args.fits):
        with np.load(path) as fit:
            sel = fit["label"] == args.label
            if not sel.any():
                raise UserWarning(f"{path} has no samples for {args.label}")
            order = np.argsort(fit["time_jd"][sel])
            if len(order) != len(mid_jd):
                raise UserWarning(f"{path} has {len(order)} samples for {args.label}, {args.ins} has {len(mid_jd)}")
            chan = int(fit["coarse_chan"])
            times, steps = step_mask(fit["p_rfi_xx"][sel][order], fit["p_rfi_yy"][sel][order],
                                     args.thresh, per_step, args.step_frac)
        by_chan[chan] = [int(i) for i in np.flatnonzero(steps)]
        flagged |= times

    flag_array = np.zeros((len(times_jd), len(freqs), npols), bool)
    flag_array[flagged] = True
    with h5py.File(args.out, "w") as h5:
        header = h5.create_group("Header")
        header["freq_array"] = freqs
        header["time_array"] = times_jd
        h5.create_group("Data")["flag_array"] = flag_array

    occ = dict(
        label=args.label, thresh=args.thresh, step_s=args.step_s, timesteps_per_step=per_step,
        ntimes=len(times_jd), nflagged=int(flagged.sum()), occupancy=float(flagged.mean()),
        flagged_timesteps=[int(i) for i in np.flatnonzero(flagged)],
        flagged_steps_by_coarse_chan={str(c): s for c, s in sorted(by_chan.items()) if s},
        fitted_coarse_chans=sorted(by_chan),
    )
    with open(args.json, "w") as f:
        json.dump(occ, f, indent=1)
    print(json.dumps(occ))


if __name__ == '__main__':
    main()
