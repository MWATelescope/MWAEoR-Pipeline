#!/usr/bin/env python

"""
Mean visibility weight of every timestep of a uvfits file.

The incoherent noise spectrum (INS) of time-differenced visibilities is only comparable between
samples when both timesteps of a difference hold the same weight. A timestep that was averaged from
fewer native samples (the first timestep of an observation, when its first native sample is
flagged) has more noise, and its difference stands out without any RFI. The INS file holds baseline
counts, not weights, so brfi_fit.py takes the weights from here and leaves such differences out.

Reads a sample of the visibilities (about `--rows` baselines-times, every `--chan-step`th channel)
straight from the random groups, without loading the file.

example:
    python templates/brfi_weights.py --uvfits birli_1472670968_8s_80kHz.uvfits --out birli_1472670968_8s_80kHz_time_weights.json
"""

import json
import shlex
import sys

import numpy as np
from astropy.io import fits


def parse_args(argv=None):
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--uvfits', required=True)
    parser.add_argument('--out', required=True, help='output json')
    parser.add_argument('--rows', type=int, default=8000, help='about how many rows (baseline-times) to read')
    parser.add_argument('--chan-step', type=int, default=8, help='read every nth fine channel')
    return parser.parse_args(argv)


def time_weights(path, rows=8000, chan_step=8):
    """(timestep centres [JD], mean weight of the unflagged visibilities, flagged fraction)."""
    with fits.open(path, memmap=True) as hdul:
        hdr = hdul[0].header
        offset = hdul[0].fileinfo()["datLoc"]
    if hdr["BITPIX"] != -32:
        raise UserWarning(f"{path}: expected 32 bit floats, BITPIX={hdr['BITPIX']}")
    gcount, pcount, naxis = hdr["GCOUNT"], hdr["PCOUNT"], hdr["NAXIS"]
    # axes of one group as numpy sees them: the last fits axis first
    axes = [(hdr.get(f"CTYPE{n}", "").strip(), hdr[f"NAXIS{n}"]) for n in range(naxis, 1, -1)]
    nelem = int(np.prod([n for _c, n in axes]))
    groups = np.memmap(path, dtype=">f4", mode="r", offset=offset, shape=(gcount, pcount + nelem))

    def row_jd(sel):
        jd = np.zeros(len(sel))
        for i in range(1, pcount + 1):
            if hdr[f"PTYPE{i}"].strip() in ("DATE", "_DATE"):
                jd += sel[:, i - 1].astype(float) * hdr.get(f"PSCAL{i}", 1.0) + hdr.get(f"PZERO{i}", 0.0)
        return jd

    index = [slice(None)]
    for ctype, n in axes:
        if ctype == "COMPLEX":
            if n < 3:
                raise UserWarning(f"{path} has no weights")
            index.append(2)
        elif ctype == "FREQ":
            index.append(slice(None, None, max(1, chan_step)))
        else:
            index.append(slice(None))

    def row_weights(sel):
        vis = sel[:, pcount:].reshape((len(sel),) + tuple(n for _c, n in axes))
        return np.asarray(vis[tuple(index)], float).reshape(len(sel), -1)

    # rows of the first timestep: with the usual layout (all baselines of a timestep together,
    # the same number in each) the same baselines are read in every timestep, so that flagged
    # tiles weigh the same in all of them
    nbl = gcount
    for a in range(0, gcount, 2048):
        jd = row_jd(groups[a:a + 2048])
        changed = np.flatnonzero(np.abs(jd - row_jd(groups[:1])[0]) > 0.01 / 86400)
        if len(changed):
            nbl = a + int(changed[0])
            break
    step = max(1, gcount // max(1, rows))
    if gcount % nbl == 0:
        # at least 500 baselines of every timestep
        step = max(1, min(step, nbl // 500))
        parts = [groups[a:a + nbl:step] for a in range(0, gcount, nbl)]
    else:
        parts = [groups[::step]]
    jd = np.concatenate([row_jd(p) for p in parts])
    weights = np.concatenate([row_weights(p) for p in parts])

    # timesteps: runs of rows less than a quarter of a timestep apart
    order = np.argsort(jd, kind="stable")
    jd, weights = jd[order], weights[order]
    steps = np.diff(jd)
    steps = steps[steps > 0.01 / 86400]
    gap = (np.median(steps) if len(steps) else 1.0) / 4
    starts = np.concatenate([[0], np.flatnonzero(np.diff(jd) > gap) + 1, [len(jd)]])
    times, mean_weight, flagged = [], [], []
    for a, b in zip(starts[:-1], starts[1:]):
        w = weights[a:b]
        good = w > 0
        times.append(float(jd[a:b].mean()))
        mean_weight.append(float(w[good].mean()) if good.any() else 0.0)
        flagged.append(float(1 - good.mean()))
    return times, mean_weight, flagged


def main():
    if len(sys.argv) > 1:
        args = parse_args()
    else:
        # is being called directly from nextflow
        args = parse_args(shlex.split("""${args}"""))
    print(vars(args))
    times, mean_weight, flagged = time_weights(args.uvfits, args.rows, args.chan_step)
    out = dict(uvfits=args.uvfits, time_jd=times, weight=mean_weight, flagged_frac=flagged)
    with open(args.out, "w") as f:
        json.dump(out, f, indent=1)
    print(json.dumps(out))


if __name__ == '__main__':
    main()
