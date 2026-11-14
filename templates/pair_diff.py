#!/usr/bin/env python

"""
Power of the differences between the paired timesteps of a calibrated uvfits, over thermal noise.

CHIPS (gridvisdiff) grids a visibility file two timesteps at a time, (0,1), (2,3).., into a sum
grid and a difference grid, and the power spectrum is the power of the sum minus the power of the
difference. The difference is meant to hold noise only. Emission that changes between the two
timesteps of a pair (a television or radio transmitter that flickers, reflected off aircraft or
the ionosphere) fills the difference grid instead, and the power spectrum of the band that holds
it goes negative. One such observation is enough to spoil a night.

For every pol and frequency band this measures, over the same pairs CHIPS uses:

    d   = v[t+1] - v[t]                             the pair difference
    s2  = <|d[c+1] - d[c]|^2> / 2                   thermal variance, from adjacent channels of d
    r1  = <|d|^2> / s2                              1 for noise
    r   = n <|mean of d over n adjacent channels|^2> / s2

`r` is the gate value: 1 for noise, and up to n times larger for a difference that is smooth
across the n channels, which is what a broadcast signal looks like and what reaches the low
k_parallel modes. Anything that changes from one channel to the next counts as noise, so quiet
observations sit a few percent above 1. Only the baselines of the power spectrum cut (--uv-min to
--uv-max wavelengths) are used, and a visibility counts only when both timesteps of the pair have
a positive weight (all n channels for `r`).

Bands are --band-mhz wide and aligned to --band-origin-mhz (7 MHz from 174 MHz: the television
channels of VHF band III); each is trimmed to a whole number of n-channel blocks from its low edge.

example:
    python templates/pair_diff.py --uvfits hyp_1472671688_noflag_30l_src8k_300it_8s_80kHz_pairs.uvfits --out pair_diff.json
"""

import json
import shlex
import sys

import numpy as np
from astropy.io import fits

STOKES = {-5: "xx", -6: "yy"}


def parse_args(argv=None):
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--uvfits', required=True)
    parser.add_argument('--out', required=True, help='output json')
    parser.add_argument('--band-mhz', type=float, default=7.0, help='width of a band')
    parser.add_argument('--band-origin-mhz', type=float, default=174.0, help='an edge of the bands')
    parser.add_argument('--avg-khz', type=float, default=640.0, help='width the difference is averaged over')
    parser.add_argument('--uv-min', type=float, default=20.0, help='shortest baseline [wavelengths]')
    parser.add_argument('--uv-max', type=float, default=50.0, help='longest baseline [wavelengths]')
    return parser.parse_args(argv)


def pair_diff(path, band_mhz=7.0, band_origin_mhz=174.0, avg_khz=640.0, uv_min=20.0, uv_max=50.0):
    with fits.open(path, memmap=True) as hdul:
        hdr = hdul[0].header
        offset = hdul[0].fileinfo()["datLoc"]
    if hdr["BITPIX"] != -32:
        raise UserWarning(f"{path}: expected 32 bit floats, BITPIX={hdr['BITPIX']}")
    gcount, pcount, naxis = hdr["GCOUNT"], hdr["PCOUNT"], hdr["NAXIS"]
    # axes of one group as numpy sees them: the last fits axis first
    fits_axes = list(range(naxis, 1, -1))
    ctypes = [hdr.get(f"CTYPE{n}", "").strip() for n in fits_axes]
    shape = tuple(hdr[f"NAXIS{n}"] for n in fits_axes)
    for want in ("FREQ", "STOKES", "COMPLEX"):
        if want not in ctypes:
            raise UserWarning(f"{path} has no {want} axis")
    if shape[ctypes.index("COMPLEX")] < 3:
        raise UserWarning(f"{path} has no weights")
    groups = np.memmap(path, dtype=">f4", mode="r", offset=offset, shape=(gcount, pcount + int(np.prod(shape))))

    def axis_values(ctype):
        n = fits_axes[ctypes.index(ctype)]
        return hdr[f"CRVAL{n}"] + (np.arange(hdr[f"NAXIS{n}"]) + 1 - hdr[f"CRPIX{n}"]) * hdr[f"CDELT{n}"]

    def param(sel, names):
        out = np.zeros(len(sel))
        for i in range(1, pcount + 1):
            if hdr[f"PTYPE{i}"].strip() in names:
                out += sel[:, i - 1].astype(float) * hdr.get(f"PSCAL{i}", 1.0) + hdr.get(f"PZERO{i}", 0.0)
        return out

    freqs = axis_values("FREQ")
    chan_hz = abs(freqs[1] - freqs[0])
    pols = {STOKES[int(s)]: i for i, s in enumerate(axis_values("STOKES")) if int(s) in STOKES}
    if not pols:
        raise UserWarning(f"{path} has neither xx nor yy")
    nav = max(2, int(round(avg_khz * 1e3 / chan_hz)))

    # baselines of the first timestep. every timestep must hold the same ones, in the same order,
    # for row a of one timestep to pair with row a of the next
    nbl = gcount
    jd0 = param(groups[:1], ("DATE", "_DATE"))[0]
    for a in range(0, gcount, 2048):
        changed = np.flatnonzero(np.abs(param(groups[a:a + 2048], ("DATE", "_DATE")) - jd0) > 0.01 / 86400)
        if len(changed):
            nbl = a + int(changed[0])
            break
    if gcount % nbl:
        raise UserWarning(f"{path}: {gcount} rows are not a whole number of timesteps of {nbl} baselines")
    nsteps = gcount // nbl
    ref_hz = hdr[f"CRVAL{fits_axes[ctypes.index('FREQ')]}"]
    uvdist = np.hypot(param(groups[:nbl], ("UU", "UU---SIN")), param(groups[:nbl], ("VV", "VV---SIN"))) * ref_hz
    cut = (uvdist >= uv_min) & (uvdist <= uv_max)

    # bands: [label, channels]
    width, origin = band_mhz * 1e6, band_origin_mhz * 1e6
    low, high = freqs.min() - chan_hz / 2, freqs.max() + chan_hz / 2
    bands = []
    for n in range(int(np.floor((low - origin) / width)), int(np.ceil((high - origin) / width))):
        lo, hi = origin + n * width, origin + (n + 1) * width
        chans = np.flatnonzero((freqs >= lo) & (freqs < hi))
        chans = chans[:(len(chans) // nav) * nav]
        if len(chans):
            edges = [round(max(lo, low) / 1e6, 1), round(min(hi, high) / 1e6, 1)]
            bands.append(["%g-%g" % tuple(edges), chans])

    # move the freq, stokes and complex axes to the end, in that order
    order = [i for i, c in enumerate(ctypes) if c not in ("FREQ", "STOKES", "COMPLEX")]
    order += [ctypes.index(c) for c in ("FREQ", "STOKES", "COMPLEX")]

    def step(t):
        rows = np.asarray(groups[t * nbl:(t + 1) * nbl][cut][:, pcount:], dtype=np.float32)
        vis = rows.reshape((len(rows),) + shape).transpose([0] + [1 + i for i in order])
        return vis.reshape(len(rows), len(freqs), shape[ctypes.index("STOKES")], -1)

    # sums and counts of |d|^2, |d[c+1] - d[c]|^2 and |block mean of d|^2
    acc = {(band, pol): np.zeros(6) for band, _c in bands for pol in pols}
    npairs = 0
    for t in range(0, nsteps - 1, 2):
        a, b = step(t), step(t + 1)
        npairs += 1
        for pol, p in pols.items():
            good = (a[:, :, p, 2] > 0) & (b[:, :, p, 2] > 0)
            diff = (b[:, :, p, 0] - a[:, :, p, 0]) + 1j * (b[:, :, p, 1] - a[:, :, p, 1])
            for band, chans in bands:
                x, g = diff[:, chans], good[:, chans]
                e, ge = x[:, 1:] - x[:, :-1], g[:, 1:] & g[:, :-1]
                xn = x.reshape(len(x), -1, nav).mean(axis=2)
                gn = g.reshape(len(g), -1, nav).all(axis=2)
                acc[(band, pol)] += [
                    float((np.abs(x[g]) ** 2).sum()), g.sum(),
                    float((np.abs(e[ge]) ** 2).sum()), ge.sum(),
                    float((np.abs(xn[gn]) ** 2).sum()), gn.sum(),
                ]

    out = dict(
        uvfits=str(path), nsteps=int(nsteps), npairs=npairs, nbls=int(cut.sum()), avg_chans=nav,
        bands=[band for band, _c in bands], r={}, r1={},
    )
    worst = None
    for pol in pols:
        out["r"][pol], out["r1"][pol] = [], []
        for band, _c in bands:
            sx, nx, se, ne, sn, nn = acc[(band, pol)]
            s2 = se / ne / 2 if ne and se else None
            r1 = float(sx / nx / s2) if s2 and nx else None
            r = float(nav * sn / nn / s2) if s2 and nn else None
            out["r1"][pol].append(r1)
            out["r"][pol].append(r)
            if r is not None and (worst is None or r > worst[0]):
                worst = (r, pol, band)
    out["max_r"], out["max_pol"], out["max_band"] = worst or (None, None, None)
    return out


def main():
    if len(sys.argv) > 1:
        args = parse_args()
    else:
        # is being called directly from nextflow
        args = parse_args(shlex.split("""${args}"""))
    print(vars(args))
    out = pair_diff(args.uvfits, args.band_mhz, args.band_origin_mhz, args.avg_khz, args.uv_min, args.uv_max)
    with open(args.out, "w") as f:
        json.dump(out, f, indent=1)
    print(json.dumps(out))


if __name__ == '__main__':
    main()
