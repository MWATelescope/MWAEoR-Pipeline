#!/usr/bin/env python

"""
Replace the fine frequency structure of per-channel calibration solutions by that of the night.

A sky calibration that solves every channel on its own puts the errors of an incomplete sky model
into the gains as frequency structure. That structure changes with the sky in the beam, from one
observation to the next, and it carries bright foregrounds from low to higher k_parallel. The
bandpass of the instrument itself (analogue chain, cable reflections) does not change over a night.

So each gain is written as (smooth) x (fine), per tile and pol:

    amplitude = polynomial(freq) x fine_amp(freq)
    phase     = line(freq)       + fine_phase(freq)

`night` takes all solutions of a night and stores the median fine structure over the observations
(per tile, pol and channel). `apply` refits the smooth part of one observation against that median
and writes a solution file in which the observation's own fine structure is replaced by the
night's. The smooth part stays per observation, so gain drifts are kept.

Only the XX and YY terms are fitted; the cross terms are copied. A phase common to all tiles in a
channel is removed first (from the XX terms, applied to the whole Jones matrix): it cancels in
every visibility, and di-cal leaves it undetermined per channel. Flagged (NaN) values stay NaN, and a
tile seen in fewer than --min-obs observations of the night keeps its own solutions.

example:
    python templates/calfit.py night --out calfit_night.npz --json calfit_night.json --solns hyp_soln_*.fits
    python templates/calfit.py apply --soln hyp_soln_1443896016_30l.fits --fit calfit_night.npz --out hyp_soln_1443896016_30l_nfit.fits
"""

import json
import shlex
import sys
import warnings

import numpy as np
from astropy.io import fits

POLS = {"xx": 0, "yy": 6}


def parse_args(argv=None):
    import argparse
    parser = argparse.ArgumentParser()
    sub = parser.add_subparsers(dest="mode", required=True)
    night = sub.add_parser("night")
    night.add_argument('--solns', nargs='+', required=True, help='hyperdrive solutions of one night')
    night.add_argument('--out', required=True, help='output npz')
    night.add_argument('--json', required=True, help='output summary')
    night.add_argument('--min-obs', type=int, default=5, help='observations a tile needs for a median')
    apply = sub.add_parser("apply")
    apply.add_argument('--soln', required=True)
    apply.add_argument('--fit', required=True, help='npz from the night mode')
    apply.add_argument('--out', required=True, help='output solutions')
    for p in (night, apply):
        p.add_argument('--amp-order', type=int, default=3, help='order of the amplitude polynomial')
        p.add_argument('--phase-order', type=int, default=1, help='order of the phase polynomial')
    return parser.parse_args(argv)


def read_soln(path, phase_order=1):
    """(gains[pol, timeblock, tile, chan] complex with the common phase removed, freqs, tile names, jones)."""
    with fits.open(path) as hdul:
        jones = np.array(hdul["SOLUTIONS"].data, dtype=float)
        freqs = np.array(hdul["CHANBLOCKS"].data["Freq"], dtype=float)
        names = [str(n).strip() for n in hdul["TILES"].data["TileName"]]
    if jones.ndim != 4 or jones.shape[-1] != 8:
        raise UserWarning(f"{path}: unexpected SOLUTIONS shape {jones.shape}")
    cplx = jones[..., 0::2] + 1j * jones[..., 1::2]  # [timeblock, tile, chan, (xx, xy, yx, yy)]
    x = freq_axis(freqs)
    with np.errstate(all="ignore"), warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        # a phase common to all tiles in a channel: first the mean xx phasor, when the tiles agree
        # well enough for it to mean something, then the median over tiles of what a smooth fit
        # per tile leaves
        unit = cplx[..., 0] / np.abs(cplx[..., 0])
        ref = np.nanmean(unit, axis=1, keepdims=True)  # [timeblock, 1, chan]
        if np.nanmedian(np.abs(ref)) > 0.5:
            ref = np.where(np.isfinite(ref) & (np.abs(ref) > 0), ref / np.abs(ref), 1.0)
            cplx = cplx * np.conj(ref)[..., None]
        for b in range(cplx.shape[0]):
            resid = np.full(cplx.shape[1:3], np.nan)
            for t in range(cplx.shape[1]):
                phase = unwrapped(np.angle(cplx[b, t, :, 0]))
                resid[t] = phase - smooth_fit(x, phase, phase_order)
            common = np.nanmedian(resid, axis=0)
            common = np.where(np.isfinite(common), common, 0.0)
            cplx[b] = cplx[b] * np.exp(-1j * common)[None, :, None]
    gains = np.stack([cplx[..., 0], cplx[..., 3]])
    return gains, freqs, names, cplx


def smooth_fit(x, y, order, clip=4.0):
    """Polynomial of y over x on the finite values, with outliers left out. NaN when too few."""
    good = np.isfinite(y)
    if good.sum() < 4 * (order + 1):
        return np.full(len(x), np.nan)
    coeffs = None
    for _ in range(3):
        coeffs = np.polynomial.polynomial.polyfit(x[good], y[good], order)
        resid = y - np.polynomial.polynomial.polyval(x, coeffs)
        sigma = 1.4826 * np.nanmedian(np.abs(resid[good] - np.nanmedian(resid[good])))
        keep = np.isfinite(y) & (np.abs(resid) <= clip * sigma) if sigma > 0 else np.isfinite(y)
        if keep.sum() < 4 * (order + 1) or keep.sum() == good.sum():
            break
        good = keep
    return np.polynomial.polynomial.polyval(x, coeffs)


def unwrapped(phase):
    """Unwrap along the last axis over the finite values only."""
    out = np.full(phase.shape, np.nan)
    good = np.isfinite(phase)
    if good.any():
        out[good] = np.unwrap(phase[good])
    return out


def split(gain, x, fine_amp, fine_phase, amp_order, phase_order):
    """(smooth amplitude, smooth phase, amplitude, unwrapped phase) of one tile and pol."""
    amp = np.abs(gain)
    phase = unwrapped(np.angle(gain))
    with np.errstate(invalid="ignore", divide="ignore"):
        smooth_amp = smooth_fit(x, amp / fine_amp, amp_order)
    smooth_phase = smooth_fit(x, phase - fine_phase, phase_order)
    return smooth_amp, smooth_phase, amp, phase


def freq_axis(freqs):
    return 2 * (freqs - freqs.mean()) / (freqs.max() - freqs.min())


def night(args):
    obs = []
    freqs = names = None
    for path in args.solns:
        gains, f, n, _c = read_soln(path, args.phase_order)
        if freqs is None:
            freqs, names = f, n
        elif len(f) != len(freqs) or not np.allclose(f, freqs) or n != names:
            raise UserWarning(f"{path}: channels or tiles differ from {args.solns[0]}")
        # every timeblock is a sample of the night
        obs.extend(gains[:, t] for t in range(gains.shape[1]))
    obs = np.array(obs)  # [sample, pol, tile, chan]
    nobs, npol, ntile, nchan = obs.shape
    x = freq_axis(freqs)
    seen = np.isfinite(obs).any(axis=3).sum(axis=0)  # [pol, tile]
    few = seen < args.min_obs
    # tiles without a median are fitted against a flat fine structure
    fine_amp = np.full((npol, ntile, nchan), np.nan)
    fine_phase = np.full((npol, ntile, nchan), np.nan)
    with np.errstate(all="ignore"), warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        # second pass: the smooth parts are refitted against the first median
        for _ in range(2):
            use_amp = np.where(np.isfinite(fine_amp), fine_amp, 1.0)
            use_phase = np.where(np.isfinite(fine_phase), fine_phase, 0.0)
            res_amp = np.full(obs.shape, np.nan)
            res_phase = np.full(obs.shape, np.nan)
            for i in range(nobs):
                for p in range(npol):
                    for t in range(ntile):
                        if not np.isfinite(obs[i, p, t]).any():
                            continue
                        sa, sp, amp, phase = split(obs[i, p, t], x, use_amp[p, t], use_phase[p, t], args.amp_order, args.phase_order)
                        res_amp[i, p, t] = amp / sa
                        res_phase[i, p, t] = phase - sp
            fine_amp = np.nanmedian(res_amp, axis=0)
            fine_phase = np.nanmedian(res_phase, axis=0)
            fine_amp[few] = np.nan
            fine_phase[few] = np.nan
        var_amp = res_amp / fine_amp[None] - 1
        var_phase = res_phase - fine_phase[None]

        def med_rms(a, axes):
            return float(np.nanmedian(np.sqrt(np.nanmean(a ** 2, axis=axes))))

        summary = dict(
            nobs=int(nobs), ntiles=int(ntile), nchans=int(nchan),
            tiles_fitted=[int((seen[p] >= args.min_obs).sum()) for p in range(npol)],
            amp_order=args.amp_order, phase_order=args.phase_order, min_obs=args.min_obs,
            # rms over channels, median over tiles: the night's fine structure, and what varies
            # from one observation to the next (which `apply` removes)
            stable_amp=[med_rms(fine_amp[p] - 1, 1) for p in range(npol)],
            stable_phase_deg=[med_rms(np.degrees(fine_phase[p]), 1) for p in range(npol)],
            variable_amp=[med_rms(var_amp[:, p], (0, 2)) for p in range(npol)],
            variable_phase_deg=[med_rms(np.degrees(var_phase[:, p]), (0, 2)) for p in range(npol)],
        )
    np.savez(
        args.out, freqs=freqs, tiles=np.array(names), fine_amp=fine_amp, fine_phase=fine_phase,
        nobs=seen, amp_order=args.amp_order, phase_order=args.phase_order,
    )
    with open(args.json, "w") as f:
        json.dump(summary, f, indent=1)
    print(json.dumps(summary))


def apply(args):
    fit = np.load(args.fit)
    gains, freqs, names, cplx = read_soln(args.soln, args.phase_order)
    if len(freqs) != len(fit["freqs"]) or not np.allclose(freqs, fit["freqs"]) or names != [str(n) for n in fit["tiles"]]:
        raise UserWarning(f"{args.soln}: channels or tiles differ from {args.fit}")
    x = freq_axis(freqs)
    fitted = kept = 0
    out = cplx.copy()  # [timeblock, tile, chan, (xx, xy, yx, yy)], common phase removed
    with np.errstate(invalid="ignore", divide="ignore"):
        for p, col in enumerate((0, 3)):
            for b in range(gains.shape[1]):
                for t in range(gains.shape[2]):
                    gain = gains[p, b, t]
                    fine_amp, fine_phase = fit["fine_amp"][p, t], fit["fine_phase"][p, t]
                    if not np.isfinite(gain).any():
                        continue
                    if not np.isfinite(fine_amp).any():
                        kept += 1
                        continue
                    sa, sp, _amp, _phase = split(gain, x, fine_amp, fine_phase, args.amp_order, args.phase_order)
                    new = sa * fine_amp * np.exp(1j * (sp + fine_phase))
                    if not np.isfinite(new[np.isfinite(gain)]).all():
                        kept += 1
                        continue
                    # flagged values stay flagged
                    out[b, t, :, col] = np.where(np.isfinite(gain), new, np.nan)
                    fitted += 1
    jones = np.empty(out.shape[:-1] + (8,))
    jones[..., 0::2] = out.real
    jones[..., 1::2] = out.imag
    with fits.open(args.soln) as hdul:
        hdul["SOLUTIONS"].data = jones.astype(hdul["SOLUTIONS"].data.dtype)
        hdul[0].header["CALFIT"] = (str(args.fit).split("/")[-1][:60], "fine structure of the night")
        hdul[0].header["CALFITN"] = (fitted, "tile pols fitted")
        hdul[0].header["CALFITK"] = (kept, "tile pols kept as solved")
        hdul.writeto(args.out, overwrite=True)
    print(json.dumps(dict(soln=args.soln, out=args.out, fitted=fitted, kept=kept)))


def main():
    if len(sys.argv) > 1:
        args = parse_args()
    else:
        # is being called directly from nextflow
        args = parse_args(shlex.split("""${args}"""))
    print(vars(args))
    if args.mode == "night":
        night(args)
    else:
        apply(args)


if __name__ == '__main__':
    main()
