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
every visibility, and di-cal leaves it undetermined per channel. The YY terms then still share a
phase per channel, that of YY against XX, which the unpolarised sky hardly constrains. It cancels
in the YY visibilities, so YY is fitted without it and it is put back when the solutions are
written.

Flagged (NaN) values stay NaN. A tile keeps its own solutions when it is seen in fewer than
--min-obs observations of the night, or when its fine structure changes over the night
(--max-variable): for such a tile the median of the night describes no single observation.

example:
    python templates/calfit.py night --out calfit_night.npz --json calfit_night.json --solns hyp_soln_*.fits
    python templates/calfit.py apply --soln hyp_soln_1443896016_30l.fits --fit calfit_night.npz --out hyp_soln_1443896016_30l_nfit.fits
"""

import json
import os
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
    apply = sub.add_parser("apply")
    apply.add_argument('--soln', required=True)
    apply.add_argument('--fit', required=True, help='npz from the night mode')
    apply.add_argument('--out', required=True, help='output solutions')
    # the same options for both, so that one set of args serves both. apply uses the orders only
    for p in (night, apply):
        p.add_argument('--min-obs', type=int, default=5, help='observations a tile needs for a median')
        p.add_argument(
            '--max-variable', type=float, default=2.0,
            help='a tile and pol keeps its own solutions when the part of its fine structure that varies over '
            'the night, averaged over --variable-chans channels, exceeds this many times the median over tiles '
            '(amplitude or phase). 0 to fit every tile',
        )
        p.add_argument('--variable-chans', type=int, default=8, help='channels averaged for --max-variable')
        p.add_argument(
            '--max-obs-variable', type=float, default=5.0,
            help='an observation is left out of the night median, and listed as an outlier, when the median over '
            'tiles of its variable amplitude (averaged over --variable-chans channels) exceeds this many times the '
            'median over observations. 0 to keep every observation',
        )
        p.add_argument('--amp-order', type=int, default=3, help='order of the amplitude polynomial')
        p.add_argument('--phase-order', type=int, default=1, help='order of the phase polynomial')
    return parser.parse_args(argv)


def common_phase(gain, x, order):
    """Phasor [timeblock, chan] of the phase that the tiles of gain[timeblock, tile, chan] share in a channel."""
    # first the mean phasor, when the tiles agree well enough for it to mean something, then the
    # median over tiles of what a smooth fit per tile leaves
    out = np.ones((gain.shape[0], gain.shape[2]), dtype=complex)
    ref = np.nanmean(gain / np.abs(gain), axis=1)
    if np.nanmedian(np.abs(ref)) > 0.5:
        out = np.where(np.isfinite(ref) & (np.abs(ref) > 0), ref / np.abs(ref), 1.0)
    gain = gain * np.conj(out)[:, None, :]
    for b in range(gain.shape[0]):
        resid = np.full(gain.shape[1:], np.nan)
        for t in range(gain.shape[1]):
            phase = unwrapped(np.angle(gain[b, t]))
            resid[t] = phase - smooth_fit(x, phase, order)
        common = np.nanmedian(resid, axis=0)
        out[b] = out[b] * np.exp(1j * np.where(np.isfinite(common), common, 0.0))
    return out


def read_soln(path, phase_order=1):
    """
    (gains[pol, timeblock, tile, chan], freqs, tile names, jones[timeblock, tile, chan, 4], yy phasor[timeblock, chan])

    The phase common to the xx terms of all tiles is removed from the gains and from the jones
    matrices. The phase common to the yy terms that is left is removed from the yy gains only.
    """
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
        cplx = cplx * np.conj(common_phase(cplx[..., 0], x, phase_order))[:, None, :, None]
        common_yy = common_phase(cplx[..., 3], x, phase_order)
    gains = np.stack([cplx[..., 0], cplx[..., 3] * np.conj(common_yy)[:, None, :]])
    return gains, freqs, names, cplx, common_yy


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
    obs, yy_phase, src = [], [], []
    freqs = names = None
    for path in args.solns:
        gains, f, n, _c, common_yy = read_soln(path, args.phase_order)
        if freqs is None:
            freqs, names = f, n
        elif len(f) != len(freqs) or not np.allclose(f, freqs) or n != names:
            raise UserWarning(f"{path}: channels or tiles differ from {args.solns[0]}")
        # every timeblock is a sample of the night
        obs.extend(gains[:, t] for t in range(gains.shape[1]))
        yy_phase.extend(np.angle(common_yy[t]) for t in range(gains.shape[1]))
        src.extend(path for _ in range(gains.shape[1]))
    obs = np.array(obs)  # [sample, pol, tile, chan]
    yy_phase = np.array(yy_phase)  # [sample, chan]
    nobs, npol, ntile, nchan = obs.shape
    x = freq_axis(freqs)
    seen = np.isfinite(obs).any(axis=3).sum(axis=0)  # [pol, tile]
    few = seen < args.min_obs
    # tiles without a median are fitted against a flat fine structure
    fine_amp = np.full((npol, ntile, nchan), np.nan)
    fine_phase = np.full((npol, ntile, nchan), np.nan)
    with np.errstate(all="ignore"), warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        def block_mean(a):
            n = max(1, args.variable_chans)
            return np.nanmean(a[..., :(nchan // n) * n].reshape(a.shape[:-1] + (nchan // n, n)), axis=-1)

        def night_median(use):
            # second pass: the smooth parts are refitted against the first median
            fa = np.full((npol, ntile, nchan), np.nan)
            fp = np.full((npol, ntile, nchan), np.nan)
            for _ in range(2):
                use_amp = np.where(np.isfinite(fa), fa, 1.0)
                use_phase = np.where(np.isfinite(fp), fp, 0.0)
                ra = np.full(obs.shape, np.nan)
                rp = np.full(obs.shape, np.nan)
                for i in range(nobs):
                    for p in range(npol):
                        for t in range(ntile):
                            if not np.isfinite(obs[i, p, t]).any():
                                continue
                            sa, sp, amp, phase = split(obs[i, p, t], x, use_amp[p, t], use_phase[p, t], args.amp_order, args.phase_order)
                            ra[i, p, t] = amp / sa
                            rp[i, p, t] = phase - sp
                fa = np.nanmedian(ra[use], axis=0)
                fp = np.nanmedian(rp[use], axis=0)
                fa[few] = np.nan
                fp[few] = np.nan
            return fa, fp, ra, rp

        # an observation whose solutions differ from the rest of the night across most tiles (a
        # different sky, beam or hardware state) would pull the median, and its data is suspect: it is
        # left out of the median and listed. Measured per sample as the median over tiles of the rms
        # over channel blocks of its variable amplitude, against the median over samples
        use = np.ones(nobs, dtype=bool)
        fine_amp, fine_phase, res_amp, res_phase = night_median(use)
        obs_var = np.nanmedian(np.sqrt(np.nanmean(block_mean(res_amp / fine_amp[None] - 1) ** 2, axis=3)), axis=2)  # [sample, pol]
        typical_obs = np.nanmedian(obs_var, axis=0)  # [pol]
        if args.max_obs_variable > 0:
            use = ~(obs_var > args.max_obs_variable * typical_obs[None]).any(axis=1)
            if not use.all() and use.sum() > 0:
                fine_amp, fine_phase, res_amp, res_phase = night_median(use)
        outlier_src = sorted({src[i] for i in np.flatnonzero(~use)})
        # an outlier is not in any statistic below
        res_amp[~use] = np.nan
        res_phase[~use] = np.nan
        nobs_used = int(use.sum())
        var_amp = res_amp / fine_amp[None] - 1
        var_phase = res_phase - fine_phase[None]

        # a tile whose fine structure changes over the night: what varies is averaged over a few
        # channels first, so that a tile which is only noisier than the others still gets the median
        def block_rms(a):
            return np.sqrt(np.nanmean(block_mean(a) ** 2, axis=(0, 3)))  # [pol, tile]

        tile_var = np.stack([block_rms(var_amp), block_rms(var_phase)])  # [(amp, phase), pol, tile]
        typical = np.nanmedian(tile_var, axis=2, keepdims=True)
        unstable = np.zeros((npol, ntile), dtype=bool)
        if args.max_variable > 0:
            unstable = (tile_var > args.max_variable * typical).any(axis=0)
        fine_amp[unstable] = np.nan
        fine_phase[unstable] = np.nan
        var_amp[:, unstable] = np.nan
        var_phase[:, unstable] = np.nan
        fitted = np.isfinite(fine_amp).any(axis=2)

        def med_rms(a, axes):
            return float(np.nanmedian(np.sqrt(np.nanmean(a ** 2, axis=axes))))

        def line_free(a):
            return np.array([row - smooth_fit(x, row, args.phase_order) for row in a])

        summary = dict(
            nobs=int(nobs), ntiles=int(ntile), nchans=int(nchan),
            # observations left out of the median (solution files), and the measure per pol: median
            # over samples, and the largest
            obs_outliers=[os.path.basename(str(o)) for o in outlier_src], nobs_used=nobs_used,
            obs_variable=[[float(typical_obs[p]), float(np.nanmax(obs_var[:, p]))] for p in range(npol)],
            max_obs_variable=args.max_obs_variable,
            tiles_fitted=[int(fitted[p].sum()) for p in range(npol)],
            # tiles that keep their own solutions: [name, variable amplitude, variable phase in
            # degrees] (both averaged over variable_chans channels), next to the median over tiles
            tiles_variable=[
                [[names[t], float(tile_var[0, p, t]), float(np.degrees(tile_var[1, p, t]))] for t in np.flatnonzero(unstable[p])]
                for p in range(npol)
            ],
            typical_variable=[[float(typical[0, p, 0]), float(np.degrees(typical[1, p, 0]))] for p in range(npol)],
            max_variable=args.max_variable, variable_chans=args.variable_chans,
            # the phase all yy terms share against xx, after a line: rms over channels of its
            # median over the night, and of what is left per observation
            yy_common_phase_deg=[
                float(np.degrees(np.nanstd(np.nanmedian(line_free(yy_phase), axis=0)))),
                med_rms(np.degrees(line_free(yy_phase) - np.nanmedian(line_free(yy_phase), axis=0)), 1),
            ],
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
    gains, freqs, names, cplx, common_yy = read_soln(args.soln, args.phase_order)
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
                    if p == 1:
                        new = new * common_yy[b]
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
