#!/usr/bin/env python

"""
Take the fine frequency structure of calibration solutions from the autocorrelations.

A sky calibration that solves every channel on its own puts the errors of an incomplete sky model
into the gains as frequency structure, which carries bright foregrounds to higher k_parallel. The
autocorrelation of a tile measures the same analogue chain, |gain|^2 times the total power the
tile receives, with a hundred times the signal to noise of a per-channel solution and with no sky
model. The total power is smooth in frequency (a zero-length baseline), so everything in the
autocorrelation spectrum beyond a short delay is the gain of the tile (after Li et al. 2019, ApJ
887:141, section 5.2).

Per tile and pol (the XX and YY terms; the cross terms are copied), with ln g = ln|g| + i phase:

    ln g  =  lowpass[ ln g_cal ]  +  highpass[ z_auto ]

- z_auto = u + i v, where u = -1/2 ln(auto) is the log amplitude from the autocorrelation
  (hyperdrive stores the inverse gain) and v is the phase that goes with it for a causal response:
  an analogue chain (filters, cable reflections) only delays the signal, so the delay transform
  of ln g is one-sided and the phase ripple follows from the amplitude ripple (a Hilbert
  transform). No reflection modes are fitted and no phase is unwrapped across a fit, which is
  where earlier implementations got phase jumps.
- g_cal is the input solution. Its part below --delay-cut (in delay) is kept: the overall level
  and slope, and the slow structure in which the total power and the response to the sky differ
  (sky spectrum, beam, receiver noise, mutual coupling). Above the cut the solution's own
  structure, which is where the sky model errors of each channel live, is replaced.

The low pass is a Gaussian in frequency (half power at --delay-cut) on what a low-order
polynomial leaves, with outlying channels left out.

A tile and pol keeps its input solution when its autocorrelation does not describe its gain
(--max-mismatch): tiles with broken dipoles see a total power with other structure than their
response towards the field.

With --common-from-cal the fine structure that all tiles share is taken from the input solutions
instead (the average over tiles, as in Li et al.), which guards against RFI and other power common
to all autocorrelations at the price of the noise of that average.

example:
    python templates/autofit.py --uvfits birli_1443896016.uvfits --soln hyp_soln_1443896016_30l.fits --out hyp_soln_1443896016_30l_auto.fits
"""

import json
import shlex
import sys
import warnings

import numpy as np

# hyperdrive Jones order XX, XY, YX, YY; uvfits pol order XX, YY, XY, YX
POLS = (("xx", 0, 0), ("yy", 3, 1))


def parse_args(argv=None):
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--uvfits', required=True, help='uncalibrated visibilities with autocorrelations')
    parser.add_argument('--soln', required=True, help='hyperdrive solutions (fits)')
    parser.add_argument('--out', required=True, help='output solutions')
    parser.add_argument('--json', default=None, help='output summary')
    parser.add_argument(
        '--delay-cut', type=float, default=0.1,
        help='delay in microseconds below which the input solution is kept and above which the '
        'autocorrelations are used (half power of a Gaussian low pass)',
    )
    parser.add_argument('--amp-order', type=int, default=3, help='order of the polynomial taken out before the low pass (log amplitude)')
    parser.add_argument('--phase-order', type=int, default=3, help='order of the polynomial taken out before the low pass (phase)')
    parser.add_argument('--clip', type=float, default=5.0, help='channels of the input solution further than this many sigma from the low pass are left out of it')
    parser.add_argument(
        '--max-mismatch', type=float, default=5.0,
        help='a tile and pol keeps its input solution when the difference between the fine structure '
        'of its solution and of its autocorrelation, averaged over --mismatch-chans channels, exceeds '
        'this many times the median over tiles. 0 to fit every tile',
    )
    parser.add_argument('--mismatch-chans', type=int, default=8, help='channels averaged for --max-mismatch')
    parser.add_argument('--common-from-cal', action='store_true', help='fine structure shared by all tiles from the input solutions')
    parser.add_argument('--no-phase', action='store_true', help='leave the phase of the input solutions as it is')
    return parser.parse_args(argv)


def read_autos(path):
    """(autos[pol(xx, yy), time, tile, chan], weights[time, tile, chan], freqs, antenna names) of a uvfits file."""
    from astropy.io import fits
    with fits.open(path, memmap=True) as hdul:
        hd = hdul[0].header
        offset = hdul[0]._data_offset
        names = [str(n).strip() for n in hdul[1].data["ANNAME"]]
        nosta = np.array(hdul[1].data["NOSTA"]).astype(int)
    pcount, gcount = hd["PCOUNT"], hd["GCOUNT"]
    shape = tuple(hd[f"NAXIS{i}"] for i in range(hd["NAXIS"], 1, -1))  # (dec, ra, freq, stokes, complex)
    nchan, npol = shape[-3], shape[-2]
    stokes = [int(round(hd["CRVAL3"] + (i + 1 - hd["CRPIX3"]) * hd["CDELT3"])) for i in range(npol)]
    if stokes[:2] != [-5, -6] or shape[-1] != 3:
        raise UserWarning(f"{path}: expected pols XX, YY first and (re, im, weight), got stokes {stokes}, {shape}")
    dtype = np.dtype([("p", ">f4", (pcount,)), ("d", ">f4", shape)])
    groups = np.memmap(path, dtype=dtype, mode="r", offset=offset, shape=(gcount,))
    ptypes = [str(hd[f"PTYPE{i + 1}"]).strip() for i in range(pcount)]
    params = np.array(groups["p"]).astype(np.float64)
    for i in range(pcount):
        params[:, i] = params[:, i] * hd.get(f"PSCAL{i + 1}", 1.0) + hd.get(f"PZERO{i + 1}", 0.0)
    if "BASELINE" in ptypes:
        bl = params[:, ptypes.index("BASELINE")]
        ant1, ant2 = (bl // 256).astype(int), (bl % 256).astype(int)
    else:
        ant1 = params[:, ptypes.index("ANTENNA1")].astype(int)
        ant2 = params[:, ptypes.index("ANTENNA2")].astype(int)
    # the date is split over two parameters; the offset of the first (PZERO) is the same for all
    # rows, so the sum orders the timesteps
    date = sum(params[:, i] - hd.get(f"PZERO{i + 1}", 0.0) for i, t in enumerate(ptypes) if t == "DATE")
    rows = np.flatnonzero(ant1 == ant2)
    if not len(rows):
        raise UserWarning(f"{path}: no autocorrelations")
    times = np.unique(date[rows])
    ntime, ntile = len(times), len(names)
    tidx = np.searchsorted(times, date[rows])
    order = np.argsort(nosta)
    aidx = order[np.searchsorted(nosta[order], ant1[rows])]
    data = np.array(groups["d"][rows]).reshape(len(rows), nchan, npol, 3)
    autos = np.full((2, ntime, ntile, nchan), np.nan)
    weights = np.zeros((ntime, ntile, nchan))
    autos[0, tidx, aidx] = data[:, :, 0, 0]
    autos[1, tidx, aidx] = data[:, :, 1, 0]
    weights[tidx, aidx] = np.minimum(data[:, :, 0, 2], data[:, :, 1, 2])
    freqs = hd["CRVAL4"] + (np.arange(nchan) + 1 - hd["CRPIX4"]) * hd["CDELT4"]
    return autos, weights, freqs, names


def auto_spectrum(autos, weights):
    """Median over the unflagged timesteps, [pol, tile, chan]. NaN where nothing is left or it is not positive."""
    with np.errstate(all="ignore"), warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        spec = np.nanmedian(np.where((weights > 0)[None], autos, np.nan), axis=1)
    spec[~(spec > 0)] = np.nan
    return spec


def polyval_fit(x, y, order, good):
    return np.polynomial.polynomial.polyval(x, np.polynomial.polynomial.polyfit(x[good], y[good], order))


def gauss_kernel(sigma_chan):
    half = int(np.ceil(4 * sigma_chan))
    k = np.exp(-0.5 * (np.arange(-half, half + 1) / sigma_chan) ** 2)
    return k / k.sum()


def lowpass(x, y, good, order, kernel, clip=0.0):
    """
    Low pass of y over the channels `good`: a polynomial of the given order plus the Gaussian
    average of what it leaves. Channels that are not good (and, with clip, outliers) carry no
    weight, so they are bridged. Returns (low pass, channels used).
    """
    good = good & np.isfinite(y)
    if good.sum() < 4 * (order + 1):
        return np.full(len(y), np.nan), good
    use = good.copy()
    low = None
    for _ in range(3):
        poly = polyval_fit(x, y, order, use)
        resid = np.where(use, y - poly, 0.0)
        norm = np.convolve(use.astype(float), kernel, mode="same")
        with np.errstate(invalid="ignore", divide="ignore"):
            low = poly + np.convolve(resid, kernel, mode="same") / norm
        low[norm < 1e-3] = np.nan
        if clip <= 0:
            break
        dev = y - low
        sigma = 1.4826 * np.nanmedian(np.abs(dev[use] - np.nanmedian(dev[use])))
        keep = good & (np.abs(dev) <= clip * sigma) if sigma > 0 else good
        if keep.sum() == use.sum() or keep.sum() < 4 * (order + 1):
            break
        use = keep
    return low, use


def causal_phase(u):
    """
    The phase v for which u + i v has a one-sided delay transform (u real, last axis frequency):
    the response of a chain that only delays. The sign is that of hyperdrive solutions (the
    inverse gain, in the uvfits convention), checked against solutions in `enhance`. u is extended
    by its mirror image at both ends to keep the ends of the band clean.
    """
    n = u.shape[-1]
    pad = n // 2
    ext = np.concatenate([u[..., pad:0:-1], u, u[..., -2:-pad - 2:-1]], axis=-1)
    m = ext.shape[-1]
    h = np.zeros(m)
    h[0] = 1
    h[m // 2 + 1:] = 2
    if m % 2 == 0:
        h[m // 2] = 1
    return np.fft.ifft(np.fft.fft(ext, axis=-1) * h, axis=-1).imag[..., pad:pad + n]


def interp_bad(y, good):
    """y with the channels that are not good filled in linearly (the ends by the nearest value)."""
    idx = np.arange(len(y))
    return np.interp(idx, idx[good], y[good])


def unwrapped(phase):
    out = np.full(phase.shape, np.nan)
    good = np.isfinite(phase)
    if good.any():
        out[good] = np.unwrap(phase[good])
    return out


def common_phase(gain, x, order=1):
    """Phasor [chan] of the phase that the tiles of gain[tile, chan] share in a channel."""
    with np.errstate(all="ignore"), warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        out = np.ones(gain.shape[1], dtype=complex)
        ref = np.nanmean(gain / np.abs(gain), axis=0)
        if np.nanmedian(np.abs(ref)) > 0.5:
            out = np.where(np.isfinite(ref) & (np.abs(ref) > 0), ref / np.abs(ref), 1.0)
        gain = gain * np.conj(out)[None, :]
        resid = np.full(gain.shape, np.nan)
        for t in range(gain.shape[0]):
            phase = unwrapped(np.angle(gain[t]))
            good = np.isfinite(phase)
            if good.sum() >= 4 * (order + 1):
                resid[t] = phase - polyval_fit(x, phase, order, good)
        common = np.nanmedian(resid, axis=0)
        return out * np.exp(1j * np.where(np.isfinite(common), common, 0.0))


def block_mean(a, n):
    n = max(1, n)
    nchan = a.shape[-1]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        return np.nanmean(a[..., :(nchan // n) * n].reshape(a.shape[:-1] + (nchan // n, n)), axis=-1)


def enhance(spec, cplx, freqs, args):
    """
    spec[pol(xx, yy), tile, chan]: autocorrelation spectra; cplx[timeblock, tile, chan, 4]: solutions.
    Returns (new solutions, summary).
    """
    ntb, ntile, nchan, _ = cplx.shape
    x = 2 * (freqs - freqs.mean()) / (freqs.max() - freqs.min())
    chan_hz = float(np.median(np.diff(freqs)))
    # Gaussian in frequency with half power at the delay cut: exp(-2 pi^2 sigma^2 tau^2) = 1/2
    sigma_chan = np.sqrt(np.log(2) / 2) / (np.pi * args.delay_cut * 1e-6) / chan_hz
    kernel = gauss_kernel(sigma_chan)
    out = cplx.copy()
    summary = dict(delay_cut_us=args.delay_cut, sigma_chan=float(sigma_chan), pols={})
    with np.errstate(all="ignore"), warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        for b in range(ntb):
            # the phase all tiles share in a channel cancels in every visibility and is left
            # undetermined by a per-channel calibration: it is taken out of the whole Jones matrix
            # (from XX). What YY then still shares (YY against XX) is taken out for the fit and
            # put back
            cx = common_phase(cplx[b, :, :, 0], x)
            out[b] = out[b] * np.conj(cx)[None, :, None]
            cy = common_phase(out[b, :, :, 3], x)
            for name, col, ipol in POLS:
                gain = out[b, :, :, col] * (np.conj(cy)[None, :] if col == 3 else 1.0)
                okc = np.isfinite(gain) & np.isfinite(spec[ipol])  # [tile, chan]
                okt = okc.sum(axis=1) >= max(16, nchan // 4)
                u = np.full((ntile, nchan), np.nan)
                v = np.full((ntile, nchan), np.nan)
                low = np.full((ntile, nchan), np.nan, dtype=complex)
                fine_cal = np.full((ntile, nchan), np.nan, dtype=complex)
                fine_auto = np.full((ntile, nchan), np.nan, dtype=complex)
                for t in np.flatnonzero(okt):
                    good = okc[t]
                    ut = interp_bad(np.where(good, -0.5 * np.log(spec[ipol, t]), 0.0), good)
                    ut_low, _ = lowpass(x, ut, np.ones(nchan, dtype=bool), args.amp_order, kernel)
                    vt = causal_phase(ut - ut_low)
                    u[t], v[t] = ut, vt
                    lg = np.log(np.abs(gain[t]))
                    ph = unwrapped(np.angle(gain[t]))
                    la, _ = lowpass(x, lg - ut, good, args.amp_order, kernel, args.clip)
                    lp, _ = lowpass(x, ph - vt, good, args.phase_order, kernel, args.clip)
                    low[t] = la + 1j * lp
                    # the fine structure of each, for the checks below
                    fine_cal[t] = (lg - ut - la) + 1j * (ph - vt - lp)  # solution minus (autos + low pass)
                    fine_auto[t] = (ut - ut_low) + 1j * vt
                    fine_cal[t][~good] = np.nan
                # what all tiles share in the difference: the fine structure common to the solutions
                # that the autos do not have (or the reverse)
                common = np.nanmedian(fine_cal[okt].real, axis=0) + 1j * np.nanmedian(fine_cal[okt].imag, axis=0)
                common = np.where(np.isfinite(common), common, 0.0)
                # tiles whose autos do not describe their gain
                mis = np.sqrt(np.nanmean(block_mean((fine_cal - common[None]).real, args.mismatch_chans) ** 2, axis=1))
                typical = np.nanmedian(mis[okt])
                keep = np.zeros(ntile, dtype=bool)
                if args.max_mismatch > 0:
                    keep = okt & (mis > args.max_mismatch * typical)
                # the sign of the causal phase: the phase ripple of the solutions has to follow it
                num = np.nansum(fine_auto[okt & ~keep].imag * (fine_cal[okt & ~keep].imag + fine_auto[okt & ~keep].imag))
                den = np.nansum(fine_auto[okt & ~keep].imag ** 2)
                phase_slope = float(num / den) if den > 0 else float("nan")
                new = u + 1j * v + low
                if args.common_from_cal:
                    new = new + common.real[None]
                if args.no_phase:
                    new = new.real + 1j * np.array([unwrapped(np.angle(gain[t])) for t in range(ntile)])
                new = np.exp(new)
                if col == 3:
                    new = new * cy[None, :]
                fit = okt & ~keep
                # flagged values stay flagged
                out[b, fit, :, col] = np.where(np.isfinite(out[b, fit, :, col]), new[fit], np.nan)
                summary["pols"].setdefault(name, []).append(dict(
                    fitted=int(fit.sum()), kept=int(keep.sum()), kept_tiles=[int(t) for t in np.flatnonzero(keep)],
                    mismatch_typical=float(typical),
                    mismatch_kept=[float(m) for m in mis[keep]],
                    phase_slope=phase_slope,
                    # rms over channels, median over tiles
                    auto_fine_amp=float(np.nanmedian(np.sqrt(np.nanmean(fine_auto[fit].real ** 2, axis=1)))),
                    cal_minus_auto_amp=float(np.nanmedian(np.sqrt(np.nanmean(fine_cal[fit].real ** 2, axis=1)))),
                    cal_minus_auto_phase_deg=float(np.degrees(np.nanmedian(np.sqrt(np.nanmean(fine_cal[fit].imag ** 2, axis=1))))),
                    common_amp=float(np.nanstd(common.real)),
                ))
    return out, summary


def main():
    if len(sys.argv) > 1:
        args = parse_args()
    else:
        # is being called directly from nextflow
        args = parse_args(shlex.split("""${args}"""))
    print(vars(args))
    from astropy.io import fits
    autos, weights, freqs, names = read_autos(args.uvfits)
    with fits.open(args.soln) as hdul:
        jones = np.array(hdul["SOLUTIONS"].data, dtype=float)
        sol_freqs = np.array(hdul["CHANBLOCKS"].data["Freq"], dtype=float)
        sol_names = [str(n).strip() for n in hdul["TILES"].data["TileName"]]
    if jones.ndim != 4 or jones.shape[-1] != 8:
        raise UserWarning(f"{args.soln}: unexpected SOLUTIONS shape {jones.shape}")
    if len(sol_freqs) != len(freqs) or not np.allclose(sol_freqs, freqs, atol=1.0) or sol_names != names:
        raise UserWarning(f"{args.soln}: channels or tiles differ from {args.uvfits}")
    cplx = jones[..., 0::2] + 1j * jones[..., 1::2]
    out, summary = enhance(auto_spectrum(autos, weights), cplx, freqs, args)
    for name in summary["pols"]:
        for s in summary["pols"][name]:
            s["kept_tiles"] = [names[t] for t in s["kept_tiles"]]
    new = np.empty(out.shape[:-1] + (8,))
    new[..., 0::2] = out.real
    new[..., 1::2] = out.imag
    fitted = sum(s["fitted"] for p in summary["pols"].values() for s in p)
    kept = sum(s["kept"] for p in summary["pols"].values() for s in p)
    with fits.open(args.soln) as hdul:
        hdul["SOLUTIONS"].data = new.astype(hdul["SOLUTIONS"].data.dtype)
        hdul[0].header["AUTOFIT"] = (str(args.uvfits).split("/")[-1][:60], "fine structure from the autos of")
        hdul[0].header["AUTOCUT"] = (args.delay_cut, "[us] delay above which the autos are used")
        hdul[0].header["AUTOFITN"] = (fitted, "tile pols fitted")
        hdul[0].header["AUTOFITK"] = (kept, "tile pols kept as solved")
        hdul.writeto(args.out, overwrite=True)
    summary.update(soln=args.soln, out=args.out, fitted=fitted, kept=kept)
    if args.json:
        with open(args.json, "w") as f:
            json.dump(summary, f, indent=1)
    print(json.dumps(summary))


if __name__ == '__main__':
    main()
