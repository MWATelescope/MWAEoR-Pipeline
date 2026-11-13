#!/usr/bin/env python

"""
Night-level bayesian-rfi fit of one coarse channel, from SSINS incoherent noise spectra.

The model is the mixture model of bayesian-rfi (https://github.com/d3v-null/bayesian-rfi, branch
feat/mwax-night-stack), without witnesses: each time-differenced sample of the channel-mean
cross-correlation INS amplitude is either clean (a smooth Legendre background plus thermal noise)
or RFI (a skewed Cauchy excess). The settings are those of the night fits behind the VS masks:
one block per run of the same pointing, non-centred Legendre coefficients, noise rescaled to
0.335, NUTS with max_tree_depth 8 and target_accept 0.6, init_to_median, seed 42.

Inputs are the per-observation files of `hyperdrive bmetrics --ins` (or `SSINS.INS.write`) for
every observation of one night. The output holds P(RFI) for each differenced sample of each
observation, for XX and YY, and is turned into a flag mask by brfi_mask.py.

Needs numpyro, jax and a bayesian-rfi checkout (`--brfi-src`, or `mixture_model` on the python
path).

example:
    python templates/brfi_fit.py --brfi-src bayesian-rfi --chan-idx 19 --out brfi_2026-09-05_ch19.npz --obs 1472670968:0:birli_1472670968_INS.h5 1472671088:0:birli_1472671088_INS.h5
"""

import json
import os
import shlex
import sys
import time

import h5py
import numpy as np

COARSE_HZ = 1.28e6
# pyuvdata polarisation numbers
POLS = {"xx": -5, "yy": -6}


def parse_args(argv=None):
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--obs', nargs="+", required=True,
                        help='one <label>:<pointing>:<ins h5> per observation of the night')
    parser.add_argument('--chan-idx', type=int, required=True,
                        help='which coarse channel of the band to fit, counted from 0 in frequency order')
    parser.add_argument('--out', required=True, help='output npz')
    parser.add_argument('--terms', type=int, default=12, help='Legendre terms of the background per block')
    parser.add_argument('--warmup', type=int, default=1000)
    parser.add_argument('--samples', type=int, default=1000)
    parser.add_argument('--chains', type=int, default=4)
    parser.add_argument('--rescale-noise', type=float, default=0.335)
    parser.add_argument('--seed', type=int, default=42)
    parser.add_argument('--brfi-src', help='a bayesian-rfi checkout, if mixture_model is not already importable')
    return parser.parse_args(argv)


def coarse_chans(freqs_hz):
    """MWA receiver coarse channel number of each fine channel."""
    return np.floor(np.asarray(freqs_hz, float) / COARSE_HZ + 0.5 + 1e-6).astype(int)


def read_ins_series(path, chan_idx):
    """
    Channel-mean INS amplitude of one coarse channel.

    Returns the centre times (JD) of the differenced samples, {pol: series}, the receiver channel
    number, the number of baselines and the number of fine channels averaged.
    """
    with h5py.File(path, "r") as h5:
        metric = np.asarray(h5["Data/metric_array"], float)
        weights = np.asarray(h5["Data/weights_array"], float)
        times = np.asarray(h5["Header/time_array"], float)
        freqs = np.asarray(h5["Header/freq_array"], float).reshape(-1)
        pols = [int(p) for p in np.asarray(h5["Header/polarization_array"]).reshape(-1)]
    cc = coarse_chans(freqs)
    uniq = np.unique(cc)
    if chan_idx >= len(uniq):
        raise UserWarning(f"{path} has {len(uniq)} coarse channels, no index {chan_idx}")
    sel = cc == uniq[chan_idx]
    series, nchans, nbls = {}, [], []
    for name, num in POLS.items():
        ip = pols.index(num)
        m = metric[:, sel, ip]
        w = weights[:, sel, ip]
        m = np.where(np.isfinite(m) & (w > 0), m, np.nan)
        counts = np.isfinite(m).sum(axis=1)
        with np.errstate(invalid="ignore"):
            series[name] = np.where(counts > 0, np.nansum(m, axis=1) / np.maximum(counts, 1), np.nan)
        if (counts > 0).any():
            nchans.append(np.median(counts[counts > 0]))
            nbls.append(np.median(w[w > 0]))
    if not nchans:
        raise UserWarning(f"{path} has no unflagged data in coarse channel {uniq[chan_idx]}")
    return times, series, int(uniq[chan_idx]), int(np.median(nbls)), int(np.median(nchans))


def night_series(obs, chan_idx):
    """Stack the observations in time order; a block is a run of observations at one pointing."""
    rows = []
    for spec in obs:
        label, pointing, path = spec.split(":", 2)
        times, series, chan, nbls, nchan = read_ins_series(path, chan_idx)
        order = np.argsort(times)
        rows.append(dict(label=label, pointing=pointing, times=times[order], chan=chan, nbls=nbls, nchan=nchan,
                         series={k: v[order] for k, v in series.items()}))
    rows.sort(key=lambda r: r["times"][0])
    chans = sorted(set(r["chan"] for r in rows))
    if len(chans) != 1:
        raise UserWarning(f"observations disagree on the coarse channel at index {chan_idx}: {chans}")
    blocks, start, n = [], 0, 0
    for i, r in enumerate(rows):
        if i > 0 and r["pointing"] != rows[i - 1]["pointing"]:
            blocks.append((start, n))
            start = n
        n += len(r["times"])
    blocks.append((start, n))
    return dict(
        label=np.concatenate([np.full(len(r["times"]), r["label"]) for r in rows]),
        time_jd=np.concatenate([r["times"] for r in rows]),
        series={pol: np.concatenate([r["series"][pol] for r in rows]) for pol in POLS},
        blocks=blocks, chan=chans[0],
        nbls=int(np.median([r["nbls"] for r in rows])), nchan=int(np.median([r["nchan"] for r in rows])),
    )


def fit_series(t, y_all, blocks, nbls, nchan, args):
    """P(RFI) of each sample of one series. Samples that are not finite get NaN."""
    import jax
    import jax.numpy as jnp
    import numpyro
    try:
        import mixture_model as mm
    except ImportError as exc:
        raise UserWarning(
            "bayesian-rfi's mixture_model is not importable: give a checkout of "
            "https://github.com/d3v-null/bayesian-rfi with --brfi-src (params.brfi_src)") from exc
    import numpyro.distributions as dist
    from numpy.polynomial.legendre import legvander
    from numpyro.diagnostics import split_gelman_rubin
    from numpyro.handlers import reparam, seed, substitute
    from numpyro.infer import MCMC, NUTS, init_to_median
    from numpyro.infer.reparam import LocScaleReparam

    L = args.terms
    # a block with no finite sample cannot be fitted: leave its samples at NaN
    blocks_fit = [(a, b) for a, b in blocks if np.isfinite(y_all[a:b]).any()]
    samples = [y_all[a:b] for a, b in blocks_fit]
    sigma_clean = mm.calculate_sigma_clean(samples, None, None, nbls=[nbls] * len(blocks_fit),
                                           nchans=[nchan] * len(blocks_fit))
    sigma_raw = np.asarray(sigma_clean)
    samples, sigma_clean, scale = mm.rescale_to_noise(samples, sigma_clean, args.rescale_noise)
    Nb, T = len(blocks_fit), max(b - a for a, b in blocks_fit)
    y_pad = np.zeros((Nb, T))
    A_pad = np.zeros((Nb, T, L))
    mask = np.zeros((Nb, T), bool)
    y_loc = np.zeros(Nb)
    y_scale = np.zeros(Nb)
    for i, ((a, b), y) in enumerate(zip(blocks_fit, samples)):
        n = b - a
        ok = np.isfinite(y)
        span = np.ptp(t[a:b])
        x = 2 * (t[a:b] - t[a:b].min()) / span - 1 if span > 0 else np.zeros(n)
        y_pad[i, :n] = np.where(ok, y, 0)
        A_pad[i, :n] = legvander(x, L - 1)
        mask[i, :n] = ok
        y_loc[i], y_scale[i] = np.median(y[ok]), np.std(y[ok])

    def model(A_pad, y_pad, mask, y_loc, y_scale, sigma_clean, rfi_scale_cutoff):
        Nn, _, Lm = A_pad.shape
        degree = Lm - 1
        w_rfi = numpyro.sample("w_rfi", dist.Beta(2, 10))
        w_clean = 1.0 - w_rfi
        loc_rfi = numpyro.sample("loc_rfi", dist.HalfNormal(0.25))
        scale_rfi = numpyro.sample("scale_rfi", dist.TransformedDistribution(
            dist.LogNormal(0, 0.5), dist.transforms.AffineTransform(loc=rfi_scale_cutoff, scale=1.0)))
        mu_beta = numpyro.sample("mu_beta", dist.Normal(jnp.zeros(degree), jnp.ones(degree)).to_event(1))
        sigma_beta = numpyro.sample("sigma_beta", dist.HalfNormal(jnp.ones(degree)).to_event(1))
        with numpyro.plate("nights", Nn, dim=-1):
            beta = numpyro.sample("beta", dist.Normal(
                jnp.concatenate([y_loc[:, None], jnp.broadcast_to(mu_beta, (Nn, degree))], -1),
                jnp.concatenate([5.0 * y_scale[:, None], jnp.broadcast_to(sigma_beta, (Nn, degree))], -1)
            ).to_event(1))
        bg = jnp.einsum("nl,ntl->nt", beta, A_pad)
        zz = y_pad - bg
        lc = jnp.log(w_clean) + dist.Normal(0.0, sigma_clean[:, None]).log_prob(zz)
        lr = jnp.log(w_rfi) + mm.SkewCauchy(loc=loc_rfi, scale=scale_rfi, a=0.99).log_prob(zz)
        numpyro.factor("log_likelihood", jnp.where(mask, jnp.logaddexp(lc, lr), 0.0))
        return lc, lr, bg

    t0 = time.perf_counter()
    margs = dict(
        A_pad=jnp.asarray(A_pad), y_pad=jnp.asarray(y_pad), mask=jnp.asarray(mask), y_loc=jnp.asarray(y_loc),
        y_scale=jnp.asarray(y_scale), sigma_clean=jnp.asarray(sigma_clean),
        rfi_scale_cutoff=2 * float(np.median(sigma_clean)))
    kern = NUTS(reparam(model, config={"beta": LocScaleReparam(centered=0)}),
                init_strategy=init_to_median(num_samples=1000), max_tree_depth=8, target_accept_prob=0.6)
    mcmc = MCMC(kern, num_warmup=args.warmup, num_samples=args.samples, num_chains=args.chains, progress_bar=False)
    mcmc.run(jax.random.key(args.seed), **margs, extra_fields=("diverging",))
    post = mcmc.get_samples(group_by_chain=True)
    ndiv = int(np.asarray(mcmc.get_extra_fields()["diverging"]).sum())
    flat = {k: v.reshape((-1,) + v.shape[2:]) for k, v in post.items()}
    names = ("w_rfi", "loc_rfi", "scale_rfi", "mu_beta", "sigma_beta", "beta")
    sel = np.linspace(0, flat["w_rfi"].shape[0] - 1, min(400, flat["w_rfi"].shape[0])).astype(int)
    P, BG = [], []
    for s in sel:
        sub = {k: v[s] for k, v in flat.items() if k in names}
        lc, lr, bg = substitute(seed(model, 0), sub)(**margs)
        P.append(np.asarray(jnp.exp(lr - jnp.logaddexp(lc, lr))))
        BG.append(np.asarray(bg))
    P, BG = np.array(P), np.array(BG).mean(0) / scale
    pm, ps = P.mean(0), P.std(0)
    out = {k: np.full(len(t), np.nan) for k in ("p_rfi", "p_rfi_std", "bg", "sigma")}
    for i, (a, b) in enumerate(blocks_fit):
        n = b - a
        out["p_rfi"][a:b] = np.where(mask[i, :n], pm[i, :n], np.nan)
        out["p_rfi_std"][a:b] = np.where(mask[i, :n], ps[i, :n], np.nan)
        out["bg"][a:b] = BG[i, :n]
        out["sigma"][a:b] = sigma_raw[i]
    rhat = {k: float(np.nanmax(split_gelman_rubin(np.asarray(v)))) for k, v in post.items() if k in names}
    info = dict(
        nchan=nchan, nbls=nbls, terms=L, warmup=args.warmup, samples=args.samples, chains=args.chains,
        wall_s=time.perf_counter() - t0, rhat=rhat, divergences=ndiv, data_scale=scale,
        post_mean={k: float(np.asarray(flat[k]).mean()) for k in ("w_rfi", "loc_rfi", "scale_rfi")},
        frac_rfi=float(np.nanmean(out["p_rfi"] >= 0.5)), mean_p=float(np.nanmean(out["p_rfi"])),
        blocks=len(blocks_fit), block_len=[int(b - a) for a, b in blocks_fit])
    return out, info


def main():
    if len(sys.argv) > 1:
        args = parse_args()
    else:
        # is being called directly from nextflow
        args = parse_args(shlex.split("""${args}"""))
    print(vars(args))
    os.environ.setdefault("JAX_PLATFORMS", "cpu")
    if args.brfi_src:
        sys.path.insert(0, args.brfi_src)
    import numpyro
    numpyro.set_host_device_count(args.chains)

    night = night_series(args.obs, args.chan_idx)
    t = night["time_jd"]
    result = dict(label=night["label"], time_jd=t, coarse_chan=night["chan"])
    info = dict(coarse_chan=night["chan"], chan_idx=args.chan_idx, nobs=len(args.obs), nsamples=len(t))
    for pol in POLS:
        out, pol_info = fit_series(t, night["series"][pol], night["blocks"], night["nbls"], night["nchan"], args)
        for k, v in out.items():
            result[f"{k}_{pol}"] = v
        result[f"ins_{pol}"] = night["series"][pol]
        info[pol] = pol_info
    result["info"] = json.dumps(info)
    np.savez(args.out, **result)
    print(json.dumps(info))


if __name__ == '__main__':
    main()
