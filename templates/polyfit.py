#!/usr/bin/env python
# this is nextflow-ified version of /astro/mwaeor/ctrott/poly_fit.py
# TODO: this deserves its own version control
from astropy.io import fits
import numpy as np
from numpy.polynomial import Polynomial
def get_unflagged_indices(n_bands, n_chan, clip_width=0):
    # Returns indices of unflagged frequency channels
    edge_flags = 2 + clip_width
    centre_flag = n_chan // 2
    channels = np.arange(n_chan)
    channels = np.delete(channels, centre_flag)
    channels = channels[edge_flags:-edge_flags]
    all_channels = []
    for n in range(n_bands):
        for c in channels:
            all_channels.append(c+(n*n_chan))
    return all_channels
def fit_hyperdrive_sols(sols, results, order=3, clip_width=0):
    # Fits polynomial to real and imaginary part of hyperdrive solutions
    # Remove flagged channels
    # Hyperdrive solutions have dimensions ((time),ant,freq,pols)
    n_bands = 24
    n_ants, n_freqs, n_pols = np.shape(sols)
    assert n_freqs % n_bands == 0
    n_chan = n_freqs // n_bands
    models_out = np.zeros_like(sols, dtype=complex)
    clipped_x = get_unflagged_indices(n_bands, n_chan, clip_width=clip_width)
    # filter any channels where result is NaN
    clipped_x = [ x for x in clipped_x if not np.isnan(results[x]) ]
    # Remove flagged tiles which are nan at first unflagged frequency and pol
    good_tiles = np.argwhere(~np.isnan(sols[:, clipped_x[0], 0]))
    freq_array = np.arange(n_freqs)
    for ant in good_tiles:
        for pol in range(n_pols):
            z_r = Polynomial.fit(clipped_x, np.real(
                sols[ant, clipped_x, pol]), deg=order)
            z_i = Polynomial.fit(clipped_x, np.imag(
                sols[ant, clipped_x, pol]), deg=order)
            models_out[ant, :, pol] = z_r(freq_array) + z_i(freq_array) * 1j
    return models_out
infile = "hyp_soln_${obsid}_${dical_name}.fits"
outfile = "hyp_soln_${obsid}_${name}.fits"
with fits.open(infile) as hdus:
    # Reads hyperdrive solutions from fits file into numpy array
    data = hdus["SOLUTIONS"].data
    for i_timeblock in range(data.shape[0]):
        sols = data[i_timeblock, :, :, ::2] + data[i_timeblock, :, :, 1::2] * 1j
        results = hdus["RESULTS"].data[i_timeblock,:]
        fit_sols = fit_hyperdrive_sols(sols, results)
        # Write fitted solutions
        hdus["SOLUTIONS"].data[i_timeblock, :, :, ::2] = np.real(fit_sols)
        hdus["SOLUTIONS"].data[i_timeblock, :, :, 1::2] = np.imag(fit_sols)
    hdus.writeto(outfile, overwrite="True")