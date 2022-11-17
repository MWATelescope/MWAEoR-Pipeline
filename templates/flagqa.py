#!/usr/bin/env python

from astropy.io import fits
from json import dump as json_dump
import numpy as np

def split_strip_filter(str):
    return list(filter(None, map(lambda tok: tok.strip(), str.split(','))))

total_occupancy = 0
total_non_preflagged_bl_occupancy = 0
num_cchans = 0
num_unflagged_cchans = 0
with \
    open('${metrics}', 'w') as out, \
    fits.open("${metafits}") as meta, \
    fits.open("${uvfits}") as uv \
:
    data = {'obsid': int('${obsid}'), 'channels': {}}
    td = meta['TILEDATA'].data
    preflagged_ants = np.unique(td[td['Flag'] > 0]['Antenna'])
    data['preflagged_ants'] = preflagged_ants.tolist()
    sky_chans = [*map(int, split_strip_filter(meta[0].header['CHANNELS']))]
    num_cchans = len(sky_chans)

    vis_hdu = uv['PRIMARY']
    baseline_array = np.int16(vis_hdu.data["BASELINE"])
    num_blts = len(baseline_array)
    num_bls = len(np.unique(baseline_array))
    data['num_baselines']=num_bls
    num_times = num_blts // num_bls
    data['num_times']=num_times
    vis_shape = vis_hdu.data.data.shape
    assert vis_shape[0] == num_blts, 'vis shape 0 != num_blts'

    num_chans = vis_hdu.header['NAXIS4']
    data['num_chans']=num_chans
    assert num_chans % num_cchans == 0, 'uvfits channels is not divisible by coarse chans'
    num_fchans = num_chans // num_cchans # number of fine chans per coarse
    assert num_fchans * num_cchans == vis_shape[3], 'vis shape 3 != NAXIS4'

    ant_2_array = baseline_array % 256 - 1
    ant_1_array = (baseline_array - ant_2_array) // 256 - 1
    antpair_array = np.stack((ant_1_array, ant_2_array), axis=1)
    preflagged_blts_mask = np.isin(antpair_array, preflagged_ants).any(axis=1)
    # unflagged_blt_idxs = np.where(np.logical_not(preflagged_blts_mask))[0]
    (unflagged_blt_idxs,) = np.where(np.logical_not(preflagged_blts_mask))
    num_unflagged_blts = len(unflagged_blt_idxs)
    num_unflagged_bls = num_unflagged_blts // num_times
    unflagged_fraction = num_unflagged_blts / num_blts

    # assumption: vis data axes 1,2 (ra/dec?) are not used
    # assumption: flags are the same for all pols, so only use pol=0
    # 4d weight array: [time, baseline, coarse channel, fine channel]
    weights = vis_hdu.data.data[unflagged_blt_idxs, 0, 0, :, 0, 2].reshape(
        (num_times, num_unflagged_bls, num_cchans, num_fchans))
    # flag count for unflagged bls by [time, coarse channel, fine channel]
    unflagged_bl_flag_count = np.sum(np.int64(weights <= 0), axis=1)
    # where unflagged_bl_flag_count is flagged for all baselines
    flagged_mask = unflagged_bl_flag_count == num_unflagged_bls
    data['flagged_timestep_idxs'] = np.where(flagged_mask.all(axis=(1,2)))[0].tolist()
    flagged_cchan_idxs = np.where(flagged_mask.all(axis=(0,2)))[0].tolist()
    data['flagged_cchan_idxs'] = flagged_cchan_idxs
    flagged_sky_chans = np.array(sky_chans)[flagged_cchan_idxs].tolist()
    data['flagged_sky_chans'] = flagged_sky_chans
    data['flagged_fchan_idxs'] = np.where(flagged_mask.all(axis=(0,1)))[0].tolist()

    # occupancy of non-preflagged baselines by coarse channel
    unflagged_bl_occupancy = np.sum(unflagged_bl_flag_count, axis=(0,2)) / (num_unflagged_blts * num_fchans)

    for chan_idx, cchan_unflagged_bl_occupancy in enumerate(unflagged_bl_occupancy):
        if chan_idx not in flagged_cchan_idxs:
            num_unflagged_cchans += 1
            total_non_preflagged_bl_occupancy += cchan_unflagged_bl_occupancy
        sky_chan = sky_chans[chan_idx]
        cchan_total_occupancy = (unflagged_fraction * cchan_unflagged_bl_occupancy) + (1 - unflagged_fraction)
        total_occupancy += cchan_total_occupancy
        data['channels'][sky_chan] = {
            'occupancy': cchan_total_occupancy,
            'non_preflagged_bl_occupancy': cchan_unflagged_bl_occupancy,
            'flagged': chan_idx in data['flagged_cchan_idxs']
        }
    if num_cchans:
        total_occupancy /= num_cchans
    data['total_occupancy'] = total_occupancy
    if num_unflagged_cchans:
        total_non_preflagged_bl_occupancy /= num_unflagged_cchans
    data['total_non_preflagged_bl_occupancy'] = total_non_preflagged_bl_occupancy
    json_dump(data, out, indent=4)