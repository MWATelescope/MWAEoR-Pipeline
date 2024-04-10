#!/usr/bin/env python

from astropy.io import fits
import astropy.units as u
from astropy.time import Time
from astropy.coordinates import EarthLocation
from astropy.coordinates.errors import UnknownSiteException
from json import dump as json_dump
import numpy as np
from math import ceil
import os
import pandas as pd
from argparse import ArgumentParser
import sys
from collections import Counter

"""
```bash
salloc --nodes=1 --mem=350G --time=8:00:00 --clusters=garrawarla --partition=workq --account=mwaeor --ntasks=20 --tmp=500G
module load singularity
mkdir /dev/shm/deleteme
cd /dev/shm/deleteme
export obsid=1087941792
birli_1087941792_2s_40kHz_norfi.uvfits
cp /astro/mwaeor/dev/nfdata/${obsid}/prep/birli_${obsid}_2s_40kHz_norfi.uvfits ${obsid}.uvfits
cp /astro/mwaeor/dev/nfdata/${obsid}/raw/${obsid}.metafits ${obsid}.metafits
singularity exec --bind `pwd` --cleanenv --home /astro/mwaeor/dev/mplhome /pawsey/mwa/singularity/mwa_qa/mwa_qa_latest.sif python \
    /pawsey/mwa/mwaeor/dev/MWAEoR-Pipeline/templates/flagqa.py \
    --metafits=${obsid}.metafits \
    --uvfits=${obsid}.uvfits \
    --obsid=${obsid}
jq -r '.' occupancy_${obsid}.json
cp ${obsid}.json /pawsey/mwa/mwaeor/dev/MWAEoR-Pipeline/test/
```
"""

def split_strip_filter(str):
    return list(filter(None, map(lambda tok: tok.strip(), str.split(','))))

def get_parser():
    parser = ArgumentParser(
        description="analyse uvfits for flag occupancy.")

    parser.add_argument('--obsid', help='obsid', default=None)
    parser.add_argument('--metafits', help='metafits file')
    parser.add_argument('--uvfits', help='source uvfits files', nargs='+')
    # parser.add_argument('--json', help='Name of output json file', default=None)
    return parser

# TODO: put all this in def main():

parser = get_parser()
if len(sys.argv) > 1:
    args = parser.parse_args()
else:
    # is being called directly from nextflow
    args = parser.parse_args([
        "--obsid=${obsid}",
        "--uvfits=${uvfits}",
        "--metafits=${metafits}",
    ])

data = {'channels': {}}
if args.obsid is not None:
    data['obsid'] = int(args.obsid)

with fits.open(args.metafits) as meta:
    mpr = meta['PRIMARY']
    sky_chans = [*map(int, split_strip_filter(mpr.header['CHANNELS']))]
    num_cchans = len(sky_chans)
    freq_res = mpr.header['FINECHAN'] * u.kHz
    bandwidth = mpr.header['BANDWDTH'] * u.MHz

    td = meta['TILEDATA'].data
    preflagged_ants = np.unique(td[td['Flag'] > 0]['Antenna'])
    data['preflagged_ants'] = preflagged_ants.tolist()
    preflagged_tilenames = np.unique(td[td['Flag'] > 0]['TileName'])
    data['preflagged_tilenames'] = preflagged_tilenames.tolist()
    metafits_inputs = pd.DataFrame(dict([
        (col, td.field(col).tolist())
        for col in ['Antenna', 'Tile', 'Rx', 'Slot', 'Pol', 'TileName', 'Flag', 'Length', 'North', 'East', 'Height', 'Gains', 'Flavors']
    ])).sort_values(['Antenna', 'Pol'])
    data['INPUTS'] = metafits_inputs.to_dict(orient='records')

cchan_bandwidth = bandwidth / num_cchans

def uvdata(data, path):
    total_occupancy = 0
    total_occupancy = 0
    total_rfi_occupancy = 0
    with fits.open(path) as uv:
        vis_hdu = uv['PRIMARY']

        # location
        # loc = None
        # if inst := vis_hdu.header.get('INSTRUME', vis_hdu.header.get('TELESCOP')):
        #     data['instrument'] = inst
        #     try:
        #         loc = EarthLocation.of_site(inst.lower())
        #     except UnknownSiteException as exc:
        #         print(exc)
        #         print(EarthLocation.get_site_names())

        # baselines
        try:
            baseline_array = np.int16(vis_hdu.data["BASELINE"])
        except (ValueError, KeyError):
            # use the ANTENNA1, ANTENNA2 columns instead
            baseline_array = np.int16(vis_hdu.data["ANTENNA1"] * 256 + vis_hdu.data["ANTENNA2"])
        unique_baselines = np.unique(baseline_array)
        num_blts = len(baseline_array)
        num_bls = len(unique_baselines)
        data['num_baselines'] = num_bls

        num_times = num_blts / num_bls

        # channels

        num_chans = vis_hdu.header['NAXIS4']
        freq_res = vis_hdu.header['CDELT4'] * u.Hz
        data['num_chans'] = num_chans

        # times
        date_cols = [col.name for col in vis_hdu.data.columns if 'DATE' in col.name]
        time_array = np.sum(np.stack([
            np.float64(vis_hdu.data[date_col])
            for date_col in date_cols
        ], axis=1), axis=1)
        unique_times = np.sort(np.unique(time_array))
        print(unique_times)

        try:
            assert num_blts == num_bls * len(unique_times), f'{num_blts=} != {num_bls=} * {num_times=}'
        except AssertionError as exc:
            print(exc)
            return data

        timestep_idx_array = np.searchsorted(unique_times, time_array)
        num_times = len(unique_times)
        data['num_times'] = num_times
        print(f"{num_times,num_bls,num_blts=}")

        # number of fine chans per coarse
        num_cchans = num_chans // ceil((cchan_bandwidth / freq_res).decompose().value)
        if num_cchans == 0:
            num_fchans = num_chans
        else:
            num_fchans = ceil((cchan_bandwidth / freq_res).decompose().value)

        vis = vis_hdu.data.data[...]
        while len(vis.shape) > 6:
            vis = vis[:, 0, ...]

        assert vis.shape[0] == num_blts, f'{vis.shape[0]=} != {num_blts=}'
        assert vis.shape[3] == num_chans, f'{vis.shape[3]=} != {num_chans=}'

        # weights / flags
        # assumption: vis data axes 1,2 (ra/dec?) are not used
        # assumption: flags are the same for all pols, so only use pol=0
        # 4d flag mask: [time, baseline, coarse channel, fine channel]
        flagged_mask = (vis_hdu.data.data[:, 0, 0, :, 0, 2] <= 0).reshape(
            (num_times, num_bls, num_cchans, num_fchans))
        print(f"{flagged_mask.shape=}")

        # flagged times
        flagged_timestep_mask = flagged_mask.all(axis=(1, 2, 3))
        (flagged_timestep_idxs,) = np.where(flagged_timestep_mask)
        (unflagged_timestep_idxs,) = np.where(np.logical_not(flagged_timestep_mask))
        data['flagged_timestep_idxs'] = flagged_timestep_idxs.tolist()
        flagged_times = unique_times[flagged_timestep_idxs]
        num_unflagged_times = num_times - len(flagged_timestep_idxs)
        non_preflagged_mask = flagged_mask[unflagged_timestep_idxs, :, :, :]

        # flagged baselines
        flagged_bl_mask = non_preflagged_mask.all(axis=(0, 2, 3))
        (flagged_bl_idxs,) = np.where(flagged_bl_mask)
        (unflagged_bl_idxs,) = np.where(np.logical_not(flagged_bl_mask))
        non_preflagged_mask = non_preflagged_mask[:, unflagged_bl_idxs, :, :]
        num_unflagged_bls = len(unflagged_bl_idxs)

        # flagged fine channels within a coarse channel
        flagged_fchan_mask = non_preflagged_mask.all(axis=(0, 1, 2))
        (flagged_fchan_idxs, ) = np.where(flagged_fchan_mask)
        (unflagged_fchan_idxs, ) = np.where(np.logical_not(flagged_fchan_mask))
        data['flagged_fchan_idxs'] = flagged_fchan_idxs.tolist()
        print(f"{flagged_fchan_idxs.tolist()=}")
        num_unflagged_fchans = len(unflagged_fchan_idxs)
        non_preflagged_mask = non_preflagged_mask[:, :, :, unflagged_fchan_idxs.tolist()]

        # flagged coarse channels
        flagged_cchan_mask = non_preflagged_mask.all(axis=(0, 1, 3))
        (flagged_cchan_idxs, ) = np.where(flagged_cchan_mask)
        (unflagged_cchan_idxs, ) = np.where(np.logical_not(flagged_cchan_mask))
        data['flagged_cchan_idxs'] = flagged_cchan_idxs.tolist()
        flagged_sky_chans = np.array(sky_chans)[flagged_cchan_idxs].tolist()
        data['flagged_sky_chans'] = flagged_sky_chans
        num_unflagged_cchans = len(unflagged_cchan_idxs)

        cchan_flagged_occupancies = np.sum(np.int64(flagged_mask), axis=(0, 1, 3)) / (num_times * num_bls * num_fchans)
        cchan_non_preflagged_occupancies = np.sum(np.int64(non_preflagged_mask), axis=(0, 1, 3)) / (num_unflagged_times * num_unflagged_bls * num_unflagged_fchans)

        for sky_chan, cchan_flagged, cchan_occupancy, cchan_non_preflagged_occupancy in \
                zip(sky_chans, flagged_cchan_mask, cchan_flagged_occupancies, cchan_non_preflagged_occupancies):
            data['channels'][sky_chan] = {
                'occupancy': cchan_occupancy,
                'rfi_occupancy': cchan_non_preflagged_occupancy,
                'flagged': sky_chan
            }
            total_occupancy += cchan_occupancy
            if not cchan_flagged:
                total_rfi_occupancy += cchan_non_preflagged_occupancy

        if num_cchans:
            total_occupancy /= num_cchans
        data['total_occupancy'] = total_occupancy
        if num_unflagged_cchans:
            total_rfi_occupancy /= num_unflagged_cchans
        data['total_rfi_occupancy'] = total_rfi_occupancy
    return data

meta = data.copy()
for path in args.uvfits:
    data = uvdata(meta.copy(), path)
    basename = os.path.basename(path)
    basename, _ = os.path.splitext(basename)
    with open(f"occupancy_{basename}.json", 'w') as out:
        json_dump(data, out, indent=4)
