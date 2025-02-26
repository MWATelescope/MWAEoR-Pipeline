#!/usr/bin/env python

import astropy
from astropy.io import fits
import astropy.units as u
from astropy.time import Time
from astropy.coordinates import EarthLocation
from astropy.coordinates.errors import UnknownSiteException
from json import dump as json_dump
import numpy as np
from math import ceil
# from mwa_qa.read_uvfits import UVfits
import sys
from os.path import abspath
from argparse import ArgumentParser
from collections import OrderedDict, defaultdict
import codecs

"""
example:

```bash
salloc --nodes=1 --mem=350G --time=8:00:00 --clusters=garrawarla --partition=workq --account=mwaeor --ntasks=20 --tmp=500G
module load singularity
mkdir -p /dev/shm/deleteme
cd /dev/shm/deleteme
export obsid=1090012424
eval cp /astro/mwaeor/dev/nfdata/\${obsid}/cal/hyp_\${obsid}_sub_30l_src4k_8s_80kHz.uvfits vis.uvfits
eval singularity exec -B \$PWD --cleanenv --home /astro/mwaeor/dev/mplhome /pawsey/mwa/singularity/mwa_qa/mwa_qa_latest.sif python \
    /pawsey/mwa/mwaeor/dev/MWAEoR-Pipeline/templates/uvmeta.py \
    --uvfits=vis.uvfits \
    --uvmeta="\${obsid}_uvmeta.json" \
    --weights
cp \${obsid}_uvmeta.json /pawsey/mwa/mwaeor/dev/MWAEoR-Pipeline/
```
"""

def make_fits_axis_array(hdu, axis):
    count = hdu.header[f"NAXIS{axis}"]
    crval = hdu.header[f"CRVAL{axis}"]
    cdelt = hdu.header[f"CDELT{axis}"]
    crpix = hdu.header[f"CRPIX{axis}"]
    return cdelt * (np.arange(count) + (1 - crpix)) + crval


def get_parser():
    parser = ArgumentParser(
        description="get basic metadata from uvfits file in json format")

    parser.add_argument('--uvfits', help='source uvfits file')
    parser.add_argument('--uvmeta', help='output json file', default=None)
    parser.add_argument('--weights', help='calculate weights', default=True, action='store_true')

    return parser


def sanitize(name):
    if type(name) is bytes:
        name = name.decode('ascii', errors='ignore')
    return name.split(chr(0))[0]


def main():

    parser = get_parser()

    if len(sys.argv) > 1:
        args = parser.parse_args()
    else:
        # is being called directly from nextflow
        args = parser.parse_args([
            "--uvfits=${vis}",
            "--uvmeta=${uvmeta}",
        ])
        # + shlex.split("\${args}")

    # uv = UVfits(args.uvfits)

    data = OrderedDict()
    # freqs = uv.freq_array.tolist()
    # mwa_loc = EarthLocation.of_site('mwa')
    # times = Time(uv.unique_times, format='jd', scale='utc',
    #              location=mwa_loc, precision=3)

    print(f"about to open {abspath(args.uvfits)}")

    with fits.open(args.uvfits) as hdus:
        vis_hdu = hdus['PRIMARY']
        for key, value in vis_hdu.header.items():
            print(f"PR {key:8} => {value}")

        if obj := vis_hdu.header.get('OBJECT'):
            data['object'] = obj
        if ra := vis_hdu.header.get('OBSRA'):
            data['ra'] = ra
        if dec := vis_hdu.header.get('OBSDEC'):
            data['dec'] = dec
        loc = None
        if inst := vis_hdu.header.get('INSTRUME', vis_hdu.header.get('TELESCOP')):
            data['instrument'] = inst
            try:
                loc = EarthLocation.of_site(inst.lower())
            except UnknownSiteException as exc:
                print(exc)
                print(EarthLocation.get_site_names())

        naxis = vis_hdu.header['NAXIS']
        weight_selection = [slice(None)] * naxis
        freq_data = None

        for axis in range(naxis):
            axis_data = OrderedDict()
            for key in ['NAXIS', 'CTYPE', 'CRVAL', 'CRPIX', 'CDELT', 'CROTA', 'CUNIT']:
                if f"{key}{axis+1}" in vis_hdu.header:
                    axis_data[key] = vis_hdu.header[f"{key}{axis+1}"]
            axis_array = None
            if all([key in axis_data for key in ['NAXIS', 'CRVAL', 'CDELT', 'CRPIX']]):
                axis_array = make_fits_axis_array(vis_hdu, axis+1).tolist()
                if axis_data['NAXIS'] > 1:
                    axis_data['_min'] = axis_array[0]
                    axis_data['_max'] = axis_array[-1]
            # axes.append(axis_data)
            if 'FREQ' in axis_data.get('CTYPE', ''):
                if freq_data is not None:
                    raise UserWarning(f"Multiple FREQ axes found. {freq_data} and {axis_data}")
                freq_data = axis_data
                freqs = axis_array
            elif 'RA' in axis_data.get('CTYPE', ''):
                if axis_data['NAXIS'] > 1:
                    print(UserWarning(f"multiple RA values found. {axis_array}"))
                data['ra'] = axis_array[0]
                weight_selection[naxis-axis] = 0
            elif 'DEC' in axis_data.get('CTYPE', ''):
                if axis_data['NAXIS'] > 1:
                    print(UserWarning(f"multiple DEC values found. {axis_array}"))
                data['dec'] = axis_array[0]
                weight_selection[naxis-axis] = 0
            elif axis > 0:
                weight_selection[naxis-axis] = 0

            print(f"{axis=}, {axis_data=}")
        # assumption: weight idx 3
        weight_selection[naxis-1] = 2
        weight_selection[0] = slice(None)
        print(f"{weight_selection=}")

        if freq_data is None:
            raise UserWarning("No FREQ axis found.")

        for grp in range(vis_hdu.header['PCOUNT']):
            grp_data = OrderedDict()
            for key in ['PTYPE', 'PSCAL', 'PZERO']:
                if f"{key}{grp+1}" in vis_hdu.header:
                    grp_data[key] = vis_hdu.header[f"{key}{grp+1}"]
            print(f"{grp=}, {grp_data=}")

        date_cols = [col.name for col in vis_hdu.data.columns if 'DATE' in col.name]
        print("reading time array")

        time_array = np.float64(vis_hdu.data[date_cols.pop()])
        if date_cols:
            time_array += np.float64(vis_hdu.data[date_cols.pop()])
        jds = np.sort(np.unique(time_array))

        if args.weights:
            jd_weights = defaultdict(lambda: 0)
            weights = vis_hdu.data['DATA'][tuple(weight_selection)].sum(axis=1)
            for jd, weight in zip(time_array, weights):
                jd_weights[jd] += weight

        antnames = []
        antnums = []

        print("reading antenna header")
        if ant_hdu := hdus['AIPS AN']:
            for key, value in ant_hdu.header.items():
                print(f"AN {key:8} => {value}")
            if dut1 := ant_hdu.header.get('UT1UTC'):
                data['dut1'] = dut1
            if loc is None:
                first_xyz = ant_hdu.data['STABXYZ'][0]
                loc = EarthLocation.from_geocentric(
                    ant_hdu.header['ARRAYX'] + first_xyz[0],
                    ant_hdu.header['ARRAYY'] + first_xyz[1],
                    ant_hdu.header['ARRAYZ'] + first_xyz[2],
                    unit='m'
                )
            antnames = ant_hdu.data['ANNAME'].tolist()
            antnums = ant_hdu.data['NOSTA'].tolist()

        times = Time(jds, format='jd', scale='utc', location=loc, precision=3)
        if dut1 := data.get('dut1'):
            print(f"updating delta_ut1_utc {times.delta_ut1_utc} -> {dut1}")
            times.delta_ut1_utc = dut1

    if data.get('ra') and data.get('dec'):
        nearest_ra = round(data['ra'])
        nearest_dec = round(data['dec'])
        print(f"nearest ra, dec: {nearest_ra}, nearest_dec")

        if nearest_ra == 0 and -30 <= nearest_dec <= -27:
            data['eorfield'] = 0
        elif nearest_ra == 0 and 340 <= nearest_dec <= 341:
            data['eorfield'] = 1

    first_mhz = round(freqs[0] / 1e6)
    print(f"first_mhz: {first_mhz}")
    if 137 <= first_mhz < 139:
        data['eorband'] = 0
    elif 166 <= first_mhz < 168:
        data['eorband'] = 1

    # lsts_rad = times.sidereal_time('apparent').radian.tolist()
    # gps_times = times.gps.tolist()
    # iso_times = times.iso.tolist()
    # jd1s, jd2s = [*zip(*map(lambda t: (t.jd1, t.jd2), times.jd.tolist()))]
    times = [
        OrderedDict(
            gps=float(f"{time.to_value('gps'):.2f}"),
            iso=time.to_value('iso', 'date_hms'),
            lst_rad=float(time.sidereal_time('apparent').radian),
            jd1=time.jd1,
            jd2=time.jd2,
        )
        for time in times
    ]
    if args.weights:
        total_weight = 0
        for i, jd in enumerate(jds):
            weight = float(jd_weights[jd])
            times[i]['weight'] = weight
            total_weight += weight
        data['total_weight'] = total_weight

    ants = [
        OrderedDict(
            nosta=antnum,
            name=sanitize(antname),
        )
        for antnum, antname in zip(antnums, antnames)
    ]

    if 'mwa' in data['instrument'].lower():
        num_lb = 0
        num_hex = 0
        for ant in ants:
            if 'lb' in ant['name'].lower():
                num_lb += 1
            if 'hex' in ant['name'].lower():
                num_hex += 1
        num_ants = len(ants)
        data['num_ants'] = num_ants
        data['num_lb'] = num_lb
        data['num_hex'] = num_hex
        if num_ants > 128:
            data['config'] = f'phase3-{num_ants}T'
        elif data['num_lb'] > 50:
            data['config'] = f'phase2b-{num_ants}T'
        elif data['num_hex'] > 50:
            data['config'] = f'phase2a-{num_ants}T'
        else:
            data['config'] = f'phase1-{num_ants}T'

    data['times'] = times
    data['freqs'] = freqs
    data['ants'] = ants

    with (open(args.uvmeta, 'w') if args.uvmeta else sys.stdout) as out_file:
        json_dump(data, out_file, indent=4)


if __name__ == '__main__':
    main()
    from os import system
    system('date -Is')
