#!/usr/bin/env python

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
from argparse import ArgumentParser
from collections import OrderedDict
import shlex
from codecs import decode

"""
example:

```bash
salloc --nodes=1 --mem=350G --time=8:00:00 --clusters=garrawarla --partition=workq --account=mwaeor --ntasks=20 --tmp=500G
module load singularity
cd /nvmetmp
export obsid=1365995952
cp /astro/mwaeor/dev/nfdata/\${obsid}/cal/hyp_\${obsid}_30l_src4k_8s_80kHz.uvfits vis.uvfits
singularity exec -B /nvmetmp --cleanenv --home /astro/mwaeor/dev/mplhome /pawsey/mwa/singularity/mwa_qa/mwa_qa_latest.sif python \
    /pawsey/mwa/mwaeor/dev/MWAEoR-Pipeline/templates/uvsplit.py \
    --uvfits=vis.uvfits \
    --split="vis_split.uvfits" \
    --times 1
cp \${obsid}_uvmeta.json /pawsey/mwa/mwaeor/dev/MWAEoR-Pipeline/
```
"""

def make_fits_axis_array(hdu, axis):
    count = hdu.header[f"NAXIS{axis}"]
    crval = hdu.header[f"CRVAL{axis}"]
    cdelt = hdu.header[f"CDELT{axis}"]
    crpix = hdu.header[f"CRPIX{axis}"]
    return cdelt * (np.arange(count) + (1 - crpix)) + crval

def update_fits_axis(hdu, axis_data, axis):
    for key in ['NAXIS', 'CTYPE', 'CRVAL', 'CRPIX', 'CDELT', 'CROTA', 'CUNIT']:
        if key in axis_data:
            hdu.header[f"{key}{axis}"] = axis_data[key]

def sanitize(name):
    if type(name) is bytes:
        name = name.decode('ascii', errors='ignore')
    return name.split(chr(0))[0]

def get_parser():
    parser = ArgumentParser(
        description="get basic metadata from uvfits file in json format")

    parser.add_argument('--uvfits', help='source uvfits file')
    parser.add_argument('--split', help='Name of split uvfits file',
                        default=None)
    parser.add_argument('--times', help='Timestep indices to split on',
                        type=int, nargs='+')
    parser.add_argument('--chans', help='channel indices to select (should be evenly spaced)',
                        type=int, nargs='+')


    return parser

def main():

    parser = get_parser()

    if len(sys.argv) > 1:
        args = parser.parse_args()
    else:
        # is being called directly from nextflow
        args = parser.parse_args([
            "--uvfits=${uvfits}",
            "--split=${split}",
        ] + shlex.split("${args}"))
    print(args)

    # uv = UVfits(args.uvfits)

    data = OrderedDict()
    # freqs = uv.freq_array.tolist()
    freqs = None
    jds = None
    time_array = None
    times = None
    mwa_loc = EarthLocation.of_site('mwa')

    with fits.open(args.uvfits) as hdus:
        vis_hdu = hdus['PRIMARY']
        print(f"{vis_hdu.header=}")
        # data['object'] = vis_hdu.header.get('OBJECT')
        # data['ra'] = vis_hdu.header.get('OBSRA')
        # data['dec'] = vis_hdu.header.get('OBSDEC')

        loc = None
        if inst := vis_hdu.header.get('INSTRUME', vis_hdu.header.get('TELESCOP')):
            data['instrument'] = inst
            try:
                loc = EarthLocation.of_site(inst.lower())
            except UnknownSiteException as exc:
                print(exc)
                print(EarthLocation.get_site_names())

        axes = []
        freq_axis = None
        freq_data = None
        naxis = vis_hdu.header['NAXIS']
        data_selection = [slice(None)] * naxis
        for axis in range(naxis):
            axis_data = OrderedDict()
            for key in ['NAXIS', 'CTYPE', 'CRVAL', 'CRPIX', 'CDELT', 'CROTA', 'CUNIT']:
                if f"{key}{axis+1}" in vis_hdu.header:
                    axis_data[key] = vis_hdu.header[f"{key}{axis+1}"]
            axis_array = None
            if all([key in axis_data for key in ['NAXIS', 'CRVAL', 'CDELT', 'CRPIX']]) and axis_data['NAXIS'] > 1:
                axis_array = make_fits_axis_array(vis_hdu, axis+1).tolist()
                axis_data['_min'] = axis_array[0]
                axis_data['_max'] = axis_array[-1]
            axes.append(axis_data)
            if 'CTYPE' in axis_data and 'FREQ' in axis_data.get('CTYPE'):
                if freq_axis is not None:
                    raise UserWarning("Multiple FREQ axes found. {freq_axis} and {axis}")
                freq_axis = axis
                freq_data = axis_data
                freqs = axis_array
            print(f"{axis=}, {axis_data=}")

        if freq_data is None:
            raise UserWarning("No FREQ axis found.")

        for grp in range(vis_hdu.header['PCOUNT']):
            grp_data = OrderedDict()
            for key in ['PTYPE', 'PSCAL', 'PZERO']:
                if f"{key}{grp+1}" in vis_hdu.header:
                    grp_data[key] = vis_hdu.header[f"{key}{grp+1}"]
            print(f"{grp=}, {grp_data=}")

        date_cols = [col.name for col in vis_hdu.data.columns if 'DATE' in col.name]
        time_array = np.float64(vis_hdu.data[date_cols.pop()])
        if date_cols:
            time_array += np.float64(vis_hdu.data[date_cols.pop()])
        jds = np.sort(np.unique(time_array))

        if ant_hdu := hdus['AIPS AN']:
            for key, value in ant_hdu.header.items():
                print(f"AN {key:8} => {value}")
            if dut1 := ant_hdu.header.get('UT1UTC'):
                data['dut1'] = dut1
            if not loc:
                first_xyz = ant_hdu.data['STABXYZ'][0]
                loc = EarthLocation.from_geocentric(
                    ant_hdu.header['ARRAYX'] + first_xyz[0],
                    ant_hdu.header['ARRAYY'] + first_xyz[1],
                    ant_hdu.header['ARRAYZ'] + first_xyz[2],
                    unit='m'
                )
            if 'ANNAME' in ant_hdu.data.columns.names:
                antnames = ant_hdu.data['ANNAME'].tolist()
                ant_hdu.data['ANNAME'] = [ sanitize(name) for name in antnames ]

        times = Time(jds, format='jd', scale='utc', location=loc, precision=3)
        if dut1 := data['dut1']:
            print(f"updating delta_ut1_utc to {dut1}")
            times.delta_ut1_utc = dut1

        # antnames = ant_hdu.data['ANNAME'].tolist()
        # antnums = ant_hdu.data['NOSTA'].tolist()

        # blt_idxs = np.arange(len(uv.time_array))
        if args.times:
            # jds = uv.unique_times
            # if args.times:
            jds = jds[sorted(args.times)]
            times = Time(jds, format='jd', scale='utc',
                location=mwa_loc, precision=3)
            print(times)

            # vis_hdu.header['']
            data_selection[0] = np.where(np.in1d(time_array, jds))[0]


        if args.chans:
            chan_idxs = np.sort(args.chans)
            freq_data['NAXIS'] = len(chan_idxs)
            freqs = np.array(freqs)[chan_idxs]
            freq_data['CRPIX'] = freq_data['NAXIS']//2 + 1
            freq_data['CRVAL'] = freqs[freq_data['CRPIX'] - 1]
            if len(freqs) > 1:
                freq_data['CDELT'] = freqs[1] - freqs[0]
                for il, ir in zip(chan_idxs[:-1], chan_idxs[1:]):
                    fl, fr = freqs[ir], freqs[il]
                    if fr - fl != freq_data['CDELT']:
                        raise UserWarning("Channels are not evenly spaced. "
                                          f"Channel {il} ({fl}) and {ir} ({fr}) "
                                          f"have delta {fr - fl} != {freq_data['CDELT']}")
            freq_data['_min'] = freqs[0]
            freq_data['_max'] = freqs[-1]
            print(f"{freq_data=}")

            update_fits_axis(vis_hdu, freq_data, freq_axis+1)

            data_selection[naxis - freq_axis] = chan_idxs

        if args.split:
            print(f"{data_selection=}")
            print(f"{vis_hdu.data.shape=}")
            print(f"{vis_hdu.data[data_selection[0]].data.shape=}")
            print(f"{vis_hdu.data.data.shape=}")
            vis_hdu.data = vis_hdu.data[data_selection[0]]
            # pars = [vis_hdu.data.par(p)[data_selection[0]] for p in range(vis_hdu.data.npar)]

            # modify vis_hdu.data with selection
            # data = vis_hdu.data.field('DATA')
            # numpy.lib.recfunctions.drop_fields(vis_hdu.data, ['DATA'])

            # for p, idxs in enumerate(data_selection[1:]):
            #     if idxs == slice(None):
            #         continue
            #     print(f"{p=}, {vis_hdu.data.par(p).shape=}, {idxs=}, {vis_hdu.data.par(p)}")
            #     vis_hdu.data.par(p)[...] = vis_hdu.data.par(p)[idxs]
            # vis_hdu.data[:] = vis_hdu.data[:].data[data_selection]
            hdus.writeto(args.split, overwrite=True)

if __name__ == '__main__':
    main()


# axis=0, axis_data=OrderedDict([('NAXIS', 0)])
# axis=1, axis_data=OrderedDict([('NAXIS', 3), ('CTYPE', 'COMPLEX'), ('CRVAL', 1.0), ('CRPIX', 1.0), ('CDELT', 1.0), ('_min', 1.0), ('_max', 3.0)])
# axis=2, axis_data=OrderedDict([('NAXIS', 4), ('CTYPE', 'STOKES'), ('CRVAL', -5), ('CRPIX', 1.0), ('CDELT', -1), ('_min', -5.0), ('_max', -8.0)])
# axis=3, axis_data=OrderedDict([('NAXIS', 768), ('CTYPE', 'FREQ'), ('CRVAL', 215695000.0), ('CRPIX', 385), ('CDELT', 40000.0), ('_min', 200335000.0), ('_max', 231015000.0)])
# axis=4, axis_data=OrderedDict([('NAXIS', 1), ('CTYPE', 'RA'), ('CRVAL', 333.607), ('CRPIX', 1), ('CDELT', 1)])
# axis=5, axis_data=OrderedDict([('NAXIS', 1), ('CTYPE', 'DEC'), ('CRVAL', -17.0267), ('CRPIX', 1), ('CDELT', 1)])
# grp=0, grp_data=OrderedDict([('PTYPE', 'UU'), ('PSCAL', 1.0), ('PZERO', 0.0)])
# grp=1, grp_data=OrderedDict([('PTYPE', 'VV'), ('PSCAL', 1.0), ('PZERO', 0.0)])
# grp=2, grp_data=OrderedDict([('PTYPE', 'WW'), ('PSCAL', 1.0), ('PZERO', 0.0)])
# grp=3, grp_data=OrderedDict([('PTYPE', 'BASELINE'), ('PSCAL', 1.0), ('PZERO', 0.0)])
# grp=4, grp_data=OrderedDict([('PTYPE', 'DATE'), ('PSCAL', 1.0), ('PZERO', 2460054.5)])