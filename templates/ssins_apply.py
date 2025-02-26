#!/usr/bin/env python

# TODO: FREQUENCY BROADCAST

import os
import json

# from SSINS import Catalog_Plot as cp, util as ssins_util
# from SSINS import SS, INS, MF, plot_lib
from matplotlib import pyplot as plt, cm
import numpy as np
import pylab
# from pyuvdata import UVData, UVFlag, utils as uvutils
import sys
import shlex
from astropy.time import Time
from astropy.io import fits
from collections import OrderedDict
import h5py

"""
singularity exec \
    --bind ${PWD} --cleanenv --home /astro/mwaeor/dev/mplhome \
    /pawsey/mwa/singularity/ssins/ssins_latest.sif python \
    /pawsey/mwa/mwaeor/dev/MWAEoR-Pipeline/templates/ssins_apply.py \
    --src /astro/mwaeor/dev/nfdata/1379177304/prep_qa/birli_1379177304_2s_40kHz__SSINS_mask.h5 \
    --dst /astro/mwaeor/dev/nfdata/1379177304/prep/birli_1379177304_2s_40kHz.uvfits
"""


def make_fits_axis_array(hdu, axis):
    count = hdu.header[f"NAXIS{axis}"]
    crval = hdu.header[f"CRVAL{axis}"]
    cdelt = hdu.header[f"CDELT{axis}"]
    crpix = hdu.header[f"CRPIX{axis}"]
    return cdelt * (np.arange(count) + (1 - crpix)) + crval


def read_h5(h5file):
    h5data = {}
    with h5py.File(h5file) as h5stream:
        h5data['freq_array'] = h5stream['Header']['freq_array'][:]
        # h5data['jd_array'] =
        jd_array = h5stream['Header']['time_array'][:]
        # weights are bool[time, channel, pol, ]
        weights = h5stream['Data']['flag_array'].astype(np.float64)
        # make weights float
        h5data['weights'] = 1.0 - weights[:, :, :]
    h5data['gps_array'] = jd_to_gps(jd_array)
    return h5data


def jd_to_gps(jd):
    return np.round(Time(jd, format='jd').gps).astype(np.int64)


def read_uvfits(uvfile):
    uvdata = {}
    with fits.open(uvfile) as hdus:
        vis_hdu = hdus['PRIMARY']
        naxis = vis_hdu.header['NAXIS']

        # build weight_selection:
        # an ndarray slice that will produce a 3d array of weights
        # when indexing
        weight_selection = [slice(None)] * naxis
        for axis in range(naxis):
            axis_data = OrderedDict()
            for key in ['NAXIS', 'CTYPE', 'CRVAL', 'CRPIX', 'CDELT', 'CROTA',
                        'CUNIT']:
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
                # freq_data = axis_data
                uvdata['freq_array'] = axis_array
            elif 'RA' in axis_data.get('CTYPE', ''):
                weight_selection[naxis-axis] = 0
            elif 'DEC' in axis_data.get('CTYPE', ''):
                weight_selection[naxis-axis] = 0
            elif 'COMPLEX' in axis_data.get('CTYPE', ''):
                weight_selection[naxis-axis] = 2
            # elif axis > 0:
            #     # get only get first pol
            #     weight_selection[naxis-axis] = 0
        uvdata['weight_selection'] = weight_selection

        date_cols = [col.name for col in vis_hdu.data.columns if 'DATE' in col.name]
        jd_array = np.float64(vis_hdu.data[date_cols.pop()])
        if date_cols:
            jd_array += np.float64(vis_hdu.data[date_cols.pop()])

        uvdata['weights'] = vis_hdu.data['DATA'][tuple(weight_selection)]
    uvdata['gps_array'] = jd_to_gps(jd_array)
    return uvdata


def set_weights(uvfile, uvdata, weights):
    new_uvfile = uvfile.replace('.uvfits', '.ssins.uvfits')
    print(f"saving to {new_uvfile}")
    with fits.open(uvfile) as hdus:
        vis_hdu = hdus['PRIMARY']
        vis_hdu.data['DATA'][tuple(uvdata['weight_selection'])] = weights
        hdus.writeto(new_uvfile, overwrite=True)


def parse_args(argv=sys.argv[1:]):
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--src', nargs="+", help='source uvfits files')
    parser.add_argument('--dst', help='destination uvfits files')
    return parser.parse_args(argv)


def concat_data_freq(datas):
    datas.sort(key=lambda x: x['freq_array'][0])
    result = {
        'freq_array': np.concatenate([d['freq_array'] for d in datas]),
        # 'jd_array': np.concatenate([d['jd_array'] for d in datas]),
        'gps_array': np.concatenate([d['gps_array'] for d in datas]),
        'weights': np.concatenate([d['weights'] for d in datas], axis=1),
    }
    assert np.all(result['freq_array'][:-1] <= result['freq_array'][1:])
    return result


def align_time(src_data, dst_data):
    src_gpss = np.sort(np.unique(src_data['gps_array']))
    src_int_time = np.median(np.diff(src_gpss))
    src_phase = np.median(np.mod(src_gpss, src_int_time) / src_int_time)
    dst_gpss = np.sort(np.unique(dst_data['gps_array']))
    dst_int_time = np.median(np.diff(dst_gpss))
    dst_phase = np.median(np.mod(dst_gpss, dst_int_time) / dst_int_time)
    assert src_int_time == dst_int_time, f"src and dst int times: {src_int_time} != {dst_int_time}"
    assert np.all(dst_gpss.round() == dst_gpss), f"dst_gpss not integer: {dst_gpss}"

    print(
        f"aligning src_gpss {len(src_gpss), src_phase}) and dst_gpss {len(dst_gpss), dst_phase}{chr(10)}"
        f"src={src_gpss[:2]}..{src_gpss[-2:]}{chr(10)}"
        f"dst={dst_gpss[:2]}..{dst_gpss[-2:]}"
    )

    if np.median(np.abs(src_phase - dst_phase) * 2).round() == 1:
        # print("alignment out of phase by 1/2")
        assert src_gpss[0] > dst_gpss[0], f"{src_gpss[0]} <= {dst_gpss[0]}"
        assert src_gpss[-1] < dst_gpss[-1], f"{src_gpss[0]} >= {dst_gpss[0]}"
        start_extra = int(((src_gpss[0] - dst_gpss[0] + src_int_time/2)/src_int_time).round())
        end_extra = int(((dst_gpss[-1] - src_gpss[-1] + src_int_time/2)/src_int_time).round())
        src_shape = src_data['weights'].shape
        ext_shape = (int(start_extra + end_extra + src_shape[0]), *src_shape[1:])
        # print(f"{start_extra=} {end_extra=} {ext_shape=}")
        # extend src to fit new shape
        ext_weights = np.zeros(ext_shape)
        ext_weights[start_extra:start_extra+src_shape[0], :, :] = src_data['weights']
        ext_gpss = np.arange(
            src_gpss[0] - start_extra * src_int_time,
            dst_gpss[-1] + src_int_time,
            src_int_time
        ).astype(np.int64)
        print(f"ext={ext_gpss[:2]}..{ext_gpss[-2:]}")
        # pairwise mean of new weights and gpss
        new_weights = np.mean(np.stack([ext_weights[:-1], ext_weights[1:]]), axis=0)
        new_gpss = np.mean(np.stack([ext_gpss[:-1], ext_gpss[1:]]), axis=0).astype(np.int64)
        print(f"new={new_gpss[:2]}..{new_gpss[-2:]}")
        src_data['weights'] = new_weights
        src_data['gps_array'] = new_gpss
        assert np.all(src_data['weights'].shape[1:] == dst_data['weights'].shape[1:]), f"{src_data['weights'].shape} != {dst_data['weights'].shape}"
        assert len(src_data['gps_array']) == len(dst_gpss), \
            f"{len(src_data['gps_array'])} != {len(dst_gpss)}"

    else:
        raise UserWarning(
            f"can't align src {len(src_gpss)} and dst {len(dst_gpss)}{chr(10)}"
            f"src={src_gpss[:2]}..{src_gpss[-2:]}{chr(10)}"
            f"dst={dst_gpss[:2]}..{dst_gpss[-2:]}"
        )


def broadcast_baseline(src_data, dst_data):
    dst_idxs = np.unique(dst_data['gps_array'], return_inverse=True)[1]
    new_weights = src_data['weights'][dst_idxs, ...]
    src_data['weights'] = new_weights
    assert np.all(src_data['weights'].shape == dst_data['weights'].shape), f"{src_data['weights'].shape} != {dst_data['weights'].shape}"


def main():
    if len(sys.argv) > 1:
        args = parse_args()
    else:
        # is being called directly from nextflow
        args = parse_args(shlex.split("""${args}"""))
    print(vars(args))
    src_datas = []
    for path in args.src:
        ext = path.split('.')[-1]
        if ext == 'uvfits':
            d = read_uvfits(path)
        elif ext == 'h5':
            d = read_h5(path)
        src_datas.append(d)
        shp = d['weights'].shape
        print(f"{path} freqs[0]={d['freq_array'][0]} weights({shp})=[{d['weights'][0,:6,0]}..{d['weights'][0,(shp[0]//2):(shp[0]//2+6),0]}]")

    if len(src_datas) == 1:
        src_data = src_datas[0]
    elif src_datas[0]['freq_array'][0] != src_datas[1]['freq_array'][0]:
        src_data = concat_data_freq(src_datas)
    else:
        raise UserWarning("src layout not supported")

    if args.dst:
        dst_data = read_uvfits(args.dst)
        dst_unique_gpss = np.unique(np.sort(dst_data['gps_array']))
        dst_data['weights'] = dst_data['weights']
        assert len(src_data['freq_array']) == len(dst_data['freq_array'])
        if len(src_data['gps_array']) != len(dst_unique_gpss):
            align_time(src_data, dst_data)
        src_unique_gpss = np.unique(np.sort(src_data['gps_array']))
        assert len(src_unique_gpss) == len(dst_unique_gpss), (
            f"src_gpss {len(src_unique_gpss)} != dst_unique_gpss {len(dst_unique_gpss)}{chr(10)}"
            f"src={src_unique_gpss[:2]}..{src_unique_gpss[-2:]}{chr(10)}"
            f"dst={dst_unique_gpss[:2]}..{dst_unique_gpss[-2:]}"
        )
        if ~np.all(dst_data['weights'].shape == src_data['weights'].shape):
            broadcast_baseline(src_data, dst_data)
        # assert np.all(dst_data['freq_array'] == src_freqs), f"{dst_data['freq_array']} != {src_freqs}" # RTS has different freqs?
        print("mixing")
        # mixed_weights = np.prod(np.stack([src_data['weights'], dst_data['weights']]), axis=0)
        set_weights(args.dst, dst_data, src_data['weights'] * dst_data['weights'])

if __name__ == '__main__':
    main()
