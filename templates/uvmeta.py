#!/usr/bin/env python

from astropy.io import fits
import astropy.units as u
from astropy.time import Time
from astropy.coordinates import EarthLocation
from json import dump as json_dump
import numpy as np
from math import ceil
from mwa_qa.read_uvfits import UVfits
import sys
from argparse import ArgumentParser


"""
example:

```bash
salloc --nodes=1 --mem=350G --time=8:00:00 --clusters=garrawarla --partition=workq --account=mwaeor --ntasks=20 --tmp=500G
module load singularity
mkdir -p /dev/shm/deleteme
cd /dev/shm/deleteme
export obsid=1090012424
cp /astro/mwaeor/dev/nfdata/${obsid}/cal/hyp_${obsid}_sub_30l_src4k_8s_80kHz.uvfits vis.uvfits
singularity exec --cleanenv --home /astro/mwaeor/dev/mplhome /pawsey/mwa/singularity/mwa_qa/mwa_qa_latest.sif python \
    /astro/mwaeor/dev/MWAEoR-Pipeline/templates/uvmeta.py \
    --uvfits=vis.uvfits \
    --output_name="${obsid}_uvmeta.json"
cp ${obsid}_uvmeta.json /astro/mwaeor/dev/MWAEoR-Pipeline/
```
"""


def get_parser():
    parser = ArgumentParser(
        description="get basic metadata from uvfits file in json format")

    parser.add_argument('--uvfits', help='source uvfits file')
    parser.add_argument(
        '--output_name', help='Name of output plot file', default=None)

    return parser


def main():

    parser = get_parser()

    if len(sys.argv) > 1:
        args = parser.parse_args()
    else:
        # is being called directly from nextflow
        args = parser.parse_args([
            "--uvfits=${vis}",
            "--output_name=${uvmeta}",
        ]  # + shlex.split("${args}")
        )

    uv = UVfits(args.uvfits)

    data = {}
    data['freqs'] = uv.freq_array.tolist()
    mwa_loc = EarthLocation.of_site('mwa')
    times = Time(uv.unique_times, format='jd', scale='utc',
                 location=mwa_loc, precision=3)

    with fits.open(args.uvfits) as hdus:
        # vis_hdu = hdus['PRIMARY']
        ant_hdu = hdus['AIPS AN']
        if dut1 := ant_hdu.header.get('UT1UTC'):
            print(f"updating delta_ut1_utc {times.delta_ut1_utc} -> {dut1}")
            times.delta_ut1_utc = dut1
    # lsts_rad = times.sidereal_time('apparent').radian.tolist()
    # gps_times = times.gps.tolist()
    # iso_times = times.iso.tolist()
    # jd1s, jd2s = [*zip(*map(lambda t: (t.jd1, t.jd2), times.jd.tolist()))]
    data['times'] = [
        {
            "gps": float(f"{time.to_value('gps'):.3f}"),
            "iso": time.to_value('iso', 'date_hms'),
            "lst_rad": time.sidereal_time('apparent').radian,
            "jd1": time.jd1,
            "jd2": time.jd2,
        }
        for time in times
    ]

    out_file = open(args.output_name, 'w') if args.output_name else sys.stdout
    json_dump(data, out_file, indent=4)


if __name__ == '__main__':
    main()
