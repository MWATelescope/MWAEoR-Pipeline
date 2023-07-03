#!/usr/bin/env python

from astropy.io import fits
from json import dump as json_dump
from collections import OrderedDict
from argparse import ArgumentParser
import pandas as pd
import numpy as np
import sys


"""
example:

```bash
salloc --nodes=1 --mem=350G --time=8:00:00 --clusters=garrawarla --partition=workq --account=mwaeor --ntasks=20 --tmp=500G
module load singularity
mkdir -p /nvmetmp/deleteme
cd /nvmetmp/deleteme
export obsid=1366000096
cp /astro/mwaeor/dev/nfdata/${obsid}/raw/${obsid}.metafits .
singularity exec --cleanenv --home /astro/mwaeor/dev/mplhome -B /nvmetmp /pawsey/mwa/singularity/mwa_qa/mwa_qa_latest.sif python \
    /pawsey/mwa/mwaeor/dev/MWAEoR-Pipeline/templates/metajson.py \
    --metafits=${obsid}.metafits \
    --json="${obsid}.json" \
    --tsv="${obsid}_tiledata.tsv"

cp ${obsid}.json /astro/mwaeor/dev/MWAEoR-Pipeline/
```
"""

def get_parser():
    parser = ArgumentParser(
        description="get basic metadata from metafits file in json format")

    parser.add_argument('--metafits', help='source metafits file')
    parser.add_argument('--json', help='Name of output json file', default=None)
    parser.add_argument('--tsv', help='Name of output tsv file', default=None)
    parser.add_argument('--txt', help='Name of output txt file', default=None)

    return parser

def main():

    parser = get_parser()

    if len(sys.argv) > 1:
        args = parser.parse_args()
    else:
        # is being called directly from nextflow
        args = parser.parse_args([
            "--metafits=${metafits}",
            "--json=${json}",
            "--tsv=${tsv}",
        ]  # + shlex.split("${args}")
        )

    with fits.open(args.metafits) as hdus:
        data = OrderedDict()
        mpr = hdus['PRIMARY']
        for key in mpr.header:
            if type(mpr.header[key]) not in [bool, int, str, float]:
                continue
            data[key] = mpr.header[key]
            if key in [ 'RECVRS', 'DELAYS', 'CHANNELS', 'CHANSEL' ]:
                data[key] = [*filter(None, map(lambda t: t.strip(), data[key].split(',')))]

        mtd = hdus['TILEDATA'].data
        metafits_inputs = pd.DataFrame(dict([
            (col, mtd.field(col).tolist())
            for col in ['Antenna', 'Tile', 'Rx', 'Slot', 'Pol', 'TileName', 'Flag', 'Length', 'North', 'East', 'Height', 'Gains', 'Flavors']
        ])).sort_values(['Antenna', 'Pol'])
        data['INPUTS'] = metafits_inputs.to_dict(orient='records')
        data['FLAGGED_INPUTS'] = [
            r['TileName'] + r['Pol']
            for r in mtd
            if r['Flag'] == 1
        ]
        data['DELAYS'] = mtd['Delays'].tolist()

    if args.tsv:
        with open(args.tsv, 'w') as out:
            metafits_inputs.to_csv(out, sep='\t', index=False)
        print(f"wrote to {args.tsv}")

    if args.json:
        with open(args.json, 'w') as out:
            json_dump(data, out, indent=4)
        print(f"wrote to {args.json}")

    if args.txt:
        with open(args.txt, 'w') as out:
            slots_per_rx = 8
            all_rxs = np.unique(mtd['Rx'])
            rx_tiles = OrderedDict((rx, [" " * 13] * slots_per_rx) for rx in all_rxs)
            for rx, slot, antenna, tile, tilename, flag in zip(
                mtd['Rx'],
                mtd['Slot'],
                mtd['Antenna'],
                mtd['Tile'],
                mtd['TileName'],
                mtd['Flag'],
            ):
                rx_tiles[rx][slot-1] = f"{antenna:>3}:{tilename:<8}{'F' if flag else ' '}"

            for rx, tiles in rx_tiles.items():
                print(f"{int(rx):>2} " + " | ".join(map(str, tiles)), file=out)
        print(f"wrote to {args.txt}")

if __name__ == "__main__":
    main()