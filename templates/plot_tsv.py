#!/usr/bin/env python

from argparse import ArgumentParser
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
import sys
import shlex
import pandas as pd
from pprint import pprint

"""
example use with singularity:

```bash
salloc --nodes=1 --mem=200G --time=2:00:00 --clusters=garrawarla --partition=workq --account=mwaeor --ntasks=10 --tmp=200G
module load singularity
mkdir /dev/shm/deleteme
cd /dev/shm/deleteme
cp /astro/mwaeor/dev/nfresults/results-cmt-eor0high-clean1284/ws_stats.tsv .
cp /astro/mwaeor/dev/nfresults/results-cmt-eor0high-clean1284/ws_stats.tsv .
singularity exec --cleanenv --home /astro/mwaeor/dev/mplhome /pawsey/mwa/singularity/mwa_qa/mwa_qa_latest.sif python \
    /pawsey/mwa/mwaeor/dev/MWAEoR-Pipeline/templates/plot_tsv.py \
    "--tsv=ws_stats.tsv" \
    "--x=OBS" \
    "--y=LST DEG" \
    "--title=Fail codes" \
    "--plot=ws_fail.png"
cp ws_fail.png /pawsey/mwa/mwaeor/dev/MWAEoR-Pipeline/test/
```
"""

def get_parser():
    parser = ArgumentParser(
        description="Plot values from a TSV file.")
    parser.add_argument('--tsv', help='tsv file to plot')
    parser.add_argument('--x', help='x column name')
    parser.add_argument('--y', help='y column name')
    parser.add_argument('--c', default=None, help='color column name')
    parser.add_argument('--plot', default=None, help='output plot name')
    parser.add_argument('--title', default=None, help='suptitle text')
    parser.add_argument('--dpi', default=300, help='dots per inch')
    parser.add_argument('--palette', default=None, help='color palette')
    parser.add_argument('--figwidth', default=30, help='figure width [inches]')
    parser.add_argument('--figheight', default=20, help='figure height [inches]')

    return parser

def main():
    parser = get_parser()
    if len(sys.argv) > 1:
        args = parser.parse_args()
    else:
        # is being called directly from nextflow with args ${args}
        args = parser.parse_args(shlex.split('${argstr}'))
    pprint(vars(args))

    df = pd.read_csv(args.tsv, sep='\t')
    pprint([*df.columns])
    for c in ["x", "y", "c"]:
        column_name = getattr(args, c)
        if column_name is None:
            continue
        if column_name not in df.columns:
            raise RuntimeError(f'column {c} not in tsv file')
    query = f'`{args.x}`.notna() and `{args.y}`.notna()'
    if args.c:
        query += f' and `{args.c}`.notna()'
    print(f"{query=}")
    df = df.query(query, engine="python")
    spkwargs = {
        # 'data': df,
        # 'x': args.x,
        # 'y': args.y,
        # 'palette': args.palette,
        # 'edgecolor': 'none',
    }
    if args.c:
        spkwargs['hue'] = args.c
    sns.scatterplot(data=df, x=args.x, y=args.y, palette=args.palette, edgecolor='none', **spkwargs)
    if args.title:
        plt.suptitle(args.title)

    if args.c:
        plt.legend(ncol=1)

    fig = plt.gcf()
    fig.set_size_inches(args.figwidth, args.figheight)

    if args.plot:
        out_name = args.plot
    else:
        out_name = args.tsv.replace('.tsv', '.png')
    plt.savefig(out_name, dpi=args.dpi, bbox_inches='tight')

if __name__ == '__main__':
    main()