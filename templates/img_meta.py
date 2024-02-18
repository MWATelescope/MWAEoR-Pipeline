#!/usr/bin/env python

from astropy.io import fits
import numpy as np
from argparse import ArgumentParser
import sys

def get_parser():
    parser = ArgumentParser(
        description="get metadata for fits image.")

    parser.add_argument('--fits', default=False, help='fits file to render')
    parser.add_argument('--csv', default=False, help='csv file to write to')
    parser.add_argument('--quantiles', nargs="+", help='quantiles to calculate',
                        default=[0.90, 0.95, 0.99, 1.0], type=float)

    return parser


def main():
    parser = get_parser()

    if len(sys.argv) > 1:
        args = parser.parse_args()
    else:
        # is being called directly from nextflow
        args = parser.parse_args([
            "--fits=${fits}",
            "--csv=${csv}",
        ])
    print(f"{args=}")

    # data = OrderedDict()

    with fits.open(args.fits) as hdus:
        header = hdus[0].header
        naxis = header['NAXIS']
        # deal with the fact that wsclean images have extra axes other than ra/dec
        # by only grabing axes with units of 'deg'
        i = [0] * naxis
        dims = []
        for a in range(1, naxis+1):
            unit = header.get(f'CUNIT{a}')
            dim = header.get(f'NAXIS{a}')
            type_ = header.get(f'CTYPE{a}')
            if unit == 'deg':
                i[naxis - a] = slice(None)
                dims.append(dim)
            elif dim > 1:
                raise ValueError(f"Unexpected dimension {dim} for axis {a}: {type_} with unit {unit}")
        try:
            print(f"indexing: {i}")
            data = hdus[0].data[i]
        except IndexError as e:
            data = hdus[0].data.reshape(dims)
        print(f"{data.shape=}")

    quantiles = np.array(args.quantiles)
    quantiles = np.concatenate((quantiles, (1.0 - quantiles)))
    quantiles.sort()
    print(f"{quantiles=}")
    values = np.quantile(data, quantiles)
    np.savetxt(
        args.csv, np.array([quantiles, values]).T,
        delimiter=',', fmt='%.6e', header='quantile,value'
    )

    # TODO: histogram?
    # import matplotlib.pyplot as plt
    # plt.style.use([astropy_mpl_style, 'dark_background'])


if __name__ == '__main__':
    main()
