#!/usr/bin/env python

from astropy.io import fits
from json import dump as json_dump
from collections import OrderedDict
import numpy as np

with \
    open('${metrics}', 'w') as out, \
    fits.open("${soln}") as hdus \
:
    data = OrderedDict()
    # Reads hyperdrive results from fits file into json
    results = hdus["RESULTS"].data[:,:]
    data["RESULTS"] = hdus["RESULTS"].data[:,:].tolist()
    data["DIPOLE_GAINS"] = hdus['TILES'].data["DipoleGains"].tolist()
    data["TOTAL_CONV"] = np.nanprod(results)
    print(repr(data))
    json_dump(data, out, indent=4)