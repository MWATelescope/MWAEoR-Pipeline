#!/usr/bin/env python
import math
import os
from itertools import product
from pathlib import Path

from PIL import Image
import numpy as np
from astropy.io import fits
import sys
from glob import glob


def make_domain_colouring(
    complex_field: np.ndarray,
    cyclic_colourmap: np.ndarray,
    gamma: float
) -> np.ndarray:
    num_colours = cyclic_colourmap.shape[0]
    radius = (complex_field.real**2+complex_field.imag**2)**(1 / 2)
    angle = np.arctan2(complex_field.imag, complex_field.real)
    turning_number = ((angle / math.tau) + 0.5) % 1
    normalised_radius = radius / radius.max()
    quantized_angle = (turning_number * num_colours).astype(int)
    angle_colours = cyclic_colourmap[quantized_angle, :]
    image = (normalised_radius[..., None]**gamma) * (angle_colours)
    return image


def open_fits_as_array(filepath):
    return fits.open(filepath)[0].data.astype(float)


if __name__ == '__main__':
    out_prefix=""
    if len(sys.argv) > 1:
        roma0 = sys.argv[-1]
    else:
        roma0 = "${roma0}"
        out_prefix = "${out_prefix}"
    stokes = {
        stoke: open_fits_as_array(glob(f"*-{stoke}-image.fits")[0])
        for stoke in ['I', 'Q', 'U', 'V']
    }
    complex_fields = {
        'I': stokes['I'] + 0j,
        'L': stokes['Q'] + stokes['U'] * 1j,
        'V': 0 + stokes['V'] * 1j,
        'IV': stokes['I'] + stokes['V'] * 1j,
    }

    cyclic_colourmap = np.load(roma0)
    gammas = [1, 1 / 2.33]

    for (field_name, gamma) in product(complex_fields, gammas):
        domain_colouring = make_domain_colouring(complex_fields[field_name], cyclic_colourmap, gamma)
        image = (domain_colouring * 256)
        save_path = str(f'{out_prefix}{field_name}-{gamma:.2f}.png')
        reshaped = image[0, 0][..., [2, 1, 0]].astype(np.uint8)
        with open(save_path, "wb") as f:
            Image.fromarray(reshaped, mode="RGB").save(f)

