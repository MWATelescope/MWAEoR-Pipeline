#!/usr/bin/env python
"""Convert a hyperdrive JSON source list to a DP3 ("makesourcedb") sky model.

Each DD direction (from the directions tsv) becomes a patch named after the source; every other
source goes into the patch "field". Flux-density lists are turned into a power law from the two
entries nearest the median list frequency (log-log slope); shapelets are written as Gaussians.
"""

import csv
import json
import math

REF_FREQ = 150e6

with open("${directions_tsv}") as f:
    directions = [row["NAME"] for row in csv.DictReader(f, delimiter="\\t")]
with open("${reduced_json}") as f:
    srclist = json.load(f)


def power_law(flux_type):
    """-> (I at REF_FREQ, spectral index) for any hyperdrive flux type."""
    if "power_law" in flux_type:
        pl = flux_type["power_law"]
        fd, si = pl["fd"], pl["si"]
        return fd["i"] * (REF_FREQ / fd["freq"]) ** si, si
    if "curved_power_law" in flux_type:
        # drop the curvature term: DP3 "LogarithmicSI" takes a polynomial in log(nu/nu0), but
        # keep it simple and comparable to hyperdrive's power law reading of DP3 models.
        pl = flux_type["curved_power_law"]
        fd, si = pl["fd"], pl["si"]
        return fd["i"] * (REF_FREQ / fd["freq"]) ** si, si
    fds = sorted(flux_type["list"], key=lambda fd: fd["freq"])
    if len(fds) == 1:
        fd = fds[0]
        si = -0.8
        return fd["i"] * (REF_FREQ / fd["freq"]) ** si, si
    # the two entries nearest the median frequency
    median = fds[len(fds) // 2]["freq"]
    fds = sorted(fds, key=lambda fd: abs(fd["freq"] - median))[:2]
    (f1, i1), (f2, i2) = [(fd["freq"], fd["i"]) for fd in sorted(fds, key=lambda fd: fd["freq"])]
    if i1 <= 0 or i2 <= 0 or f1 == f2:
        si = -0.8
        return i1 * (REF_FREQ / f1) ** si, si
    si = math.log(i2 / i1) / math.log(f2 / f1)
    return i1 * (REF_FREQ / f1) ** si, si


def deg_to_hms(ra):
    h = ra / 15.0
    hh = int(h)
    mm = int((h - hh) * 60)
    ss = ((h - hh) * 60 - mm) * 60
    return "%02d:%02d:%09.6f" % (hh, mm, ss)


def deg_to_dms(dec):
    sign = "-" if dec < 0 else "+"
    d = abs(dec)
    dd = int(d)
    mm = int((d - dd) * 60)
    ss = ((d - dd) * 60 - mm) * 60
    return "%s%02d.%02d.%08.5f" % (sign, dd, mm, ss)


def clean(name):
    return "".join(c if c.isalnum() or c in "+-_." else "_" for c in name)


rows = []
patches = {}
for name, comps in srclist.items():
    patch = clean(name) if name in directions else "field"
    if patch not in patches:
        patches[patch] = (comps[0]["ra"], comps[0]["dec"])
    for k, comp in enumerate(comps):
        flux, si = power_law(comp["flux_type"])
        comp_type = comp["comp_type"]
        if comp_type == "point":
            shape = ("POINT", "", "", "")
        else:
            kind, geom = next(iter(comp_type.items()))
            # hyperdrive sizes are FWHM arcsec, PA degrees east of north: same as DP3
            shape = ("GAUSSIAN", "%.4f" % geom["maj"], "%.4f" % geom["min"], "%.3f" % geom["pa"])
        rows.append((
            "%s_%d" % (clean(name), k), shape[0], patch,
            deg_to_hms(comp["ra"]), deg_to_dms(comp["dec"]),
            "%.6g" % flux, "[%.4f]" % si, "true", "%g" % REF_FREQ,
            shape[1], shape[2], shape[3],
        ))

with open("${skymodel}", "w") as out:
    out.write(
        "Format = Name, Type, Patch, Ra, Dec, I, SpectralIndex, LogarithmicSI, "
        "ReferenceFrequency='%g', MajorAxis, MinorAxis, Orientation\\n" % REF_FREQ
    )
    # directions first, then the field patch
    for patch in [clean(d) for d in directions if clean(d) in patches] + (["field"] if "field" in patches else []):
        ra, dec = patches[patch]
        out.write(", , %s, %s, %s\\n" % (patch, deg_to_hms(ra), deg_to_dms(dec)))
    for row in rows:
        out.write(", ".join(row) + "\\n")

print("wrote %d components in %d patches (%d directions) to ${skymodel}" % (len(rows), len(patches), len(directions)))
