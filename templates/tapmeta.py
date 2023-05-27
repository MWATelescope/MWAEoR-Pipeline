#!/usr/bin/env python3
from pyvo.dal import TAPService
import json

tap = TAPService("http://vo.mwatelescope.org/mwa_asvo/tap")
table = tap.search("SELECT obs_id,${obs_fields.join(',')} FROM mwa.observation WHERE obs_id = ${obsid}").table
row = table[0]

with open("${tapmeta}", "w") as f:
    json.dump(dict(row), f, indent=4, default=str)