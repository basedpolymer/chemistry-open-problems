#!/usr/bin/env python3
"""Regenerate the integer mesh certificates listed in certificates/claims.json
(deterministic relaxations with run_dfm.py, then export on the 2^-44 grid).
usage: make_certificates.py [NAME ...]     Runtime: minutes.  Single thread."""
import os, sys, json, shutil
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import dfm  # noqa
import numpy as np
import pickle
from run_dfm import run_one
from certify_upper_bound import mesh_to_cert

ROOT = os.path.join(HERE, "..")
claims = json.load(open(os.path.join(ROOT, "certificates", "claims.json")))
only = set(sys.argv[1:])
tmp = os.path.join(ROOT, "work", "cert_runs")
for c in claims:
    if only and c["name"] not in only:
        continue
    run_one(c["name"], c["box"], c["level"], tmp)
    st = os.path.join(tmp, c["name"].replace(":", "_") + "__" + c["box"], "state.pkl")
    M = pickle.load(open(st, "rb"))
    assert M.level == c["level"]
    np.savez_compressed(os.path.join(ROOT, "certificates", "meshes", c["file"]), **mesh_to_cert(M))
    print("exported", c["file"])
