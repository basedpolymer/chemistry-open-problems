#!/usr/bin/env python3
"""Check every integer mesh certificate in certificates/meshes against the
claims in certificates/claims.json (exact integer/rational arithmetic).
Each claim: the explicit partition stored in the file is embedded (exact cone
test) and its DFM energy is <= F_upper (exact), and F_upper < threshold.
Prints PASS/FAIL per claim; exit status 1 if any FAIL.
Runtime: ~1-3 minutes in total.  Environment: Python 3, numpy (file loading)."""
import os, sys, json
from fractions import Fraction
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from certify_upper_bound import certify, load

ROOT = os.path.join(HERE, "..")
claims = json.load(open(os.path.join(ROOT, "certificates", "claims.json")))
bad = 0
for c in claims:
    path = os.path.join(ROOT, "certificates", "meshes", c["file"])
    if not os.path.exists(path):
        print("FAIL", c["name"], "missing", c["file"]); bad += 1; continue
    ok, res = certify(load(path), verbose=False)
    thr = Fraction(c["threshold"])
    passed = ok and res["F_upper"] < thr
    bad += (not passed)
    print("%s %-8s box=%-5s L%d ntri=%7d  inf F <= %.9f  < %s  (%s)" % ("PASS" if passed else "FAIL", c["name"], c["box"],
          c["level"], res["ntri"], res["F_upper"], c["threshold"], c["meaning"]))
sys.exit(1 if bad else 0)
