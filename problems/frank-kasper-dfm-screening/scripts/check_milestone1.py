#!/usr/bin/env python3
"""Milestone 1: reproduction of the published DFM values of Reddy et al.,
PNAS 115 (2018) 10233 (arXiv:1806.01726v1, Figs. 7-18), with curved cells.

Claims checked (PASS/FAIL each; exit 1 on any FAIL):
 (a) Flat 'Initial Voronoi' rows (A, I) reproduced to the 6 printed decimals
     (tolerance 1.5e-6) with our inputs for BCC, A15, C15, H, P, delta,
     psigma, M  -> confirms conventions, coordinates and boxes.
 (b) Equal-volume minimal area (Kelvin problem): BCC (Kelvin foam) 1.097251 and
     A15 (Weaire-Phelan) 1.093541 reproduced by Richardson extrapolation of
     levels 1-3 (recomputed here), tolerance 2e-6.
 (c) Unconstrained DFM value F (fixed paper box), from the stored level-0..4
     runs (data/dfm_fixedbox_paper_structures.csv, produced by run_dfm.py +
     analyze.py): |F_extrap - F_published| <= 2e-6 for A15, C15, Z, C14, delta.
 (d) Exact certificates (scripts/check_certificates.py): for BCC, P, M, H, mu a
     certified upper bound strictly below the published value.
 Reported but not claimed: sigma, psigma (box ratio ambiguous: the paper's own
 CVP row for sigma matches a/c ~ 1.93, the caption says 1.9).
usage: check_milestone1.py [--skip-b] ; runtime ~2-5 min (claim b), single thread.
"""
import os
import sys
import csv
import subprocess

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "..")
sys.path.insert(0, HERE)
import dfm  # noqa: E402
import numpy as np  # noqa: E402
from dfm import registry, topology, mesh, energy, constrained  # noqa: E402
from run_dfm import initial_mesh  # noqa: E402

pub = {}
for r in csv.reader(open(os.path.join(ROOT, "data", "published_reddy2018.csv"))):
    if not r or r[0].startswith("#"):
        continue
    pub[(r[0], r[4])] = r
NAMEMAP = {"BCC": "BCC", "A15": "A15", "C15": "C15", "H": "H", "P": "P", "delta": "delta", "psigma": "psigma",
           "M": "M", "Z": "Z", "C14": "C14", "sigma": "sigma", "mu": "mu"}
nfail = 0


def report(ok, msg):
    global nfail
    nfail += (not ok)
    print(("PASS " if ok else "FAIL ") + msg, flush=True)


# (a) flat Voronoi rows
print("(a) flat Voronoi partitions of the paper's inputs")
for nm in ["BCC", "A15", "C15", "H", "P", "delta", "psigma", "M", "Z", "C14", "sigma", "mu"]:
    st = registry.get(nm)
    M = initial_mesh(st)
    d = energy.dfm_values(M)
    r = pub[(NAMEMAP[nm], "InitialVoronoi")]
    Ap, Ip = float(r[5]), float(r[6])
    ok = abs(d["A"] - Ap) <= 1.5e-6 and abs(d["I"] - Ip) <= 1.5e-6
    claimed = nm in ("BCC", "A15", "C15", "H", "P", "delta", "psigma", "M")
    line = "%-7s A=%.6f (pub %.6f)  I=%.6f (pub %.6f)" % (nm, d["A"], Ap, d["I"], Ip)
    if claimed:
        report(ok, line)
    else:
        print("info " + line + ("  [matches]" if ok else "  [inputs differ from the paper's; not claimed]"))

# (b) Kelvin / Weaire-Phelan
if "--skip-b" not in sys.argv:
    print("(b) equal-volume minimal area (Kelvin problem), levels 1-3 + Richardson")
    for nm, target in (("BCC", 1.097251), ("A15", 1.093541)):
        M = initial_mesh(registry.get(nm))
        vals = []
        for lev in range(4):
            if lev:
                M = mesh.refine(M)
            vals.append(constrained.solve(M, objective="area")["A"])
        R = vals[3] - (vals[2] - vals[3]) / 3.0
        report(abs(R - target) <= 2e-6, "%s minimal equal-volume area: levels %s -> %.7f (pub %.6f)"
               % (nm, ["%.7f" % v for v in vals], R, target))

# (c) unconstrained DFM, fixed box
print("(c) unconstrained DFM value, fixed paper box (stored runs, levels 0-4)")
rows = {r["name"]: r for r in csv.DictReader(open(os.path.join(ROOT, "data", "dfm_fixedbox_paper_structures.csv")))}
for nm in ["A15", "C15", "Z", "C14", "delta", "sigma", "psigma", "H", "P", "M", "mu", "BCC"]:
    r = rows[nm]
    Fe, err = float(r["F_extrap"]), float(r["err"])
    Fp = float(pub[(NAMEMAP[nm], "Unconstrained")][9])
    line = "%-7s F_L4=%s  F_extrap=%.7f +- %.1e  published %.6f  diff %+.2e" % (nm, r["F_K"], Fe, err, Fp, Fe - Fp)
    if nm in ("A15", "C15", "Z", "C14", "delta"):
        report(abs(Fe - Fp) <= 2e-6, line)
    else:
        print("info " + line)

# (d) certificates
print("(d) exact certificates")
rc = subprocess.call([sys.executable, os.path.join(HERE, "check_certificates.py")])
report(rc == 0, "all certificates in certificates/claims.json")
print("\n%d FAIL" % nfail)
sys.exit(1 if nfail else 0)
