#!/usr/bin/env python3
"""Second, independently written exact check of the networks in certificates/solutions/: each one is a feasible
point of the ORIGINAL transshipment MILP of Letsios-Kouyialis-Misener (2018) (their Pyomo model
lib/exact_methods/transshipment_mip_solver.py), with their own big-M parameters U_ij.

usage: python3 check_network_lkm18.py [instance ...]      (default: every *.json in certificates/solutions/)

This script shares no code with hen_instance.py / hen_exact.py: it has its own .dat and .bgm parsers.
For each network (rational heats q_ijt from the JSON file) it sets y_ij = 1 iff sum_t q_ijt > 0 and
r_{i,0} = 0, r_{i,t+1} = r_{i,t} + QH_it - sum_j q_ijt, and verifies in exact rational arithmetic
(fractions.Fraction) every constraint of the LKM18 transshipment model:
  * q_ijt >= 0;  r_{i,t} >= 0 for t = 0..k;  sum_i r_{i,0} = 0 and sum_i r_{i,k} = 0;
  * hot conservation  sum_j q_ijt + r_{i,t+1} = QH_it + r_{i,t}  (holds by construction of r);
  * cold conservation sum_i q_ijt = QC_jt;
  * matched streams   sum_t q_ijt <= U_ij y_ij  for the three big-M files of the LKM18 repository
    (data/bigM_parameters/<set>/{simple,gundersen,greedy}/<instance>.bgm);
  * the objective  matches = sum_ij y_ij  equals the number of matches stated in the JSON file and in
    certificates/solutions/claimed_upper_bounds.json.
Hence the optimum of the LKM18 MILP (with any of the three big-M files) is <= the stated number of matches.
PASS/FAIL per network and per claim; exit code 1 on any failure.
Data: clone github.com/dimletsios/min_matches_heuristics into certificates/sources/min_matches_heuristics
(or set HEN_REPO).  Environment: Python >= 3.10, standard library only.  Runtime: about 10 s for 22 networks.
"""
import glob
import json
import os
import re
import sys
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.environ.get("HEN_REPO", os.path.join(HERE, "..", "certificates", "sources", "min_matches_heuristics"))
SOL = os.path.join(HERE, "..", "certificates", "solutions")


def read_dat(path):
    text = open(path).read()
    n = int(re.search(r"^n=(\d+)", text, re.M).group(1))
    m = int(re.search(r"^m=(\d+)", text, re.M).group(1))
    k = int(re.search(r"^k=(\d+)", text, re.M).group(1))
    QH = {}
    QC = {}
    for kind, idx, body in re.findall(r"^Q([HC])\[(\d+)\]:(.*)$", text, re.M):
        tgt = QH if kind == "H" else QC
        tok = body.split()
        for a in range(0, len(tok), 2):
            t = int(tok[a].lstrip("T"))
            key = (int(idx), t)
            if key in tgt:
                raise ValueError("duplicate entry")
            tgt[key] = Fraction(tok[a + 1])
    return n, m, k, QH, QC


def read_bgm(path):
    U = {}
    for i, j, v in re.findall(r"U\[(\d+),(\d+)\]=\s*(\S+)", open(path).read()):
        U[int(i), int(j)] = Fraction(v)
    return U


def find_one(pattern):
    hits = glob.glob(pattern)
    if len(hits) != 1:
        raise FileNotFoundError(pattern)
    return hits[0]


def check(name, claimed):
    out = []
    dat = find_one(os.path.join(REPO, "data", "mip_instances", "*", name + ".dat"))
    tset = os.path.basename(os.path.dirname(dat))
    n, m, k, QH, QC = read_dat(dat)
    d = json.load(open(os.path.join(SOL, name + ".json")))
    q = {}
    for i, j, t, v in d["q"]:
        if (i, j, t) in q:
            raise ValueError("duplicate q")
        q[i, j, t] = Fraction(v)
    ok_range = all(0 <= i < n and 0 <= j < m and 0 <= t < k for (i, j, t) in q)
    out.append(("indices in range", ok_range))
    out.append(("q >= 0", all(v >= 0 for v in q.values())))
    tot = {}
    for (i, j, t), v in q.items():
        tot[i, j] = tot.get((i, j), Fraction(0)) + v
    y = {(i, j): (1 if tot.get((i, j), 0) > 0 else 0) for i in range(n) for j in range(m)}
    # residuals r_{i,0..k} of the LKM18 model
    rmin_ok, rend_ok = True, True
    for i in range(n):
        r = Fraction(0)
        for t in range(k):
            r = r + QH.get((i, t), Fraction(0)) - sum(q.get((i, j, t), Fraction(0)) for j in range(m))
            if r < 0:
                rmin_ok = False
        if r != 0:
            rend_ok = False
    out.append(("r_{i,t} >= 0 for all i, t", rmin_ok))
    out.append(("sum_i r_{i,k} = 0 (and r_{i,0} = 0)", rend_ok))
    cold_ok = all(sum(q.get((i, j, t), Fraction(0)) for i in range(n)) == QC.get((j, t), Fraction(0))
                  for j in range(m) for t in range(k))
    out.append(("cold conservation", cold_ok))
    for variant in ("simple", "gundersen", "greedy"):
        U = read_bgm(os.path.join(REPO, "data", "bigM_parameters", tset, variant, name + ".bgm"))
        ok = len(U) == n * m and all(tot.get((i, j), Fraction(0)) <= U[i, j] * y[i, j]
                                     for i in range(n) for j in range(m))
        out.append((f"sum_t q_ijt <= U_ij y_ij ({variant} big-M)", ok))
    nm = sum(y.values())
    listed = {tuple(p) for p in d["matches"]}
    out.append((f"matches = sum y = {nm} = stated {d['number_of_matches']} = claimed {claimed}",
                nm == d["number_of_matches"] == claimed and listed == {p for p, v in y.items() if v}))
    return out


def main():
    claimed = json.load(open(os.path.join(SOL, "claimed_upper_bounds.json")))
    names = sys.argv[1:] or sorted(claimed)
    bad = 0
    for name in names:
        res = check(name, claimed[name])
        good = all(ok for _, ok in res)
        bad += not good
        for claim, ok in res:
            if not ok:
                print(f"FAIL {name}: {claim}")
        print(("PASS" if good else "FAIL") + f" {name}: feasible point of the LKM18 transshipment MILP with "
              f"{claimed[name]} matches (all {len(res)} checks)")
    if bad:
        print(f"SOME FAIL ({bad} networks)")
        sys.exit(1)
    print(f"ALL PASS ({len(names)} networks)")


if __name__ == "__main__":
    main()
