#!/usr/bin/env python3
"""Side findings on 37sp-yfyv (exact integer max-flow), PASS/FAIL per claim.
 (a) the 36-match set returned as optimal by floating-point SCIP is infeasible: exact max-flow deficit = 1/10;
 (b) none of the six CPLEX/Gurobi match lists of the Letsios et al. repository for 37sp-yfyv can be reduced
     below 37 matches by deleting pairs (greedy deletion with exact feasibility tests; each list has 37 or 38 pairs,
     and the greedy result has 37 for every list).  Note: (b) concerns these lists only; a different exact 36-match
     network exists (certificates/solutions/37sp-yfyv.json).
Environment: Python >= 3.10, standard library.  Runtime ~10 s.
"""
import os, sys, json, glob, re
from fractions import Fraction
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from hen_instance import parse_dat
from hen_exact import flow_for_matches, prune_matches
DATA = os.path.join(HERE, "..", "certificates", "sources", "min_matches_heuristics", "data")
inst = parse_dat(os.path.join(DATA, "mip_instances", "furman_sahinidis", "37sp-yfyv.dat"))
ok_all = True
d = json.load(open(os.path.join(HERE, "..", "certificates", "side", "37sp-yfyv_scip_float_36.json")))
ok, deficit = flow_for_matches(inst, [tuple(p) for p in d["matches"]])
a = (not ok) and deficit == Fraction(1, 10) and len(d["matches"]) == 36
print(("PASS" if a else "FAIL"), "(a) float-SCIP 36-match set infeasible, exact deficit", deficit if not ok else 0)
ok_all &= a
for sol in sorted(glob.glob(os.path.join(DATA, "mip_solutions", "furman_sahinidis", "*", "37sp-yfyv_*.sol"))):
    M = [tuple(map(int, x)) for x in re.findall(r"^y\[(\d+),(\d+)\]=1", open(sol).read(), re.M)]
    P = prune_matches(inst, M)
    b = len(P) == 37
    ok_all &= b
    print(("PASS" if b else "FAIL"), "(b)", "/".join(sol.split("/")[-2:]), f"{len(M)} listed pairs -> {len(P)} after exact greedy pruning")
print("ALL PASS" if ok_all else "SOME FAIL"); sys.exit(0 if ok_all else 1)
