#!/usr/bin/env python3
"""Exact feasibility check of every solution file certificates/solutions/<instance>.json.

For each file: the instance is re-read from the Letsios-Kouyialis-Misener data file (exact decimals), and
check_solution() verifies in rational arithmetic (fractions.Fraction):
  * q_ijt >= 0 and positive only on declared matches; each declared match carries positive heat;
  * every cold balance sum_i q_ijt = QC[j,t];
  * every hot residual r_it = sum_{s<=t}(QH[i,s] - sum_j q_ijs) >= 0, and r_{i,k-1} = 0
    (temperature feasibility: heat is only sent to the same or colder intervals);
  * the number of matches equals the claimed value in the table below.
Prints PASS/FAIL per instance; exit code 1 if any FAIL.

Environment: Python >= 3.10, standard library only.  Runtime: about 10 s in total.
Data: certificates/sources/min_matches_heuristics (git clone of github.com/dimletsios/min_matches_heuristics,
commit 2019-02-15; set HEN_DATA to override the data directory).
"""
import os, sys, glob, json, hashlib
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from hen_instance import parse_dat, consistency_report
from hen_exact import load_solution, check_solution

DATA = os.environ.get("HEN_DATA", os.path.join(HERE, "..", "certificates", "sources", "min_matches_heuristics", "data",
                                               "mip_instances"))
SOLDIR = os.path.join(HERE, "..", "certificates", "solutions")

# claimed number of matches of each solution file (upper bounds)
CLAIMED = json.load(open(os.path.join(SOLDIR, "claimed_upper_bounds.json")))


def main():
    allok = True
    for path in sorted(glob.glob(os.path.join(SOLDIR, "*.json"))):
        if path.endswith("claimed_upper_bounds.json"):
            continue
        d, matches, q = load_solution(path)
        name = d["instance"]
        cands = glob.glob(os.path.join(DATA, "*", name + ".dat"))
        if len(cands) != 1:
            print(f"FAIL {name}: data file not found")
            allok = False
            continue
        inst = parse_dat(cands[0])
        rep = consistency_report(inst)
        errs = []
        if not (rep["balanced"] and rep["cascade_nonneg"] and rep["R_matches_file"]):
            errs.append("instance data inconsistent")
        errs += check_solution(inst, matches, q)
        if len(matches) != d["number_of_matches"]:
            errs.append("number_of_matches field wrong")
        claim = CLAIMED.get(name)
        if claim is None or len(matches) != claim:
            errs.append(f"claimed {claim}, file has {len(matches)} matches")
        sha = hashlib.sha256(open(cands[0], "rb").read()).hexdigest()[:16]
        status = "PASS" if not errs else "FAIL"
        allok &= not errs
        print(f"{status} {name}: feasible network with {len(matches)} matches (exact rational check; data sha256 {sha}...)"
              + ("" if not errs else "  errors: " + "; ".join(errs[:5])))
    print("ALL PASS" if allok else "SOME FAIL")
    sys.exit(0 if allok else 1)


if __name__ == "__main__":
    main()
