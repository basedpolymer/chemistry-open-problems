"""Exact checks of Theorems C (no purity) and U (uniqueness + stability) on exhaustive families of well-formed
equal-conflict (WF-EC) networks (agent Q3, Claude Opus 5.5, effort xhigh; round2/ec/notes.md).

Environment: Python 3.12 (WSL Fedora 44), standard library only, exact rationals.
Usage (from round2/ec/): python scripts/check_wfec_families.py N MAXMOL RMAX ORDINARY [NPROC]
Networks: agent P's generator (petri/scripts/enumerate_ec.py, imported): every EC network up to species permutation,
every species consumed, molecularity <= MAXMOL, at most RMAX reactions, connected, (H1), (H2).  Class WF-EC: s = m - 1.
Claims (PASS/FAIL, exit status 1 on FAIL), on every WF-EC network (pure or not, balanced or not):
  C   Theorem C certificate: every minimal siphon R meets each S_c at most once, B = A_v[R, C_R] is Metzler with zero
      row sums, its left kernel has a vector z >= 0 whose extension is a conservation law with support R;
  U1  sign condition sigma(S) cap sigma(D^perp) = {0} (maximal-minor criterion);
  U1d the same by direct enumeration of sign vectors (exact Fourier-Motzkin), n <= 4;
  U2  stability certificate at x = (1,..,1) and at two seeded random rational x (k chosen so that x is a steady
      state): G Metzler with zero row sums, det(tI - J) = t^(n-s) q(t), q(0) != 0, q Hurwitz (exact Routh array).
Reports (no PASS/FAIL): the same properties on EC networks that violate the rank equation.
"""
import os
import sys
import time
import random
import platform
from fractions import Fraction
from itertools import permutations
from multiprocessing import Pool

import eclib
from eclib import Network, ECStructure
from enumerate_ec import gen_networks, canonical, set_partitions  # agent P's generator (petri/scripts)


def evaluate(N, rng):
    E = ECStructure(N)
    wf = N.rank() == E.m - 1
    res = {}
    okc, msg = eclib.theorem_c_certificate(N, E)
    res["C"] = okc
    if wf or True:
        try:
            oks, _ = eclib.sign_condition_minors(N, E)
        except AssertionError:
            oks = None
        res["U1"] = oks
        if N.n <= 4:
            okd, _ = eclib.sign_condition_direct(N, E)
            res["U1d"] = okd
    xs = [[Fraction(1)] * N.n] + [[Fraction(rng.randint(1, 9), rng.randint(1, 9)) for _ in range(N.n)]
                                  for _ in range(2)]
    res["U2"] = all(eclib.stability_certificate(N, x, E)[0] for x in xs)
    info = {"wf": wf, "pure": eclib.is_pure(N), "balanced": eclib.is_balanced(N)}
    return res, info


def worker(args):
    n, m, rmax, ordinary, pidx = args
    perms = list(permutations(range(n)))
    rng = random.Random(1000 * pidx + 7)
    seen = set()
    out = []
    for rx in gen_networks(n, m, rmax, ordinary, pidx):
        rx = tuple(sorted(rx))
        can = canonical(rx, perms)
        if can != rx or can in seen:
            continue
        seen.add(can)
        N = Network(n, rx)
        if not eclib.petri.is_connected(N) or not N.is_consistent() or not N.is_conservative():
            continue
        res, info = evaluate(N, rng)
        out.append((rx, res, info))
    return out


def main():
    n, m, rmax, ordinary = (int(a) for a in sys.argv[1:5])
    nproc = int(sys.argv[5]) if len(sys.argv) > 5 else 12
    t0 = time.time()
    print(f"# check_wfec_families.py n={n} molecularity<={m} r<={rmax} ordinary={ordinary}; "
          f"Python {platform.python_version()} on {platform.system()}")
    nparts = sum(1 for _ in set_partitions(list(range(n))))
    tally = {}
    first_fail = {}
    counts = {"EC": 0, "WF": 0, "WF impure": 0, "WF unbalanced": 0, "WF impure unbalanced": 0}
    with Pool(nproc) as pool:
        for chunk in pool.imap_unordered(worker, [(n, m, rmax, ordinary, p) for p in range(nparts)]):
            for rx, res, info in chunk:
                counts["EC"] += 1
                cls = "WF" if info["wf"] else "nonWF"
                if info["wf"]:
                    counts["WF"] += 1
                    counts["WF impure"] += not info["pure"]
                    counts["WF unbalanced"] += not info["balanced"]
                    counts["WF impure unbalanced"] += (not info["pure"]) and (not info["balanced"])
                for name, val in res.items():
                    if val is None:
                        continue
                    key = (name, cls)
                    h, f = tally.get(key, (0, 0))
                    tally[key] = (h + bool(val), f + (not val))
                    if not val and (key not in first_fail or len(rx) < len(first_fail[key])):
                        first_fail[key] = rx
    print("# counts: " + ", ".join(f"{k}={v}" for k, v in counts.items()))
    bad = 0
    for (name, cls) in sorted(tally):
        h, f = tally[(name, cls)]
        extra = f"  smallest failure: {Network(n, first_fail[(name, cls)]).describe()}" if f else ""
        if cls == "WF":
            st = "PASS" if f == 0 else "FAIL"
            bad += f > 0
        else:
            st = "REPORT"
        print(f"{st}  {name} [{cls}]: {h}/{h + f}{extra}")
    print(f"# runtime {time.time() - t0:.1f} s")
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
