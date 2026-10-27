"""Q3-K8s: K8 in the strong form (both networks with (H1)-(H3) and rows and columns of Gamma pairwise linearly
independent).  Same question as the first-round k8.py (its helpers are imported unchanged), restricted to the strong
class, with a memory-light canonical-key formulation.

Environment: Python 3.12 (WSL Fedora 44), standard library only; records of the first round (~/gacE-data) and the
round-2 molecularity-3 samples (~/gacE2-data).  Usage (from round2/computation/): python scripts/k8_strong.py [NPROC]
Side N: two strong networks A != B on the labelled species set {0..n-1} (n <= 6) with the same twin classes and the
same minimal siphons but different P(N*).  For each network let K be the lexicographically least relabelling of
(twin classes, minimal siphons) and O the set of relabelled P(N*) over all relabellings reaching K (an orbit of the
automorphism group of K).  A labelled witness pair exists iff two networks have the same K and the union of their O
has at least two elements.  Side D: the same on the reaction set (doubling classes, minimal siphons of N*, P(N)),
r <= 6.  A search: always exits 0; prints the number of witness keys and the smallest witness pair.
"""
import os
import sys
import time
import platform
from itertools import permutations
from multiprocessing import Pool

import e2lib
from e2lib import CNet
import k8  # first-round helpers (twin_partition, permset)


def skey(fam):
    return tuple(sorted(tuple(sorted(b)) for b in fam))


def canon(args):
    m, tw, MS, P = args
    best = None
    orbit = set()
    for pi in permutations(range(m)):
        key = (skey(k8.permset(tw, pi)), skey(k8.permset(MS, pi)))
        if best is None or key < best:
            best = key
            orbit = {skey(k8.permset(P, pi))}
        elif key == best:
            orbit.add(skey(k8.permset(P, pi)))
    return (m, best), frozenset(orbit)


def main():
    t0 = time.time()
    nproc = int(sys.argv[1]) if len(sys.argv) > 1 else 4
    print(f"# k8_strong.py  Python {platform.python_version()} on {platform.system()}; nproc={nproc}")
    files = e2lib.E_FILES + [os.path.join(e2lib.DATA_E2, f) for f in ("sample_5_3_5_8_s3.pkl.gz",
                                                                        "sample_6_3_5_8_s4.pkl.gz")]
    recs = [r for r in e2lib.load(files) if e2lib.H(r) and r["rows_pi"] and r["cols_pi"]]
    print(f"# {len(recs)} strong networks ((H1)-(H3), rows and columns of Gamma pairwise linearly independent)")
    for side in ("N", "D"):
        jobs, idx = [], []
        for i, r in enumerate(recs):
            n, rx = r["n"], r["rx"]
            if side == "N":
                m, rows = n, [tuple(a[t] for a, b in rx) for t in range(n)]
                MS, P = r["N"]["MS"], r["D"]["P"]
            else:
                m, rows = r["r"], [a for a, b in rx]
                MS, P = r["D"]["MS"], r["N"]["P"]
            if m > 6:
                continue
            jobs.append((m, k8.twin_partition(rows), frozenset(MS), frozenset(P)))
            idx.append(i)
        with Pool(nproc) as pool:
            out = pool.map(canon, jobs, chunksize=50)
        groups = {}
        for i, (K, O) in zip(idx, out):
            groups.setdefault(K, []).append((i, O))
        wit = []
        for K, lst in groups.items():
            if len(lst) < 2:
                continue
            allO = set()
            for i, O in lst:
                allO |= O
            if len(allO) >= 2:
                Od = dict(lst)
                ids = sorted(Od)
                pair = min(((recs[a]["n"] + recs[a]["r"] + recs[b]["n"] + recs[b]["r"], a, b)
                            for a in ids for b in ids if a < b and len(Od[a] | Od[b]) >= 2))
                wit.append((pair, K))
        print(f"# side {side}: {len(jobs)} networks, {len(groups)} canonical keys, {len(wit)} witness keys")
        if wit:
            wit.sort()
            (size, a, b), K = wit[0]
            for lab, i in (("first", a), ("second", b)):
                r = recs[i]
                N = CNet(r["n"], r["rx"])
                print(f"#   {lab}: N = {N.describe()}   N* = {N.dual().describe()}   WR={r['N']['WR']}")
            print(f"#   canonical key (up to relabelling): ground set {K[0]}, classes {K[1][0]}, minimal siphons "
                  f"{K[1][1]}")
    print(f"# runtime {time.time() - t0:.1f} s")


if __name__ == "__main__":
    main()
