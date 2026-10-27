"""Exhaustive enumeration of ALL small networks that are conservative and consistent ((H1) and (H2)), with records.

Environment: Python 3.12 (WSL Fedora 44), standard library only, exact rational arithmetic (crn.py via ecore.py).
Usage: python enum_all.py N_SPECIES MAX_MOLECULARITY MAX_REACTIONS [--nproc K] [--out PATH]

Networks: sets of at most R distinct reactions y -> y' (y != y', nonzero complexes of molecularity <= m), every
species occurring, one representative per species-permutation orbit (ecore.canonical), as in enumerate_small.py of
result 17 (same family definition, independent code).  A record (ecore.record) is kept for every network with (H1)
and (H2); (H3) is a flag of the record.  Prints the counts and PASS/FAIL for:
  C1  representatives pairwise distinct;
  C2  when the family is one of those of enumerate_small.py (result 17, certificates/enumeration_*.log), the numbers
      of networks with (H1)-(H3), with (H1)-(H3)+WR, and with (H1),(H2) but equal rows, agree with those logs.
"""
import sys
import time
import gzip
import pickle
import platform
from itertools import combinations
from multiprocessing import Pool

from ecore import CNet, canonical, complexes, record

ARGS = None


def worker(first):
    n, m, R = ARGS
    C = complexes(n, m)
    RX = sorted((a, b) for a in C for b in C if a != b)
    rest = RX[first + 1:]
    out = []
    for k in range(R):
        for tail in combinations(rest, k):
            rx = (RX[first],) + tail
            if not all(any(a[i] or b[i] for a, b in rx) for i in range(n)):
                continue
            ok = True
            for i in range(n):  # quick necessary test for (H1)+(H2): each species both gained and lost, or neither
                cons = any(b[i] < a[i] for a, b in rx)
                prod_ = any(b[i] > a[i] for a, b in rx)
                if cons != prod_:
                    ok = False
                    break
            if not ok:
                continue
            if canonical(n, rx) != rx:
                continue
            N = CNet(n, rx)
            if not (N.is_consistent() and N.is_conservative()):
                continue
            out.append(record(rx, n))
    return out


def main():
    global ARGS
    argv = sys.argv[1:]
    n, m, R = int(argv[0]), int(argv[1]), int(argv[2])
    nproc = int(argv[argv.index("--nproc") + 1]) if "--nproc" in argv else 8
    out = argv[argv.index("--out") + 1] if "--out" in argv else None
    ARGS = (n, m, R)
    t0 = time.time()
    print(f"# enum_all.py n={n} m={m} R={R} nproc={nproc}; Python {platform.python_version()} on {platform.system()}")
    C = complexes(n, m)
    nR = len(C) * (len(C) - 1)
    recs = []
    with Pool(nproc) as pool:
        for part in pool.imap_unordered(worker, range(nR)):
            recs.extend(part)
    recs.sort(key=lambda r: r["rx"])
    for r in recs:
        r["family"] = ("all", n, m, R)
    c1 = len({r["rx"] for r in recs}) == len(recs)
    nH = sum(r["N"]["rows_distinct"] for r in recs)
    print(f"# {len(recs)} networks with (H1),(H2); {nH} also with (H3); "
          f"{sum(r['N']['WR'] for r in recs if r['N']['rows_distinct'])} of these weakly reversible")
    if out:
        with gzip.open(out, "wb") as fh:
            pickle.dump(recs, fh, protocol=5)
        print(f"# records written to {out}")
    print(f"{'PASS' if c1 else 'FAIL'} C1: representatives pairwise distinct")
    expected = {(2, 2, 5): (100, 16, 0), (3, 2, 4): (3881, 112, 2), (3, 3, 3): (11332, 123, 77),
                (4, 1, 12): (84, 84, 0), (4, 2, 3): (116, 4, 35)}
    got = (nH, sum(r["N"]["WR"] for r in recs if r["N"]["rows_distinct"]), len(recs) - nH)
    c2 = True
    if (n, m, R) in expected:
        c2 = got == expected[(n, m, R)]
        print(f"{'PASS' if c2 else 'FAIL'} C2: (H, H+WR, H1H2 with equal rows) = {got}, result 17 logs: "
              f"{expected[(n, m, R)]}")
    c1 = c1 and c2
    print(f"# runtime {time.time() - t0:.1f} s")
    sys.exit(0 if c1 else 1)


if __name__ == "__main__":
    main()
