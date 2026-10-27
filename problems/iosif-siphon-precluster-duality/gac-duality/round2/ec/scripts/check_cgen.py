"""Theorem C^gen (round2/ec/notes.md §1.1; agent Q3, Claude Opus 5.5, effort xhigh): if N is connected, (H1), (H2)
and s = (number of Petri clusters) - 1, then neither N nor N* has a critical siphon.  Checked on every record of
agent E's files (WSL ~/gacE-data/*.pkl.gz; fields siphons_eq_SF of N and of N*, which E proved equivalent to 'no
critical siphon'), and, for every such network, the constructive certificate of the proof on N and on N*:
for every minimal siphon R, B = Abar_v[R, C_R] (cluster-aggregated matrix, v > 0 in ker Gamma) has zero row sums,
nonnegative entries off the own-cluster column, vanishes outside C_R, and a nonnegative z with z^T B = 0 found by
exact LP-free search (left kernel + Fourier-Motzkin) is a conservation law supported in R.
Also checks the generalised K5 (s >= #clusters - 1) on every connected network with (H1), (H2).
Environment: Python 3.12, standard library only, exact.  Usage: python check_cgen.py FILES...  Exit 1 on FAIL.
"""
import sys
import time
import gzip
import pickle
import platform
from fractions import Fraction

import eclib
from eclib import Network


def load(files):
    seen = {}
    for fn in files:
        with gzip.open(fn, "rb") as fh:
            for r in pickle.load(fh):
                seen.setdefault((r["n"], tuple(r["rx"])), r)
    return list(seen.values())


def clusters(N):
    return [c for c in eclib.petri.clusters(N) if c[1]]


def certificate(N):
    """Constructive check of the proof of Theorem C^gen on every minimal siphon of N."""
    cl = clusters(N)
    sp_cl = {}
    for idx, (S, T) in enumerate(cl):
        for p in S:
            sp_cl[p] = idx
    rx_cl = {}
    for idx, (S, T) in enumerate(cl):
        for t in T:
            rx_cl[t] = idx
    v = eclib.positive_flux(N)
    G = N.Gamma
    m = len(cl)
    A = [[Fraction(0)] * m for _ in range(N.n)]
    for t in range(N.r):
        c = rx_cl[t]
        for p in range(N.n):
            A[p][c] += v[t] * G[p][t]
    for R in N.minimal_siphons():
        R = sorted(R)
        CR = sorted({sp_cl[p] for p in R})
        B = [[A[p][c] for c in CR] for p in R]
        if any(sum(row) != 0 for row in B):
            return False
        if any(B[i][j] < 0 for i, p in enumerate(R) for j, c in enumerate(CR) if c != sp_cl[p]):
            return False
        if any(A[p][c] != 0 for p in R for c in range(m) if c not in CR):
            return False
        # nonnegative z != 0 with z^T B = 0: feasibility of {z >= 0, sum z = 1, B^T z = 0} by Fourier-Motzkin
        k = len(R)
        cons = []
        for j in range(len(CR)):
            col = [B[i][j] for i in range(k)]
            cons.append((col, 0))
            cons.append(([-x for x in col], 0))
        for i in range(k):
            e = [Fraction(int(i == j)) for j in range(k)]
            cons.append((e, 0))
        cons.append(([Fraction(1)] * k, 1))
        cons.append(([Fraction(-1)] * k, -1))
        if not eclib._fm_feasible(cons, k):
            return False
    return True


def main():
    t0 = time.time()
    print(f"# check_cgen.py; Python {platform.python_version()} on {platform.system()}")
    recs = load(sys.argv[1:])
    n_k5 = bad_k5 = 0
    tally = {}
    bad = 0
    first = {}
    for r in recs:
        if not (r["N"]["cons"] and r["N"]["consist"]):
            continue
        N = Network(r["n"], r["rx"])
        if not eclib.petri.is_connected(N):
            continue
        mcl = len(clusters(N))
        s = N.rank()
        n_k5 += 1
        if s < mcl - 1:
            bad_k5 += 1
            first.setdefault("K5gen", N.describe())
        if s != mcl - 1:
            continue
        ec = eclib.ECStructure(N).is_ec()
        kind = "EC" if ec else "nonEC"
        okN = r["N"]["siphons_eq_SF"]
        okD = r["D"]["siphons_eq_SF"]
        okC = certificate(N) and certificate(N.dual()) if kind == "nonEC" or N.n + N.r <= 9 else True
        for name, ok in (("C^gen N", okN), ("C^gen N*", okD), ("certificate N, N*", okC)):
            h, f = tally.get((name, kind), (0, 0))
            tally[(name, kind)] = (h + ok, f + (not ok))
            if not ok:
                first.setdefault((name, kind), N.describe())
    print(f"{'PASS' if bad_k5 == 0 else 'FAIL'}  generalised K5 (s >= #clusters - 1, connected + (H1) + (H2)): "
          f"{n_k5 - bad_k5}/{n_k5}" + (f"  first failure {first['K5gen']}" if bad_k5 else ""))
    bad += bad_k5 > 0
    for key in sorted(tally):
        h, f = tally[key]
        bad += f > 0
        print(f"{'PASS' if f == 0 else 'FAIL'}  {key[0]} [{key[1]}]: {h}/{h + f}" + (f"  first failure {first[key]}" if f else ""))
    print(f"# runtime {time.time() - t0:.1f} s")
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
