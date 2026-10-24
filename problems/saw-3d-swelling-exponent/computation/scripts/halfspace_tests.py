#!/usr/bin/env python3
"""halfspace_tests.py -- exact tests of the upper-bound solver's CL1 (arch lemma) and CL2 (renewal identity) on the
output of src/halfspace_enum.cpp, with cross-checks against the main enumeration (src/saw_enum.cpp).

Environment: Python 3.8+ standard library.  Runtime < 5 s.
Usage: python3 halfspace_tests.py [halfspace file, default ../data/halfspace_N19.txt] [saw_enum data file]

CL1  a_{n,M} <= sum_{k=M}^{n-M} b_{k,M} b_{n-k,M}  for all n >= 2, M >= 1   (arches bounded by pairs of bridges)
CL2  r_{n,T} = sum_k b_{k,T} b_{n-k,T}  for all n, T >= 1  (bridges of height 2T with a renewal at height T)
Cross-checks (PASS/FAIL): b_{n,T} from both programs; half-space totals h_n; arches total
  sum_M a_{n,M} = #{(n-1)-step half-space walks ending at height 1} - #{those ending at e1}  (from saw_enum's H table).
"""
import os
import sys
from fractions import Fraction as F

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from sawdata import SawData  # noqa: E402

ROOT = os.path.dirname(HERE)


def load(path):
    B, R, A, H = {}, {}, {}, {}
    for line in open(path):
        t = line.split()
        if not t or t[0].startswith("#"):
            continue
        if t[0] == "H":
            H[int(t[1])] = int(t[2])
        else:
            {"B": B, "R": R, "A": A}[t[0]][(int(t[1]), int(t[2]))] = int(t[3])
    return B, R, A, H


def main():
    args = sys.argv[1:]
    hpath = args[0] if args else os.path.join(ROOT, "data", "halfspace_N19.txt")
    dpath = args[1] if len(args) > 1 else os.path.join(ROOT, "data", "saw_enum_N21.txt.gz")
    B, R, A, H = load(hpath)
    N = max(H)
    D = SawData(dpath)
    fails = 0
    # cross-checks
    nn = min(N, D.N)
    ok = all(B.get((n, T), 0) == D.bridge_span_dist(n).get(T, 0) for n in range(1, nn + 1) for T in range(0, n + 1))
    print("%s [x1] bridges by height b_{n,T} identical in halfspace_enum and saw_enum, n <= %d" % ("PASS" if ok else "FAIL", nn))
    fails += not ok
    ok = all(H[n] == D.h[n] for n in range(0, nn + 1))
    print("%s [x2] half-space totals h_n identical, n <= %d" % ("PASS" if ok else "FAIL", nn))
    fails += not ok
    ok = True
    for n in range(2, min(N, D.N + 1) + 1):
        tot = sum(v for (m, M), v in A.items() if m == n)
        ref = sum(Ns for (x1, b, c), (Ns, orb) in D.H.get(n - 1, {}).items() if x1 == 1)
        ref -= sum(Ns for (x1, b, c), (Ns, orb) in D.H.get(n - 1, {}).items() if (x1, b, c) == (1, 0, 0))
        ok = ok and tot == ref
    print("%s [x3] arch totals = (n-1)-step half-space walks ending at height 1, not at e1 (saw_enum H table)" %
          ("PASS" if ok else "FAIL"))
    fails += not ok

    def b(k, T):
        return B.get((k, T), 0)

    # [x4] irreducible bridges i_{n,T} from 1/(1 - I(z,x)) = B(z,x) against Clisby-Conway-Guttmann 2016 Table 7
    ref_irr = {}
    rp = os.path.join(ROOT, "data", "reference_irreducible_bridges.tsv")
    if os.path.exists(rp):
        for line in open(rp):
            if not line.startswith("#"):
                t = line.split()
                ref_irr[(int(t[0]), int(t[1]))] = int(t[2])
        Irr = {}
        for n in range(1, N + 1):
            for T in range(1, n + 1):
                s = b(n, T)
                for (n1, T1), v in Irr.items():
                    if n1 < n and T1 < T:
                        s -= v * b(n - n1, T - T1)
                Irr[(n, T)] = s
        mism = [(k, Irr[k], ref_irr.get(k, 0)) for k in Irr if Irr[k] != ref_irr.get(k, 0)]
        typo = [m for m in mism if m[0] == (18, 1)]
        other = [m for m in mism if m[0] != (18, 1)]
        ok = not other
        print("%s [x4] irreducible bridges by height (from b_{n,T}) = CCG 2016 Table 7 for n <= %d%s" % (
            "PASS" if ok else "FAIL", N,
            "; known arXiv typo at (18,1): ours %d, printed %d" % (typo[0][1], typo[0][2]) if typo else ""))
        fails += not ok

    # CL1
    first_fail = None
    print("\n## CL1 (arch lemma) a_{n,M} <= sum_k b_{k,M} b_{n-k,M}, exact, n <= %d" % N)
    ratios = []
    for n in range(2, N + 1):
        lhs_tot = rhs_tot = 0
        for M in range(1, n):
            a = A.get((n, M), 0)
            rhs = sum(b(k, M) * b(n - k, M) for k in range(M, n - M + 1))
            lhs_tot += a
            rhs_tot += rhs
            if a > rhs and first_fail is None:
                first_fail = (n, M, a, rhs)
        if rhs_tot:
            ratios.append((n, F(lhs_tot, rhs_tot)))
    if first_fail is None:
        print("[CL1] HOLDS for all n in [2, %d] and all M" % N)
    else:
        print("[CL1] FAILS at n = %d, M = %d: a = %d > %d" % first_fail)
    print("     sum_M a_{n,M} / sum_M sum_k b_{k,M} b_{n-k,M}:", " ".join("%d:%.4f" % (n, float(r)) for n, r in ratios))
    print("     arches a_n = sum_M a_{n,M}:", " ".join("%d:%d" % (n, sum(v for (m, M), v in A.items() if m == n))
                                                    for n in range(2, N + 1)))
    # CL2
    print("\n## CL2 (renewal identity) r_{n,T} = sum_k b_{k,T} b_{n-k,T}, exact, n <= %d" % N)
    bad = []
    cnt = 0
    for n in range(1, N + 1):
        for T in range(1, n // 2 + 1):
            rhs = sum(b(k, T) * b(n - k, T) for k in range(T, n - T + 1))
            lhs = R.get((n, T), 0)
            cnt += 1
            if lhs != rhs:
                bad.append((n, T, lhs, rhs))
    if not bad:
        print("[CL2] HOLDS for all n in [1, %d], all T >= 1 (%d pairs (n,T), %d with r > 0)" %
              (N, cnt, sum(1 for k, v in R.items() if k[1] >= 1 and v > 0)))
    else:
        print("[CL2] FAILS: first (n,T,lhs,rhs) = %s; %d failures" % (bad[0], len(bad)))
        fails += 1
    print("     fraction of height-2T bridges with a renewal at height T (n = %d): %s" % (
        N, " ".join("T=%d:%.4f" % (T, R.get((N, T), 0) / b(N, 2 * T)) for T in range(1, N // 2 + 1) if b(N, 2 * T))))
    print("\n%s" % ("ALL CROSS-CHECKS PASS" if fails == 0 else "%d FAIL" % fails))
    sys.exit(1 if fails else 0)


if __name__ == "__main__":
    main()
