#!/usr/bin/env python3
"""lemma_testbed.py -- evaluate candidate finite-n statements on exact SAW data on Z^3.

Environment: Python 3.8+ standard library only (fractions, math).  Runtime < 1 min.
Usage: python3 lemma_testbed.py [data file, default ../data/saw_enum_N21.txt.gz]

Inputs: the exact enumeration (src/saw_enum.cpp output: endpoint distributions, bridges, half-space walks,
intermediate-point distributions, n <= N) and the published series c_n, P_n = sum |omega_n|^2 for n <= 36
(data/reference_series.tsv; Schram-Barkema-Bisseling 2011 Table I, agreeing with our enumeration for n <= N).

Every verdict of the form "HOLDS for n <= N" / "FAILS at n = ..." is an exact statement (integer / rational
arithmetic).  Lines tagged NUMERICS are floating-point summaries (local exponents etc.) and prove nothing.
Exit code 0 unless a sanity check of the data fails.
"""
import math
import os
import sys
from fractions import Fraction as F

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from sawdata import SawData, read_reference  # noqa: E402

ROOT = os.path.dirname(HERE)


def verdict(name, first_fail, nmax, detail=""):
    if first_fail is None:
        print("[%s] HOLDS for all tested n <= %d%s" % (name, nmax, ("; " + detail) if detail else ""))
    else:
        print("[%s] FAILS at n = %s (exact)%s" % (name, first_fail, ("; " + detail) if detail else ""))


def main():
    path = sys.argv[1] if len(sys.argv) > 1 else os.path.join(ROOT, "data", "saw_enum_N21.txt.gz")
    D = SawData(path)
    N = D.N
    ref = {int(r[0]): r for r in read_reference(os.path.join(ROOT, "data", "reference_series.tsv"))}
    R = max(ref)
    c = {n: int(ref[n][1]) for n in ref}
    P = {n: int(ref[n][2]) for n in ref}
    # sanity: own data agrees with the reference
    for n in range(N + 1):
        assert D.c[n] == c[n] and D.S(n) == P[n], n
    E = {n: F(P[n], c[n]) for n in range(R + 1)}  # E_n |omega_n|^2, exact, n <= 36
    print("# lemma testbed: own enumeration n <= %d (%s); published c_n, P_n for n <= %d" % (N, os.path.basename(path), R))
    print()
    print("## A. Mean-square end-to-end distance (exact, n <= %d from published series)" % R)

    # A1  E_n >= n
    ff = next((n for n in range(1, R + 1) if E[n] < n), None)
    mins = min((E[n] / n, n) for n in range(1, R + 1))
    verdict("A1 E_n|w_n|^2 >= n", ff, R, "min_n E_n/n = %.6f at n = %d" % (float(mins[0]), mins[1]))
    # A1b E_n >= n + something: report E_n - n
    # A2 E_n/n nondecreasing
    ff = next((n for n in range(1, R) if E[n + 1] / (n + 1) < E[n] / n), None)
    verdict("A2 E_n/n nondecreasing in n", ff, R)
    # A2b E_n nondecreasing
    ff = next((n for n in range(0, R) if E[n + 1] < E[n]), None)
    verdict("A2b E_n nondecreasing in n", ff, R)
    # A2c E_n / n^2 nonincreasing
    ff = next((n for n in range(1, R) if E[n + 1] / (n + 1) ** 2 > E[n] / n ** 2), None)
    verdict("A2c E_n/n^2 nonincreasing", ff, R)
    # A2d E_{n+1} - E_n nondecreasing (convexity of n -> E_n)
    ff = next((n for n in range(1, R) if E[n + 1] - E[n] < E[n] - E[n - 1]), None)
    verdict("A2d n -> E_n convex (E_{n+1}-E_n nondecreasing)", ff, R)
    # A2e increments E_{n+1}-E_n >= 1
    ff = next((n for n in range(0, R) if E[n + 1] - E[n] < 1), None)
    verdict("A2e E_{n+1} - E_n >= 1", ff, R)
    # A3 superadditivity
    ff = None
    for s in range(2, R + 1):
        for m in range(1, s // 2 + 1):
            if E[s] < E[m] + E[s - m]:
                ff = "%d (= %d + %d)" % (s, m, s - m)
                break
        if ff:
            break
    verdict("A3 superadditivity E_{n+m} >= E_n + E_m", ff, R)
    # A3b strict superadditive margin: min over n+m<=R of E_{n+m} - E_n - E_m
    mg = min((E[a + b] - E[a] - E[b], a, b) for a in range(1, R) for b in range(a, R - a + 1))
    print("     min_{n<=m, n+m<=%d} (E_{n+m} - E_n - E_m) = %.6f at (n,m) = (%d,%d)" % (R, float(mg[0]), mg[1], mg[2]))
    # A3c E_{2n} >= 2 E_n  and stronger E_{2n} >= 2^{2*0.5} ... report ratio E_{2n}/E_n
    rat = ["%d:%.4f" % (n, float(E[2 * n] / E[n])) for n in range(1, R // 2 + 1)]
    print("     NUMERICS E_{2n}/E_n (= 2^{2 nu_eff}):", " ".join(rat))
    # A4 Madras 2014 bound
    ff = next((n for n in range(1, R + 1) if E[n] < n ** (4.0 / 9) / 36), None)
    verdict("A4 Madras 2014 Prop 1.1 (p=2,d=3): E_n >= n^{4/9}/36 (float comparison, margin huge)", ff, R)
    # A5 log-concavity of c_n (c_{n+1}/c_n nonincreasing) and of P_n
    ff = next((n for n in range(1, R) if c[n] * c[n] < c[n - 1] * c[n + 1]), None)
    verdict("A5 c_n^2 >= c_{n-1} c_{n+1} (c_{n+1}/c_n nonincreasing)", ff, R)
    ff = next((n for n in range(2, R) if P[n] * P[n] < P[n - 1] * P[n + 1]), None)
    verdict("A5b P_n^2 >= P_{n-1} P_{n+1}", ff, R)
    # A5c: c_{n+2}/c_n nonincreasing
    ff = next((n for n in range(1, R - 2) if F(c[n + 3], c[n + 1]) > F(c[n + 2], c[n])), None)
    verdict("A5c c_{n+2}/c_n nonincreasing", ff, R)
    # A6 local exponent (numerics)
    print("     NUMERICS nu_eff(n) = log(E_{n+1}/E_{n-1}) / (2 log((n+1)/(n-1))):")
    line = []
    for n in range(2, R):
        v = math.log(E[n + 1] / E[n - 1]) / (2 * math.log((n + 1) / (n - 1)))
        line.append("%d:%.4f" % (n, v))
    print("       " + " ".join(line))

    print()
    print("## B. Endpoint distribution (exact, n <= %d, own enumeration)" % N)
    # B1 argmax and max probability
    rows = []
    for n in range(1, N + 1):
        best, arg = D.max_point(n)
        rows.append((n, best, arg, F(best, D.c[n])))
    print("     n  argmax-class  max_x P_n(w_n=x)")
    for n, best, arg, p in rows:
        print("     %2d  %-14s %.6e" % (n, ",".join("(%d,%d,%d)" % a for a in arg), float(p)))
    # B1a: max attained at nearest neighbour for odd n
    ff = next((n for n, b, arg, p in rows if n % 2 == 1 and arg != [(1, 0, 0)]), None)
    verdict("B1a odd n: c_n(x) maximised (uniquely, up to symmetry) at x = e1", ff, N)
    ff = next((n for n, b, arg, p in rows if n % 2 == 0 and n >= 2 and arg != [(1, 1, 0)]), None)
    verdict("B1b even n >= 2: c_n(x) maximised (uniquely, up to symmetry) at x = e1+e2", ff, N)
    # decay: max_x P_n is nonincreasing in n within each parity
    pm = {n: p for n, b, a, p in rows}
    ff = next((n for n in range(1, N - 1) if pm[n + 2] > pm[n]), None)
    verdict("B1c max_x P_n(w_n=x) nonincreasing in n within each parity class", ff, N)
    print("     NUMERICS local exponent of max_x P_n (same parity): log(p_{n-2}/p_n)/log(n/(n-2)):")
    print("       " + " ".join("%d:%.3f" % (n, math.log(pm[n - 2] / pm[n]) / math.log(n / (n - 2)))
                               for n in range(4, N + 1)))
    # B2 P_n(w_n = x) for fixed small x, nonincreasing in n (same parity)
    pts = [(1, 0, 0), (2, 0, 0), (1, 1, 0), (3, 0, 0), (2, 1, 0), (1, 1, 1), (4, 0, 0), (2, 2, 0), (3, 1, 0),
           (2, 1, 1), (5, 0, 0), (3, 2, 0), (1, 1, 1)]
    for x in pts[:10]:
        par = sum(x) % 2
        ns = [n for n in range(1, N + 1) if n % 2 == par and n >= sum(x)]
        pr = {n: F(D.cnx(n, x), D.c[n]) for n in ns}
        ff = next((n for n in ns if n + 2 in pr and pr[n] > 0 and pr[n + 2] > pr[n]), None)
        nz = [n for n in ns if pr[n] > 0]
        loc = ""
        if len(nz) >= 3:
            n1, n0 = nz[-1], nz[-2]
            loc = "NUMERICS local exponent at n=%d: %.3f" % (n1, math.log(pr[n0] / pr[n1]) / math.log(n1 / n0))
        # nonincreasing only from the first n where the point is reached with positive probability + 2 (rise first)
        first_dec = None
        ups = [n for n in nz if n + 2 in pr and pr[n + 2] > pr[n]]
        print("[B2 x=%s] P_n(w_n=x) over n=%d..%d (step 2): increases at n -> n+2 for n in %s; %s" %
              (x, nz[0], nz[-1], ups if ups else "none", loc))
    # P_n(|w_n| = 1) exact
    print("     P_n(|w_n|=1) = 6 c_n(e1)/c_n (odd n):", " ".join("%d:%.5e" % (n, 6 * D.cnx(n, (1, 0, 0)) / D.c[n])
                                                         for n in range(1, N + 1, 2)))
    # B3 monotonicity along axes / diagonals / coordinatewise
    fa = fd = fd3 = fco = fco2 = None
    for n in range(1, N + 1):
        for k in range(0, n + 1):
            if (k - n) % 2 == 0 and k + 2 <= n and fa is None:
                if D.cnx(n, (k + 2, 0, 0)) > D.cnx(n, (k, 0, 0)) and not (k == 0):
                    fa = "%d (k=%d)" % (n, k)
        for k in range(0, n):
            if n % 2 == 0 and 2 * (k + 1) <= n and k >= 1 and fd is None:
                if D.cnx(n, (k + 1, k + 1, 0)) > D.cnx(n, (k, k, 0)):
                    fd = "%d (k=%d)" % (n, k)
            if (3 * k - n) % 2 == 0 and 3 * (k + 2) <= n and k >= 1 and fd3 is None:
                if D.cnx(n, (k + 2, k + 2, k + 2)) > D.cnx(n, (k, k, k)):
                    fd3 = "%d (k=%d)" % (n, k)
        # coordinatewise: c_n(x) >= c_n(x + 2 e_i) for x_i >= 0, x != 0 (x=0 excluded: c_n(0)=0)
        for (a, b, cc) in list(D.A[n].keys()) + [(0, 0, 0)]:
            x = (a, b, cc)
            if x == (0, 0, 0):
                continue
            for i in range(3):
                y = list(x)
                y[i] += 2
                if D.cnx(n, y) > D.cnx(n, x) and fco is None:
                    fco = "%d (x=%s -> %s)" % (n, x, tuple(y))
                # move +1 in two coordinates (parity preserving): x -> x + e_i + e_j, i != j
                for j in range(i + 1, 3):
                    z = list(x)
                    z[i] += 1
                    z[j] += 1
                    if D.cnx(n, z) > D.cnx(n, x) and fco2 is None:
                        fco2 = "%d (x=%s -> %s)" % (n, x, tuple(z))
    verdict("B3a axis monotonicity: c_n(k e1) >= c_n((k+2) e1), k >= 1", fa, N)
    verdict("B3b diagonal: c_n(k,k,0) >= c_n(k+1,k+1,0), k >= 1", fd, N)
    verdict("B3c body diagonal: c_n(k,k,k) >= c_n(k+2,k+2,k+2), k >= 1", fd3, N)
    verdict("B3d coordinatewise: c_n(x) >= c_n(x + 2 e_i) for x in the closed positive octant, x != 0", fco, N)
    verdict("B3e c_n(x) >= c_n(x + e_i + e_j) (i != j) for x in the closed positive octant, x != 0", fco2, N)
    # B3f: |x| larger => fewer walks (Euclidean norm), same parity
    ff = None
    for n in range(1, N + 1):
        pts_n = sorted(((a * a + b * b + cc * cc), Ns // orb, (a, b, cc)) for (a, b, cc), (Ns, orb) in D.A[n].items())
        for i in range(len(pts_n)):
            for j in range(len(pts_n)):
                if pts_n[i][0] < pts_n[j][0] and pts_n[i][1] < pts_n[j][1]:
                    ff = "%d (|x|^2=%d %s has %d < %d for |y|^2=%d %s)" % (
                        n, pts_n[i][0], pts_n[i][2], pts_n[i][1], pts_n[j][1], pts_n[j][0], pts_n[j][2])
                    break
            if ff:
                break
        if ff:
            break
    verdict("B3f radial monotonicity: |x| < |y| (same parity) => c_n(x) >= c_n(y)", ff, N)
    # B3g: l1-radial: |x|_1 < |y|_1 => c_n(x) >= c_n(y)
    ff = None
    for n in range(1, N + 1):
        L = sorted(((a + b + cc), Ns // orb, (a, b, cc)) for (a, b, cc), (Ns, orb) in D.A[n].items())
        for i in range(len(L)):
            for j in range(len(L)):
                if L[i][0] < L[j][0] and L[i][1] < L[j][1]:
                    ff = "%d (%s: %d < %s: %d)" % (n, L[i][2], L[i][1], L[j][2], L[j][1])
                    break
            if ff:
                break
        if ff:
            break
    verdict("B3g l1-radial monotonicity: |x|_1 < |y|_1 => c_n(x) >= c_n(y)", ff, N)

    # B4 hyperplane
    print()
    print("## C. Hyperplane x_1 = 0, bridges, half-space walks")
    q = {n: F(D.ending_on_plane(n), D.c[n]) for n in range(N + 1)}
    print("     P_n(x_1(w_n)=0) exact values:", " ".join("%d:%.5f" % (n, float(q[n])) for n in range(1, N + 1)))
    ff = next((n for n in range(1, N - 1) if q[n + 2] > q[n]), None)
    verdict("C1 P_n(x_1(w_n) = 0) nonincreasing in n within parity", ff, N)
    print("     NUMERICS local exponent of P_n(x_1=0) (same parity):",
          " ".join("%d:%.3f" % (n, math.log(q[n - 2] / q[n]) / math.log(n / (n - 2))) for n in range(4, N + 1)))
    print("     NUMERICS sqrt(E_n) * P_n(x_1=0):", " ".join("%d:%.4f" % (n, math.sqrt(E[n]) * float(q[n]))
                                                         for n in range(1, N + 1)))
    bc = {n: F(D.b[n], D.c[n]) for n in range(N + 1)}
    print("     b_n/c_n:", " ".join("%d:%.5f" % (n, float(bc[n])) for n in range(1, N + 1)))
    ff = next((n for n in range(1, N) if bc[n + 1] > bc[n]), None)
    verdict("C2 b_n/c_n nonincreasing", ff, N)
    print("     NUMERICS local exponent of b_n/c_n: log(r_{n-1}/r_{n+1})/log((n+1)/(n-1)):",
          " ".join("%d:%.3f" % (n, math.log(bc[n - 1] / bc[n + 1]) / math.log((n + 1) / (n - 1))) for n in range(2, N)))
    # supermultiplicativity of bridges (proved, sanity)
    ff = next(("%d+%d" % (m, k) for m in range(1, N) for k in range(1, N - m + 1) if D.b[m] * D.b[k] > D.b[m + k]), None)
    verdict("C3 b_m b_n <= b_{m+n} (PROVED in the literature; sanity check)", ff, N)
    hc = {n: F(D.h[n], D.c[n]) for n in range(N + 1)}
    print("     h_n/c_n:", " ".join("%d:%.5f" % (n, float(hc[n])) for n in range(1, N + 1)))
    Eb = {n: F(D.S_bridge(n), D.b[n]) for n in range(1, N + 1)}
    ff = next((n for n in range(1, N + 1) if Eb[n] < E[n]), None)
    verdict("C4 E[|w_n|^2 | bridge] >= E_n|w_n|^2", ff, N,
            "ratio at n=%d: %.4f" % (N, float(Eb[N] / E[N])))
    Eh = {n: F(D.S_half(n), D.h[n]) for n in range(1, N + 1)}
    ff = next((n for n in range(1, N + 1) if Eh[n] < E[n]), None)
    verdict("C5 E[|w_n|^2 | half-space] >= E_n|w_n|^2", ff, N, "ratio at n=%d: %.4f" % (N, float(Eh[N] / E[N])))
    ff = next((n for n in range(1, N + 1) if Eb[n] < Eh[n]), None)
    verdict("C6 E[|w_n|^2 | bridge] >= E[|w_n|^2 | half-space]", ff, N)
    # mean span of bridges
    sp = {n: F(sum(k * v for k, v in D.bridge_span_dist(n).items()), D.b[n]) for n in range(1, N + 1)}
    ff = next((n for n in range(1, N + 1) if sp[n] ** 2 * 3 < E[n]), None)
    verdict("C7 3 (E[span | bridge])^2 >= E_n|w_n|^2", ff, N)
    ff = next((n for n in range(1, N) if sp[n + 1] < sp[n]), None)
    verdict("C8 mean bridge span nondecreasing in n", ff, N)

    print()
    print("## D. Intermediate points (exact, n <= %d)" % N)
    Ek = {(k, n): F(D.S_mid(k, n), D.c[n]) for n in range(1, N + 1) for k in range(0, n + 1)}
    # D1 E_n|w_k|^2 >= E_k|w_k|^2
    ff = None
    worst = None
    for n in range(2, N + 1):
        for k in range(1, n):
            r = Ek[(k, n)] / E[k]
            if worst is None or r < worst[0]:
                worst = (r, k, n)
            if Ek[(k, n)] < E[k] and ff is None:
                ff = "%d (k=%d)" % (n, k)
    verdict("D1 E_n|w_k|^2 >= E_k|w_k|^2 for 1 <= k < n (initial segments of longer walks more extended)", ff, N,
            "min ratio %.6f at (k,n)=(%d,%d)" % (float(worst[0]), worst[1], worst[2]))
    # D2 monotone in n for fixed k
    ff = next(("%d (k=%d)" % (n, k) for n in range(2, N) for k in range(1, n + 1) if Ek[(k, n + 1)] < Ek[(k, n)]), None)
    verdict("D2 n -> E_n|w_k|^2 nondecreasing for fixed k (n >= k)", ff, N)
    # D3 E_n|w_k|^2 >= k
    ff = next(("%d (k=%d)" % (n, k) for n in range(1, N + 1) for k in range(1, n + 1) if Ek[(k, n)] < k), None)
    verdict("D3 E_n|w_k|^2 >= k for all k <= n", ff, N)
    # D4 k -> E_n|w_k|^2 nondecreasing
    ff = next(("%d (k=%d)" % (n, k) for n in range(1, N + 1) for k in range(0, n) if Ek[(k + 1, n)] < Ek[(k, n)]), None)
    verdict("D4 k -> E_n|w_k|^2 nondecreasing in k for fixed n", ff, N)
    # D5 E_n|w_n|^2 >= E_n|w_k|^2 + E_n|w_n - w_k|^2 ?  (needs w_n - w_k: not available in general) -- the special
    # case k = n/2 via reversal: E_n|w_n - w_{n-k}|^2 = E_n|w_k|^2.  So test E_n|w_n|^2 >= E_n|w_k|^2 + E_n|w_{n-k}|^2.
    ff = next(("%d (k=%d)" % (n, k) for n in range(2, N + 1) for k in range(1, n) if E[n] < Ek[(k, n)] + Ek[(n - k, n)]),
              None)
    verdict("D5 E_n|w_n|^2 >= E_n|w_k|^2 + E_n|w_{n-k}|^2 for 1 <= k < n",
            ff, N)
    print("     note D5: by reversal E_n|w_n - w_k|^2 = E_n|w_{n-k}|^2, so D5 <=> E_n[w_k . (w_n - w_k)] >= 0 for all k.")
    # D6 midpoint sup probability (DGHM Prop 1.3: <= C n^{-1/2})
    mp = {n: F(D.max_point_mid(n // 2, n), D.c[n]) for n in range(2, N + 1)}
    print("     sup_x P_n(w_{floor(n/2)} = x) * n^{1/2}:", " ".join("%d:%.4f" % (n, float(mp[n]) * math.sqrt(n))
                                                             for n in range(2, N + 1)))
    print("     NUMERICS sup_x P_n(w_{n/2}=x) * n^{3 nu}, nu=0.5876:",
          " ".join("%d:%.4f" % (n, float(mp[n]) * n ** (3 * 0.587597)) for n in range(2, N + 1, 2)))
    print()
    print("## E. NUMERICS: amplitudes")
    print("     E_n/n^{2 nu} (nu = 0.587597):", " ".join("%d:%.4f" % (n, float(E[n]) / n ** (2 * 0.587597))
                                                  for n in range(1, R + 1)))


if __name__ == "__main__":
    main()
