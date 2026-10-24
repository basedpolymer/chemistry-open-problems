#!/usr/bin/env python3
"""lemma_testbed.py -- evaluate candidate finite-n statements on exact SAW data on Z^3.

Environment: Python 3.8+ standard library only (fractions, math).  Runtime < 1 min.
Usage: python3 lemma_testbed.py [data file, default ../data/saw_enum_N21.txt.gz] [gram file, optional]

Inputs: the exact enumeration (src/saw_enum.cpp output: endpoint distributions, bridges, half-space walks,
intermediate-point distributions, n <= N), the published series c_n, P_n = sum |omega_n|^2 for n <= 36
(data/reference_series.tsv; Schram-Barkema-Bisseling 2011 Table I, agreeing with our enumeration for n <= N),
and optionally the step-correlation (Gram) data of src/saw_gram.cpp.

Verdicts.  For a statement S(n) tested on a range of n:
  HOLDS for all n in [a, b]                         -- exact, every tested n
  FAILS at n = n1 (exact) ...; failing n: [...]     -- exact counterexamples, with the full list of failing n,
                                                       so that "fails only for small n" is visible
Lines tagged NUMERICS are floating-point summaries (local exponents etc.): they prove nothing.
"""
import math
import os
import sys
from fractions import Fraction as F

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from sawdata import SawData, read_reference  # noqa: E402

ROOT = os.path.dirname(HERE)


def scan(name, ns, pred):
    """pred(n) -> None if S(n) holds, else a short string describing the counterexample at n."""
    ns = list(ns)
    fails = [(n, pred(n)) for n in ns]
    fails = [(n, d) for n, d in fails if d is not None]
    if not fails:
        print("[%s] HOLDS for all n in [%d, %d]" % (name, ns[0], ns[-1]))
        return
    fl = [n for n, _ in fails]
    tail = ""
    if fl[-1] < ns[-1]:
        tail = "; HOLDS for %d < n <= %d" % (fl[-1], ns[-1])
    print("[%s] FAILS at n = %d (exact): %s; failing n (%d of %d tested): %s%s" %
          (name, fl[0], fails[0][1], len(fl), len(ns), compress(fl), tail))


def compress(lst):
    out, i = [], 0
    while i < len(lst):
        j = i
        step = 1
        if j + 1 < len(lst):
            step = lst[j + 1] - lst[j]
        while j + 1 < len(lst) and lst[j + 1] - lst[j] == step and step in (1, 2):
            j += 1
        if j - i >= 2:
            out.append("%d..%d%s" % (lst[i], lst[j], "" if step == 1 else " step 2"))
            i = j + 1
        else:
            out.append(str(lst[i]))
            i += 1
    return "[" + ", ".join(out) + "]"


def fl(x):
    return float(x)


def main():
    args = sys.argv[1:]
    path = args[0] if args else os.path.join(ROOT, "data", "saw_enum_N21.txt.gz")
    gram_path = args[1] if len(args) > 1 else None
    D = SawData(path)
    N = D.N
    ref = {int(r[0]): r for r in read_reference(os.path.join(ROOT, "data", "reference_series.tsv"))}
    R = max(ref)
    c = {n: int(ref[n][1]) for n in ref}
    P = {n: int(ref[n][2]) for n in ref}
    for n in range(N + 1):  # sanity: own data agrees with the published series
        assert D.c[n] == c[n] and D.S(n) == P[n], n
    E = {n: F(P[n], c[n]) for n in range(R + 1)}  # E_n |omega_n|^2, exact
    print("# lemma testbed: own enumeration n <= %d (%s); published c_n, P_n for n <= %d" % (N, os.path.basename(path), R))
    print()
    print("## A. Mean-square end-to-end distance E_n = E_n|w_n|^2 (exact; n <= %d from the published series)" % R)
    rng = range(1, R + 1)
    scan("A1 E_n >= n", rng, lambda n: None if E[n] >= n else "E_n/n = %.5f" % fl(E[n] / n))
    print("     E_n/n:", " ".join("%d:%.4f" % (n, fl(E[n] / n)) for n in rng))
    scan("A2 E_{n+1}/(n+1) >= E_n/n", range(1, R), lambda n: None if E[n + 1] / (n + 1) >= E[n] / n else "")
    scan("A2b E_{n+1} >= E_n", range(0, R), lambda n: None if E[n + 1] >= E[n] else "")
    scan("A2c E_{n+1}/(n+1)^2 <= E_n/n^2", range(1, R), lambda n: None if E[n + 1] / (n + 1) ** 2 <= E[n] / n ** 2 else "")
    scan("A2d convexity E_{n+1}-E_n >= E_n-E_{n-1}", range(1, R),
         lambda n: None if E[n + 1] - E[n] >= E[n] - E[n - 1] else "%.4f < %.4f" % (fl(E[n + 1] - E[n]), fl(E[n] - E[n - 1])))
    scan("A2e convexity within parity E_{n+2}-E_n >= E_n-E_{n-2}", range(2, R - 1),
         lambda n: None if E[n + 2] - E[n] >= E[n] - E[n - 2] else "")
    scan("A2f E_{n+1} - E_n >= 1", range(0, R), lambda n: None if E[n + 1] - E[n] >= 1 else "")
    scan("A2g E_{n+2} - E_n >= 2 E_n/n (i.e. E_n/n increasing at step 2 by...)", range(1, R - 1),
         lambda n: None if E[n + 2] - E[n] >= 2 * E[n] / n else "")

    def sup(s):
        for m in range(1, s // 2 + 1):
            if E[s] < E[m] + E[s - m]:
                return "E_%d < E_%d + E_%d" % (s, m, s - m)
        return None
    scan("A3 superadditivity E_{m+k} >= E_m + E_k (all m+k = n)", range(2, R + 1), sup)
    mg = min((E[a + b] - E[a] - E[b], a, b) for a in range(1, R) for b in range(a, R - a + 1))
    print("     min over m <= k, m+k <= %d of E_{m+k} - E_m - E_k = %s = %.5f at (m,k) = (%d,%d)" %
          (R, mg[0], fl(mg[0]), mg[1], mg[2]))

    def supmult(s):  # E_{m+k}/(m+k)^{2a} type: test E_{2n} >= 2^{2*0.5} E_n trivially implied by A3; test 2^{1.1}
        return None
    for alpha in (1.0, 1.1, 1.15, 1.17, 1.175):
        scan("A3b E_{2n} >= 2^{%.3f} E_n (float)" % alpha, range(1, R // 2 + 1),
             lambda n, a=alpha: None if fl(E[2 * n]) >= 2 ** a * fl(E[n]) else "ratio %.4f" % fl(E[2 * n] / E[n]))
    print("     NUMERICS E_{2n}/E_n:", " ".join("%d:%.4f" % (n, fl(E[2 * n] / E[n])) for n in range(1, R // 2 + 1)))
    scan("A4 Madras 2014 Prop 1.1 (p=2, d=3): E_n >= n^{4/9}/36 (float; margin huge)", rng,
         lambda n: None if fl(E[n]) >= n ** (4.0 / 9) / 36 else "")
    scan("A5 log-concavity c_n^2 >= c_{n-1} c_{n+1}", range(1, R),
         lambda n: None if c[n] ** 2 >= c[n - 1] * c[n + 1] else "")
    scan("A5b c_{n+2}/c_n nonincreasing", range(1, R - 2),
         lambda n: None if F(c[n + 3], c[n + 1]) <= F(c[n + 2], c[n]) else "")
    scan("A5c P_n^2 >= P_{n-1} P_{n+1}", range(2, R), lambda n: None if P[n] ** 2 >= P[n - 1] * P[n + 1] else "")
    scan("A5d c_{n+1}/c_n >= mu_lower = 4.572 (float)", range(0, R),
         lambda n: None if c[n + 1] / c[n] >= 4.572 else "%.4f" % (c[n + 1] / c[n]))
    print("     NUMERICS nu_eff(n) = log(E_{n+1}/E_{n-1}) / (2 log((n+1)/(n-1))):")
    print("       " + " ".join("%d:%.4f" % (n, math.log(E[n + 1] / E[n - 1]) / (2 * math.log((n + 1) / (n - 1))))
                               for n in range(2, R)))

    print()
    print("## B. Endpoint distribution c_n(x) (exact, n <= %d, own enumeration)" % N)
    rows = {}
    for n in range(1, N + 1):
        best, arg = D.max_point(n)
        rows[n] = (best, arg, F(best, D.c[n]))
    print("     n  argmax class(es)   max_x P_n(w_n = x)   |argmax|/sqrt(E_n)")
    for n in range(1, N + 1):
        best, arg, p = rows[n]
        r = math.sqrt(sum(v * v for v in arg[0]) / fl(E[n]))
        print("     %2d  %-18s %.6e         %.3f" % (n, ",".join("(%d,%d,%d)" % a for a in arg), fl(p), r))
    scan("B1a odd n: max_x c_n(x) attained at x = e1 only", range(1, N + 1, 2),
         lambda n: None if rows[n][1] == [(1, 0, 0)] else "argmax %s" % rows[n][1])
    scan("B1b even n: max_x c_n(x) attained at |x|_1 = 2 only", range(2, N + 1, 2),
         lambda n: None if all(sum(a) == 2 for a in rows[n][1]) else "argmax %s" % rows[n][1])
    scan("B1c max_x P_n(w_n = x) >= max_x P_{n+2}(w_{n+2} = x)", range(1, N - 1),
         lambda n: None if rows[n + 2][2] <= rows[n][2] else "")
    scan("B1d max_x P_n(w_n = x) >= max_x P_{n+1}(w_{n+1} = x)", range(1, N),
         lambda n: None if rows[n + 1][2] <= rows[n][2] else "")
    pm = {n: rows[n][2] for n in rows}
    print("     NUMERICS local exponent of max_x P_n (same parity) log(p_{n-2}/p_n)/log(n/(n-2)):")
    print("       " + " ".join("%d:%.3f" % (n, math.log(pm[n - 2] / pm[n]) / math.log(n / (n - 2))) for n in range(4, N + 1)))
    print("     NUMERICS n^{3 nu} max_x P_n (nu = 0.5876):", " ".join("%d:%.4f" % (n, fl(pm[n]) * n ** (3 * 0.587597))
                                                                 for n in range(1, N + 1)))
    # fixed points
    pts = [(1, 0, 0), (2, 0, 0), (1, 1, 0), (3, 0, 0), (2, 1, 0), (1, 1, 1), (4, 0, 0), (2, 2, 0), (3, 1, 0), (2, 1, 1)]
    for x in pts:
        par = sum(x) % 2
        ns = [n for n in range(sum(x), N + 1) if n % 2 == par]
        pr = {n: F(D.cnx(n, x), D.c[n]) for n in ns}
        n0 = next(n for n in ns if pr[n] > 0)
        scan("B2 P_n(w_n=%s) >= P_{n+2}(w_{n+2}=%s)" % (str(x).replace(" ", ""), str(x).replace(" ", "")),
             [n for n in ns if n >= n0 and n + 2 <= N], lambda n: None if pr[n + 2] <= pr[n] else "")
        last = [n for n in ns if pr[n] > 0][-2:]
        print("     NUMERICS local exponent at n=%d: %.3f;  n^{1+gamma-alpha}-scaling check n^1.92 P_n: %s" %
              (last[1], math.log(pr[last[0]] / pr[last[1]]) / math.log(last[1] / last[0]),
               " ".join("%d:%.4f" % (n, fl(pr[n]) * n ** 1.92) for n in ns[-4:])))
    # line bound (lower-bound attempts.md A1): P_n(w_n on the x1-axis), and on a ray {t e1 : t >= 1}
    line = {n: F(sum(D.cnx(n, (t, 0, 0)) for t in range(-n, n + 1)), D.c[n]) for n in range(1, N + 1)}
    ray = {n: F(sum(D.cnx(n, (t, 0, 0)) for t in range(1, n + 1)), D.c[n]) for n in range(1, N + 1)}
    print("     P_n(w_n in the x1-axis line):", " ".join("%d:%.4e" % (n, fl(line[n])) for n in range(1, N + 1)))
    print("     NUMERICS local exponent of P_n(w_n in line) (same parity):",
          " ".join("%d:%.3f" % (n, math.log(line[n - 2] / line[n]) / math.log(n / (n - 2))) for n in range(4, N + 1)))
    scan("B2b P_n(w_n in x1-axis) >= P_{n+2}(w_{n+2} in x1-axis)", range(1, N - 1),
         lambda n: None if line[n + 2] <= line[n] else "")
    print("     (ray {t e1, t >= 1}: P = line/2 exactly by symmetry; lower-bound A1 needs a decay exponent b > 2/3)")
    assert all(2 * ray[n] == line[n] for n in ray)
    print("     P_n(|w_n| = 1) = 6 c_n(e1)/c_n (odd n):", " ".join("%d:%.5e" % (n, 6 * D.cnx(n, (1, 0, 0)) / D.c[n])
                                                               for n in range(1, N + 1, 2)))

    def axis(n):
        for k in range(1, n - 1):
            if (k - n) % 2 == 0 and D.cnx(n, (k + 2, 0, 0)) > D.cnx(n, (k, 0, 0)):
                return "c_%d(%d e1) = %d < c_%d(%d e1) = %d" % (n, k, D.cnx(n, (k, 0, 0)), n, k + 2, D.cnx(n, (k + 2, 0, 0)))
        return None
    scan("B3a axis: c_n(k e1) >= c_n((k+2) e1) for all k >= 1", range(1, N + 1), axis)

    def axis_far(n):  # monotone beyond the maximum?
        vals = [(k, D.cnx(n, (k, 0, 0))) for k in range(n % 2 if n % 2 else 2, n + 1, 2)]
        kmax = max(vals, key=lambda t: t[1])[0]
        for (k1, v1), (k2, v2) in zip(vals, vals[1:]):
            if k1 >= kmax and v2 > v1:
                return "k=%d" % k1
        return None
    scan("B3a' axis, unimodal: k -> c_n(k e1) nonincreasing beyond its maximum", range(1, N + 1), axis_far)

    def diag(n):
        if n % 2:
            return None
        for k in range(1, n // 2):
            if D.cnx(n, (k + 1, k + 1, 0)) > D.cnx(n, (k, k, 0)):
                return "k=%d" % k
        return None
    scan("B3b plane diagonal: c_n(k,k,0) >= c_n(k+1,k+1,0), k >= 1 (even n)", range(2, N + 1, 2), diag)

    def diag3(n):
        for k in range(1, n):
            if 3 * (k + 2) <= n and (3 * k - n) % 2 == 0 and D.cnx(n, (k + 2,) * 3) > D.cnx(n, (k,) * 3):
                return "k=%d" % k
        return None
    scan("B3c body diagonal: c_n(k,k,k) >= c_n(k+2,k+2,k+2), k >= 1", range(1, N + 1), diag3)

    def coord(n, far=False):
        for (a, b, cc) in D.A[n]:
            x = (a, b, cc)
            for i in range(3):
                y = list(x)
                y[i] += 2
                if D.cnx(n, y) > D.cnx(n, x):
                    if far and x[i] * x[i] * 3 < fl(E[n]):
                        continue
                    return "c(%s)=%d < c(%s)=%d" % (x, D.cnx(n, x), tuple(y), D.cnx(n, y))
        return None
    scan("B3d coordinatewise: c_n(x) >= c_n(x + 2 e_i) for x in the closed positive octant, x != 0",
         range(1, N + 1), coord)
    scan("B3d' same, restricted to x_i >= sqrt(E_n/3)", range(1, N + 1), lambda n: coord(n, True))

    def radial(n, l1=False):
        L = sorted((((a + b + cc) if l1 else (a * a + b * b + cc * cc)), Ns // orb, (a, b, cc))
                   for (a, b, cc), (Ns, orb) in D.A[n].items())
        for i in range(len(L)):
            for j in range(len(L)):
                if L[i][0] < L[j][0] and L[i][1] < L[j][1]:
                    return "c(%s)=%d < c(%s)=%d" % (L[i][2], L[i][1], L[j][2], L[j][1])
        return None
    scan("B3f radial: |x| < |y| => c_n(x) >= c_n(y)", range(1, N + 1), radial)
    scan("B3g l1-radial: |x|_1 < |y|_1 => c_n(x) >= c_n(y)", range(1, N + 1), lambda n: radial(n, True))

    def radial_far(n):
        r0 = 1.0 * fl(E[n])
        L = sorted(((a * a + b * b + cc * cc), Ns // orb, (a, b, cc)) for (a, b, cc), (Ns, orb) in D.A[n].items())
        L = [t for t in L if t[0] >= r0]
        for i in range(len(L)):
            for j in range(len(L)):
                if L[i][0] < L[j][0] and L[i][1] < L[j][1]:
                    return "c(%s)=%d < c(%s)=%d" % (L[i][2], L[i][1], L[j][2], L[j][1])
        return None
    scan("B3h radial beyond the rms radius: E_n <= |x|^2 < |y|^2 => c_n(x) >= c_n(y)", range(1, N + 1), radial_far)

    print()
    print("## C. Hyperplane x_1 = 0, bridges, half-space walks (exact, n <= %d)" % N)
    q = {n: F(D.ending_on_plane(n), D.c[n]) for n in range(N + 1)}
    print("     P_n(x_1(w_n) = 0):", " ".join("%d:%.5f" % (n, fl(q[n])) for n in range(1, N + 1)))
    scan("C1 P_n(x_1(w_n)=0) >= P_{n+2}(x_1(w_{n+2})=0)", range(1, N - 1), lambda n: None if q[n + 2] <= q[n] else "")
    scan("C1b P_n(x_1(w_n)=0) >= P_{n+1}(x_1(w_{n+1})=0)", range(1, N), lambda n: None if q[n + 1] <= q[n] else "")
    print("     NUMERICS local exponent of P_n(x_1 = 0) (same parity):",
          " ".join("%d:%.3f" % (n, math.log(q[n - 2] / q[n]) / math.log(n / (n - 2))) for n in range(4, N + 1)))
    print("     NUMERICS E_n^{1/2} P_n(x_1 = 0):", " ".join("%d:%.4f" % (n, math.sqrt(E[n]) * fl(q[n])) for n in range(1, N + 1)))
    bc = {n: F(D.b[n], D.c[n]) for n in range(N + 1)}
    print("     b_n/c_n:", " ".join("%d:%.5f" % (n, fl(bc[n])) for n in range(1, N + 1)))
    scan("C2 b_n/c_n >= b_{n+1}/c_{n+1}", range(1, N), lambda n: None if bc[n + 1] <= bc[n] else "")
    print("     NUMERICS local exponent of b_n/c_n:",
          " ".join("%d:%.3f" % (n, math.log(bc[n - 1] / bc[n + 1]) / math.log((n + 1) / (n - 1))) for n in range(2, N)))
    scan("C3 b_m b_k <= b_{m+k} (PROVED, sanity)", range(2, N + 1),
         lambda n: next(("m=%d" % m for m in range(1, n) if D.b[m] * D.b[n - m] > D.b[n]), None))
    scan("C3b b_n^2 >= b_{n-1} b_{n+1} (log-concave bridges)", range(1, N),
         lambda n: None if D.b[n] ** 2 >= D.b[n - 1] * D.b[n + 1] else "")
    hc = {n: F(D.h[n], D.c[n]) for n in range(N + 1)}
    print("     h_n/c_n:", " ".join("%d:%.5f" % (n, fl(hc[n])) for n in range(1, N + 1)))
    Eb = {n: F(D.S_bridge(n), D.b[n]) for n in range(1, N + 1)}
    Eh = {n: F(D.S_half(n), D.h[n]) for n in range(1, N + 1)}
    scan("C4 E[|w_n|^2 | bridge] >= E_n", range(1, N + 1), lambda n: None if Eb[n] >= E[n] else "")
    scan("C5 E[|w_n|^2 | half-space] >= E_n", range(1, N + 1), lambda n: None if Eh[n] >= E[n] else "")
    scan("C6 E[|w_n|^2 | bridge] >= E[|w_n|^2 | half-space]", range(1, N + 1), lambda n: None if Eb[n] >= Eh[n] else "")
    print("     E[|w_n|^2|bridge]/E_n:", " ".join("%d:%.4f" % (n, fl(Eb[n] / E[n])) for n in range(1, N + 1)))
    sp = {n: F(sum(k * v for k, v in D.bridge_span_dist(n).items()), D.b[n]) for n in range(1, N + 1)}
    scan("C7 E[span | bridge] nondecreasing", range(1, N), lambda n: None if sp[n + 1] >= sp[n] else "")
    print("     E[span|bridge]/sqrt(E_n):", " ".join("%d:%.4f" % (n, fl(sp[n]) / math.sqrt(E[n])) for n in range(1, N + 1)))

    print()
    print("## D. Intermediate points (exact, n <= %d): E_n|w_k|^2 for k < n" % N)
    Ek = {(k, n): F(D.S_mid(k, n), D.c[n]) for n in range(1, N + 1) for k in range(0, n + 1)}
    scan("D1 E_n|w_k|^2 >= E_k|w_k|^2 for all 1 <= k < n", range(2, N + 1),
         lambda n: next(("k=%d" % k for k in range(1, n) if Ek[(k, n)] < E[k]), None))
    worst = min((Ek[(k, n)] / E[k], k, n) for n in range(2, N + 1) for k in range(1, n))
    print("     min ratio E_n|w_k|^2 / E_k|w_k|^2 = %.6f at (k,n) = (%d,%d)" % (fl(worst[0]), worst[1], worst[2]))
    scan("D2 E_{n+1}|w_k|^2 >= E_n|w_k|^2 for all k <= n", range(1, N),
         lambda n: next(("k=%d" % k for k in range(1, n + 1) if Ek[(k, n + 1)] < Ek[(k, n)]), None))
    scan("D2b E_{n+2}|w_k|^2 >= E_n|w_k|^2 for all k <= n", range(1, N - 1),
         lambda n: next(("k=%d" % k for k in range(1, n + 1) if Ek[(k, n + 2)] < Ek[(k, n)]), None))
    scan("D3 E_n|w_k|^2 >= k for all k <= n", range(1, N + 1),
         lambda n: next(("k=%d" % k for k in range(1, n + 1) if Ek[(k, n)] < k), None))
    scan("D4 E_n|w_{k+1}|^2 >= E_n|w_k|^2 for all k < n", range(1, N + 1),
         lambda n: next(("k=%d" % k for k in range(0, n) if Ek[(k + 1, n)] < Ek[(k, n)]), None))
    scan("D5 E_n|w_n|^2 >= E_n|w_k|^2 + E_n|w_{n-k}|^2, i.e. E_n[w_k.(w_n - w_k)] >= 0, all 1 <= k < n",
         range(2, N + 1), lambda n: next(("k=%d" % k for k in range(1, n) if E[n] < Ek[(k, n)] + Ek[(n - k, n)]), None))
    print("     note D5: by reversal (w_n - w_{n-i})_i is again uniform, so E_n|w_n - w_k|^2 = E_n|w_{n-k}|^2.")
    # lower-bound CL3 margins: Delta_1(n) = E_n|w_{n-1}|^2 - E_{n-1}, Delta_2(n) = E_n[w_n . w_1] - 1.
    # By reversal (w_n - w_{n-i})_i, E_n[w_n . w_1] = E_n[w_n . s_n] = E_n[w_{n-1} . s_n] + 1 and
    # |w_n|^2 = |w_{n-1}|^2 + 2 w_{n-1}.s_n + 1, so Delta_2(n) = (E_n - E_n|w_{n-1}|^2 - 1)/2 (exact, from the M table).
    d1 = {n: Ek[(n - 1, n)] - E[n - 1] for n in range(2, N + 1)}
    d2 = {n: (E[n] - Ek[(n - 1, n)] - 1) / 2 for n in range(2, N + 1)}
    print("     lower-bound CL3 margins (exact; Delta_1 + 2 Delta_2 = E_n - E_{n-1} - 1):")
    for n in range(2, N + 1):
        print("       n=%2d  Delta_1 = %-34s = %.6f   Delta_2 = %-34s = %.6f" % (n, d1[n], fl(d1[n]), d2[n], fl(d2[n])))
    scan("D6 (I1) Delta_1(n) = E_n|w_{n-1}|^2 - E_{n-1}|w_{n-1}|^2 >= 0", range(2, N + 1),
         lambda n: None if d1[n] >= 0 else "%s" % d1[n])
    scan("D7 (I2) Delta_2(n) = E_n[w_n . w_1] - 1 >= 0", range(2, N + 1), lambda n: None if d2[n] >= 0 else "%s" % d2[n])
    scan("D6b Delta_1(n) nondecreasing", range(2, N), lambda n: None if d1[n + 1] >= d1[n] else "")
    scan("D6c Delta_1(n) nondecreasing within parity", range(2, N - 1), lambda n: None if d1[n + 2] >= d1[n] else "")
    scan("D7b Delta_2(n) nondecreasing", range(2, N), lambda n: None if d2[n + 1] >= d2[n] else "")
    print("     NUMERICS Delta_2(n)/n^{2 nu - 1}:", " ".join("%d:%.4f" % (n, fl(d2[n]) / n ** (2 * 0.587597 - 1))
                                                       for n in range(2, N + 1)))
    print("     E_n|w_{n/2}|^2 / E_{n/2}:", " ".join("%d:%.4f" % (n, fl(Ek[(n // 2, n)] / E[n // 2])) for n in range(2, N + 1, 2)))
    mp = {n: F(D.max_point_mid(n // 2, n), D.c[n]) for n in range(2, N + 1)}
    print("     sup_x P_n(w_{floor(n/2)} = x) * n^{1/2} (DGHM 2016 Prop 1.3: bounded):",
          " ".join("%d:%.4f" % (n, fl(mp[n]) * math.sqrt(n)) for n in range(2, N + 1)))
    print("     NUMERICS sup_x P_n(w_{n/2} = x) * n^{3 nu}:",
          " ".join("%d:%.4f" % (n, fl(mp[n]) * n ** (3 * 0.587597)) for n in range(2, N + 1)))

    print()
    print("## F. Sup-norm and large deviations of the endpoint (exact, n <= %d)" % N)

    def dist_points(n):
        for (a, b, cc), (Ns, orb) in D.A[n].items():
            yield (a, b, cc), Ns
    Einf = {n: F(sum(Ns * a * a for (a, b, cc), Ns in dist_points(n)), D.c[n]) for n in range(1, N + 1)}
    scan("F1 E_n ||w_n||_inf^2 >= n^{2/3} (would be the 'nu >= 1/d' rate for the endpoint)", range(1, N + 1),
         lambda n: None if Einf[n] ** 3 >= n * n else "%.4f" % fl(Einf[n]))
    print("     E_n||w_n||_inf^2 / n^{2/3}:", " ".join("%d:%.4f" % (n, fl(Einf[n]) / n ** (2 / 3)) for n in range(1, N + 1)))
    print("     E_n||w_n||_inf^2 / E_n:", " ".join("%d:%.4f" % (n, fl(Einf[n] / E[n])) for n in range(1, N + 1)))
    scan("F2 E_n||w_n||_inf^2 / n nondecreasing", range(1, N),
         lambda n: None if Einf[n + 1] / (n + 1) >= Einf[n] / n else "")
    print("     NUMERICS large deviations: -(1/n) log P_n(|w_n|_1 >= v n), at the largest n of each parity:")
    for v in (0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9):
        vals = []
        for n in (N - 3, N - 2, N - 1, N):
            pr = F(sum(Ns for (a, b, cc), Ns in dist_points(n) if a + b + cc >= v * n), D.c[n])
            vals.append("%d:%.4f" % (n, -math.log(pr) / n if pr > 0 else float("inf")))
        print("       v=%.1f  %s" % (v, " ".join(vals)))
    print("     NUMERICS -(1/n) log P_n(|w_n|_2 >= v n):")
    for v in (0.3, 0.4, 0.5, 0.6, 0.7):
        vals = []
        for n in (N - 3, N - 2, N - 1, N):
            pr = F(sum(Ns for (a, b, cc), Ns in dist_points(n) if a * a + b * b + cc * cc >= v * v * n * n), D.c[n])
            vals.append("%d:%.4f" % (n, -math.log(pr) / n if pr > 0 else float("inf")))
        print("       v=%.1f  %s" % (v, " ".join(vals)))

    if gram_path:
        gram_section(gram_path)
    es = os.path.join(ROOT, "data", "end_steps_N18.txt")
    if os.path.exists(es):
        print()
        print("## G'. First/last step correlation from src/end_steps.cpp (independent program), data/end_steps_N18.txt")
        rows = {}
        for line in open(es):
            t = line.split()
            if t and t[0] == "E":
                v = list(map(int, t[2:]))
                rows[int(t[1])] = (6 * (v[0] - v[1]), 6 * sum(v))
        for n in sorted(rows):
            assert rows[n][1] == c[n]
        scan("G1' E_n[s_1 . s_n] >= 0", sorted(rows), lambda n: None if rows[n][0] >= 0 else
             "sum over n-step SAWs of s_1.s_n = %d, E_n[s_1.s_n] = %s" % (rows[n][0], F(rows[n][0], rows[n][1])))
        print("     E_n[s_1 . s_n]:", " ".join("%d:%.4e" % (n, rows[n][0] / rows[n][1]) for n in sorted(rows)))

    print()
    print("## E. NUMERICS: amplitudes")
    print("     E_n/n^{2 nu} (nu = 0.587597):", " ".join("%d:%.4f" % (n, fl(E[n]) / n ** (2 * 0.587597))
                                                  for n in range(1, R + 1)))


def gram_section(path):
    """step-correlation data: lines 'G n i j value' = sum over n-step SAWs of w_i . w_j (i <= j <= n)."""
    G = {}
    cn = {}
    for line in open(path):
        t = line.split()
        if not t or t[0].startswith("#"):
            continue
        if t[0] == "G":
            G[(int(t[1]), int(t[2]), int(t[3]))] = int(t[4])
        elif t[0] == "T":
            cn[int(t[1])] = int(t[2])
    Ng = max(cn)

    def g(n, i, j):
        if i == 0 or j == 0:
            return 0
        if i > j:
            i, j = j, i
        return G[(n, i, j)]

    def corr(n, i, j):  # sum over walks of s_i . s_j, s_i = w_i - w_{i-1}
        return g(n, i, j) - g(n, i - 1, j) - g(n, i, j - 1) + g(n, i - 1, j - 1)
    print()
    print("## G. Step correlations E_n[s_i . s_j] (exact, n <= %d, src/saw_gram.cpp)" % Ng)
    # sanity: diagonal = c_n (|s_i| = 1)
    assert all(corr(n, i, i) == cn[n] for n in range(1, Ng + 1) for i in range(1, n + 1))

    def pos(n):
        for i in range(1, n + 1):
            for j in range(i + 1, n + 1):
                if corr(n, i, j) < 0:
                    return "E_%d[s_%d . s_%d] = %s" % (n, i, j, F(corr(n, i, j), cn[n]))
        return None
    scan("G1 E_n[s_i . s_j] >= 0 for all 1 <= i < j <= n", range(1, Ng + 1), pos)
    neg = [(n, i, j) for n in range(1, Ng + 1) for i in range(1, n + 1) for j in range(i + 1, n + 1) if corr(n, i, j) < 0]
    print("     number of negative step correlations per n:",
          " ".join("%d:%d" % (n, sum(1 for t in neg if t[0] == n)) for n in range(1, Ng + 1)))
    if neg:
        print("     negative (n,i,j) for the largest n: %s" % [t[1:] for t in neg if t[0] == Ng][:40])
    print("     min_{i<j} E_n[s_i . s_j] and its location (i,j):")
    for n in range(3, Ng + 1):
        v, i, j = min((F(corr(n, i, j), cn[n]), i, j) for i in range(1, n + 1) for j in range(i + 1, n + 1))
        print("       n=%2d: %.3e at (%d,%d)" % (n, float(v), i, j))

    def rowsum(n):  # sum_j E[s_i.s_j] >= 0 for each i  (=> E|w_n|^2 >= ... )
        for i in range(1, n + 1):
            if sum(corr(n, i, j) for j in range(1, n + 1)) < cn[n]:
                return "i=%d" % i
        return None
    scan("G2 E_n[s_i . w_n] >= 1 for all i (row sums of the correlation matrix >= diagonal)", range(1, Ng + 1), rowsum)

    def wk(n):  # E_n[w_k . s_{k+1}] >= 0
        for k in range(1, n):
            if g(n, k, k + 1) - g(n, k, k) < 0:
                return "k=%d" % k
        return None
    scan("G3 E_n[w_k . s_{k+1}] >= 0 for all k (the next step points outward on average)", range(1, Ng + 1), wk)

    def wkn(n):  # E_n[w_k . (w_n - w_k)] >= 0
        for k in range(1, n):
            if g(n, k, n) - g(n, k, k) < 0:
                return "k=%d" % k
        return None
    scan("G4 E_n[w_k . (w_n - w_k)] >= 0 for all k (same as D5)", range(1, Ng + 1), wkn)

    def mono(n):  # j -> E[s_1 . s_j] nonincreasing
        v = [corr(n, 1, j) for j in range(2, n + 1)]
        for a, b in zip(v, v[1:]):
            if b > a:
                return "E[s_1.s_j] increases somewhere"
        return None
    scan("G5 j -> E_n[s_1 . s_j] nonincreasing", range(3, Ng + 1), mono)
    n = Ng
    print("     E_%d[s_1 . s_j], j = 2..%d:" % (n, n), " ".join("%.5f" % (corr(n, 1, j) / cn[n]) for j in range(2, n + 1)))
    m = n // 2
    print("     E_%d[s_%d . s_j], j = 1..%d:" % (n, m, n), " ".join("%.5f" % (corr(n, m, j) / cn[n]) for j in range(1, n + 1)))


if __name__ == "__main__":
    main()
