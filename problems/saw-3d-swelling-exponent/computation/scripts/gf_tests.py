#!/usr/bin/env python3
"""gf_tests.py -- generating-function tests on the exact data (requests of the solvers' candidate-lemmas.md files).

Environment: Python 3.8+ standard library (fractions, math).  Runtime < 1 min.
Usage: python3 gf_tests.py [data file, default ../data/saw_enum_N19.txt.gz]

Sections
  U4  (upper-bound CL4)  B_{T+1}(z) <= B_T(z) for z in {1/10, 3/20, 9/50, 1/5}: RIGOROUS enclosures of
      B_T(z) = sum_n b_{n,T} z^n (b_{n,T} = n-step bridges of height T, b_{0,0} = 1) in exact rational arithmetic:
        lower = sum_{n<=N} b_{n,T} z^n,
        upper = lower + tail,  tail = sum_{n>N} c_n z^n bounded through submultiplicativity c_{a+b} <= c_a c_b:
        sum_{n>=N} c_n z^n <= rho/(1-rho) * sum_{r<N} c_r z^r,  rho = c_N z^N < 1  (write n = qN + r, q >= 1).
      Only our own enumerated numbers (n <= N) are used, so the verdicts are certificates (exact arithmetic, no
      published input).  A verdict "VERIFIED" means upper(B_{T+1}) < lower(B_T); "REFUTED" means
      lower(B_{T+1}) > upper(B_T); otherwise "UNDECIDED".
  U3  (upper-bound CL3)  NUMERICS: truncated sums B_T^{(N)}(z_c) and a ratio-method tail extrapolation;
      z_c = 1/4.684039931 (Clisby 2013).  Fit of B_T(z_c) ~ C T^{-a}.
  U5  (upper-bound CL5)  NUMERICS: phi_n(v) = -(1/n) log(#{bridges with x1(w_n) >= v n} / b_n).
  L4  (lower-bound CL4)  NUMERICS: S(M) = sum_{m<=M} m p_m mu^{-m} (polygons p_m, CLS 2007 for m <= 32), local
      exponents of the increments; RIGOROUS enclosures of G_z(x) = sum_n c_n(x) z^n at z in {1/10, 3/20, 9/50, 1/5}
      and the location of max_x G_z(x).
"""
import math
import os
import sys
from fractions import Fraction as F

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from sawdata import SawData, read_reference  # noqa: E402

ROOT = os.path.dirname(HERE)
MU = 4.684039931
ZC = 1 / MU
ZS = [F(1, 10), F(3, 20), F(9, 50), F(1, 5)]


def tail_bound(c, N, z):
    """rigorous upper bound on sum_{n>N} c_n z^n using c_n, n <= N (exact rationals)."""
    rho = c[N] * z ** N
    if rho >= 1:
        return None
    chi = sum(c[r] * z ** r for r in range(N))
    total_ge_N = rho / (1 - rho) * chi  # >= sum_{n >= N} c_n z^n
    return total_ge_N - c[N] * z ** N  # remove the exact n = N term


def main():
    path = sys.argv[1] if len(sys.argv) > 1 else os.path.join(ROOT, "data", "saw_enum_N19.txt.gz")
    D = SawData(path)
    N = D.N
    c = D.c
    span = {n: D.bridge_span_dist(n) for n in range(1, N + 1)}
    Tmax = N
    # b[n][T]
    b = {0: {0: 1}}
    for n in range(1, N + 1):
        b[n] = span[n]
    print("# gf_tests on %s (exact data n <= %d)" % (os.path.basename(path), N))

    print("\n## U4. B_{T+1}(z) <= B_T(z): rigorous enclosures (exact rationals, own data n <= %d only)" % N)
    for z in ZS:
        tb = tail_bound(c, N, z)
        if tb is None:
            print("  z = %s: rho = c_N z^N >= 1, no rigorous tail" % z)
            continue
        lo = {T: sum(b[n].get(T, 0) * z ** n for n in range(0, N + 1)) for T in range(0, Tmax + 1)}
        up = {T: lo[T] + tb for T in lo}
        res = []
        first_undecided = None
        for T in range(0, Tmax):
            if up[T + 1] < lo[T]:
                res.append("V")
            elif lo[T + 1] > up[T]:
                res.append("R")
            else:
                res.append("?")
                if first_undecided is None:
                    first_undecided = T
        verified = [T for T in range(Tmax) if res[T] == "V"]
        refuted = [T for T in range(Tmax) if res[T] == "R"]
        print("  z = %-5s (%.3f): tail <= %.3e; VERIFIED for T in %s; REFUTED for T in %s; UNDECIDED from T = %s"
              % (z, float(z), float(tb), compress(verified), refuted if refuted else "none", first_undecided))
        print("      B_T(z) enclosures T=0..6: " + "  ".join("[%.5f, %.5f]" % (float(lo[T]), float(up[T]))
                                                          for T in range(0, 7)))

    print("\n## U4'. Same test with CITED inputs for a sharper tail (exact rationals):")
    print("  tail_{n>N} <= sum_{N<n<=28} b_n z^n + sum_{n>28} (mu_+ z)^n, with b_n from Clisby-Conway-Guttmann 2016")
    print("  Table 6 (data/reference_bridges.tsv; equal to our b_n for n <= N), b_n <= mu^n (BDGS lectures,")
    print("  arXiv:1206.2092, eq. (2.4)) and mu <= mu_+ = 4.7388 (Ponitz-Tittmann, Electron. J. Combin. 7 (2000) R21,")
    print("  doi:10.37236/1499, Table 2: 4.7387; quoted as 4.73871 by Kim-Pinna 2025).")
    rb = {int(r[0]): r for r in read_reference(os.path.join(ROOT, "data", "reference_bridges.tsv"))}
    bref = {n: int(rb[n][1]) for n in rb if rb[n][1] != "-"}
    assert all(bref[n] == D.b[n] for n in range(1, N + 1))
    mup = F(47388, 10000)
    for z, use_ccg in [(z, u) for z in ZS for u in (False, True)]:
        q = mup * z
        if q >= 1:
            continue
        if use_ccg:  # published b_n for N < n <= 28 (not re-verified by us beyond N)
            tb = sum(bref[n] * z ** n for n in range(N + 1, 29)) + q ** 29 / (1 - q)
        else:  # only b_n <= mu^n beyond our own data
            tb = q ** (N + 1) / (1 - q)
        lo = {T: sum(b[n].get(T, 0) * z ** n for n in range(0, N + 1)) for T in range(0, Tmax + 1)}
        verified = [T for T in range(0, Tmax) if lo[T + 1] + tb < lo[T]]
        refuted = [T for T in range(0, Tmax) if lo[T + 1] > lo[T] + tb]
        und = next((T for T in range(Tmax) if T not in verified and T not in refuted), None)
        print("  z = %-5s %-22s: tail <= %.3e; VERIFIED for T in %s; REFUTED for T in %s; first UNDECIDED T = %s"
              % (z, "(with CCG b_n, n<=28)" if use_ccg else "(own b_n + mu^n)", float(tb), compress(verified),
                 refuted if refuted else "none", und))

    print("\n## U3. NUMERICS: B_T(z_c), z_c = 1/%.9f" % MU)
    est = {}
    for T in range(1, 9):
        terms = [(n, b[n].get(T, 0) * ZC ** n) for n in range(T, N + 1) if b[n].get(T, 0) > 0]
        trunc = sum(t for _, t in terms)
        # ratio-method tail: q = geometric mean of the last two term ratios
        if len(terms) >= 4:
            q1 = terms[-1][1] / terms[-2][1]
            q2 = terms[-2][1] / terms[-3][1]
            q = math.sqrt(q1 * q2)
            tail = terms[-1][1] * q / (1 - q) if q < 1 else float("inf")
        else:
            q, tail = float("nan"), float("nan")
        est[T] = (trunc, trunc + tail, q)
        print("  T=%d  B_T^(N=%d)(z_c) = %.5f   last term ratio q = %.4f   extrapolated %.5f   (tail/total %.2f)"
              % (T, N, trunc, q, trunc + tail, tail / (trunc + tail) if trunc + tail > 0 else float("nan")))
    Ts = [T for T in range(1, 5) if math.isfinite(est[T][1])]
    if len(Ts) >= 2:
        xs = [math.log(T) for T in Ts]
        ys = [math.log(est[T][1]) for T in Ts]
        mx, my = sum(xs) / len(xs), sum(ys) / len(ys)
        slope = sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / sum((x - mx) ** 2 for x in xs)
        print("  log-log slope of the extrapolated B_T(z_c), T = %s: a ~ %.3f (prediction 0.661; T this small is"
              " far from asymptotic, and the extrapolation is uncontrolled for T >= 3)" % (Ts, -slope))

    print("\n## U5. NUMERICS: bridge large deviations phi_n(v) = -(1/n) log P(x1(w_n) >= v n | bridge)")
    for v in (0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9):
        row = []
        for n in range(N - 5, N + 1):
            num = sum(cnt for s, cnt in span[n].items() if s >= v * n)
            row.append("%d:%.4f" % (n, -math.log(num / D.b[n]) / n if num else float("inf")))
        print("  v=%.1f  %s" % (v, " ".join(row)))

    print("\n## L4. NUMERICS: S(M) = sum_{m<=M} m p_m mu^{-m}, polygons from CLS (m <= 32), mu = %.9f" % MU)
    ref = {int(r[0]): r for r in read_reference(os.path.join(ROOT, "data", "reference_series.tsv"))}
    p = {m: int(ref[m][3]) for m in ref if ref[m][3] not in ("-",) and m >= 4 and m % 2 == 0}
    # own polygons for consistency
    for m in p:
        if m - 1 <= N and m >= 4:
            assert 3 * D.cnx(m - 1, (1, 0, 0)) == m * p[m]
    S = 0.0
    prev = None
    rows = []
    for m in sorted(p):
        inc = m * p[m] / MU ** m
        S += inc
        loc = "" if prev is None else "%.3f" % (math.log(prev[1] / inc) / math.log(m / prev[0]))
        rows.append("m=%d S=%.6f inc=%.3e loc.exp=%s" % (m, S, inc, loc))
        prev = (m, inc)
    print("  " + "\n  ".join(rows))
    print("  (local exponent of the increments m p_m mu^{-m}; convergence needs > 1; expected 3 - alpha - 1 = 1.76)")
    print("  G_{z_c}(e1) = (mu/3) S(inf) would follow if finite; partial value (mu/3) S(32) = %.5f" % (MU / 3 * S))

    print("\n## L4b. max_x G_z(x), G_z(x) = sum_n c_n(x) z^n: rigorous enclosures (own data n <= %d)" % N)
    classes = sorted({k for n in D.A for k in D.A[n]})
    for z in ZS:
        tb = tail_bound(c, N, z)
        if tb is None:
            continue
        G = {x: sum(D.cnx(n, x) * z ** n for n in range(0, N + 1)) for x in classes if x != (0, 0, 0)}
        best = max(G, key=lambda x: G[x])
        # rigorous: is G_z(best) - tail > G_z(y) for all other y?  (tail bounds each G_z(y) - truncation)
        rivals = [y for y in G if y != best and G[y] + tb >= G[best]]
        print("  z = %-5s: argmax of the truncated G_z is %s (lower bound %.6f); classes not excluded rigorously: %s"
              % (z, best, float(G[best]), rivals[:8] if rivals else "none (argmax certified)"))


def compress(lst):
    if not lst:
        return "[]"
    out, i = [], 0
    while i < len(lst):
        j = i
        while j + 1 < len(lst) and lst[j + 1] == lst[j] + 1:
            j += 1
        out.append(str(lst[i]) if i == j else "%d..%d" % (lst[i], lst[j]))
        i = j + 1
    return "[" + ", ".join(out) + "]"


if __name__ == "__main__":
    main()
