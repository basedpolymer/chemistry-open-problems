#!/usr/bin/env python3
"""
verify_symbolic.py -- exact computer-algebra check of the Jaschke-Skogestad marginal-cost formula.

Environment : WSL Fedora, Python 3 + python-flint 0.9 (fmpz_mpoly).  Single thread.
Runtime     : ~65 s with the defaults, peak RSS ~220 MB (see certificates/verify_symbolic.log).
              (The raw finite-wh model at M = 6 has 19 indeterminates and exceeds 4 GB; it is covered by L1 below.)
Usage       : python -I verify_symbolic.py [Mmax] [Minf]
              Mmax = largest M for the raw model with finite hot heat capacities (default 5),
              Minf = largest M for constant hot temperatures, wh = inf (default 8).

Model (Jaschke & Skogestad 2014, Sec. 3.1, Assumptions 1-3), one branch j with M exchangers in series,
shifted temperatures theta = T - T0 (theta_0 = 0):
    Q_i = w (theta_i - theta_{i-1})                                  (cold-side energy balance, eq. 9)
    Q_i = wh_i (thh_i - tho_i)                                       (hot-side energy balance, eq. 10)
    Q_i = UA_i [ (thh_i - theta_i) + (tho_i - theta_{i-1}) ] / 2     (arithmetic mean temperature difference, eq. 12)
All of w, UA_i, wh_i, thh_i are independent indeterminates (fully generic, nothing is eliminated by hand):
the 3x3 linear system of each exchanger is solved by Cramer's rule in the field Q(w, UA, wh, thh).

Claims checked (each prints PASS/FAIL; exit code 1 if any FAIL):
  S0(M)    the solution satisfies the three model equations of every exchanger;
  S1(M,i)  dQ_i/dw == a_i / 2   for all i <= M, where (CORRECTED recursion, cumulative sum)
               a_i = c_i (theta_i + theta_{i-1} - sum_{l<i} a_l),   c_i = (theta_i - theta_{i-1})/(thh_i - theta_{i-1});
  S2(M)    the paper's conjectured recursion  a^P_i = c_i (theta_i + theta_{i-1} - a^P_{i-1})  is identical to the
           corrected one for i <= 2 (symbolic identity) and DIFFERS for every i >= 3; also the paper's
           T_J = sum_i a^P_i differs from 2*d(sum Q_i)/dw, and a^P_3 - a_3 = c_3 a_1.  Non-vanishing of these
           rational functions is certified by exact evaluation at one rational point where no denominator vanishes
           (a rational function with a nonzero value is not identically zero); dQ_i/dw is evaluated from its
           symbolic expression, independently of S1;
  S3(M)    the same holds with all hot heat capacities infinite (constant hot temperatures, the paper's
           "proved" M_j = 3 case): corrected formula exact, paper's formula wrong for i >= 3;
  S4       the Maple output printed in the first line of the paper's eq. (A.12) equals the corrected a_3,
           and differs from the "simplified" last line of (A.12) by exactly c_3 * a_1;
  L1       one exchanger with generic inlet temperature: the finite-wh solution equals the constant-hot-temperature
           solution with UA' = (1/UA + 1/(2 wh))^-1 (cold side depends on (UA, wh) only through k = 2/UA + 1/wh).
           Since UA' does not depend on w, a branch with finite hot capacities has the same theta_i(w), Q_i(w) as
           the wh = inf branch with conductances UA'_i, so S1/S3 for wh = inf (M <= Minf) imply them for finite wh.
"""
import sys
import time
import flint
from fractions import Fraction

Mmax = int(sys.argv[1]) if len(sys.argv) > 1 else 5
FAILED = []


def report(name, ok):
    print(f"{'PASS' if ok else 'FAIL'}  {name}", flush=True)
    if not ok:
        FAILED.append(name)


class RF:
    """Rational function num/den over Z[vars], kept reduced (gcd) with a normalised denominator."""
    __slots__ = ("n", "d")

    def __init__(self, n, d=None):
        if d is None:
            d = n.context().from_dict({}) + 1
        g = n.gcd(d)
        if not g.is_one():
            n = n / g
            d = d / g
        if d.leading_coefficient() < 0:
            n, d = -n, -d
        self.n, self.d = n, d

    def __add__(self, o):
        o = lift(o, self)
        return RF(self.n * o.d + o.n * self.d, self.d * o.d)

    __radd__ = __add__

    def __neg__(self):
        return RF(-self.n, self.d)

    def __sub__(self, o):
        return self + (-lift(o, self))

    def __rsub__(self, o):
        return lift(o, self) - self

    def __mul__(self, o):
        o = lift(o, self)
        return RF(self.n * o.n, self.d * o.d)

    __rmul__ = __mul__

    def __truediv__(self, o):
        o = lift(o, self)
        if o.n.is_zero():
            raise ZeroDivisionError
        return RF(self.n * o.d, self.d * o.n)

    def deriv(self, k):
        return RF(self.n.derivative(k) * self.d - self.n * self.d.derivative(k), self.d * self.d)

    def is_zero(self):
        return self.n.is_zero()


def lift(x, ref):
    if isinstance(x, RF):
        return x
    ctx = ref.n.context()
    return RF(ctx.from_dict({}) + int(x))


def solve_branch(M, infinite_wh=False):
    names = ["w"] + [f"UA{i}" for i in range(1, M + 1)] + [f"thh{i}" for i in range(1, M + 1)]
    if not infinite_wh:
        names += [f"wh{i}" for i in range(1, M + 1)]
    ctx = flint.fmpz_mpoly_ctx.get(tuple(names), "lex")
    g = ctx.gens()
    one = RF(ctx.from_dict({}) + 1)
    zero = one - one
    w = RF(g[0])
    UA = [RF(g[i]) for i in range(1, M + 1)]
    TH = [RF(g[M + i]) for i in range(1, M + 1)]
    WH = [RF(g[2 * M + i]) for i in range(1, M + 1)] if not infinite_wh else None

    theta = [zero]
    tho_l, Q = [], []
    for i in range(M):
        tp = theta[-1]
        if infinite_wh:
            # constant hot temperature: tho = thh, unknowns (t, q):
            #   q - w t = -w tp ;  q + (UA/2) t = UA (thh - tp/2)
            # => t (w + UA/2) = UA thh - UA tp/2 + w tp
            t = (UA[i] * TH[i] - UA[i] * tp / 2 + w * tp) / (w + UA[i] / 2)
            q = w * (t - tp)
            tho = TH[i]
        else:
            # unknowns x = (t, tho, q); rows:  -w t + 0 tho + q = -w tp
            #                                   0 t + wh tho + q = wh thh
            #                            (UA/2) t - (UA/2) tho + q = UA thh/2 - UA tp/2
            A = [[-w, zero, one], [zero, WH[i], one], [UA[i] / 2, -UA[i] / 2, one]]
            b = [-w * tp, WH[i] * TH[i], UA[i] * TH[i] / 2 - UA[i] * tp / 2]

            def det3(m):
                return (m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
                        - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
                        + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]))

            D = det3(A)
            sol = []
            for col in range(3):
                Ac = [[(b[r] if c == col else A[r][c]) for c in range(3)] for r in range(3)]
                sol.append(det3(Ac) / D)
            t, tho, q = sol
        theta.append(t)
        tho_l.append(tho)
        Q.append(q)
    return ctx, w, UA, TH, WH, theta, tho_l, Q


def check(M, infinite_wh=False):
    tag = f"M={M}{', wh=inf' if infinite_wh else ''}"
    t0 = time.time()
    ctx, w, UA, TH, WH, theta, tho, Q = solve_branch(M, infinite_wh)
    # S0: model equations
    ok = True
    for i in range(M):
        r1 = Q[i] - w * (theta[i + 1] - theta[i])
        r3 = Q[i] - UA[i] * ((TH[i] - theta[i + 1]) + (tho[i] - theta[i])) / 2
        ok &= r1.is_zero() and r3.is_zero()
        if not infinite_wh:
            r2 = Q[i] - WH[i] * (TH[i] - tho[i])
            ok &= r2.is_zero()
    report(f"S0 [{tag}] model equations satisfied", ok)
    # derivatives (variable index 0 is w)
    dQ = [q.deriv(0) for q in Q]
    # corrected recursion, computed symbolically from temperatures only
    a_cor = []
    cum = theta[0] - theta[0]
    for i in range(1, M + 1):
        c = (theta[i] - theta[i - 1]) / (TH[i - 1] - theta[i - 1])
        ai = c * (theta[i] + theta[i - 1] - cum)
        a_cor.append(ai)
        cum = cum + ai
    for i in range(M):
        report(f"S1 [{tag}] dQ_{i+1}/dw == a_{i+1}/2 (corrected, cumulative recursion)",
               (dQ[i] - a_cor[i] / 2).is_zero())
    # S2/S3: the paper's recursion a^P_i = c_i (theta_i + theta_{i-1} - a^P_{i-1}).
    #  - for i = 1, 2 it is compared with a_cor symbolically (identity);
    #  - for i >= 3, a rational function is shown to be NONZERO by exact evaluation at one rational point
    #    (a nonzero value at a point where all denominators are nonzero proves the function is not identically 0).
    s2ok = True
    a_pap = []
    for i in range(1, min(M, 2) + 1):
        c = (theta[i] - theta[i - 1]) / (TH[i - 1] - theta[i - 1])
        prev = a_pap[-1] if a_pap else theta[0] - theta[0]
        a_pap.append(c * (theta[i] + theta[i - 1] - prev))
        s2ok &= (a_pap[-1] - a_cor[i - 1]).is_zero()
    if M >= 3:
        nv = len(ctx.gens())
        # deterministic integer point: w, UA_i, thh_i (increasing), wh_i
        pt = [37] + [50 + 13 * i for i in range(M)] + [40 + 17 * i for i in range(M)]
        if not infinite_wh:
            pt += [90 + 29 * i for i in range(M)]
        assert len(pt) == nv

        def ev(rf):
            d = int(rf.d(*pt))
            assert d != 0, "denominator vanishes at evaluation point"
            return Fraction(int(rf.n(*pt)), d)

        th_v = [ev(t) for t in theta]
        thh_v = [Fraction(x) for x in pt[1 + M:1 + 2 * M]]
        dQ_v = [ev(d) for d in dQ]
        ap, ac, cumv = [], [], Fraction(0)
        for i in range(1, M + 1):
            assert thh_v[i - 1] != th_v[i - 1]
            c = (th_v[i] - th_v[i - 1]) / (thh_v[i - 1] - th_v[i - 1])
            ac.append(c * (th_v[i] + th_v[i - 1] - cumv))
            cumv += ac[-1]
            ap.append(c * (th_v[i] + th_v[i - 1] - (ap[-1] if ap else 0)))
        for i in range(M):
            s2ok &= (ap[i] == ac[i]) if i < 2 else (ap[i] != ac[i])
            s2ok &= (dQ_v[i] == ac[i] / 2)          # consistency with S1 at the point
        s2ok &= sum(ap) != 2 * sum(dQ_v)
        s2ok &= (ap[2] - ac[2]) == ac[0] * (th_v[3] - th_v[2]) / (thh_v[2] - th_v[2])   # d_3 = c_3 a_1
    name = "S3" if infinite_wh else "S2"
    report(f"{name} [{tag}] paper recursion == corrected for i<=2 (identity); != for every i>=3 and"
           f" paper T_J != 2 d(sum Q)/dw (exact evaluation at a rational point); d_3 == c_3 a_1"
           if M >= 3 else f"{name} [{tag}] paper recursion == corrected (identity, M<=2)", s2ok)
    print(f"      ({tag}: {time.time()-t0:.1f} s)", flush=True)


def check_A12():
    ctx = flint.fmpz_mpoly_ctx.get(("t1", "t2", "t3", "h1", "h2", "h3"), "lex")
    t1, t2, t3, h1, h2, h3 = [RF(x) for x in ctx.gens()]
    # numerator and denominator exactly as printed in the first line of eq. (A.12), p. 47
    N = (-t1*t1*t2*h1 + t1*t1*t2*h2 + t1*t1*t3*h1 - t1*t1*t3*h2 - t1*t1*t2*t2 + t3*t3*h1*h2
         + t1*t2*t2*h1 - t2*t2*h1*h2 - t1*t3*t3*h1 + t1*t1*t2*t3 - t2*t2*t3*h1 + t2*t2*t2*h1)
    Dn = h1 * (-t2*h2 + h2*h3 + t1*t2 - t1*h3)
    line1 = N / Dn
    a1 = t1 * t1 / h1
    a2 = (t2 - t1) * (t2 + t1 - a1) / (h2 - t1)
    c3 = (t3 - t2) / (h3 - t2)
    a3_paper = c3 * (t3 + t2 - a2)          # last line of (A.12)
    a3_corr = c3 * (t3 + t2 - a1 - a2)       # corrected recursion
    report("S4 paper (A.12) line 1 (Maple output) == corrected a_3", (line1 - a3_corr).is_zero())
    report("S4 paper (A.12) last line - line 1 == c_3 * a_1 (nonzero)",
           (a3_paper - line1 - c3 * a1).is_zero() and not (c3 * a1).is_zero())


def check_L1():
    """L1: one exchanger with GENERIC inlet temperature tin: the finite-wh solution of (9),(10),(12) equals the
    constant-hot-temperature solution with UA' = (1/UA + 1/(2 wh))^{-1}, i.e. the exchanger acts on the cold
    stream through k = 2/UA + 1/wh only.  Hence every M-exchanger branch with finite hot capacities has exactly the
    same theta_i(w), Q_i(w) (and therefore the same w-derivatives, UA' being independent of w) as the branch with
    infinite hot capacities and conductances UA'_i.  This transfers S1/S3 for wh = inf to finite wh."""
    ctx = flint.fmpz_mpoly_ctx.get(("w", "UA", "wh", "thh", "tin"), "lex")
    w, UA, WH, TH, tin = [RF(x) for x in ctx.gens()]
    one = RF(ctx.from_dict({}) + 1)
    zero = one - one
    # finite wh: Cramer on the raw 3x3 system (same rows as in solve_branch)
    A = [[-w, zero, one], [zero, WH, one], [UA / 2, -UA / 2, one]]
    b = [-w * tin, WH * TH, UA * TH / 2 - UA * tin / 2]

    def det3(m):
        return (m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
                - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
                + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]))

    D = det3(A)
    t_fin = det3([[(b[r] if c == 0 else A[r][c]) for c in range(3)] for r in range(3)]) / D
    UAp = one / (one / UA + one / (2 * WH))
    t_inf = (UAp * TH - UAp * tin / 2 + w * tin) / (w + UAp / 2)
    k = 2 * one / UA + one / WH
    t_k = tin + 2 * (TH - tin) / (one + w * k)
    report("L1 finite-wh exchanger == constant-hot-temperature exchanger with UA' = (1/UA+1/(2wh))^-1 "
           "== tin + 2 (thh - tin)/(1 + w k)", (t_fin - t_inf).is_zero() and (t_fin - t_k).is_zero())


if __name__ == "__main__":
    # Mmax  : largest M checked with the raw finite-wh model (19 indeterminates at M = 6 needs > 4 GB; default 5)
    # Minf  : largest M checked with constant hot temperatures (default 8)
    Minf = int(sys.argv[2]) if len(sys.argv) > 2 else 8
    T = time.time()
    check_A12()
    check_L1()
    for M in range(1, Mmax + 1):
        check(M, infinite_wh=False)
    for M in range(1, Minf + 1):
        check(M, infinite_wh=True)
    print(f"total {time.time()-T:.1f} s")
    if FAILED:
        print(f"FAIL: {len(FAILED)} claim(s) failed")
        sys.exit(1)
    print("ALL PASS")
