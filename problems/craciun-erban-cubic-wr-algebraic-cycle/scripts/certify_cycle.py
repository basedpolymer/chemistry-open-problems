#!/usr/bin/env python3
"""Rigorous certificate: a cubic weakly reversible mass-action system with an algebraic limit cycle.

Answers Craciun & Erban (J. Math. Biol. 2025, arXiv:2406.05057): "it remains an open question whether a cubic
weakly reversible system can have an algebraic limit cycle" (W^a(3) >= 0 in their Table 2).

Two networks are certified: N10 (10 reactions, integer rate constants, minimal support found by MILP) and N17
(17 reactions, found during exploration). For each one, the checks are (exact rational arithmetic with sympy,
rigorous ball arithmetic with python-flint/Arb):
  1. All complexes have molecularity <= 3 (so the network lies in M_3), and the reaction graph is strongly
     connected (so it is weakly reversible, hence in W_3); all rate constants are positive.
  2. The circle h = (x-2)^2 + (y-2)^2 - 1 = 0 (contained in the open positive quadrant) is invariant:
     f h_x + g h_y = s h exactly, with an explicit polynomial cofactor s.
  3. The tangential speed P(theta) = F . (-sin theta, cos theta) on x = 2 + cos theta, y = 2 + sin theta is
     strictly positive (certified two ways: Arb enclosures on 2^12 subintervals, and an exact resultant with no real
     root in [1, 3]); so the circle carries no equilibrium and is a periodic orbit, run counter-clockwise.
  4. The characteristic exponent  Lambda = oint s dt = int_0^{2 pi} s(theta)/P(theta) d theta  is enclosed by
     rigorous Arb integration and is nonzero; the same holds for oint div F dt (which equals Lambda).
     Hence the periodic orbit is hyperbolic, thus isolated: an algebraic limit cycle (stable iff Lambda < 0).
Runtime ~ 1 min. Environment: Python 3.12, sympy 1.14, python-flint 0.9. Exits 1 on any failure.
"""
import sys
import sympy as sp
from flint import arb, acb, ctx, fmpq

ctx.prec = 128
x, y = sp.symbols("x y")
R = sp.Rational
# (source, product, rate) with complexes as (#X, #Y)
# N10: minimal-support network (10 reactions, integer rate constants); N17: network found during exploration.
NETS = {
    "N10": [((0, 0), (0, 2), 434), ((0, 0), (1, 1), 448), ((0, 1), (0, 0), 3374), ((0, 2), (0, 3), 1728),
            ((0, 3), (0, 1), 194), ((0, 3), (2, 1), 11), ((1, 1), (1, 2), 328), ((1, 2), (0, 3), 250),
            ((2, 1), (3, 0), 216), ((3, 0), (1, 2), 28)],
    "N17": [((0, 0), (1, 0), R(801, 1000000)), ((0, 1), (0, 0), R(104261883, 56000000)),
            ((0, 2), (0, 3), R(14893029, 14000000)), ((0, 3), (0, 0), R(14847677, 168000000)),
            ((0, 3), (0, 1), R(211, 1000000)), ((0, 3), (0, 2), R(213, 500000)),
            ((1, 0), (2, 0), R(1, 560000)), ((1, 0), (3, 0), R(115777, 250000)),
            ((1, 1), (0, 1), R(37358609, 28000000)), ((1, 1), (2, 0), R(748283, 560000)),
            ((1, 2), (0, 3), R(213551, 400000)), ((2, 0), (3, 0), R(7497783, 5600000)),
            ((2, 1), (3, 0), R(7477489, 56000000)), ((3, 0), (0, 3), R(11166451, 84000000)),
            ((3, 0), (1, 1), R(59813, 56000000)), ((3, 0), (1, 2), R(107, 500000)),
            ((3, 0), (2, 1), R(53, 1000000))],
}
fails = 0


def check(c, m):
    global fails
    print(("PASS " if c else "FAIL ") + m)
    fails += (not c)


def name(c):
    t = ([f"{c[0]}X" if c[0] > 1 else "X"] if c[0] else []) + ([f"{c[1]}Y" if c[1] > 1 else "Y"] if c[1] else [])
    return " + ".join(t) or "0"


def to_arb_poly(expr):
    poly = sp.Poly(expr, x, y)
    return [(i, j, fmpq(int(c.p), int(c.q))) for (i, j), c in poly.terms()]


def ev(coeffs, X, Y):
    tot = X * 0
    for i, j, c in coeffs:
        tot += c * X ** i * Y ** j
    return tot



def certify(label, NET):
    print(f"===== {label}: {len(NET)} reactions =====")
    print("Network:")
    for s_, p_, k_ in NET:
        print(f"   {name(s_):8s} -> {name(p_):8s}   k = {k_}")

    # 1. class W_3
    cx = {c for s_, p_, _ in NET for c in (s_, p_)}
    check(all(k > 0 for *_, k in NET), f"all {len(NET)} rate constants are positive rationals")
    check(max(sum(c) for c in cx) == 3, f"{len(cx)} complexes, maximal molecularity 3 (network in M_3)")
    adj = {c: set() for c in cx}
    for s_, p_, _ in NET:
        adj[s_].add(p_)


    def reach(u):
        seen, st = {u}, [u]
        while st:
            a = st.pop()
            for b in adj[a]:
                if b not in seen:
                    seen.add(b); st.append(b)
        return seen


    check(all(reach(u) == cx for u in cx), "reaction graph strongly connected => weakly reversible (network in W_3)")

    # 2. invariance
    f = sp.expand(sum(k * (p_[0] - s_[0]) * x ** s_[0] * y ** s_[1] for s_, p_, k in NET))
    g = sp.expand(sum(k * (p_[1] - s_[1]) * x ** s_[0] * y ** s_[1] for s_, p_, k in NET))
    h = (x - 2) ** 2 + (y - 2) ** 2 - 1
    num = sp.expand(f * sp.diff(h, x) + g * sp.diff(h, y))
    s, r = sp.div(num, h, x, y)
    check(r == 0, "f h_x + g h_y = s h exactly (circle (x-2)^2+(y-2)^2=1 invariant)")
    print("   f =", f)
    print("   g =", g)
    print("   cofactor s =", sp.expand(s))
    divF = sp.expand(sp.diff(f, x) + sp.diff(g, y))

    # 3. tangential speed has no zero on the circle
    P = sp.expand(-f * (y - 2) + g * (x - 2))
    res = sp.Poly(sp.resultant(P, h, y), x)
    check(res.count_roots(1, 3) == 0, f"resultant_y(P, h) (degree {res.degree()}) has no real root in [1, 3]")
    check(P.subs({x: 3, y: 2}) > 0, f"P(3, 2) = {P.subs({x: 3, y: 2})} > 0, so P > 0 on the whole circle")
    Pc, sc, dc = to_arb_poly(P), to_arb_poly(sp.expand(s)), to_arb_poly(divF)

    # interval check of P > 0 (independent of the resultant)
    N = 4096
    two_pi = 2 * arb.pi()
    minP = None
    ok = True
    for i in range(N):
        th = (arb(i) + arb("0 +/- 0.5") + arb(0.5)) * two_pi / N      # ball covering [i, i+1] * 2pi/N
        val = ev(Pc, 2 + th.cos(), 2 + th.sin())
        if not (val > 0):
            ok = False
        lo = val.lower()
        minP = lo if minP is None or lo < minP else minP
    check(ok, f"Arb: P(theta) > 0 on all {N} subintervals of [0, 2 pi] (certified lower bound {float(minP):.4e})")


    # 4. rigorous integrals
    def integrand(coeffs):
        def F(z, analytic):
            X, Y = 2 + z.cos(), 2 + z.sin()
            return ev(coeffs, X, Y) / ev(Pc, X, Y)
        return F


    one = [(0, 0, fmpq(1))]
    period = acb.integral(integrand(one), 0, 2 * acb.pi(), rel_tol=arb(2) ** -60)
    lam = acb.integral(integrand(sc), 0, 2 * acb.pi(), rel_tol=arb(2) ** -60)
    lamdiv = acb.integral(integrand(dc), 0, 2 * acb.pi(), rel_tol=arb(2) ** -60)
    print("   period          T      =", period.real.str(15, radius=True))
    print("   oint s dt       Lambda =", lam.real.str(15, radius=True))
    print("   oint div F dt          =", lamdiv.real.str(15, radius=True))
    check(period.real > 0 and abs(period.imag) < 1e-20, "period enclosure is a positive real number")
    stab = "orbitally asymptotically stable" if lam.real < 0 else "unstable (repelling)"
    check(lam.real < 0 or lam.real > 0, f"Lambda = oint s dt != 0 (certified): hyperbolic limit cycle, {stab}")
    check(not lamdiv.real.contains(0) and (lam.real - lamdiv.real).contains(0),
          "oint div F dt != 0 and its enclosure is consistent with Lambda")


for label, NET in NETS.items():
    certify(label, NET)
print("ALL PASS" if not fails else f"{fails} FAILURES")
sys.exit(1 if fails else 0)
