#!/usr/bin/env python3
"""Classification of all bimolecular (2,4,2) networks up to dynamical equivalence (and X <-> Y), with the
reason why each class does or does not admit two asymptotically stable equilibria in the closed quadrant.

The reasons refer to the results of paper.tex (Lemma lem:few, Propositions prop:bdry, prop:axispos,
prop:four, Corollary cor:pos, Theorem thm:index):
  T   dynamically trivial: no kappa > 0 with Gamma kappa = 0, hence no positive equilibrium; by
      Proposition prop:bdry every pair of AS equilibria contains a positive one.
  S2  at most two distinct sources: no isolated positive equilibrium.
  S3  three distinct sources: at most one isolated positive equilibrium; excluded because a pair of AS
      equilibria must be (origin, positive) and then Theorem thm:index needs a second positive equilibrium.
  M4  four distinct sources, Phi monotone or constant on I (Q has no zero in I; certificate
      positive_equilibria.py): at most one positive equilibrium; excluded as for S3.
  F-in  Table 1 network with an inflow reaction 0 -> C: the origin is not an equilibrium, so no pair.
  F-un  Table 1 network without inflow whose origin is unstable for every kappa (checked symbolically here):
      (i) an axis is invariant and the coordinate grows linearly on it, or (ii) J(0) is irreducible with
      det J(0) < 0 for all kappa (hyperbolic saddle whose unstable eigenvector points into the open quadrant).
  YES one of the 10 orange networks (certificate certify_ten.py).

Checks (PASS/FAIL): totals 27195 networks / 5716 classes; every class gets exactly one reason; the YES classes
are exactly the 10 of BBH Table 1; the F-un classes are exactly 5 and each passes (i) or (ii); 30 = 15 + 5 + 10.
Writes ../certificates/classes.csv.

Runtime: ~1 min.  Environment: Python 3.12, sympy 1.14 (WSL Fedora 44).  Exit status 1 on any FAIL.
"""
import os
import sys
import csv
import sympy as sp
from networks import *
from positive_equilibria import analyse as analyse4, sign_changes_on_I, t as tsym

FAILS = []


def check(cond, msg):
    print(("PASS " if cond else "FAIL ") + msg)
    if not cond:
        FAILS.append(msg)


def collinear(pts):
    pts = list(pts)
    a = pts[0]
    return all(cross((p[0] - a[0], p[1] - a[1]), (q[0] - a[0], q[1] - a[1])) == 0 for p in pts for q in pts)


def origin_never_AS(net):
    """symbolic check, valid for all kappa > 0; returns a string describing the argument or None"""
    ks = sp.symbols('k1:5', positive=True)
    x, y = sp.symbols('x y')
    F = sp.Matrix([0, 0])
    for kj, r in zip(ks, net):
        F += kj * x**r[0][0] * y**r[0][1] * sp.Matrix(vec(r))
    F = sp.expand(F)
    if F.subs({x: 0, y: 0}) != sp.Matrix([0, 0]):
        return None
    # (i) invariant x-axis with linear growth along it
    if sp.expand(F[1].subs(y, 0)) == 0:
        lin = sp.Poly(F[0].subs(y, 0), x).coeff_monomial(x)
        if lin != 0 and all(c > 0 for c in sp.Poly(lin, *ks).coeffs()):
            return "x-axis invariant, xdot = (%s) x + O(x^2) on it" % lin
    if sp.expand(F[0].subs(x, 0)) == 0:
        lin = sp.Poly(F[1].subs(x, 0), y).coeff_monomial(y)
        if lin != 0 and all(c > 0 for c in sp.Poly(lin, *ks).coeffs()):
            return "y-axis invariant, ydot = (%s) y + O(y^2) on it" % lin
    # (ii) irreducible Metzler J(0) with det < 0 for all kappa
    J0 = F.jacobian([x, y]).subs({x: 0, y: 0})
    d = sp.expand(J0.det())
    a12, a21 = J0[0, 1], J0[1, 0]
    pos = lambda e: e != 0 and all(c > 0 for c in sp.Poly(e, *ks).coeffs())
    if d != 0 and all(c < 0 for c in sp.Poly(d, *ks).coeffs()) and pos(a12) and pos(a21):
        return "J(0) irreducible, det J(0) = %s < 0" % d
    return None


def main():
    nets = all_networks()
    cl = classes()
    check(len(nets) == 27195, "27195 rank-2 bimolecular networks with 4 reactions (got %d)" % len(nets))
    check(len(cl) == 5716, "5716 classes up to dynamical equivalence and X<->Y (got %d)" % len(cl))
    orange = set(canon_class(parse_net(s)) for s in ORANGE10)
    table1 = set(canon_class(parse_net(s)) for s in TABLE1)
    rows = []
    counts = {}
    fun = []
    t1nets = {canon_class(parse_net(s_)): parse_net(s_) for s_ in TABLE1}
    for key, mem in sorted(cl.items(), key=lambda kv: net_str(kv[1][0])):
        if key in t1nets:   # use the form printed in BBH Table 1 as representative
            mem = [t1nets[key]] + [m for m in mem if m != t1nets[key]]
        net = mem[0]
        srcs = sorted(set(r[0] for r in net))
        nontriv = cone_canon([vec(r) for r in net]) == ('plane',)
        if not nontriv:
            reason = 'T'
        elif len(srcs) <= 2:
            reason = 'S2'
        elif len(srcs) == 3:
            reason = 'S3'
        else:
            d = analyse4(net)
            Q = d['Q']
            nz = False
            if Q != 0:
                PQ = sp.Poly(Q, tsym)
                nz = PQ.count_roots(0, 1) - (PQ.eval(0) == 0) - (PQ.eval(1) == 0) > 0
            if not nz:
                reason = 'M4'
            else:
                check(key in table1, "four-source class with a fold is in Table 1: " + net_str(net))
                if (0, 0) in srcs:
                    reason = 'F-in'
                elif key in orange:
                    reason = 'YES'
                else:
                    arg = origin_never_AS(net)
                    check(arg is not None, "origin never AS for " + net_str(net) + ": " + str(arg))
                    reason = 'F-un'
                    fun.append((net_str(net), arg))
        counts[reason] = counts.get(reason, 0) + 1
        rows.append([net_str(net), len(mem), " ; ".join(net_str(m) for m in mem[1:]), len(srcs),
                     int(nontriv), reason])
    print("reasons:", counts)
    for f in fun:
        print("   F-un:", f[0], "|", f[1])
    check(counts.get('YES', 0) == 10, "exactly 10 classes admit two AS equilibria")
    check(counts.get('F-in', 0) == 15 and counts.get('F-un', 0) == 5, "Table 1: 15 with inflow, 5 with unstable origin")
    check(sum(counts.values()) == 5716, "every class classified once")
    yes = set(canon_class(parse_net(r[0])) for r in rows if r[5] == 'YES')
    check(yes == orange, "the YES classes are the 10 orange networks of BBH Table 1")
    out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "certificates", "classes.csv")
    with open(out, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["representative", "n_networks_in_class", "other_members", "n_sources", "dyn_nontrivial", "reason"])
        w.writerows(rows)
    print("written certificates/classes.csv (%d rows)" % len(rows))
    print("SUMMARY:", "ALL PASS" if not FAILS else "FAILURES: %d" % len(FAILS))
    sys.exit(1 if FAILS else 0)


if __name__ == "__main__":
    main()
