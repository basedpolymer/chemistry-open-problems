"""Verifier V3 (Claude Opus 5.5, effort xhigh), 27 October 2026: exact re-checks of the explicit examples used as
counterexamples or certificates by agents G, E, Q1, E2, D2.  Independent of the authors' scripts (only v3lib.py,
written here, and sympy are used).  Every claim prints PASS/FAIL; exit status 1 on any FAIL.

Environment: Python 3.12 (WSL Fedora 44), sympy 1.14; exact rational arithmetic throughout.
Usage (from verification/V3/scripts/): python v3_examples.py
"""
import sys
import time
import platform
import random
from fractions import Fraction as F
import sympy as sp

import v3lib as L

FAILS = 0


def claim(name, ok, msg=""):
    global FAILS
    print(("PASS" if ok else "FAIL") + ": " + name + ((" -- " + msg) if msg else ""))
    if not ok:
        FAILS += 1


def net(n, txt_reactions):
    """txt_reactions: list of (dict src, dict prod) with 1-based species indices."""
    R = []
    for s, p in txt_reactions:
        R.append((tuple(s.get(i + 1, 0) for i in range(n)), tuple(p.get(i + 1, 0) for i in range(n))))
    return R


def ode(n, R, kap, x):
    f = [0] * n
    for j, (s, p) in enumerate(R):
        rate = kap[j]
        for i in range(n):
            rate *= x[i] ** s[i]
        for i in range(n):
            f[i] += rate * (p[i] - s[i])
    return [sp.expand(e) for e in f]


def same_up_to_relabel(n, R1, R2):
    return L.canonical(n, R1) == L.canonical(n, R2)


def H123(n, R):
    return L.conservative(n, R), L.consistent(n, R), L.H3(n, R)


def drain_species(n, R):
    """E's Prop. E.4 hypotheses: species i, reaction k0 with source d0*e_i consuming X_i, every reaction producing
    X_i has X_i-coefficient > d0 in its source.  Returns list of (i, d0)."""
    out = []
    for i in range(n):
        for s, p in R:
            if s[i] >= 1 and all(s[l] == 0 for l in range(n) if l != i) and p[i] < s[i]:
                d0 = s[i]
                if all(s2[i] > d0 for s2, p2 in R if p2[i] > s2[i]):
                    out.append((i, d0))
    return sorted(set(out))


def is_drain_set(n, R, Z, w, d):
    """Q1's Theorem D hypotheses (D1), (D2) for species set Z (list), weights w, degrees d (dicts)."""
    def dV(s, p):
        return sum(F(w[i]) * (p[i] - s[i]) for i in Z)
    for i in Z:
        ok = any(s[i] == d[i] and all(s[l] == 0 for l in range(n) if l != i) and dV(s, p) < 0 for s, p in R)
        if not ok:
            return False
    for s, p in R:
        if dV(s, p) > 0 and not sum(F(s[i], d[i]) for i in Z) > 1:
            return False
    return True


def main():
    t0 = time.time()
    print("# v3_examples.py  Python %s on %s; sympy %s; exact" % (platform.python_version(), platform.system(),
                                                                   sp.__version__))
    # =========================================================== agent G
    print("## agent G")
    M = net(2, [({1: 3}, {1: 2, 2: 1}), ({1: 2, 2: 1}, {1: 3})])
    nD, D = L.dual(2, M)
    Dstated = [((3, 2), (2, 3)), ((0, 1), (1, 0))]
    claim("G3.dual  D = M* = 3U+2V -> 2U+3V, V -> U", D == Dstated)
    claim("G3.Gamma both [[-1,1],[1,-1]]", L.Gamma(2, M) == [[-1, 1], [1, -1]] == L.Gamma(2, D))
    claim("G3.H  (H1)-(H3) for M and for D", H123(2, M) == (True, True, True) == H123(2, D))
    claim("G3.counts  (m,l,s,delta): M (2,1,1,0), D (4,2,1,1)", L.counts(2, M) == (2, 1, 1, 0) and
          L.counts(2, D) == (4, 2, 1, 1))
    claim("G3.WR  M reversible/WR, D not WR", L.is_wr(M) and not L.is_wr(D))
    # D not endotactic: w=(0,1): both reactions essential, max source height 2 at 3U+2V whose reaction increases V
    w = (0, 1)
    ess = [(s, p) for s, p in D if sum(w[i] * (p[i] - s[i]) for i in range(2)) != 0]
    hmax = max(sum(w[i] * s[i] for i in range(2)) for s, p in ess)
    bad = [(s, p) for s, p in ess if sum(w[i] * s[i] for i in range(2)) == hmax and
           sum(w[i] * (p[i] - s[i]) for i in range(2)) > 0]
    claim("G3.endo  w=(0,1) witnesses that D is not endotactic", len(ess) == 2 and hmax == 2 and len(bad) == 1)
    u, v, T, k1, k2 = sp.symbols('u v T k1 k2', positive=True)
    fD = ode(2, D, [1, 1], [u, v])
    claim("G3.ode  unit rates: v' = u^3 v^2 - v", sp.expand(fD[1] - (u**3 * v**2 - v)) == 0)
    h = (2 - v)**3 * v
    claim("G3.h'  h'(v) = (2-v)^2 (2-4v)", sp.expand(sp.diff(h, v) - (2 - v)**2 * (2 - 4 * v)) == 0)
    claim("G3.h(1/8) = 3375/4096 and 1 - h(1/8) = 721/4096",
          h.subs(v, sp.Rational(1, 8)) == sp.Rational(3375, 4096) and 1 - sp.Rational(3375, 4096) == sp.Rational(721, 4096))
    claim("G3.eq  (1,1) is an equilibrium of the class T=2", [e.subs({u: 1, v: 1}) for e in fD] == [0, 0])
    roots = sp.Poly(sp.expand(h - 1), v).real_roots()
    inside = [r for r in roots if 0 < r < 2]
    claim("E.G3  class u+v=2 has exactly two positive equilibria (v=1 and one root in (1/8,1/2))",
          len(inside) == 2 and sp.Integer(1) in inside and
          any(sp.Rational(1, 8) < r < sp.Rational(1, 2) for r in inside), str([sp.N(r, 8) for r in inside]))
    fDk = ode(2, D, [k1, k2], [u, v])
    claim("E.G3b  all rates: v' = v (k1 u^3 v - k2)", sp.expand(fDk[1] - v * (k1 * u**3 * v - k2)) == 0)
    # bound: for 0<v<=vstar=min(T, k2/(2 k1 T^3)), u<=T: k1 u^3 v <= k1 T^3 vstar <= k2/2
    claim("E.G3b  k1 T^3 * (k2/(2 k1 T^3)) = k2/2 (so v' <= -(k2/2) v below v*)",
          sp.simplify(k1 * T**3 * (k2 / (2 * k1 * T**3)) - k2 / 2) == 0)
    a, b = sp.symbols('a b', positive=True)
    fM = ode(2, M, [k1, k2], [a, b])
    claim("G3.M  a' = a^2 (k2 b - k1 a) for M (single positive equilibrium on each segment, global attraction)",
          sp.expand(fM[0] - a**2 * (k2 * b - k1 * a)) == 0)
    # 0 <-> A and A <-> B
    Z0 = [((0,), (1,)), ((1,), (0,))]
    n0, Z0d = L.dual(1, Z0)
    claim("G2.1(2)  dual of 0<->A is V -> U (conservative)", Z0d == [((0, 1), (1, 0))] and L.conservative(n0, Z0d)
          and not L.conservative(1, Z0))
    AB = net(2, [({1: 1}, {2: 1}), ({2: 1}, {1: 1})])
    claim("G2.1  A<->B: same Gamma as M, dual counts (2,1,1,0)", L.Gamma(2, AB) == L.Gamma(2, M) and
          L.counts(*L.dual(2, AB)) == (2, 1, 1, 0))
    # Prop 4.1, five-species example
    R5 = []
    for mult in (1, 2):
        for i in range(5):
            s = [0] * 5
            p = [0] * 5
            s[i] += mult
            s[(i + 2) % 5] += mult
            p[(i + 1) % 5] += mult
            p[(i + 3) % 5] += mult
            R5.append((tuple(s), tuple(p)))
    G5 = L.Gamma(5, R5)
    claim("G4.1a  5-species: 10 reactions, (m,l,s,delta) = (10,2,4,4), q = 1, WR, H3",
          len(R5) == 10 and L.counts(5, R5) == (10, 2, 4, 4) and L.is_wr(R5) and L.H3(5, R5))
    n5d, R5d = L.dual(5, R5)
    claim("G4.1b  dual is WR, one linkage class, deficiency 0 (a single 5-cycle)",
          L.is_wr(R5d) and L.counts(n5d, R5d)[1] == 1 and L.counts(n5d, R5d)[3] == 0 and
          len(L.complexes(R5d)) == 5)
    spans = []
    for lc in L.linkage_classes(R5):
        cols = [[p[i] - s[i] for i in range(5)] for s, p in R5 if s in lc]
        spans.append(L.rank(cols))
    claim("G4.1c  each linkage class spans S (rank 4)", spans == [4, 4])
    Sig = (0, 1, 2)
    claim("G4.1d  {X1,X2,X3} is a siphon containing no semiflow support",
          L.is_siphon(5, R5, Sig) and not L.semipositive_in_kernel_rows(G5, Sig))
    claim("G4.1e  realizable at (0,0,0,2,3) in the class of (1,1,1,1,1); codim rank Gamma[Sigma,:] = 3",
          sum((0, 0, 0, 2, 3)) == 5 and L.rank([G5[i] for i in Sig]) == 3)
    J0 = [j for j, (s, p) in enumerate(R5) if all(s[i] == 0 for i in Sig)]
    S0dim = 4 - L.rank([G5[i] for i in Sig])
    claim("G4.1f  no blind reaction (J0 empty), dim S0 = 1 (Wiuf (a) fails)", J0 == [] and S0dim == 1)
    # randomized corroboration of strong endotacticity (GMS Def. 3.14): w integer, not orthogonal to S
    rnd = random.Random(20261027)
    okSE = True
    for _ in range(3000):
        w = [rnd.randint(-6, 6) for _ in range(5)]
        if len(set(w)) == 1:
            continue  # orthogonal to S (S = sum-zero space)
        hts = [sum(w[i] * s[i] for i in range(5)) for s, p in R5]
        dl = [sum(w[i] * (p[i] - s[i]) for i in range(5)) for s, p in R5]
        ess = [j for j in range(10) if dl[j] != 0]
        hm = max(hts[j] for j in ess)
        endo = all(dl[j] < 0 for j in ess if hts[j] == hm)
        Hm = max(hts)
        strong = any(dl[j] < 0 and hts[j] == Hm for j in range(10))
        if not (endo and strong):
            okSE = False
            break
    claim("G4.1g  strongly endotactic test (GMS Def 3.14) passes for 3000 random integer w (corroboration only)",
          okSE)
    # four-species counterexample to dropping primality
    R4 = net(4, [({1: 1}, {2: 1}), ({2: 1}, {3: 1}), ({3: 1}, {4: 1}), ({4: 1}, {1: 1}),
                 ({1: 1, 3: 1}, {2: 1, 4: 1}), ({2: 1, 4: 1}, {1: 1, 3: 1})])
    q4 = 4 - L.rank(L.Gamma(4, R4))
    n4d, R4d = L.dual(4, R4)
    w = (3, 2, 1, 2)
    hts = [sum(w[i] * s[i] for i in range(4)) for s, p in R4]
    dl = [sum(w[i] * (p[i] - s[i]) for i in range(4)) for s, p in R4]
    Hm = max(hts)
    claim("G4.1np  A->B->C->D->A, A+C<->B+D: q=1, H3, WR, dual WR; w=(3,2,1,2) not orthogonal to S, no strictly "
          "inward reaction at the maximal reactant height (not strongly endotactic)",
          q4 == 1 and L.H3(4, R4) and L.is_wr(R4) and L.is_wr(R4d) and len(set(w)) > 1 and Hm == 4 and
          not any(dl[j] < 0 and hts[j] == Hm for j in range(6)), "heights %s, deltas %s" % (hts, dl))
    # Section 5: H
    H = net(2, [({1: 2}, {2: 2}), ({2: 2}, {1: 2}), ({1: 2}, {1: 1, 2: 1}), ({1: 1, 2: 1}, {1: 2})])
    nHd, Hd = L.dual(2, H)
    claim("G5.dual  H* = 2U1+2U3+U4 -> 2U2+U3+2U4, 2U2+U4 -> 2U1+U3",
          Hd == [((2, 0, 2, 1), (0, 2, 1, 2)), ((0, 2, 0, 1), (2, 0, 1, 0))])
    claim("G5.H  H is WR with (H1)-(H3), minimal siphons = [{A,B}] containing the support of a+b",
          L.is_wr(H) and H123(2, H) == (True, True, True) and L.minimal_siphons(2, H) == [frozenset({0, 1})]
          and L.semipositive_in_kernel_rows(L.Gamma(2, H), (0, 1)))
    GHd = L.Gamma(nHd, Hd)
    claim("G5.crit  {U4} is a siphon of H*, row 4 of Gamma* = (1,-1), H* has (H1)-(H3)",
          L.is_siphon(nHd, Hd, (3,)) and GHd[3] == [1, -1] and H123(nHd, Hd) == (True, True, True))
    # Example 8.1
    E8 = net(3, [({1: 3, 3: 2}, {1: 2, 2: 3}), ({1: 2, 2: 3}, {2: 2, 3: 3}), ({2: 2, 3: 3}, {1: 3, 3: 2})])
    G8 = L.Gamma(3, E8)
    claim("G8.1a  Gamma = [[-1,-2,3],[3,-1,-2],[-2,3,-1]], rank 2, delta 0, flux (1,1,1), law a+b+c",
          G8 == [[-1, -2, 3], [3, -1, -2], [-2, 3, -1]] and L.counts(3, E8) == (3, 1, 2, 0) and
          all(sum(r) == 0 for r in G8) and all(sum(G8[i][j] for i in range(3)) == 0 for j in range(3)))
    Ye8 = L.Ye(3, E8)
    claim("G8.1b  support condition (Ye)_ij > 0 iff Gamma_ij < 0", all((Ye8[i][j] > 0) == (G8[i][j] < 0)
                                                                       for i in range(3) for j in range(3)))
    n8d, E8d = L.dual(3, E8)
    claim("G8.1c  dual WR with deficiency 0", L.is_wr(E8d) and L.counts(n8d, E8d)[3] == 0)
    alt = net(3, [({1: 10, 3: 2}, {1: 9, 2: 3}), ({1: 2, 2: 1}, {3: 3}), ({2: 10, 3: 1}, {1: 3, 2: 8})])
    claim("G8.1d  alternative realisation has the same Gamma", L.Gamma(3, alt) == G8)
    xa, xb, xc = sp.symbols('xa xb xc', positive=True)
    g = [xa**10 * xc**2, xa**2 * xb, xb**10 * xc]
    Gm = sp.Matrix([[sp.diff(gi, xx).subs({xa: 1, xb: 1, xc: 1}) for xx in (xa, xb, xc)] for gi in g])
    claim("G8.1e  G = Dg(1,1,1) = [[10,0,2],[2,1,0],[0,10,1]]", Gm == sp.Matrix([[10, 0, 2], [2, 1, 0], [0, 10, 1]]))
    J = sp.Matrix(G8) * Gm
    B = sp.Matrix([[1, 0], [0, 1], [-1, -1]])
    # J maps S into S; coordinates: J*B = B*JS, solve JS = (B^T B)^{-1} B^T J B
    JS = (B.T * B).inv() * B.T * J * B
    claim("G8.1f  J_S = [[-15,27],[24,-25]], det = -273 (positive eigenvalue)",
          J * B == B * JS and JS == sp.Matrix([[-15, 27], [24, -25]]) and JS.det() == -273, str(JS.tolist()))

    # =========================================================== agent E
    print("## agent E")
    NE = net(4, [({4: 1}, {3: 3}), ({3: 3}, {1: 1, 2: 2}), ({1: 1, 2: 2}, {2: 3}), ({2: 3}, {4: 1})])
    claim("E.ex  X4 -> 3X3 -> X1+2X2 -> 3X2 -> X4: (m,l,s,delta) = (4,1,3,0), WR, (H1)-(H3)",
          L.counts(4, NE) == (4, 1, 3, 0) and L.is_wr(NE) and H123(4, NE) == (True, True, True))
    nEd, NEd = L.dual(4, NE)
    stated = net(4, [({4: 1}, {2: 1}), ({3: 3, 4: 2}, {2: 2, 4: 3}), ({2: 3}, {1: 3}), ({1: 1}, {3: 1})])
    claim("E.ex  dual equals k4 -> k2, 3k3+2k4 -> 2k2+3k4, 3k2 -> 3k1, k1 -> k3 up to relabelling",
          same_up_to_relabel(4, NEd, stated))
    claim("E.ex  dual satisfies (H1)-(H3)", H123(nEd, NEd) == (True, True, True))
    ds = drain_species(nEd, NEd)
    claim("E.ex  the dual has a drain species (E.4 hypotheses), d0 = 1", len(ds) == 1 and ds[0][1] == 1, str(ds))
    claim("E.ex  N itself has no drain species (consistent with WR)", drain_species(4, NE) == [])
    # E.4 proof inequality on the dual: x' <= x^{d0}(C x - c0): production terms have exponent >= d0+1 in x_i
    i0, d0 = ds[0]
    prod_ok = all(s[i0] >= d0 + 1 for s, p in NEd if p[i0] > s[i0])
    claim("E.4  every producing reaction of the drain species has order >= d0+1", prod_ok)

    # =========================================================== agent Q1
    print("## agent Q1")
    k1, k2, k3, k4 = sp.symbols('k1 k2 k3 k4', positive=True)
    x1, x2, x3, x4 = sp.symbols('x1 x2 x3 x4', positive=True)
    X = net(2, [({2: 2}, {1: 1, 2: 1}), ({1: 1, 2: 1}, {2: 2})])
    nXd, Xd = L.dual(2, X)
    claim("Q1.X  N = 2X2 <-> X1+X2: WR, (H1)-(H3), delta 0; dual = k2 -> k1, 2k1+k2 -> k1+2k2",
          L.is_wr(X) and H123(2, X) == (True, True, True) and L.counts(2, X)[3] == 0 and
          Xd == [((0, 1), (1, 0)), ((2, 1), (1, 2))])
    f = ode(2, Xd, [k1, k2], [x1, x2])
    claim("Q1.X  k2' = k2 (kappa2 k1^2 - kappa1)", sp.expand(f[1] - x2 * (k2 * x1**2 - k1)) == 0)
    vX, ftX = L.verdict_rank_one(nXd, Xd)
    claim("Q1.X  face classification: NONE, face {k2} is ESC", vX == 'NONE' and
          any(Z == (1,) and t == 'ESC' for s, Z, t in ftX), str(ftX))
    claim("Q1.X  dual has no drain species", drain_species(nXd, Xd) == [])
    # drain sets of the dual of X: Z subsets of {k1,k2}; (D1) needs pure-source reaction d_i k_i for each i in Z.
    # Only pure source is k2 (reaction k2 -> k1), so Z = {k2}: (D2): reaction 2k1+k2 -> k1+2k2 increases k2 and
    # has sum y_i/d_i = 1/d (d=1): 1, not > 1.  So no drain set.
    claim("Q1.X  dual has no drain set ({k2} with d=1 fails (D2): ratio 1)",
          not is_drain_set(nXd, Xd, [1], {1: 1}, {1: 1}))
    # Example C
    C = net(3, [({3: 2}, {2: 1, 3: 1}), ({2: 1, 3: 1}, {1: 1, 3: 1}), ({1: 1, 3: 1}, {3: 2})])
    nCd, Cd = L.dual(3, C)
    claim("Q1.C  N = 2X3 -> X2+X3 -> X1+X3 -> 2X3: WR, (H1)-(H3), (m,l,s,d) = (3,1,2,0); dual as stated",
          L.is_wr(C) and H123(3, C) == (True, True, True) and L.counts(3, C) == (3, 1, 2, 0) and
          Cd == [((0, 0, 1), (0, 1, 0)), ((0, 1, 0), (1, 0, 0)), ((2, 1, 1), (1, 1, 2))])
    f = ode(3, Cd, [k1, k2, k3], [x1, x2, x3])
    claim("Q1.C  V = k2 + 2k3: V' = -kappa1 k3 - kappa2 k2 + 2 kappa3 k1^2 k2 k3",
          sp.expand(f[1] + 2 * f[2] - (-k1 * x3 - k2 * x2 + 2 * k3 * x1**2 * x2 * x3)) == 0)
    claim("Q1.C  {k2,k3} drain set with w=(1,2), d=(1,1); no drain species",
          is_drain_set(nCd, Cd, [1, 2], {1: 1, 2: 2}, {1: 1, 2: 1}) and drain_species(nCd, Cd) == [])
    # Example P
    P = net(3, [({1: 2, 3: 1}, {1: 3, 2: 2}), ({1: 3, 2: 2}, {1: 2, 3: 1}),
                ({2: 2, 3: 1}, {1: 1, 2: 4}), ({1: 1, 2: 4}, {2: 2, 3: 1})])
    nPd, Pd = L.dual(3, P)
    vP, ftP = L.verdict_rank_one(nPd, Pd)
    esc = [(s, Z) for s, Z, t in ftP if t == 'ESC']
    single_ok = all(t == 'IN' for s, Z, t in ftP if len(Z) == 1)
    claim("Q1.P  WR, s=1, (H1)-(H3); dual: all singleton faces IN, the face {k1,k3} (forward reactions) is ESC, "
          "verdict NONE", L.is_wr(P) and H123(3, P) == (True, True, True) and L.rank(L.Gamma(3, P)) == 1 and
          vP == 'NONE' and single_ok and [Z for s, Z in esc] == [(0, 2)], str(ftP))
    # Example S
    S = net(3, [({1: 1, 2: 1}, {2: 2, 3: 2}), ({1: 1, 2: 1, 3: 2}, {1: 2})])
    nSd, Sd = L.dual(3, S)
    claim("Q1.S  (H1)-(H3), s=1, not WR; dual = k1+k2 -> 2k2, k1+k2 -> 2k1, 2k2 -> 2k1",
          H123(3, S) == (True, True, True) and L.rank(L.Gamma(3, S)) == 1 and not L.is_wr(S) and
          Sd == [((1, 1), (0, 2)), ((1, 1), (2, 0)), ((0, 2), (2, 0))])
    f = ode(2, Sd, [k1, k2, k3], [x1, x2])
    claim("Q1.S  k1' = k2 ((kappa2 - kappa1) k1 + 2 kappa3 k2)",
          sp.expand(f[0] - x2 * ((k2 - k1) * x1 + 2 * k3 * x2)) == 0)
    vS, ftS = L.verdict_rank_one(nSd, Sd)
    doms = [(s, Z) for s, Z, t in ftS if t == 'DOM']
    claim("Q1.S  dual: DOM faces only at one end, all others IN (Theorem S1(3) case)", vS == 'DOM-only' and
          len({s for s, Z in doms}) == 1, str(ftS))
    # Example B
    Bn = net(2, [({1: 1, 2: 1}, {2: 2}), ({1: 1, 2: 1}, {1: 2})])
    nBd, Bd = L.dual(2, Bn)
    f = ode(2, Bd, [k1, k2], [x1, x2])
    vB, ftB = L.verdict_rank_one(nBd, Bd)
    claim("Q1.B  k1' = (kappa2 - kappa1) k1 k2; DOM faces at both ends", sp.expand(f[0] - (k2 - k1) * x1 * x2) == 0
          and vB == 'DOM-only' and len({s for s, Z, t in ftB if t == 'DOM'}) == 2, str(ftB))
    # Example Q
    Qn = net(3, [({1: 1}, {2: 1}), ({1: 1}, {3: 1}), ({2: 1, 3: 1}, {1: 2})])
    nQd, Qd = L.dual(3, Qn)
    f = ode(3, Qd, [k1, k2, k3], [x1, x2, x3])
    claim("Q1.Q  pure, (H1)-(H3), s=2; dual = k1+k2 -> 2k3, k3 -> k1, k3 -> k2; (k2-k1)' = (kappa3-kappa2) k3",
          H123(3, Qn) == (True, True, True) and L.rank(L.Gamma(3, Qn)) == 2 and
          Qd == [((1, 1, 0), (0, 0, 2)), ((0, 0, 1), (1, 0, 0)), ((0, 0, 1), (0, 1, 0))] and
          sp.expand(f[1] - f[0] - (k3 - k2) * x3) == 0)
    # C-pure needs (H1): N = 0 <-> A is pure, WR, consistent, s = 1, but not conservative; N* = k2 -> k1 loses k2
    Zp = [((0,), (1,)), ((1,), (0,))]
    nZp, Zpd = L.dual(1, Zp)
    vZp, ftZp = L.verdict_rank_one(nZp, Zpd)
    claim("V3.Cpure  0 <-> A: pure, WR, (H2), s=1, not (H1); dual k2 -> k1 is conservative and NONE (OUT face): "
          "C-pure needs (H1)", L.is_wr(Zp) and L.consistent(1, Zp) and not L.conservative(1, Zp) and
          Zpd == [((0, 1), (1, 0))] and L.conservative(nZp, Zpd) and vZp == 'NONE' and
          any(t == 'OUT' for _, _, t in ftZp), str(ftZp))
    # V3's own example: all faces IN although the network is not W_I-endotactic (Zhang-Gao-Dochain Def. 3.1)
    Zx = net(2, [({1: 1, 2: 1}, {1: 2}), ({1: 2}, {1: 1, 2: 1}), ({1: 3, 2: 1}, {1: 4})])
    vZ, ftZ = L.verdict_rank_one(2, Zx)
    # siphons and w_I-endotacticity
    sips = [Zs for Zs in ((0,), (1,), (0, 1)) if L.is_siphon(2, Zx, Zs)]

    def w_endotactic(R, w):
        ess = [(s, p) for s, p in R if sum(w[i] * (p[i] - s[i]) for i in range(len(w))) != 0]
        if not ess:
            return True
        hs = [sum(w[i] * s[i] for i in range(len(w))) for s, p in ess]
        ds_ = [sum(w[i] * (p[i] - s[i]) for i in range(len(w))) for s, p in ess]
        return all(d < 0 for h_, d in zip(hs, ds_) if h_ == max(hs)) and \
            all(d > 0 for h_, d in zip(hs, ds_) if h_ == min(hs))
    wi = all(w_endotactic(Zx, [1 if i in Zs else 0 for i in range(2)]) for Zs in sips)
    a_, b_ = sp.symbols('a b', positive=True)
    f = ode(2, Zx, [k1, k2, k3], [a_, b_])
    claim("V3.ZGD  A+B -> 2A, 2A -> A+B, 3A+B -> 4A: conservative rank one, every face IN (S1: persistent for "
          "all rates), siphons %s, but NOT W_I-endotactic (ZGD Thm 3.10 does not apply)" % sips,
          vZ == 'ALL' and not wi and L.conservative(2, Zx) and
          sp.expand(f[0] - a_ * (k1 * b_ - k2 * a_ + k3 * a_**2 * b_)) == 0, str(ftZ))

    # =========================================================== agent E2
    print("## agent E2")
    E27 = net(3, [({3: 2}, {2: 1}), ({2: 1}, {1: 1, 3: 1}), ({1: 1, 3: 1}, {3: 2})])
    n27, E27d = L.dual(3, E27)
    claim("E2.7ex  N = 2X3 -> X2 -> X1+X3 -> 2X3; dual = k3 -> k2, k2 -> k1, 2k1+k3 -> k2+2k3",
          E27d == [((0, 0, 1), (0, 1, 0)), ((0, 1, 0), (1, 0, 0)), ((2, 0, 1), (0, 1, 2))] and
          H123(3, E27) == (True, True, True) and L.is_wr(E27))
    f = ode(3, E27d, [k1, k2, k3], [x1, x2, x3])
    Tt = sp.symbols('T', positive=True)
    # vertex (T,0,0) in coordinates (x2, x3), x1 = T - x2 - x3
    sub = {x1: Tt - x2 - x3}
    Jv = sp.Matrix([[sp.diff(f[i].subs(sub), xx).subs({x2: 0, x3: 0}) for xx in (x2, x3)] for i in (1, 2)])
    ev = sorted(Jv.eigenvals().keys(), key=str)
    claim("E2.7ex  linearisation at (T,0,0): eigenvalues -kappa2 and kappa3 T^2 - kappa1",
          set(Jv.eigenvals().keys()) == {-k2, k3 * Tt**2 - k1}, str(Jv.tolist()))
    # Q1-6 refutation
    Q6 = net(4, [({4: 1}, {2: 1, 3: 1}), ({2: 1, 3: 1}, {4: 1}), ({3: 3}, {1: 1, 4: 2}), ({1: 1, 4: 2}, {3: 3})])
    pure = all(not (s[i] > 0 and p[i] > 0) for s, p in Q6 for i in range(4))
    nq6, Q6d = L.dual(4, Q6)
    statedQ6 = net(4, [({4: 1}, {2: 1}), ({3: 1}, {1: 1}), ({2: 3, 3: 1}, {1: 1, 4: 3}), ({1: 1, 4: 2}, {2: 2, 3: 1})])
    claim("E2.Q1-6  N = X4 <-> X2+X3, 3X3 <-> X1+2X4: pure, WR, (H1)-(H3), s=2, delta 0; dual as stated (relabelled)",
          pure and L.is_wr(Q6) and H123(4, Q6) == (True, True, True) and L.counts(4, Q6) == (4, 2, 2, 0) and
          same_up_to_relabel(4, Q6d, statedQ6))
    # find a drain set {two species} in my labelling: try all pairs with weights in 1..4
    found = []
    for i in range(4):
        for j in range(i + 1, 4):
            for wi_ in range(1, 5):
                for wj in range(1, 5):
                    if is_drain_set(nq6, Q6d, [i, j], {i: wi_, j: wj}, {i: 1, j: 1}):
                        found.append(((i, j), (wi_, wj)))
    claim("E2.Q1-6  the dual has a drain set of two species with d = (1,1) (Q1-6 refuted)", len(found) > 0,
          str(found[:3]))
    # generic SOME example
    SO = net(4, [({1: 1, 3: 1, 4: 1}, {1: 2, 2: 1}), ({1: 1, 2: 2, 3: 1, 4: 1}, {2: 1, 3: 2, 4: 2}),
                 ({1: 1, 2: 1, 4: 2}, {3: 1, 4: 3})])
    nSO, SOd = L.dual(4, SO)
    statedSO = net(3, [({1: 1, 2: 1, 3: 1}, {1: 2}), ({2: 2, 3: 1}, {1: 1, 2: 1}), ({1: 1, 2: 1}, {2: 2, 3: 1}),
                       ({1: 1, 2: 1, 3: 2}, {2: 2, 3: 3})])
    vSO, ftSO = L.verdict_rank_one(nSO, SOd)
    claim("E2.S  generic SOME example: s=1, (H1), (H2), not WR; dual as stated; only non-IN face is DOM",
          L.rank(L.Gamma(4, SO)) == 1 and L.conservative(4, SO) and L.consistent(4, SO) and not L.is_wr(SO) and
          SOd == statedSO and vSO == 'DOM-only', str([x for x in ftSO if x[2] != 'IN']))
    # leading coefficient at face {k2 = 0}: k_A, k_B, k_C, k_D = rates of dual reactions 1..4
    # exact expansion: classes are k0 + t*b; compute along the face via the S1 lemma: L = sum_{j in J} c_j kappa_j ...
    g, c = L.rank_one_data(nSO, SOd)
    kap = sp.symbols('kA kB kC kD', positive=True)
    y1, y3 = sp.symbols('y1 y3', positive=True)
    Zf = (1,)
    act = [j for j in range(4) if c[j] != 0]
    order = {j: sum(SOd[j][0][i] for i in Zf) for j in act}
    mu = min(order.values())
    Jm = [j for j in act if order[j] == mu]
    Lc = sum(sp.Rational(c[j].numerator, c[j].denominator) * kap[j] * abs(g[1])**SOd[j][0][1] *
             y1**SOd[j][0][0] * y3**SOd[j][0][2] for j in Jm)
    s1 = Lc.subs({kap[0]: 1, kap[1]: 1, kap[2]: sp.Rational(3, 2), kap[3]: sp.Rational(3, 2)})
    s2 = Lc.subs({kap[0]: 3, kap[1]: 3, kap[2]: 1, kap[3]: 1})
    # sigma: the face {k2} is on the side where g_2 has sign; inward iff sigma*c_j<0; persistent iff sigma*L<0 always
    sigma = 1 if g[1] < 0 else -1
    q1 = sp.Poly(sp.expand(sigma * s1 / y1), y3)
    disc1 = q1.discriminant()
    claim("E2.S  face {k2}: sigma*L = y1*(quadratic in y3); kappa=(1,1,3/2,3/2): negative for all y3>0 "
          "(leading coeff < 0, discriminant < 0); kappa=(3,3,1,1): positive at y3=1",
          sp.expand(sigma * s1).has(y1) and q1.LC() < 0 and disc1 < 0 and (sigma * s2).subs({y1: 1, y3: 1}) > 0,
          "sigma*L(kappa1) = %s ; sigma*L(kappa2) = %s" % (sp.factor(sigma * s1), sp.factor(sigma * s2)))

    # =========================================================== agent D2
    print("## agent D2")
    U = net(3, [({3: 3}, {1: 2, 2: 1}), ({2: 1, 3: 2}, {3: 3}), ({1: 2, 2: 1}, {2: 1, 3: 2})])
    nUd, Ud = L.dual(3, U)
    claim("D2.U  N = 3X3 -> 2X1+X2 (R1), X2+2X3 -> 3X3 (R2), 2X1+X2 -> X2+2X3 (R3): one 3-cycle, WR, (H1)-(H3), s=2;"
          " dual = 2k3 -> 2k1, k2+k3 -> k1+k3, 3k1+2k2 -> 3k2+2k3",
          L.is_wr(U) and H123(3, U) == (True, True, True) and L.rank(L.Gamma(3, U)) == 2 and
          Ud == [((0, 0, 2), (2, 0, 0)), ((0, 1, 1), (1, 0, 1)), ((3, 2, 0), (0, 3, 2))])
    f = ode(3, Ud, [k1, k2, k3], [x1, x2, x3])
    claim("D2.U  k2' = -kappa2 k2 k3 + kappa3 k1^3 k2^2, k3' = -2 kappa1 k3^2 + 2 kappa3 k1^3 k2^2",
          sp.expand(f[1] - (-k2 * x2 * x3 + k3 * x1**3 * x2**2)) == 0 and
          sp.expand(f[2] - (-2 * k1 * x3**2 + 2 * k3 * x1**3 * x2**2)) == 0)
    K = k3 * x1**3
    z = x2 / x3
    zdot = sp.simplify((f[1] * x3 - x2 * f[2]) / x3**2)
    claim("D2.U  z = a/b: z' = a[(2kappa1 - kappa2) + K z (1 - 2z)] (exact identity)",
          sp.simplify(zdot - x2 * ((2 * k1 - k2) + K * z * (1 - 2 * z))) == 0)
    uu = x3 / x2
    udot = sp.simplify((f[2] * x2 - x3 * f[1]) / x2**2)
    claim("D2.U  u = b/a: u' = a[(kappa2 - 2kappa1)u^2 - K u + 2K], a' = a^2 (K - kappa2 u)",
          sp.simplify(udot - x2 * ((k2 - 2 * k1) * uu**2 - K * uu + 2 * K)) == 0 and
          sp.simplify(f[1] - x2**2 * (K - k2 * uu)) == 0)
    # kappa2 > 2 kappa1 regime bound: b' = b^2(-2 kappa1 + 2 K z^2) < 0 if K z0^2 < kappa1 with z0 <= sqrt(kappa1/(2 kappa3 T^3))
    claim("D2.U  b' = b^2 (-2 kappa1 + 2 K z^2) (identity)", sp.simplify(f[2] - x3**2 * (-2 * k1 + 2 * K * z**2)) == 0)
    # CM example
    CM = net(3, [({3: 2}, {2: 1, 3: 1}), ({2: 1, 3: 1}, {1: 1, 2: 1}), ({1: 1, 2: 1}, {3: 2})])
    nCM, CMd = L.dual(3, CM)
    claim("D2.CM  N = 2X3 -> X2+X3 -> X1+X2 -> 2X3; dual = k3 -> k2, k2+k3 -> k1+k2, 2k1+k2 -> k1+2k3",
          CMd == [((0, 0, 1), (0, 1, 0)), ((0, 1, 1), (1, 1, 0)), ((2, 1, 0), (1, 0, 2))])
    f = ode(3, CMd, [k1, k2, k3], [x1, x2, x3])
    # chart at (0,T,0): a = x1, b = x3, x2 = T - a - b
    A_, B_ = sp.symbols('A B')
    fa = sp.expand(f[0].subs({x1: A_, x3: B_, x2: Tt - A_ - B_}))
    fb = sp.expand(f[2].subs({x1: A_, x3: B_, x2: Tt - A_ - B_}))
    claim("D2.CM  a' = kappa2 T b - kappa3 T a^2 + h.o.t., b' = -(kappa1 + kappa2 T) b + 2 kappa3 T a^2 + h.o.t.",
          sp.Poly(fa, A_, B_).coeff_monomial(B_) == k2 * Tt and sp.Poly(fa, A_, B_).coeff_monomial(A_**2) == -k3 * Tt
          and sp.Poly(fa, A_, B_).coeff_monomial(A_) == 0 and
          sp.Poly(fb, A_, B_).coeff_monomial(B_) == -(k1 + k2 * Tt) and
          sp.Poly(fb, A_, B_).coeff_monomial(A_**2) == 2 * k3 * Tt and sp.Poly(fb, A_, B_).coeff_monomial(A_) == 0)
    # center manifold b = h2 a^2 + ..., reduced a' = c2 a^2 + ...
    h2 = sp.symbols('h2')
    eqh = sp.expand(fb.subs(B_, h2 * A_**2))
    h2sol = sp.solve(sp.Poly(eqh, A_).coeff_monomial(A_**2), h2)[0]
    c2 = sp.simplify(sp.Poly(sp.expand(fa.subs(B_, h2sol * A_**2)), A_).coeff_monomial(A_**2))
    claim("D2.CM  reduced coefficient on the centre manifold = kappa3 T (kappa2 T - kappa1)/(kappa1 + kappa2 T)",
          sp.simplify(c2 - k3 * Tt * (k2 * Tt - k1) / (k1 + k2 * Tt)) == 0, str(c2))
    # L example
    LL = net(3, [({3: 2}, {2: 2, 3: 1}), ({2: 2, 3: 1}, {1: 2, 3: 1}), ({1: 2, 3: 1}, {3: 2})])
    nLL, LLd = L.dual(3, LL)
    f = ode(3, LLd, [k1, k2, k3], [x1, x2, x3])
    claim("D2.L  N = 2X3 -> 2X2+X3 -> 2X1+X3 -> 2X3: V = k2 + 2k3, V' = -2kappa1 k3^2 - 2kappa2 k2^2 + "
          "2kappa3 k1^2 k2 k3", sp.expand(f[1] + 2 * f[2] - (-2 * k1 * x3**2 - 2 * k2 * x2**2 +
                                                          2 * k3 * x1**2 * x2 * x3)) == 0)
    # negative definiteness for kappa3^2 T^4 < 4 kappa1 kappa2 with x1 <= T: quadratic form in (k2,k3)
    claim("D2.L  quadratic form [[-2k2, k3 T^2],[k3 T^2, -2k1]] negative definite iff kappa3^2 T^4 < 4 kappa1 kappa2",
          sp.expand(sp.Matrix([[-2 * k2, k3 * Tt**2], [k3 * Tt**2, -2 * k1]]).det() - (4 * k1 * k2 - k3**2 * Tt**4)) == 0)
    # ALL example
    AL = net(3, [({3: 2}, {2: 1}), ({2: 1}, {1: 2, 3: 1}), ({1: 2, 3: 1}, {3: 2})])
    nAL, ALd = L.dual(3, AL)
    f = ode(3, ALd, [k1, k2, k3], [x1, x2, x3])
    claim("D2.ALL  N = 2X3 -> X2 -> 2X1+X3 -> 2X3; dual = 2k3 -> 2k2, k2 -> k1, 2k1+k3 -> k2+2k3; "
          "k3' = k3 (kappa3 k1^2 - 2 kappa1 k3)",
          ALd == [((0, 0, 2), (0, 2, 0)), ((0, 1, 0), (1, 0, 0)), ((2, 0, 1), (0, 1, 2))] and
          sp.expand(f[2] - x3 * (k3 * x1**2 - 2 * k1 * x3)) == 0)
    sips1 = [Zs for Zs in ((0,), (1,), (2,)) if L.is_siphon(3, ALd, Zs)]
    sips2 = [Zs for Zs in ((0, 1), (0, 2), (1, 2)) if L.is_siphon(3, ALd, Zs)]
    claim("D2.ALL  edge (facet) zero sets that are siphons: only {k3} (edges k1=0, k2=0 not siphons); vertex "
          "siphons: only {k2,k3} = the vertex (T,0,0)", sips1 == [(2,)] and sips2 == [(1, 2)], str((sips1, sips2)))
    # E example (zero edge flow, mu = 1)
    EE = net(3, [({3: 2}, {2: 2, 3: 1}), ({2: 2, 3: 1}, {1: 1, 2: 1, 3: 1}), ({1: 1, 2: 1, 3: 1}, {3: 2})])
    nEE, EEd = L.dual(3, EE)
    ft = [x for x in L.facet_types_any_rank(3, EE) if x[0] == (2,)]
    claim("D2.E  N = 2X3 -> 2X2+X3 -> X1+X2+X3 -> 2X3: facet {k3} of the dual is ESC and a siphon; every dual "
          "reaction has k3 in its source (zero edge flow), minimal load 1",
          ft and ft[0][1] == 'ESC' and L.is_siphon(nEE, EEd, (2,)) and all(s[2] >= 1 for s, p in EEd if
                                                                             any(p[i] != s[i] for i in range(3)))
          and min(s[2] for s, p in EEd if p[2] != s[2]) == 1, str(ft))
    print("# runtime %.1f s" % (time.time() - t0))
    print("ALL CLAIMS PASS" if FAILS == 0 else "%d FAIL" % FAILS)
    sys.exit(1 if FAILS else 0)


if __name__ == "__main__":
    main()
