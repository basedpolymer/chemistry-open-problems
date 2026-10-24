#!/usr/bin/env python3
"""Exact / ball-arithmetic checks of the geometric facts used in result.md.

G1 (genericity is a.e.): the 15 orientation determinants det(v_b-v_a, v_c-v_a, v_d-v_a) of the T135
   action-angle hexagon are real-analytic on int(P) x T^3; each is certified NONZERO (python-flint arb,
   200 bits) at the explicit point d = (1, 9/10, 4/5), theta = (3/10, 1/2, 7/10), so each zero set is null.
G2 (equilateral): at that point all six edges have length 1 (ball contains 1, radius < 1e-50) and the
   diagonals |v1v3|, |v3v5|, |v5v1| equal d1, d2, d3.
G3 (symmetries, Section 2.4): at 200 random rational points, the hexagon alpha(d2,d3,d1,th2,th3,th1) is
   congruent to the relabelled hexagon (v3,v4,v5,v6,v1,v2) of alpha(d,th) (all 15 pairwise distances and all
   15 orientation determinants agree), and alpha(d3,d2,d1,th3,th2,th1) is congruent by an orientation-
   REVERSING isometry to (v1,v6,v5,v4,v3,v2) (distances agree, determinants change sign); and
   alpha(d,-th) is the mirror image of alpha(d,th).  (ball arithmetic, all differences certified < 1e-40)
G4 (coordinates): the (d1,s,t) chart d2 = d1(1+s)/2, d3 = d1((1-s)/2+ts) maps [0,2]x[0,1]^2 onto the
   chamber {2 >= d1 >= d2 >= d3 >= 0, d1 <= d2+d3} with Jacobian d1^2 s/2 (sympy, exact), the chamber has
   volume 2/3 = Vol(P)/6, and 4 r3^2-(1-r2^2+r3^2)^2 = t s (1-s^2(1-t)^2)(2+ts) (sympy, exact).
Environment: Python 3.13, python-flint 0.9, sympy 1.14.  Runtime: ~5 s.
"""
import sys, itertools, random
import flint
from flint import arb
import sympy as sp

flint.ctx.prec = 200

def hexagon(d1, d2, d3, t1, t2, t3):
    D = (2*d1**2*d2**2 + 2*d1**2*d3**2 + 2*d2**2*d3**2 - d1**4 - d2**4 - d3**4).sqrt()
    v1 = (arb(0), arb(0), arb(0)); v3 = (d1, arb(0), arb(0))
    v5 = ((d1**2 - d2**2 + d3**2) / (2*d1), D / (2*d1), arb(0))
    V = [v1, None, v3, None, v5, None]
    for k, (a, b, dd, t) in enumerate([(v1, v3, d1, t1), (v3, v5, d2, t2), (v5, v1, d3, t3)]):
        mx = (a[0] + b[0]) / 2; my = (a[1] + b[1]) / 2
        nx = -(b[1] - a[1]) / dd; ny = (b[0] - a[0]) / dd
        h = (1 - dd*dd/4).sqrt()
        V[2*k+1] = (mx + h*t.cos()*nx, my + h*t.cos()*ny, h*t.sin())
    return V

def det(V, a, b, c, d):
    u = [V[b][i] - V[a][i] for i in range(3)]; v = [V[c][i] - V[a][i] for i in range(3)]; w = [V[d][i] - V[a][i] for i in range(3)]
    return u[0]*(v[1]*w[2]-v[2]*w[1]) - u[1]*(v[0]*w[2]-v[2]*w[0]) + u[2]*(v[0]*w[1]-v[1]*w[0])

def dist2(V, a, b):
    return sum((V[a][i] - V[b][i])**2 for i in range(3))

def small(x, eps=arb("1e-40")):
    return abs(x) < eps

def main():
    ok = True
    def chk(name, cond):
        nonlocal ok
        print(("PASS" if cond else "FAIL"), name); ok = ok and cond
    Q = lambda p, q: arb(p) / q
    V = hexagon(arb(1), Q(9, 10), Q(4, 5), Q(3, 10), Q(1, 2), Q(7, 10))
    quads = list(itertools.combinations(range(6), 4))
    nz = [det(V, *q) for q in quads]
    chk("G1 all 15 orientation determinants nonzero at the test point", all(x > 0 or x < 0 for x in nz))
    edges = [dist2(V, i, (i+1) % 6) for i in range(6)]
    chk("G2 six unit edges", all(small(e - 1, arb("1e-50")) for e in edges))
    chk("G2 diagonals d1,d2,d3", small(dist2(V, 0, 2) - 1) and small(dist2(V, 2, 4) - Q(81, 100)) and small(dist2(V, 4, 0) - Q(16, 25)))
    rnd = random.Random(7)
    good_c = good_t = good_m = True
    for _ in range(200):
        while True:
            d = [Q(rnd.randint(1, 199), 100) for _ in range(3)]
            if all(d[i] < d[(i+1) % 3] + d[(i+2) % 3] for i in range(3)):
                break
        th = [Q(rnd.randint(1, 628), 100) for _ in range(3)]
        H = hexagon(*d, *th)
        # cyclic
        Hc = hexagon(d[1], d[2], d[0], th[1], th[2], th[0])
        lab = [2, 3, 4, 5, 0, 1]
        for (a, b) in itertools.combinations(range(6), 2):
            if not small(dist2(Hc, a, b) - dist2(H, lab[a], lab[b])): good_c = False
        for q in quads:
            if not small(det(Hc, *q) - det(H, *[lab[i] for i in q])): good_c = False
        # transposition (orientation reversing)
        Ht = hexagon(d[2], d[1], d[0], th[2], th[1], th[0])
        lab = [0, 5, 4, 3, 2, 1]
        for (a, b) in itertools.combinations(range(6), 2):
            if not small(dist2(Ht, a, b) - dist2(H, lab[a], lab[b])): good_t = False
        for q in quads:
            if not small(det(Ht, *q) + det(H, *[lab[i] for i in q])): good_t = False
        # mirror
        Hm = hexagon(*d, -th[0], -th[1], -th[2])
        for i in range(6):
            if not (small(Hm[i][0] - H[i][0]) and small(Hm[i][1] - H[i][1]) and small(Hm[i][2] + H[i][2])): good_m = False
    chk("G3 cyclic relabelling = coordinate permutation (proper congruence)", good_c)
    chk("G3 reversal relabelling = coordinate transposition (improper congruence)", good_t)
    chk("G3 theta -> -theta is the mirror image z -> -z", good_m)
    # G4 exact symbolic checks
    d1, s, t = sp.symbols("d1 s t", positive=True)
    r2 = (1 + s) / 2; r3 = (1 - s) / 2 + t * s
    J = sp.Matrix([d1, d1 * r2, d1 * r3]).jacobian([d1, s, t]).det()
    chk("G4 Jacobian = d1^2 s / 2", sp.simplify(J - d1**2 * s / 2) == 0)
    chk("G4 chamber volume = 2/3", sp.integrate(d1**2 * s / 2, (d1, 0, 2), (s, 0, 1), (t, 0, 1)) == sp.Rational(2, 3))
    chk("G4 factorisation of 4 r3^2 - (1 - r2^2 + r3^2)^2", sp.expand(4*r3**2 - (1 - r2**2 + r3**2)**2 - t*s*(1 - s**2*(1 - t)**2)*(2 + t*s)) == 0)
    # chamber description: d1 - d2 = d1(1-s)/2 >= 0, d2 - d3 = d1 s (1-t) >= 0, d2 + d3 - d1 = d1 t s >= 0
    chk("G4 d1-d2, d2-d3, d2+d3-d1 = d1(1-s)/2, d1 s(1-t), d1 t s", sp.simplify(d1 - d1*r2 - d1*(1-s)/2) == 0 and sp.simplify(d1*r2 - d1*r3 - d1*s*(1-t)) == 0 and sp.simplify(d1*r2 + d1*r3 - d1 - d1*t*s) == 0)
    sys.exit(0 if ok else 1)

if __name__ == "__main__":
    main()
