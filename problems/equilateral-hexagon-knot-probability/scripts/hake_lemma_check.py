#!/usr/bin/env python3
"""Certified counterexamples to Hake's Lemmas 3.9 / 3.12 (arXiv:1810.11882v1), used in her proof of Theorem 4.1.

Hake, Lemma 3.9 (J = (1,1)) and Lemma 3.12 (J = (-1,1)), same wording: "Suppose d1, d2, and d3 are distinct
and let di > dj, dk. If J(H) = (+-1, 1) then theta_i in (0, pi) and theta_j, theta_k in (0, pi/2). Moreover, if
di > sqrt(dj^2 + dk^2) then theta_i in (pi/2, pi)."  Her coordinates are the T135 action-angle coordinates of
Section 3.1, which coincide with ours (result.md, Section 1.1; theta = 0 folds the flap onto the central
triangle, theta = pi gives the planar convex hexagon).

For each explicit rational parameter point below this script certifies with python-flint ball arithmetic
(300 bits): the six edges have length 1; all 15 orientation determinants are nonzero (generic); curl = +1
(z(v2) > 0, i.e. sin theta1 > 0); the xy-projection is a regular diagram (all signs certified) whose Fox
3-colouring count is 9, so the hexagon is KNOTTED (a trefoil); and the stated inequalities on d and theta.
Example 1: d1 > d2 > d3, knotted with curl +1, but theta2 > pi/2      (contradicts 'theta_j in (0, pi/2)').
Example 2: d1^2 > d2^2 + d3^2, knotted with curl +1, but theta1 < pi/2 (contradicts the 'Moreover' part).
Whatever the chirality convention, J(H) = (+-1, 1) for a knotted hexagon with curl +1, so each example
contradicts Lemma 3.9 or Lemma 3.12.  Hake's Theorem 4.1 integrates only over the theta-region allowed by these
lemmas, so its proof is incomplete (the stated bound (14-3pi)/192 is nevertheless true: it is far above the
certified upper bound of result.md).
Environment: Python 3.13, python-flint 0.9.  Runtime: < 1 s.
"""
import sys, itertools
from fractions import Fraction as Fr
import flint
from flint import arb

flint.ctx.prec = 300

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

def sgn(x):
    if x > 0: return 1
    if x < 0: return -1
    raise ValueError("sign not certified")

def det3(V, a, b, c, d):
    u = [V[b][i] - V[a][i] for i in range(3)]; v = [V[c][i] - V[a][i] for i in range(3)]; w = [V[d][i] - V[a][i] for i in range(3)]
    return u[0]*(v[1]*w[2]-v[2]*w[1]) - u[1]*(v[0]*w[2]-v[2]*w[0]) + u[2]*(v[0]*w[1]-v[1]*w[0])

def orient2(p, q, r):
    return (q[0]-p[0])*(r[1]-p[1]) - (q[1]-p[1])*(r[0]-p[0])

def fox3(V):
    """Fox 3-colouring count of the xy-diagram, all geometric signs certified (raises otherwise)."""
    cr = []
    for i in range(6):
        sgn(orient2(V[i], V[(i+1) % 6], V[(i+2) % 6]))       # consecutive edges not collinear in projection
    for i in range(6):
        for j in range(i+2, 6):
            if i == 0 and j == 5: continue
            p, q, r, s = V[i], V[(i+1) % 6], V[j], V[(j+1) % 6]
            o1, o2, o3, o4 = orient2(p, q, r), orient2(p, q, s), orient2(r, s, p), orient2(r, s, q)
            a1, a2, a3, a4 = sgn(o1), sgn(o2), sgn(o3), sgn(o4)
            if a1 == a2 or a3 == a4: continue
            t = o3 / (o3 - o4); u = o1 / (o1 - o2)
            zi = p[2] + t*(q[2]-p[2]); zj = r[2] + u*(s[2]-r[2])
            if sgn(zi - zj) > 0: cr.append((i + t, j + u))     # (over parameter, under parameter)
            else: cr.append((j + u, i + t))
    m = len(cr)
    if m == 0: return 3
    allp = [x for c in cr for x in c]
    for a, b in itertools.combinations(allp, 2):
        sgn(a - b)                                             # all crossing parameters certified distinct
    mids = lambda x: float(x.mid())
    us = sorted(mids(c[1]) for c in cr)
    def arc(x): return (sum(1 for u in us if u < x) - 1) % m
    M = [[0]*m for _ in range(m)]
    for k, (ov, un) in enumerate(cr):
        x = mids(un); idx = sum(1 for u in us if u < x)
        for a in (arc(mids(ov)), idx % m, (idx - 1) % m):
            M[k][a] = (M[k][a] + 2) % 3
    rank = 0
    for col in range(m):
        piv = next((r for r in range(rank, m) if M[r][col] % 3), None)
        if piv is None: continue
        M[rank], M[piv] = M[piv], M[rank]
        inv = 1 if M[rank][col] == 1 else 2
        M[rank] = [(x * inv) % 3 for x in M[rank]]
        for r in range(m):
            if r != rank and M[r][col]:
                f = M[r][col]; M[r] = [(M[r][c] - f*M[rank][c]) % 3 for c in range(m)]
        rank += 1
    return 3 ** (m - rank)

def certify(name, d, th, extra):
    d1, d2, d3 = [arb(x.numerator) / x.denominator for x in d]
    t = [arb(x.numerator) / x.denominator for x in th]
    V = hexagon(d1, d2, d3, *t)
    ok = True
    edges = [sum((V[i][c] - V[(i+1) % 6][c])**2 for c in range(3)) for i in range(6)]
    ok &= all(abs(e - 1) < arb("1e-80") for e in edges)
    ok &= all(det3(V, *q) != 0 and (det3(V, *q) > 0 or det3(V, *q) < 0) for q in itertools.combinations(range(6), 4))
    ok &= V[1][2] > 0                                          # curl +1
    f = fox3(V)
    ok &= (f == 9)
    ok &= extra(d, t)
    print(("PASS" if ok else "FAIL"), name, "d =", [str(x) for x in d], "theta =", [str(x) for x in th], "Fox-3 count =", f)
    return ok

def main():
    pi = arb.pi()
    ok = True
    # Example 1: d1 > d2 > d3 (distinct), theta2 > pi/2
    ok &= certify("Example 1 (theta2 > pi/2 although d1 is the largest diagonal)",
                  (Fr(9086, 10000), Fr(8830, 10000), Fr(6928, 10000)), (Fr(451, 10000), Fr(20783, 10000), Fr(148, 10000)),
                  lambda d, t: d[0] > d[1] > d[2] and d[0] < d[1] + d[2] and t[1] > pi / 2)
    ok &= certify("Example 2 (d1^2 > d2^2 + d3^2 but theta1 < pi/2)",
                  (Fr(16039, 10000), Fr(12848, 10000), Fr(7104, 10000)), (Fr(2226, 10000), Fr(6245, 10000), Fr(4309, 10000)),
                  lambda d, t: d[0] > d[1] > d[2] and d[0] < d[1] + d[2] and d[0]**2 > d[1]**2 + d[2]**2 and t[0] < pi / 2)
    sys.exit(0 if ok else 1)
if __name__ == "__main__":
    main()
