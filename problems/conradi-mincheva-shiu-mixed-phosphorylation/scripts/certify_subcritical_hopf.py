#!/usr/bin/env python3
"""Certificate: subcritical Hopf bifurcation in the mixed-mechanism phosphorylation network.

Conradi, Mincheva, Shiu, Bull. Math. Biol. 81 (2019) 1829-1852, network (1), ODEs (2).

Checks, for the data in certificates/hopf_point.json (default: integer rate constants
k = (1,1,1,1,300,1,3,1,1,100), x1 = 1/3, x6 = 1/4, bifurcation parameter mu = x2):

  C1  the monomial map phi of the source (Theorem/Proposition 3.2) gives steady states (exact, symbolic in x1,x2,x6);
  C2  the reduced 6x6 Jacobian A(mu) (coordinates z = (x3,x4,x5,x7,x8,x9)) has characteristic polynomial
      p(lam) with lam^3 p(lam) = char. poly. of the full 9x9 Jacobian (exact, symbolic in mu);
  C3  h(mu) = det H5(mu) changes sign on the stored rational interval [a, b] (exact rationals): mu0 in [a,b];
  C4  on [a,b]: det H1..det H4 > 0, b6 > 0, h'(mu) != 0 (ball arithmetic)  [Yang / Liu criterion];
  C5  omega0^2 is identified as a root of b1 s^2 - b3 s + b5 = 0, the cofactor quartic is Hurwitz, so
      the spectrum at mu0 is {+-i omega0} plus four eigenvalues with negative real part (ball arithmetic);
  C6  first Lyapunov coefficient l1(mu0) > 0 (Kuznetsov projection formula, ball arithmetic);
  C7  crossing speed Re lambda'(mu0) > 0 (ball arithmetic);
  C8  the totals T(mu) = W x*(mu) are the stated affine functions of mu (exact), with enclosures at mu0.

Every claim prints PASS or FAIL; the exit code is non-zero if any claim fails.

Environment: Python 3.12, sympy 1.14, python-flint 0.9 (arb/acb ball arithmetic), WSL Fedora.
Runtime: about 10 seconds.
Usage: python certify_subcritical_hopf.py [path/to/hopf_point.json]
"""
import json
import os
import sys
import time

import sympy as sp
from flint import arb, acb, acb_mat, fmpq, fmpq_poly, ctx

t_start = time.time()
ctx.prec = 600
HERE = os.path.dirname(os.path.abspath(__file__))
path = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, '..', 'certificates', 'hopf_point.json')
data = json.load(open(path))
k = [sp.Rational(s) for s in data['k']]
X1 = sp.Rational(data['x1']); X6 = sp.Rational(data['x6'])
a_q = sp.Rational(data['mu_interval'][0]); b_q = sp.Rational(data['mu_interval'][1])
FAILS = []


def claim(name, ok, detail=''):
    print('%-4s %s  %s' % ('PASS' if ok else 'FAIL', name, detail))
    if not ok:
        FAILS.append(name)


# ---------------------------------------------------------------- network
# species: x1 S0, x2 K, x3 S0K, x4 S1K, x5 S2, x6 P, x7 S2P, x8 S1, x9 S1P  (source, Table 3)
REACTIONS = [((0, 1), (2,)),   # k1  S0 + K -> S0K
             ((2,), (0, 1)),   # k2  S0K -> S0 + K
             ((2,), (3,)),     # k3  S0K -> S1K
             ((3,), (4, 1)),   # k4  S1K -> S2 + K
             ((4, 5), (6,)),   # k5  S2 + P -> S2P
             ((6,), (4, 5)),   # k6  S2P -> S2 + P
             ((6,), (7, 5)),   # k7  S2P -> S1 + P
             ((7, 5), (8,)),   # k8  S1 + P -> S1P
             ((8,), (7, 5)),   # k9  S1P -> S1 + P
             ((8,), (0, 5))]   # k10 S1P -> S0 + P
NS, NR = 9, 10
Nm = sp.zeros(NS, NR)
for j, (re, pr) in enumerate(REACTIONS):
    for i in re:
        Nm[i, j] -= 1
    for i in pr:
        Nm[i, j] += 1
Wm = sp.Matrix([[0, 1, 1, 1, 0, 0, 0, 0, 0],     # Ktot
                [0, 0, 0, 0, 0, 1, 1, 0, 1],     # Ptot
                [1, 0, 1, 1, 1, 0, 1, 1, 1]])    # Stot
ZIDX = [2, 3, 4, 6, 7, 8]
E = sp.zeros(9, 6)
for j, i in enumerate(ZIDX):
    E[i, j] = 1
E[0, :] = sp.Matrix([[-1] * 6]); E[1, :] = sp.Matrix([[-1, -1, 0, 0, 0, 0]]); E[5, :] = sp.Matrix([[0, 0, 0, -1, 0, -1]])
Psel = sp.zeros(6, 9)
for j, i in enumerate(ZIDX):
    Psel[j, i] = 1
claim('C0 conservation laws W N = 0, W E = 0, P E = I, rank N = 6',
      Wm * Nm == sp.zeros(3, NR) and Wm * E == sp.zeros(3, 6) and Psel * E == sp.eye(6) and Nm.rank() == 6)

def g(kk, x):
    r = []
    for j, (re, pr) in enumerate(REACTIONS):
        v = kk[j]
        for i in re:
            v = v * x[i]
        r.append(v)
    return Nm * sp.Matrix(r)


def dg(kk, x):
    D = sp.zeros(NR, NS)
    for j, (re, pr) in enumerate(REACTIONS):
        for a_ in re:
            v = kk[j]
            for b_ in re:
                if b_ != a_:
                    v = v * x[b_]
            D[j, a_] += v
    return Nm * D


def phi(kk, x1, x2, x6):
    k1, k2, k3, k4, k5, k6, k7, k8, k9, k10 = kk
    return sp.Matrix([x1, x2, k1 / (k2 + k3) * x1 * x2, k1 * k3 / ((k2 + k3) * k4) * x1 * x2,
                      k1 * k3 * (k6 + k7) / ((k2 + k3) * k5 * k7) * x1 * x2 / x6, x6,
                      k1 * k3 / ((k2 + k3) * k7) * x1 * x2, k1 * k3 * (k9 + k10) / ((k2 + k3) * k8 * k10) * x1 * x2 / x6,
                      k1 * k3 / ((k2 + k3) * k10) * x1 * x2])

# C1: parametrization of steady states, with symbolic rate constants and symbolic x1, x2, x6
ksym = sp.symbols('k1:11', positive=True)
s1, s2, s6 = sp.symbols('s1 s2 s6', positive=True)
claim('C1 g(phi(x1,x2,x6)) = 0 identically (symbolic k, x)',
      all(sp.simplify(e) == 0 for e in g(ksym, phi(ksym, s1, s2, s6))))

mu, lam = sp.symbols('mu lam')
xs = phi(k, X1, mu, X6).applyfunc(sp.expand)
A = (Psel * dg(k, xs) * E).applyfunc(sp.expand)
p = sp.Poly((lam * sp.eye(6) - A).det(method='berkowitz'), lam)
p9 = sp.Poly((lam * sp.eye(9) - dg(k, xs)).det(method='berkowitz'), lam)
claim('C2 char.poly(9x9 Jacobian) = lam^3 * char.poly(reduced A(mu))',
      sp.expand(p9.as_expr() - lam ** 3 * p.as_expr()) == 0 and p.degree() == 6 and p.LC() == 1,
      'A(mu) is affine in mu: %s' % all(sp.Poly(e, mu).degree() <= 1 for e in A))
bco = [sp.expand(c) for c in p.all_coeffs()]          # [1, b1, ..., b6] in Q[mu]


def hurwitz(b):
    n = len(b) - 1
    def c(m):
        return b[m] if 0 <= m <= n else 0
    return [sp.expand(sp.Matrix(i, i, lambda r, s: c(2 * (r + 1) - (s + 1))).det(method='berkowitz')) for i in range(1, n)]

Hd = hurwitz(bco)          # det H1..det H5


def qpoly(e):
    cs = sp.Poly(e, mu).all_coeffs()[::-1]
    return fmpq_poly([fmpq(int(sp.Rational(c).p), int(sp.Rational(c).q)) for c in cs])

h = qpoly(Hd[4])
a_f = fmpq(int(a_q.p), int(a_q.q)); b_f = fmpq(int(b_q.p), int(b_q.q))
claim('C3 h = det H5 changes sign on [a,b] (exact rationals)', a_f < b_f and h(a_f) * h(b_f) < 0,
      'deg h = %d, b - a = %.3e' % (h.degree(), float(b_q - a_q)))
M = arb(a_f).union(arb(b_f))


def ev_arb(e, X):
    cs = sp.Poly(e, mu).all_coeffs()
    r = arb(0)
    for c in cs:
        r = r * X + arb(fmpq(int(sp.Rational(c).p), int(sp.Rational(c).q)))
    return r

bv = [ev_arb(c, M) for c in bco[1:]]
Hv = [ev_arb(e, M) for e in Hd]
dh = ev_arb(sp.diff(Hd[4], mu), M)
claim('C4 det H1..H4 > 0, b6 > 0, h\'(mu) != 0 on [a,b]',
      all(Hv[i] > 0 for i in range(4)) and bv[5] > 0 and (dh > 0 or dh < 0),
      'h\' in %s' % dh.str(5))

b1, b2, b3, b4, b5, b6 = bv
disc = b3 * b3 - 4 * b1 * b5
roots = [(b3 + disc.sqrt()) / (2 * b1), (b3 - disc.sqrt()) / (2 * b1)]
reeq = [-s ** 3 + b2 * s ** 2 - b4 * s + b6 for s in roots]
good = [i for i in range(2) if reeq[i].contains(0)]
bad = [i for i in range(2) if (reeq[i] > 0 or reeq[i] < 0)]
ok5 = disc > 0 and len(good) == 1 and len(bad) == 1
if ok5:
    s0 = roots[good[0]]
    ok5 = s0 > 0
    om = s0.sqrt()
    c1 = b1; c2 = b2 - s0; c3 = b3 - s0 * c1; c4 = b4 - s0 * c2
    ok5 = ok5 and all([c1 > 0, c1 * c2 - c3 > 0, c3 * (c1 * c2 - c3) - c1 * c1 * c4 > 0, c4 > 0])
    ok5 = ok5 and (b5 - s0 * c3).contains(0) and (b6 - s0 * c4).contains(0)
claim('C5 spectrum at mu0: simple pair +-i*omega0, other 4 eigenvalues in Re < 0', ok5,
      'omega0 in %s' % (om.str(15) if ok5 else '?'))

# ---------------------------------------------------------------- eigenvectors and l1 (ball arithmetic)
Ab = acb_mat([[acb(ev_arb(A[i, j], M)) for j in range(6)] for i in range(6)])
I6 = acb_mat([[acb(1 if i == j else 0) for j in range(6)] for i in range(6)])
Mq = Ab - acb(0, om) * I6
jq, rq, jw, rw = data['normalization']   # component fixed to 1 and dropped equation (right / left vector)


def nullvec(Mat, j, r):
    rows = [i for i in range(6) if i != r]
    cols = [c for c in range(6) if c != j]
    S = acb_mat([[Mat[i, c] for c in cols] for i in rows])
    sol = S.solve(acb_mat([[-Mat[i, j]] for i in rows]))
    it = iter(range(5))
    return [acb(1) if c == j else sol[next(it), 0] for c in range(6)]

q = nullvec(Mq, jq, rq)                # (A - i om) q = 0, exact eigenvector enclosed
w = nullvec(Mq.transpose(), jw, rw)    # w^T (A - i om) = 0
wq = sum((w[i] * q[i] for i in range(6)), acb(0))
pbar = [wi / wq for wi in w]           # pbar^T q = 1
PN = Psel * Nm
Eb = [[int(E[i, j]) for j in range(6)] for i in range(9)]
BIM = [(arb(fmpq(int(k[j].p), int(k[j].q))), re[0], re[1], [int(PN[i, j]) for i in range(6)])
       for j, (re, pr) in enumerate(REACTIONS) if len(re) == 2]


def Bf(u, v):
    """B(u,v) = P D^2 g (Eu, Ev); the vector field is quadratic, so the third derivative C vanishes."""
    Eu = [sum((u[j] * Eb[i][j] for j in range(6)), acb(0)) for i in range(9)]
    Ev = [sum((v[j] * Eb[i][j] for j in range(6)), acb(0)) for i in range(9)]
    out = [acb(0)] * 6
    for kj, i1, i2, col in BIM:
        s = acb(kj) * (Eu[i1] * Ev[i2] + Eu[i2] * Ev[i1])
        out = [out[i] + s * col[i] for i in range(6)]
    return out


def solve(Mat, v):
    sol = Mat.solve(acb_mat([[x] for x in v]))
    return [sol[i, 0] for i in range(6)]

qb = [x.conjugate() for x in q]
y1 = solve(Ab, Bf(q, qb))                         # A^{-1} B(q, qbar)
y2 = solve(acb(0, 2 * om) * I6 - Ab, Bf(q, q))     # (2 i om - A)^{-1} B(q, q)
t1 = Bf(q, y1); t2 = Bf(qb, y2)
val = sum((pbar[i] * (-2 * t1[i] + t2[i]) for i in range(6)), acb(0))
l1 = val.real / (2 * om)
claim('C6 first Lyapunov coefficient l1(mu0) > 0', l1 > 0,
      'l1 in %s (normalisation q_%d = 1)' % (l1.str(12), jq + 1))
dA = acb_mat([[acb(arb(fmpq(int(sp.Rational(sp.diff(A[i, j], mu)).p), int(sp.Rational(sp.diff(A[i, j], mu)).q))))
               for j in range(6)] for i in range(6)])
dAq = dA * acb_mat([[x] for x in q])
lamp = sum((pbar[i] * dAq[i, 0] for i in range(6)), acb(0))
claim("C7 crossing speed Re lambda'(mu0) > 0", lamp.real > 0, 'in %s' % lamp.real.str(12))

T = (Wm * xs).applyfunc(sp.expand)
Tstated = [sp.sympify(s_, locals={'mu': mu}) for s_ in data['totals_of_mu']]
claim('C8 totals T(mu) = W phi(x1, mu, x6) as stated', all(sp.expand(T[i] - Tstated[i]) == 0 for i in range(3)),
      'T(mu0) in (%s, %s, %s)' % tuple(ev_arb(T[i], M).str(15) for i in range(3)))
print('mu0 in %s' % M.str(25))
print('runtime %.1f s' % (time.time() - t_start))
if FAILS:
    print('FAILED:', FAILS)
    sys.exit(1)
print('ALL CLAIMS PASS')
