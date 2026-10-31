"""
modal_lib.py -- ball-arithmetic calculus for D_n-invariant polynomials f(p,q) on the regular n-gon P_n
(apothem 1, side 0 = {x = 1, |y| <= t}, t = tan(pi/n), R = sec(pi/n)), with p = x^2 + y^2 and
q = Re((x+iy)^n).  Scaled monomials  m_{a,k} = (p/R^2)^a (q/R^n)^k.

All integrals are exact: for a homogeneous function F_d of degree d,
   <F_d> = (2/(d+2)) * (1/t) int_0^t F_d(1,tau) dtau            (average over P_n, cone coordinates)
and the integrand is a polynomial in tau, integrated exactly by Gauss-Legendre with ball enclosures of
the nodes and weights (arb.legendre_p_root).

Environment: Python 3.12, python-flint 0.9.0.  Floating point is used only to choose coefficients.
"""
import math
import numpy as np
from flint import arb, arb_mat, acb, ctx


def ball(x):
    """exact arb from a Python float (dyadic number)"""
    return arb(float(x))


def ub(x):
    return arb(x.upper())


def lb(x):
    return arb(x.lower())


class Gauss:
    """Gauss-Legendre nodes on [0,1], weights summing to 1."""
    cache = {}

    def __init__(self, m):
        key = (m, ctx.prec)
        if key in Gauss.cache:
            self.x, self.w = Gauss.cache[key]
            return
        xs, ws = [], []
        half = (m + 1) // 2
        for i in range(half):
            r, w = arb.legendre_p_root(m, i, weight=True)
            xs.append((1 + r) / 2)
            ws.append(w / 2)
            if not (m % 2 == 1 and i == half - 1):
                xs.append((1 - r) / 2)
                ws.append(w / 2)
        self.x, self.w = xs, ws
        Gauss.cache[key] = (xs, ws)


class Geometry:
    def __init__(self, n):
        self.n = n
        self.t = (arb.pi() / n).tan()
        self.R = 1 / (arb.pi() / n).cos()
        self.R2 = self.R * self.R
        self.Rn = self.R ** n
        self.CP = self.R / arb.pi()              # Payne-Weinberger on the half-sector (diameter R)
        self.nu = arb.pi() ** 2 / self.R2        # lower bound for the first invariant eigenvalue
        self.tr = (2 * self.CP ** 2 + 2 * self.R * self.CP).sqrt()   # trace constant


class Moments:
    """Table I[a][k] = < phat^a qhat^k >_P for 0 <= a <= amax, 0 <= k <= kmax."""

    def __init__(self, G, amax, kmax):
        n = G.n
        self.G = G
        self.amax, self.kmax = amax, kmax
        dmax = 2 * amax + n * kmax
        m = dmax // 2 + 2
        gs = Gauss(m)
        self.m = m
        P, Q = [], []
        for x in gs.x:
            tau = G.t * x
            P.append((1 + tau * tau) / G.R2)
            Q.append((acb(1, tau) ** n).real / G.Rn)
        VP = [[arb(1)] * m]
        for a in range(amax):
            VP.append([u * v for u, v in zip(VP[-1], P)])
        VQ = [[wi for wi in gs.w]]
        for k in range(kmax):
            VQ.append([u * v for u, v in zip(VQ[-1], Q)])
        Mx = arb_mat(amax + 1, m, [e for row in VP for e in row]) * \
            arb_mat(m, kmax + 1, [VQ[k][j] for j in range(m) for k in range(kmax + 1)])
        self.I = [[2 * Mx[a, k] / (2 * a + n * k + 2) for k in range(kmax + 1)] for a in range(amax + 1)]

    def __call__(self, a, k):
        return self.I[a][k]


# ------- polynomials as dicts {(a,k): coefficient} in the scaled monomials -------------------------

def padd(f, g, cg=1):
    r = dict(f)
    for mon, c in g.items():
        r[mon] = r.get(mon, arb(0)) + c * cg
    return r


def pscale(f, s):
    return {m: c * s for m, c in f.items()}


def pmul(f, g):
    r = {}
    for (a, k), c in f.items():
        for (b, l), d in g.items():
            mon = (a + b, k + l)
            r[mon] = r.get(mon, arb(0)) + c * d
    return r


def harmonic_polys(G, K):
    """h_k = Re((w/R)^(nk)) for 0 <= k <= K in the scaled monomials: h_{k+1} = 2 qhat h_k - phat^n h_{k-1}."""
    n = G.n
    hs = [{(0, 0): arb(1)}, {(0, 1): arb(1)}]
    for k in range(1, K):
        nxt = {}
        for (a, kk), c in hs[k].items():
            nxt[(a, kk + 1)] = nxt.get((a, kk + 1), arb(0)) + 2 * c
        for (a, kk), c in hs[k - 1].items():
            nxt[(a + n, kk)] = nxt.get((a + n, kk), arb(0)) - c
        hs.append(nxt)
    return hs[:K + 1]


def lap_poly(G, f):
    """Laplacian of a polynomial: Lap m_{a,k} = R^-2 [4a(a+nk) m_{a-1,k} + n^2 k(k-1) m_{a+n-1,k-2}]"""
    n = G.n
    r = {}
    iR2 = 1 / G.R2
    for (a, k), c in f.items():
        if a > 0:
            mon = (a - 1, k)
            r[mon] = r.get(mon, arb(0)) + c * iR2 * (4 * a * (a + n * k))
        if k > 1:
            mon = (a + n - 1, k - 2)
            r[mon] = r.get(mon, arb(0)) + c * iR2 * (n * n * k * (k - 1))
    return r


def mean_poly(M, f):
    s = arb(0)
    for (a, k), c in f.items():
        s += c * M(a, k)
    return s


def inner_poly(M, f, g):
    """<f g> by the moment table"""
    s = arb(0)
    for (a, k), c in f.items():
        for (b, l), d in g.items():
            s += c * d * M(a + b, k + l)
    return s


def cholesky(A):
    """Lower Cholesky factor of a symmetric arb matrix given as a list of lists (ball arithmetic)."""
    N = len(A)
    L = [[arb(0)] * N for _ in range(N)]
    for j in range(N):
        s = A[j][j]
        Lj = L[j]
        for k in range(j):
            s -= Lj[k] * Lj[k]
        if not s > 0:
            raise ArithmeticError(f"Cholesky failed at {j}: {s}")
        d = s.sqrt()
        Lj[j] = d
        for i in range(j + 1, N):
            s = A[i][j]
            Li = L[i]
            for k in range(j):
                s -= Li[k] * Lj[k]
            Li[j] = s / d
    return L


def side_nodes(G, m):
    """Gauss nodes on the half side: returns taus, weights (summing to 1 on [0,t])."""
    gs = Gauss(m)
    return [G.t * x for x in gs.x], gs.w


def tomat(L):
    N = len(L)
    M = len(L[0])
    return arb_mat(N, M, [e for r in L for e in r])
