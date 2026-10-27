#!/usr/bin/env python3
"""
s1_direct.py -- numerical exploration for E2, structure S1: reactant complexes {0, E, C+D, A+B} (species A..E).
f = N_E (x_E - c1) + N_CD (x_C x_D - c2) + N_AB (x_A x_B - c3); columns random within the bimolecular cones;
N_0 = -(c1 N_E + c2 N_CD + c3 N_AB) must be >= 0.  In theta-coordinates x = x0 + B theta the steady states of a
class are the intersections of the plane x_E = c1 with two quadrics; the class Jacobian is
M = D phi(x) B (Metzler).  Looks for classes with >= 2 Hurwitz steady states.
Usage: python s1_direct.py seed seconds
"""
import numpy as np, sys, time
from numpy.polynomial import polynomial as Pl

rng = np.random.default_rng(int(sys.argv[1]))
TMAX = float(sys.argv[2])
A_, B_, C_, D_, E_ = range(5)


def rand_cols():
    while True:
        NE = rng.exponential(1, 5) * (rng.random(5) < 0.6)
        NE[E_] = -rng.exponential(2)
        NCD = rng.exponential(1, 5) * (rng.random(5) < 0.6)
        NCD[C_] = -rng.exponential(1) * (rng.random() < 0.8)
        NCD[D_] = -rng.exponential(1) * (rng.random() < 0.8)
        if NCD.sum() > 0:
            NCD[C_] -= NCD.sum() * rng.uniform(1, 2)
        NAB = rng.exponential(1, 5) * (rng.random(5) < 0.6)
        NAB[A_] = -rng.exponential(1) * (rng.random() < 0.8)
        NAB[B_] = -rng.exponential(1) * (rng.random() < 0.8)
        if NAB.sum() > 0:
            NAB[A_] -= NAB.sum() * rng.uniform(1, 2)
        B = np.column_stack([NE, NCD, NAB])
        if np.linalg.matrix_rank(B) < 3:
            continue
        c = np.exp(rng.uniform(-2, 2, 3))
        N0 = -B @ c
        if np.all(N0 >= -1e-12) and np.any(N0 > 1e-9):
            return B, c


def solve_class(B, c, x0):
    # plane x_E = c1: B[E] . theta = c1 - x0_E
    bE = B[E_]
    # particular solution and basis of plane
    th0 = bE * (c[0] - x0[E_]) / (bE @ bE)
    Q, _ = np.linalg.qr(np.column_stack([bE, rng.standard_normal(3), rng.standard_normal(3)]))
    P1, P2 = Q[:, 1], Q[:, 2]
    L = lambda i: (x0[i] + B[i] @ th0, B[i] @ P1, B[i] @ P2)  # affine in (s,t)

    def quad(i, j, cc):
        a0, a1, a2 = L(i)
        b0, b1, b2 = L(j)
        # (a0 + a1 s + a2 t)(b0 + b1 s + b2 t) - cc = A t^2 + (Bc0 + Bc1 s) t + (C0 + C1 s + C2 s^2)
        return (a2 * b2, np.array([a0 * b2 + a2 * b0, a1 * b2 + a2 * b1]),
                np.array([a0 * b0 - cc, a0 * b1 + a1 * b0, a1 * b1]))
    a1, b1, c1 = quad(C_, D_, c[1])
    a2, b2, c2 = quad(A_, B_, c[2])
    u = Pl.polysub(a1 * c2, a2 * c1)
    v = Pl.polysub(a1 * b2, a2 * b1)
    w = Pl.polysub(Pl.polymul(b1, c2), Pl.polymul(b2, c1))
    res = Pl.polysub(Pl.polymul(u, u), Pl.polymul(v, w))
    roots = Pl.polyroots(res) if np.any(np.abs(res) > 1e-14) else []
    out = []
    for s in roots:
        if abs(s.imag) > 1e-8 * (1 + abs(s)):
            continue
        s = s.real
        vv = Pl.polyval(s, v)
        if abs(vv) < 1e-12:
            continue
        t = -Pl.polyval(s, u) / vv
        th = th0 + s * P1 + t * P2
        x = x0 + B @ th
        if np.all(x > 1e-10):
            phi = np.array([x[E_] - c[0], x[C_] * x[D_] - c[1], x[A_] * x[B_] - c[2]])
            if np.max(np.abs(phi)) < 1e-7 * (1 + np.max(x) ** 2):
                out.append(x)
    # dedupe
    ded = []
    for x in out:
        if all(np.max(np.abs(x - y)) > 1e-7 * (1 + np.max(np.abs(x))) for y in ded):
            ded.append(x)
    return ded


def Mmat(B, x):
    Dphi = np.zeros((3, 5))
    Dphi[0, E_] = 1
    Dphi[1, C_] = x[D_]; Dphi[1, D_] = x[C_]
    Dphi[2, A_] = x[B_]; Dphi[2, B_] = x[A_]
    return Dphi @ B


t0 = time.time()
tried = 0
multi = 0
hits = 0
hist = {}
while time.time() - t0 < TMAX:
    B, c = rand_cols()
    tried += 1
    for rep in range(20):
        x0 = np.exp(rng.uniform(-3, 3, 5))
        sts = solve_class(B, c, x0)
        k = len(sts)
        if k < 2:
            continue
        stab = []
        for x in sts:
            ev = np.linalg.eigvals(Mmat(B, x))
            stab.append(np.max(ev.real) < -1e-9 * np.max(np.abs(ev)))
        key = (k, sum(stab))
        hist[key] = hist.get(key, 0) + 1
        if sum(stab) >= 2:
            hits += 1
            print("HIT", np.round(B, 4).tolist(), np.round(c, 4).tolist(), np.round(x0, 4).tolist(),
                  [np.round(x, 4).tolist() for x in sts], stab, flush=True)
print("tried", tried, "hist (n_steady, n_stable):", hist, "hits", hits)
