# Exploration script (floating point; not a certificate). Environment: WSL Fedora, Python 3.12, numpy, scipy, numba.
# Run from scripts/exploration/ (it imports the sibling modules). Typical runtime: seconds to minutes.
"""Floating-point model of the mixed-mechanism network (exploration only, not a certificate).

Species order (Conradi-Mincheva-Shiu 2019, Table 3):
x1 S0, x2 K, x3 S0K, x4 S1K, x5 S2, x6 P, x7 S2P, x8 S1, x9 S1P.
"""
import numpy as np

# reactions: (reactant indices, product indices) 0-based
REACTIONS = [
    ((0, 1), (2,)),      # k1: S0 + K -> S0K
    ((2,), (0, 1)),      # k2
    ((2,), (3,)),        # k3: S0K -> S1K
    ((3,), (4, 1)),      # k4: S1K -> S2 + K
    ((4, 5), (6,)),      # k5: S2 + P -> S2P
    ((6,), (4, 5)),      # k6
    ((6,), (7, 5)),      # k7: S2P -> S1 + P
    ((7, 5), (8,)),      # k8: S1 + P -> S1P
    ((8,), (7, 5)),      # k9
    ((8,), (0, 5)),      # k10: S1P -> S0 + P
]
NS, NR = 9, 10
N = np.zeros((NS, NR))
for j, (re, pr) in enumerate(REACTIONS):
    for i in re:
        N[i, j] -= 1
    for i in pr:
        N[i, j] += 1
U, s, Vt = np.linalg.svd(N)
V = U[:, :6]          # orthonormal basis of stoichiometric subspace S
W = np.array([[0, 1, 1, 1, 0, 0, 0, 0, 0],
              [0, 0, 0, 0, 0, 1, 1, 0, 1],
              [1, 0, 1, 1, 1, 0, 1, 1, 1]], float)
assert np.allclose(W @ N, 0)


def steady_state(k, x1, x2, x6):
    k1, k2, k3, k4, k5, k6, k7, k8, k9, k10 = k
    J = k1 * k3 / (k2 + k3) * x1 * x2
    x3 = k1 / (k2 + k3) * x1 * x2
    x4 = J / k4
    x5 = (k6 + k7) * J / (k7 * k5 * x6)
    x7 = J / k7
    x8 = (k9 + k10) * J / (k10 * k8 * x6)
    x9 = J / k10
    return np.array([x1, x2, x3, x4, x5, x6, x7, x8, x9])


def rates(k, x):
    out = np.empty(NR)
    for j, (re, pr) in enumerate(REACTIONS):
        v = k[j]
        for i in re:
            v = v * x[i]
        out[j] = v
    return out


def g(k, x):
    return N @ rates(k, x)


def jac(k, x):
    D = np.zeros((NR, NS))
    for j, (re, pr) in enumerate(REACTIONS):
        for a in re:
            v = k[j]
            for b in re:
                if b != a:
                    v = v * x[b]
            D[j, a] += v
    return N @ D


def hess_bilinear(k, u, v):
    """D^2 g(u, v) (independent of x: quadratic field)."""
    r = np.zeros(NR, dtype=complex)
    for j, (re, pr) in enumerate(REACTIONS):
        if len(re) == 2:
            a, b = re
            r[j] = k[j] * (u[a] * v[b] + u[b] * v[a])
    return N @ r


def reduced_A(k, x):
    return V.T @ jac(k, x) @ V


def B6(k, u, v):
    return V.T @ hess_bilinear(k, V @ u, V @ v)


def first_lyapunov(k, x):
    """Kuznetsov's projection formula (C = 0 for a quadratic field).
    Returns (l1, omega, eigenvalues). Uses the pair with smallest |Re|."""
    A = reduced_A(k, x)
    ev, R = np.linalg.eig(A)
    idx = np.argmin(np.abs(ev.real) + 1e-300 * (ev.imag <= 0) + (ev.imag <= 0) * 1e300)
    lam = ev[idx]
    omega = lam.imag
    sc = np.abs(ev).max()
    if not (omega > 1e-9 * sc and abs(lam.real) < 1e-6 * omega):
        return np.nan, 0.0, ev
    q = R[:, idx]
    q = q / np.linalg.norm(q)
    evl, L = np.linalg.eig(A.T)
    jdx = np.argmin(np.abs(evl - np.conj(lam)))
    p = L[:, jdx]
    # normalize <p,q> = conj(p).q = 1 with A^T p = conj(lam) p
    p = p / np.conj(np.vdot(p, q))
    Bqq = B6(k, q, q)
    Bqqb = B6(k, q, np.conj(q))
    s1 = np.linalg.solve(A, Bqqb)
    s2 = np.linalg.solve(2j * omega * np.eye(6) - A, Bqq)
    val = np.vdot(p, -2 * B6(k, q, s1) + B6(k, np.conj(q), s2))
    return val.real / (2 * omega), omega, ev


def maxre(k, x):
    return np.max(np.linalg.eigvals(reduced_A(k, x)).real)
