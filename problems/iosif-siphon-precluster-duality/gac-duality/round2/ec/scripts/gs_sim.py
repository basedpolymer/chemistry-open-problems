"""Numerical evidence (floating point, NOT a proof) for Theorem GS and the open question GS* (round2/ec/notes.md):
global convergence of WF-EC networks N and of their duals N*, in particular N* with unbalanced joins and s >= 3.

Agent Q3, Claude Opus 5.5 (effort xhigh).  Environment: Python 3.12, numpy, scipy (WSL Fedora 44).
Usage (from round2/ec/): python scripts/gs_sim.py NSPECIES NNETS SEED [NPROC]
Networks: random (seeded) WF-EC networks on NSPECIES species with the unbalanced source 2X1+X2 and the other
species as single-species sources a*Xi (a in {1,2}); products random nonzero complexes of molecularity <= 3; r in
[m, m+3]; kept if WF-EC (exact test).  For each network: random rate constants k (N) and kappa (N*), three random
initial points in one class of N and three in one class of N*; integration with LSODA up to t = 2000.
Reported: for N, max over the run of the increase of M(t) = max_c x^{y_c}/u_c(k) and of the decrease of
mu(t) = min_c ... (Theorem GS: both <= 0 up to rounding); for both sides, the spread of the three end points and
the norm of the vector field at the end points.
"""
import sys
import time
import platform
import warnings
from itertools import product
from multiprocessing import Pool

import numpy as np
from scipy.integrate import solve_ivp

import eclib
from eclib import Network, ECStructure

warnings.filterwarnings("ignore")


def random_wfec(n, rng):
    cpx = [c for c in product(range(4), repeat=n) if 0 < sum(c) <= 3]
    for _ in range(200000):
        srcs = [tuple([2, 1] + [0] * (n - 2))]
        for i in range(2, n):
            v = [0] * n
            v[i] = int(rng.integers(1, 3))
            srcs.append(tuple(v))
        m = len(srcs)
        r = int(rng.integers(m, m + 4))
        rx = set()
        for c in range(m):
            b = cpx[int(rng.integers(len(cpx)))]
            if b != srcs[c]:
                rx.add((srcs[c], b))
        while len(rx) < r:
            c = int(rng.integers(m))
            b = cpx[int(rng.integers(len(cpx)))]
            if b != srcs[c]:
                rx.add((srcs[c], b))
        N = Network(n, sorted(rx))
        if len({a for a, b in N.reactions}) == m and eclib.is_wf_ec(N):
            return N
    return None


def simulate(F, x0, T=2000.0):
    sol = solve_ivp(lambda t, x: F(np.maximum(x, 0)), (0, T), x0, method="LSODA", rtol=1e-10, atol=1e-13,
                    dense_output=False)
    return sol.t, sol.y


def run(args):
    n, seed = args
    rng = np.random.default_rng(seed)
    N = random_wfec(n, rng)
    if N is None:
        return None
    E = ECStructure(N)
    G = np.array(N.Gamma, float)
    Ye = np.array(N.Ye, float)
    Ys = np.array([list(y) for y in E.sources], float)          # m x n
    k = np.exp(rng.uniform(-1.5, 1.5, N.r))
    A = np.zeros((N.n, E.m))
    for t, (a, b) in enumerate(N.reactions):
        A[:, E.cl_of_reaction[t]] += k[t] * (np.array(b) - np.array(a))
    u = np.linalg.svd(A)[2][-1]
    u = u / u.sum()
    F = lambda x: G @ (k * np.exp(Ye.T @ np.log(np.maximum(x, 1e-300))))
    x_ref = np.exp(rng.uniform(-1, 1, N.n))
    Sb = np.linalg.svd(G)[0][:, :N.rank()]
    ends, worstM, worstmu = [], 0.0, 0.0
    for _ in range(3):
        x0 = x_ref + Sb @ rng.normal(0, 0.3, Sb.shape[1])
        if np.any(x0 <= 0):
            x0 = x_ref.copy()
        t, X = simulate(F, x0)
        psi = np.exp(Ys @ np.log(np.maximum(X, 1e-300))) / u[:, None]
        M, mu = psi.max(axis=0), psi.min(axis=0)
        worstM = max(worstM, float(np.max(np.diff(M) / M[:-1])) if len(M) > 1 else 0.0)
        worstmu = max(worstmu, float(np.max(-np.diff(mu) / mu[:-1])) if len(mu) > 1 else 0.0)
        ends.append(X[:, -1])
    spreadN = max(np.max(np.abs(e - ends[0])) for e in ends)
    resN = max(np.max(np.abs(F(e))) for e in ends)
    # dual
    kap = np.exp(rng.uniform(-1.5, 1.5, N.n))
    Fd = lambda q: G.T @ (kap * np.exp(Ye @ np.log(np.maximum(q, 1e-300))))
    q_ref = np.exp(rng.uniform(-1, 1, N.r))
    Sd = np.linalg.svd(G.T)[0][:, :N.rank()]
    endsD = []
    for _ in range(3):
        q0 = q_ref + Sd @ rng.normal(0, 0.3, Sd.shape[1])
        if np.any(q0 <= 0):
            q0 = q_ref.copy()
        t, Q = simulate(Fd, q0)
        endsD.append(Q[:, -1])
    spreadD = max(np.max(np.abs(e - endsD[0])) for e in endsD)
    resD = max(np.max(np.abs(Fd(e))) for e in endsD)
    return (N.describe(), N.rank(), eclib.is_balanced(N), worstM, worstmu, spreadN, resN, spreadD, resD)


def main():
    n, nnets, seed = (int(a) for a in sys.argv[1:4])
    nproc = int(sys.argv[4]) if len(sys.argv) > 4 else 6
    t0 = time.time()
    print(f"# gs_sim.py n={n} nets={nnets} seed={seed}; Python {platform.python_version()}, numpy {np.__version__}")
    with Pool(nproc) as pool:
        res = [r for r in pool.map(run, [(n, seed + i) for i in range(nnets)]) if r]
    bysd = {}
    for desc, s, bal, wM, wmu, spN, rN, spD, rD in res:
        bysd[s] = bysd.get(s, 0) + 1
        flag = "  <-- CHECK" if (wM > 1e-6 or wmu > 1e-6 or spN > 1e-5 or rN > 1e-7 or spD > 1e-5 or rD > 1e-7) else ""
        print(f"NET s={s} balanced={bal} M-increase={wM:.1e} mu-decrease={wmu:.1e} N: spread={spN:.1e} "
              f"|f|={rN:.1e}  N*: spread={spD:.1e} |f|={rD:.1e}  {desc}{flag}")
    print(f"# networks by s: {dict(sorted(bysd.items()))}; runtime {time.time() - t0:.1f} s")


if __name__ == "__main__":
    main()
