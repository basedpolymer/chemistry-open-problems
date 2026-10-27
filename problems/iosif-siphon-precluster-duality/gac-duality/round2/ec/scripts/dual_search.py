"""Numerical exploration (NOT a proof, floating point): multistationarity and instability of N* when N is a
well-formed equal-conflict (WF-EC) network with unbalanced joins (question Q3-Udual of round2/ec/candidates.md).

Agent Q3, Claude Opus 5.5 (effort xhigh).  Environment: Python 3.12 with numpy and scipy (WSL Fedora 44).
Usage (from round2/ec/): python scripts/dual_search.py NSPECIES RMAX MAXMOL NKAPPA NCLASS NSTART SEED [NPROC]
Networks N: exhaustive, up to relabelling of the remaining species: species X1..Xn, one unbalanced source
2X1 + X2 (n = 3: the other source a*X3; n = 4: the other sources are a block partition of {X3, X4} with arbitrary
coefficients), products arbitrary nonzero complexes of molecularity <= MAXMOL different from the source, at most
RMAX reactions; kept if WF-EC (exact test with eclib).  For each such N: NKAPPA log-uniform kappa in [1/10, 10],
NCLASS random classes of N* (k0 log-uniform), NSTART random starts of a Newton-type solver on the system
"k in class, Gamma^T omega(k) = 0" in log coordinates; distinct positive solutions are counted (relative tolerance
1e-6) and the eigenvalues of the Jacobian of N* restricted to its stoichiometric subspace are computed.
Output: counts; every (N, kappa, class) with >= 2 solutions or an eigenvalue with positive real part is printed
with all data needed for an exact re-check.
"""
import sys
import time
import platform
from itertools import combinations, product
from multiprocessing import Pool

import numpy as np
import warnings
warnings.filterwarnings("ignore")
np.seterr(all="ignore")
from scipy.optimize import root

import eclib
from eclib import Network


def complexes(n, m):
    return [c for c in product(range(m + 1), repeat=n) if 0 < sum(c) <= m]


def source_sets(n):
    """Unbalanced source 2X1+X2 plus sources on the remaining species (all block partitions, coefficients 1..2)."""
    rest = list(range(2, n))
    y0 = tuple([2, 1] + [0] * (n - 2))
    out = []
    if not rest:
        return [[y0]]
    # partitions of rest
    def parts(items):
        if not items:
            yield []
            return
        f, r = items[0], items[1:]
        for p in parts(r):
            for k in range(len(p)):
                yield p[:k] + [[f] + p[k]] + p[k + 1:]
            yield [[f]] + p
    for p in parts(rest):
        choices = []
        for block in p:
            opts = []
            for coeffs in product(range(1, 3), repeat=len(block)):
                v = [0] * n
                for i, c in zip(block, coeffs):
                    v[i] = c
                opts.append(tuple(v))
            choices.append(opts)
        for srcs in product(*choices):
            out.append([y0] + list(srcs))
    return out


def gen(n, rmax, maxmol):
    cpx = complexes(n, maxmol)
    for srcs in source_sets(n):
        k = len(srcs)
        prods = [[c for c in cpx if c != s] for s in srcs]

        def rec(idx, budget, acc):
            if idx == k:
                yield list(acc)
                return
            rem = k - idx - 1
            for size in range(1, budget - rem + 1):
                for ps in combinations(prods[idx], size):
                    yield from rec(idx + 1, budget - size, acc + [(srcs[idx], p) for p in ps])
        yield from rec(0, rmax, [])


def solve_dual(N, kappa, k0, nstart, rng):
    D = N.dual()                                   # species = reactions of N (r), reactions = species of N (n)
    G = np.array(N.Gamma, dtype=float)             # n x r ;  Gamma(N*) = G^T (r x n)
    Ye = np.array(N.Ye, dtype=float)               # n x r ; rate of rho_p = kappa_p prod_t k_t^{Ye[p,t]}
    r = N.r
    U, sv, Vt = np.linalg.svd(G)
    s = int((sv > 1e-9).sum())
    kerG = Vt[s:]                                  # rows: basis of ker Gamma = conservation laws of N*
    Sbasis = U[:, :s]                              # columns: basis of S = im Gamma (in R^n)

    def omega(k):
        return kappa * np.exp(Ye @ np.log(k))

    def F(q):
        k = np.exp(q)
        return np.concatenate([kerG @ (k - k0), Sbasis.T @ omega(k)])

    sols = []
    for _ in range(nstart):
        q0 = np.log(k0) + rng.normal(0, 1.5, r)
        try:
            res = root(F, q0, method="hybr", options={"xtol": 1e-13, "maxfev": 4000})
        except Exception:
            continue
        if not res.success:
            continue
        q = res.x
        if not np.all(np.isfinite(q)) or np.max(np.abs(F(q))) > 1e-9 * (1 + np.max(np.abs(kappa))):
            continue
        k = np.exp(q)
        if any(np.max(np.abs(k - k1) / (np.abs(k1) + 1e-12)) < 1e-6 for k1 in sols):
            continue
        sols.append(k)
    eigs = []
    for k in sols:
        om = omega(k)
        # J* = Gamma^T diag(om) Ye diag(1/k)   (r x r)
        J = G.T @ np.diag(om) @ Ye @ np.diag(1.0 / k)
        # restrict to im Gamma^T (dimension s): basis Q of row space of G
        Q = Vt[:s].T                                # r x s orthonormal basis of im Gamma^T
        Jr = Q.T @ J @ Q
        eigs.append(np.linalg.eigvals(Jr))
    return sols, eigs


def worker(args):
    rx, n, nk, ncl, nst, seed = args
    N = Network(n, rx)
    rng = np.random.default_rng(seed)
    hits = []
    worst = -np.inf
    nmax = 0
    for _ in range(nk):
        kappa = np.exp(rng.uniform(np.log(0.1), np.log(10), N.n))
        for _ in range(ncl):
            k0 = np.exp(rng.uniform(np.log(0.1), np.log(10), N.r))
            sols, eigs = solve_dual(N, kappa, k0, nst, rng)
            nmax = max(nmax, len(sols))
            for e in eigs:
                worst = max(worst, float(np.max(e.real)))
            if len(sols) >= 2 or any(np.max(e.real) > 1e-8 for e in eigs):
                hits.append((list(kappa), list(k0), [list(k) for k in sols], [list(e) for e in eigs]))
            if len(sols) == 0:
                hits.append(("NO SOLUTION FOUND", list(kappa), list(k0)))
    return rx, nmax, worst, hits


def main():
    n, rmax, maxmol, nk, ncl, nst, seed = (int(a) for a in sys.argv[1:8])
    nproc = int(sys.argv[8]) if len(sys.argv) > 8 else 8
    t0 = time.time()
    print(f"# dual_search.py n={n} r<={rmax} mol<={maxmol} NKAPPA={nk} NCLASS={ncl} NSTART={nst} seed={seed}; "
          f"Python {platform.python_version()}, numpy {np.__version__} on {platform.system()}")
    nets = []
    seen = set()
    for rx in gen(n, rmax, maxmol):
        key = tuple(sorted(rx))
        if key in seen:
            continue
        seen.add(key)
        N = Network(n, key)
        if eclib.is_wf_ec(N) and not eclib.is_balanced(N):
            nets.append(key)
    print(f"# unbalanced WF-EC networks generated: {len(nets)} ({time.time() - t0:.1f} s)")
    tasks = [(rx, n, nk, ncl, nst, seed + 7919 * i) for i, rx in enumerate(nets)]
    nmulti = nunstable = nnone = 0
    hist = {}
    with Pool(nproc) as pool:
        for rx, nmax, worst, hits in pool.imap_unordered(worker, tasks):
            hist[nmax] = hist.get(nmax, 0) + 1
            N = Network(n, rx)
            for h in hits:
                if h[0] == "NO SOLUTION FOUND":
                    nnone += 1
                    print(f"NOSOL  {N.describe()}  kappa={h[1]} k0={h[2]}")
                    continue
                kappa, k0, sols, eigs = h
                if len(sols) >= 2:
                    nmulti += 1
                if any(max(x.real for x in e) > 1e-8 for e in eigs):
                    nunstable += 1
                print(f"HIT  {N.describe()}  kappa={kappa}  k0={k0}  sols={sols}  maxRe={[max(x.real for x in e) for e in eigs]}")
    print(f"# histogram of max number of solutions found per network: {dict(sorted(hist.items()))}")
    print(f"# (N, kappa, class) with >=2 solutions: {nmulti}; with an unstable solution: {nunstable}; "
          f"with no solution found: {nnone}")
    print(f"# runtime {time.time() - t0:.1f} s")


if __name__ == "__main__":
    main()
