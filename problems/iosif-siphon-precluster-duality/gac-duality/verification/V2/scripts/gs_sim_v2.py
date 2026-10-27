"""Verifier V2: numerical red team of Theorem GS (round2/ec/notes.md §3) -- floating point, NOT a proof.

For random WF-EC networks (own generator, impure and unbalanced allowed, n = 4..6) and random k, integrates the
mass-action ODE (scipy solve_ivp, Radau, rtol 1e-10, atol 1e-13) from (i) random interior points, (ii) points with
some coordinates at 1e-6 (near the boundary), (iii) points with two clusters tied for the maximum of psi_c/u_c,
and records: the largest relative increase of M(t) = max_c x^{y_c}/u_c and decrease of mu(t) = min_c ..., on a dense
grid; the spread of the end points of trajectories in one class; the residual |f| at the end; min_t min_p x_p(t).
Also the dual N* (GS* open for s >= 3) from three points of one class.

Environment: WSL Fedora 44, Python 3.12, numpy 2.5, scipy 1.18.  Usage: python gs_sim_v2.py SEED COUNT
"""
import sys, os, time, random
import numpy as np
from scipy.integrate import solve_ivp
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from v2lib import *  # noqa
from itertools import product

T0 = time.time()
seed = int(sys.argv[1]) if len(sys.argv) > 1 else 1
COUNT = int(sys.argv[2]) if len(sys.argv) > 2 else 20
rng = random.Random(seed)
nprng = np.random.default_rng(seed)


def gen_wfec(n, m):
    while True:
        lab = [rng.randrange(m) for _ in range(n)]
        if len(set(lab)) == m:
            break
    B = [[p for p in range(n) if lab[p] == c] for c in range(m)]
    srcs = []
    for blk in B:
        y = [0] * n
        for p in blk:
            y[p] = rng.choice([1, 1, 2, 3])
        srcs.append(tuple(y))
    rows = n - m + 1
    while True:
        Z = [[rng.choice([0, 0, 1, 1, 2]) for _ in range(n)] for _ in range(rows)]
        if all(any(Z[i][p] for i in range(rows)) for p in range(n)) and rank(Z) == rows:
            break
    cand = list(product(range(4), repeat=n))
    zv = lambda y: tuple(sum(Z[i][p] * y[p] for p in range(n)) for i in range(rows))
    N = []
    for y in srcs:
        prods = [c for c in cand if c != y and zv(c) == zv(y)]
        if not prods:
            return None
        for _ in range(rng.randint(1, 3)):
            N.append((y, rng.choice(prods)))
    N = sorted(set(N))
    return N if is_WFEC(N) else None


def mk_rhs(N, k):
    Ye_ = np.array(Ye(N), dtype=float)
    G = np.array(Gamma(N), dtype=float)
    kk = np.array([float(t) for t in k])

    def f(t, x):
        xp = np.maximum(x, 0.0)
        w = kk * np.prod(xp[:, None] ** Ye_, axis=0)
        return G @ w

    def jac(t, x):
        xp = np.maximum(x, 1e-300)
        w = kk * np.prod(xp[:, None] ** Ye_, axis=0)
        D = (Ye_ * w[None, :]) / xp[:, None]  # dw_t/dx_p = Ye[p,t] w_t / x_p
        return G @ D.T
    return f, jac


worst_M = 0.0
worst_mu = 0.0
worst_spread = 0.0
worst_res = 0.0
worst_spread_dual = 0.0
worst_res_dual = 0.0
minx_global = np.inf
done = 0
tries = 0
while done < COUNT and tries < 5000 * COUNT:
    tries += 1
    n = rng.randint(4, 6)
    m = rng.randint(2, min(4, n))
    N = gen_wfec(n, m)
    if N is None:
        continue
    done += 1
    r = len(N)
    srcs, ct = ec_structure(N)
    k = [F(rng.randint(1, 20), rng.randint(1, 20)) for _ in range(r)]
    u = np.array([float(t) for t in u_of_k(N, k)])
    S = np.array(srcs, dtype=float)
    f, jac = mk_rhs(N, k)
    s = rank(Gamma(N))
    Sbasis = np.array([[float(t) for t in row] for row in row_basis(transpose(Gamma(N)))])
    x0 = nprng.uniform(0.2, 3.0, n)
    starts = [x0]
    # same class: move along S, keep positivity
    for _ in range(2):
        for _try in range(100):
            d = nprng.normal(size=Sbasis.shape[0]) @ Sbasis
            lam = nprng.uniform(0.1, 1.0)
            # scale so that the point stays positive
            neg = d < 0
            tmax = np.min(-x0[neg] / d[neg]) if neg.any() else 1.0
            x1 = x0 + lam * min(tmax, 1.0) * d * 0.99
            if (x1 > 0).all():
                starts.append(x1)
                break
    # near-boundary point in the same class: push towards a face
    d = nprng.normal(size=Sbasis.shape[0]) @ Sbasis
    neg = d < 0
    if neg.any():
        tmax = np.min(-x0[neg] / d[neg])
        xb = x0 + tmax * d * (1 - 1e-6)
        if (xb > 0).all():
            starts.append(xb)
    ends = []
    for xs in starts:
        sol = solve_ivp(f, (0, 400.0), xs, method="Radau", jac=jac, rtol=1e-10, atol=1e-13, dense_output=True)
        tt = np.concatenate([np.linspace(0, 1, 400), np.linspace(1, 400, 1600)])
        X = sol.sol(tt)
        minx_global = min(minx_global, X.min() / max(1.0, xs.max()))
        psi = np.exp(S @ np.log(np.maximum(X, 1e-300)))  # m x len(tt)
        ratio = psi / u[:, None]
        M = ratio.max(axis=0)
        mu = ratio.min(axis=0)
        worst_M = max(worst_M, np.max((M[1:] - M[:-1]) / M[:-1]))
        worst_mu = max(worst_mu, np.max((mu[:-1] - mu[1:]) / mu[:-1]))
        ends.append(sol.y[:, -1])
        worst_res = max(worst_res, np.abs(f(0, sol.y[:, -1])).max() / max(1.0, np.abs(sol.y[:, -1]).max()))
    E = np.array(ends)
    worst_spread = max(worst_spread, np.abs(E - E[0]).max())
    # dual
    Dn = dual(N)
    kap = [F(rng.randint(1, 20), rng.randint(1, 20)) for _ in range(n)]
    fd, jd = mk_rhs(Dn, kap)
    Sd = np.array([[float(t) for t in row] for row in row_basis(transpose(Gamma(Dn)))])
    k0 = nprng.uniform(0.2, 3.0, r)
    dstarts = [k0]
    for _ in range(2):
        for _try in range(100):
            d = nprng.normal(size=Sd.shape[0]) @ Sd
            neg = d < 0
            tmax = np.min(-k0[neg] / d[neg]) if neg.any() else 1.0
            k1 = k0 + nprng.uniform(0.1, 1.0) * min(tmax, 1.0) * d * 0.99
            if (k1 > 0).all():
                dstarts.append(k1)
                break
    dends = []
    for ks in dstarts:
        sol = solve_ivp(fd, (0, 400.0), ks, method="Radau", jac=jd, rtol=1e-10, atol=1e-13)
        dends.append(sol.y[:, -1])
        worst_res_dual = max(worst_res_dual, np.abs(fd(0, sol.y[:, -1])).max() / max(1.0, np.abs(sol.y[:, -1]).max()))
    DE = np.array(dends)
    worst_spread_dual = max(worst_spread_dual, np.abs(DE - DE[0]).max())
print(f"seed {seed}: {done} WF-EC networks (n = 4..6), tries {tries}")
print(f"N : max relative increase of M(t) on the grid {worst_M:.2e}; max relative decrease of mu(t) {worst_mu:.2e}")
print(f"N : max spread of end points within a class {worst_spread:.2e}; max scaled residual |f| at t = 400: {worst_res:.2e}")
print(f"N : min over trajectories of min_p x_p(t) / max(1, max x0) = {minx_global:.2e} (near-boundary starts at 1e-6 scale)")
print(f"N*: max spread of end points within a class {worst_spread_dual:.2e}; max scaled residual {worst_res_dual:.2e}")
print(f"runtime {time.time() - T0:.1f} s")
