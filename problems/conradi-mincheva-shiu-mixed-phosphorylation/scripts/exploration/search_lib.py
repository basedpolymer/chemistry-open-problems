# Exploration script (floating point; not a certificate). Environment: WSL Fedora, Python 3.12, numpy, scipy, numba.
# Run from scripts/exploration/ (it imports the sibling modules). Typical runtime: seconds to minutes.
"""Random search over rate constants and steady states for Hopf points with l1 > 0 (floating point)."""
import sys, time
import numpy as np
from scipy.optimize import brentq
from model import *



def ss_batch(k, x1, x2, x6):
    k1, k2, k3, k4, k5, k6, k7, k8, k9, k10 = k
    J = k1 * k3 / (k2 + k3) * x1 * x2
    one = np.ones_like(x2)
    return np.stack([x1 * one, x2, k1 / (k2 + k3) * x1 * x2, J / k4, (k6 + k7) * J / (k7 * k5 * x6),
                     x6 * one, J / k7, (k9 + k10) * J / (k10 * k8 * x6), J / k10], axis=-1)


def jac_batch(k, X):
    m = X.shape[0]
    D = np.zeros((m, NR, NS))
    for j, (re, pr) in enumerate(REACTIONS):
        if len(re) == 1:
            D[:, j, re[0]] = k[j]
        else:
            a, b = re
            D[:, j, a] = k[j] * X[:, b]
            D[:, j, b] = k[j] * X[:, a]
    return N @ D


def maxre_batch(k, X):
    A = V.T @ jac_batch(k, X) @ V
    return np.linalg.eigvals(A).real.max(axis=1)


if __name__ == '__main__':
    seed = int(sys.argv[1]); nsamp = int(sys.argv[2]); span = float(sys.argv[3])
    rng = np.random.default_rng(seed)
    grid = np.logspace(-5, 6, 160)
    out = []
    t0 = time.time()
    nh = 0; npos = 0
    for it in range(nsamp):
        k = np.exp(rng.uniform(-span, span, 10))
        if rng.random() < 0.5:
            k[[1, 5, 8]] = k[1]  # equal dissociation constants half the time
        x1 = np.exp(rng.uniform(-6, 6)); x6 = np.exp(rng.uniform(-6, 6))
        X = ss_batch(k, x1, grid, x6)
        mr = maxre_batch(k, X)
        sc = np.where(np.sign(mr[:-1]) != np.sign(mr[1:]))[0]
        for i in sc:
            f = lambda lx2: maxre(k, steady_state(k, x1, np.exp(lx2), x6))
            try:
                lx2 = brentq(f, np.log(grid[i]), np.log(grid[i + 1]), xtol=1e-13)
            except Exception:
                continue
            x = steady_state(k, x1, np.exp(lx2), x6)
            l1, om, ev = first_lyapunov(k, x)
            evs = np.sort(ev.real)
            # normalise l1 by scale: dimensionless sign is what matters
            nh += 1
            if l1 > 0:
                npos += 1
            out.append(np.concatenate([k, [x1, np.exp(lx2), x6, l1, om, evs[-3]]]))
    np.save('search_%d.npy' % seed, np.array(out))
    print('seed', seed, 'samples', nsamp, 'hopf points', nh, 'l1>0:', npos, 'time %.1f' % (time.time() - t0))
