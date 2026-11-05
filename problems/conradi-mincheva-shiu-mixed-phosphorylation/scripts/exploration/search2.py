# Exploration script (floating point; not a certificate). Environment: WSL Fedora, Python 3.12, numpy, scipy, numba.
# Run from scripts/exploration/ (it imports the sibling modules). Typical runtime: seconds to minutes.
"""Biased random search: rate constants in the region of Theorem 4.2 (k10 large, k5/k8 large),
with unequal dissociation constants allowed. Records dimensionless L = l1*|x|^2/omega."""
import sys, time
import numpy as np
from scipy.optimize import brentq
from model import *
from search_lib import ss_batch, maxre_batch  # noqa (search.py guarded below)

seed = int(sys.argv[1]); nsamp = int(sys.argv[2]); span = float(sys.argv[3])
rng = np.random.default_rng(seed)
grid = np.logspace(-6, 7, 200)
out = []
t0 = time.time(); nh = 0; npos = 0
for it in range(nsamp):
    k = np.exp(rng.uniform(-span, span, 10))
    if rng.random() < 0.3:
        k[[1, 5, 8]] = k[1]
    k3, k4, k7 = k[2], k[3], k[6]
    k[9] = (k3 + k4) * (k3 + k7) * (k4 + k7) / (k3 * k4) * np.exp(rng.uniform(0, 7))
    k[4] = k[7] * np.exp(rng.uniform(0, 10))
    x1 = np.exp(rng.uniform(-8, 8)); x6 = np.exp(rng.uniform(-8, 8))
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
        if om <= 0:
            continue
        L = l1 * np.dot(x, x) / om
        nh += 1
        if l1 > 0:
            npos += 1
        out.append(np.concatenate([k, [x1, np.exp(lx2), x6, l1, om, L]]))
np.save('search2_%d.npy' % seed, np.array(out))
print('seed', seed, 'samples', nsamp, 'hopf points', nh, 'l1>0:', npos, 'time %.1f' % (time.time() - t0))
