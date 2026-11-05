# Exploration script (floating point; not a certificate). Environment: WSL Fedora, Python 3.12, numpy, scipy, numba.
# Run from scripts/exploration/ (it imports the sibling modules). Typical runtime: seconds to minutes.
"""Count eigenvalues with positive real part at random steady states (float exploration)."""
import sys
import numpy as np
from model import *
from search_lib import ss_batch, maxre_batch
from collections import Counter

rng = np.random.default_rng(int(sys.argv[1])); n = int(sys.argv[2])
cnt = Counter(); ex = {}
for it in range(n):
    k = np.exp(rng.uniform(-4, 4, 10))
    if rng.random() < 0.7:
        k3, k4, k7 = k[2], k[3], k[6]
        k[9] = (k3 + k4) * (k3 + k7) * (k4 + k7) / (k3 * k4) * np.exp(rng.uniform(0, 8))
        k[4] = k[7] * np.exp(rng.uniform(0, 12))
    x1, x2, x6 = np.exp(rng.uniform(-8, 8, 3))
    x = steady_state(k, x1, x2, x6)
    ev = np.linalg.eigvals(reduced_A(k, x))
    m = int((ev.real > 1e-9 * np.abs(ev).max()).sum())
    cnt[m] += 1
    if m >= 3 and m not in ex:
        ex[m] = (k, x1, x2, x6, ev)
print(dict(cnt))
for m, v in ex.items():
    print(m, v)
