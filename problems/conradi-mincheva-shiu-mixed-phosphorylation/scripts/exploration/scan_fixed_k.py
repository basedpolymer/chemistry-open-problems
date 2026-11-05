# Exploration script (floating point; not a certificate). Environment: WSL Fedora, Python 3.12, numpy, scipy, numba.
# Run from scripts/exploration/ (it imports the sibling modules). Typical runtime: seconds to minutes.
"""For fixed rate constants, scan the Hopf surface over a grid of (x1, x6) and report the sign of l1."""
import sys
import numpy as np
from scipy.optimize import brentq
from model import *
from search_lib import ss_batch, maxre_batch

k = np.array([float(s) for s in sys.argv[1].split(',')])
n = int(sys.argv[2]) if len(sys.argv) > 2 else 60
grid = np.logspace(-7, 8, 300)
res = []
for x1 in np.logspace(-4, 4, n):
    for x6 in np.logspace(-5, 4, n):
        X = ss_batch(k, x1, grid, x6)
        mr = maxre_batch(k, X)
        sc = np.where(np.sign(mr[:-1]) != np.sign(mr[1:]))[0]
        for i in sc:
            f = lambda lx2: maxre(k, steady_state(k, x1, np.exp(lx2), x6))
            try:
                lx2 = brentq(f, np.log(grid[i]), np.log(grid[i + 1]), xtol=1e-13)
            except ValueError:
                continue
            x = steady_state(k, x1, np.exp(lx2), x6)
            l1, om, ev = first_lyapunov(k, x)
            if om > 0:
                res.append([x1, np.exp(lx2), x6, l1, om, l1 * np.dot(x, x) / om] + list(W @ x))
res = np.array(res)
print('Hopf points:', len(res), ' l1>0:', (res[:, 3] > 0).sum())
np.set_printoptions(precision=5, linewidth=220)
if (res[:, 3] > 0).any():
    pos = res[res[:, 3] > 0]
    print('positive (x1,x2,x6,l1,om,L,Ktot,Ptot,Stot):')
    print(pos[np.argsort(-pos[:, 5])][:15])
np.save('scan_%s.npy' % sys.argv[1].replace(',', '_')[:60], res)
