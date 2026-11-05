# Exploration script (floating point; not a certificate). Environment: WSL Fedora, Python 3.12, numpy, scipy, numba.
# Run from scripts/exploration/ (it imports the sibling modules). Typical runtime: seconds to minutes.
"""usage: python hopf_float.py 'k1,..,k10' x1 x6  -> all Hopf points in x2 with l1 (float)."""
import sys
import numpy as np
from fractions import Fraction
from scipy.optimize import brentq
from model import *
from search_lib import ss_batch, maxre_batch

k = np.array([float(Fraction(s)) for s in sys.argv[1].split(',')])
x1 = float(Fraction(sys.argv[2])); x6 = float(Fraction(sys.argv[3]))
grid = np.logspace(-8, 9, 2000)
mr = maxre_batch(k, ss_batch(k, x1, grid, x6))
for i in np.where(np.sign(mr[:-1]) != np.sign(mr[1:]))[0]:
    f = lambda lx2: maxre(k, steady_state(k, x1, np.exp(lx2), x6))
    try:
        lx2 = brentq(f, np.log(grid[i]), np.log(grid[i + 1]), xtol=1e-14)
    except ValueError:
        continue
    x = steady_state(k, x1, np.exp(lx2), x6)
    l1, om, ev = first_lyapunov(k, x)
    d = 1e-6
    dre = (maxre(k, steady_state(k, x1, np.exp(lx2) * (1 + d), x6)) - maxre(k, steady_state(k, x1, np.exp(lx2) * (1 - d), x6))) / (2 * d * np.exp(lx2))
    print('x2 = %.12g  omega = %.6g  l1 = %.6g  L = %.4g  dRe/dx2 = %.4g  totals = %s' % (np.exp(lx2), om, l1, l1 * x.dot(x) / om, dre, W @ x))
