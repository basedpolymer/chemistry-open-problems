# Exploration script (floating point; not a certificate). Environment: WSL Fedora, Python 3.12, numpy, scipy, numba.
# Run from scripts/exploration/ (it imports the sibling modules). Typical runtime: seconds to minutes.
"""For each orbit at mu, find crossings of z_j = c_j (steady-state value) and report transversality."""
import sys, json
import numpy as np
from fractions import Fraction
from scipy.integrate import solve_ivp
from fastmodel import g9
from model import steady_state

k = np.array([float(Fraction(s)) for s in sys.argv[1].split(',')])
x1 = float(Fraction(sys.argv[2])); x6 = float(Fraction(sys.argv[3])); mu = float(Fraction(sys.argv[4]))
orb = json.load(open(sys.argv[5]))
xs = steady_state(k, x1, mu, x6)
ZIDX = [2, 3, 4, 6, 7, 8]
out = {}
for label, o in orb.items():
    x0 = np.array(o['x0']); T = o['T']
    s = solve_ivp(lambda t, x: g9(x, k), [0, T], x0, method='DOP853', rtol=1e-13, atol=1e-14, dense_output=True)
    ts = np.linspace(0, T, 20001); X = s.sol(ts).T
    print(label, 'closure error', np.linalg.norm(X[-1] - X[0]))
    for jj, i in enumerate(ZIDX):
        d = X[:, i] - xs[i]
        cr = np.where((d[:-1] < 0) & (d[1:] >= 0))[0]
        cr2 = np.where((d[:-1] > 0) & (d[1:] <= 0))[0]
        if len(cr) == 1:
            tc = ts[cr[0]]
            from scipy.optimize import brentq
            tc = brentq(lambda t: s.sol(t)[i] - xs[i], ts[cr[0]], ts[cr[0] + 1], xtol=1e-15)
            xc = s.sol(tc); f = g9(xc, k)
            print('  z%d (x%d): up-crossings %d down %d  x_i range [%.5g, %.5g] (ss %.5g)  f_i/|f| at crossing %.3f' % (
                jj, i + 1, len(cr), len(cr2), X[:, i].min(), X[:, i].max(), xs[i], f[i] / np.linalg.norm(f)))
            out.setdefault(label, {})[str(jj)] = {'x': xc.tolist(), 'c': xs[i]}
json.dump(out, open('section_points.json', 'w'), indent=1)
