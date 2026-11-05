# Exploration script (floating point; not a certificate). Environment: WSL Fedora, Python 3.12, numpy, scipy, numba.
# Run from scripts/exploration/ (it imports the sibling modules). Typical runtime: seconds to minutes.
"""Locate the unstable and stable periodic orbits at fixed mu (float), in z = (x3,x4,x5,x7,x8,x9) coordinates.
Output: points on a coordinate section, period, multipliers. Used to seed the CAPD proof."""
import sys, json
import numpy as np
from fractions import Fraction
from scipy.integrate import solve_ivp
from fastmodel import g9, J9, var_rhs
from po_cont2 import PO, E, P

k = [float(Fraction(s)) for s in sys.argv[1].split(',')]
x1 = float(Fraction(sys.argv[2])); x6 = float(Fraction(sys.argv[3]))
mu = float(Fraction(sys.argv[4]))
branch = np.load(sys.argv[5])          # continuation output with u = (z0, T, mu, mults)
po = PO(k, x1, x6, rtol=1e-13, atol=1e-14)
res = {}
for label, sel in [('unstable', lambda r: r[8] > 1 + 1e-9), ('stable', lambda r: r[8] <= 1 + 1e-9 and r[9] < 1 - 1e-9)]:
    rows = [r for r in branch if sel(r)]
    if not rows or min(abs(r[7]-mu) for r in rows) > 0.6:
        continue
    r = min(rows, key=lambda r: abs(r[7] - mu))
    u = r[:8].copy()
    # move to the requested mu: Newton with mu fixed
    u[7] = mu
    zref = u[:6].copy(); nref = P @ g9(po.xs(mu) + E @ u[:6], po.k); nref /= np.linalg.norm(nref)
    u, J, Mz, ok = po.newton(u, zref, nref, tol=1e-13)
    xs = po.xs(mu)
    x0 = xs + E @ u[:6]
    mults = np.linalg.eigvals(Mz)
    print(label, 'ok', ok, 'T=%.12f' % u[6], 'multipliers', np.sort_complex(mults)[::-1][:6])
    print('   x0 =', x0)
    res[label] = {'x0': x0.tolist(), 'T': u[6], 'mults': [abs(m) for m in mults]}
json.dump(res, open('orbits_mu_%s.json' % sys.argv[4].replace('/', '_'), 'w'), indent=1)
