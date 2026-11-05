# Exploration script (floating point; not a certificate). Environment: WSL Fedora, Python 3.12, numpy, scipy, numba.
# Run from scripts/exploration/ (it imports the sibling modules). Typical runtime: seconds to minutes.
"""Independent check of l1 sign by integrating the 9-dim ODE at the Hopf point and slightly on the stable side."""
import sys
import numpy as np
from scipy.integrate import solve_ivp
from model import *

a = np.load(sys.argv[1]); row = int(sys.argv[2])
r = a[a[:, 13] > 0][row] if sys.argv[3] == 'pos' else a[a[:, 13] < 0][row]
k = r[:10]; x1, x2, x6 = r[10:13]
print('k =', k, '\nx1,x2,x6 =', x1, x2, x6, ' l1 =', r[13], 'omega =', r[14])
x = steady_state(k, x1, x2, x6)
A = reduced_A(k, x)
ev, R = np.linalg.eig(A)
print('eigs at Hopf:', np.sort_complex(ev))
i = np.argmax(ev.imag)
q = V @ R[:, i]; q = np.real(q) / np.linalg.norm(np.real(q))
T = 2 * np.pi / ev[i].imag


def rhs(t, y):
    return g(k, y)


def amp_series(x0, xs, ncyc):
    """integrate and return the distance to xs at successive returns (sampled once per period)"""
    sol = solve_ivp(rhs, [0, ncyc * T], x0, method='LSODA', rtol=1e-11, atol=1e-14,
                    dense_output=True)
    ts = np.linspace(0, ncyc * T, ncyc * 40 + 1)
    d = np.linalg.norm(sol.sol(ts).T - xs, axis=1)
    return np.array([d[j * 40:(j + 1) * 40].max() for j in range(ncyc)])

scale = np.linalg.norm(x)
for eps in [1e-3, 3e-3, 1e-2, 3e-2]:
    x0 = x + eps * scale * q
    if (x0 <= 0).any():
        print('eps', eps, 'leaves orthant'); continue
    am = amp_series(x0, x, 30)
    print('at Hopf  eps=%g  amplitude ratio after 30 periods: %.6f' % (eps, am[-1] / am[0]))
