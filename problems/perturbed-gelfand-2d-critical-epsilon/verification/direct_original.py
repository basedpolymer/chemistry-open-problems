"""Verifier's NON-RIGOROUS direct check on the original problem (no Du-Lou reduction).

Delta u + lambda exp(u/(1+eps u)) = 0 in the unit disk, u = 0 on the boundary.
Shooting in x = log(s^2/4):  w_x = P, P_x = -e^x f(w), w(-inf) = alpha;  v = w_alpha solves the
variational system.  lambda = 4 exp(X) at the first zero X of w, and sign(lambda_alpha) = sign(v(X))
because lambda_alpha = -4 e^X v(X) / P(X) with P(X) < 0.
For each eps we scan alpha on a fine grid, count sign changes of v(X), refine each fold with brentq,
and report the transformed amplitude a = eps + eps^2 alpha of each fold.
Python 3 + numpy + scipy, one thread, about one minute.  Floating point only: evidence, not proof.
"""
import os
os.environ.setdefault('OMP_NUM_THREADS', '1')
os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
import math
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq


def shoot(alpha, eps):
    def rhs(x, y):
        w, P, v, Q = y
        D = 1.0 + eps * w
        f = math.exp(w / D) if w > -1.0 / eps * 0.5 else 0.0
        fw = f / (D * D)
        ex = math.exp(x)
        return [P, -ex * f, Q, -ex * fw * v]

    def ev(x, y):
        return y[0]
    ev.terminal = True
    ev.direction = -1
    # start where e^x f(alpha) is negligible: f <= e^{1/eps}
    x0 = -1.0 / eps - 40.0
    sol = solve_ivp(rhs, (x0, 60.0), [alpha, 0.0, 1.0, 0.0], method='DOP853',
                    rtol=1e-12, atol=1e-14, events=ev)
    if not sol.t_events[0].size:
        raise RuntimeError('no zero')
    X = sol.t_events[0][0]
    w, P, v, Q = sol.y_events[0][0]
    return X, v, P


def lam(alpha, eps):
    return 4.0 * math.exp(shoot(alpha, eps)[0])


def folds(eps, amax):
    grid = np.concatenate([np.geomspace(1e-3, 1.0, 60, endpoint=False),
                           np.arange(1.0, amax, 0.01 if amax < 50 else (0.05 if amax < 200 else 0.2))])
    vals = [shoot(a, eps)[1] for a in grid]
    out = []
    for a1, a2, v1, v2 in zip(grid[:-1], grid[1:], vals[:-1], vals[1:]):
        if v1 * v2 < 0:
            r = brentq(lambda al: shoot(al, eps)[1], a1, a2, xtol=1e-12)
            out.append((r, 'max' if v1 > 0 else 'min', lam(r, eps)))
    return out, vals[0] > 0, vals[-1] > 0


def main():
    cases = [(0.05, 1500.0), (0.15, 200.0), (0.23, 60.0), (0.2420, 30.0), (0.2422, 30.0), (0.245, 30.0)]
    for eps, amax in cases:
        fl, first_pos, last_pos = folds(eps, amax)
        print(f'eps = {eps}: {len(fl)} fold(s); lambda_alpha > 0 at alpha=1e-3: {first_pos}, '
              f'at alpha={amax:g}: {last_pos}')
        for r, kind, L in fl:
            a = eps + eps * eps * r
            print(f'    alpha = {r:.10f}  ({kind} of lambda, lambda = {L:.10f})   a = eps+eps^2 alpha = {a:.8f}')
    # second fold for small eps in the transformed amplitude
    print('second fold, transformed amplitude, small eps:')
    for eps in [0.08, 0.04, 0.02, 0.01, 0.005]:
        amin = (1.4 - eps) / eps**2
        amax_ = (1.7 - eps) / eps**2
        grid = np.linspace(amin, amax_, 121)
        vals = [shoot(a, eps)[1] for a in grid]
        roots = [brentq(lambda al: shoot(al, eps)[1], a1, a2, xtol=1e-9)
                 for a1, a2, v1, v2 in zip(grid[:-1], grid[1:], vals[:-1], vals[1:]) if v1 * v2 < 0]
        print(f'    eps = {eps}: folds with a in (1.4,1.7): ' +
              ', '.join(f'alpha = {r:.4f}, a = {eps + eps*eps*r:.6f}' for r in roots))


if __name__ == '__main__':
    main()
