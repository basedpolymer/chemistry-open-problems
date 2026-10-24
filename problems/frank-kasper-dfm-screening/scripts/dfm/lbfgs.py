"""Plain L-BFGS (two-loop recursion) with Armijo backtracking, for smooth
objectives with many variables.  Deterministic."""
import numpy as np


def lbfgs(fg, x0, m=20, maxiter=10000, gtol=1e-13, stall_window=100, stall_tol=1e-14, callback=None, noise=4e-15):
    x = x0.copy()
    f, g = fg(x)
    S, Y, RHO = [], [], []
    hist = [f]
    nfev = 1
    nback = 0
    it = 0
    msg = "maxiter"
    for it in range(1, maxiter + 1):
        # two-loop recursion
        q = g.copy()
        alpha = []
        for s, y, rho in zip(reversed(S), reversed(Y), reversed(RHO)):
            a = rho * (s @ q)
            alpha.append(a)
            q -= a * y
        if S:
            gamma = (S[-1] @ Y[-1]) / (Y[-1] @ Y[-1])
        else:
            gamma = 1e-3 / max(np.abs(g).max(), 1e-300)
        r = gamma * q
        for (s, y, rho), a in zip(zip(S, Y, RHO), reversed(alpha)):
            b = rho * (y @ r)
            r += (a - b) * s
        d = -r
        gd = g @ d
        if gd >= 0:  # not a descent direction: reset
            S, Y, RHO = [], [], []
            d = -gamma * g if gamma > 0 else -g
            gd = g @ d
        if -gd < noise * abs(f):
            msg = "noise floor"
            break
        t = 1.0
        accepted = False
        for ls in range(30):
            xn = x + t * d
            fn, gn = fg(xn)
            nfev += 1
            if fn <= f + 1e-4 * t * gd:
                accepted = True
                break
            t *= 0.5
            nback += 1
        if not accepted:
            msg = "line search failed"
            break
        s = xn - x
        y = gn - g
        sy = s @ y
        x, f, g = xn, fn, gn
        if sy > 1e-30 * (np.linalg.norm(s) * np.linalg.norm(y) + 1e-300):
            S.append(s); Y.append(y); RHO.append(1.0 / sy)
            if len(S) > m:
                S.pop(0); Y.pop(0); RHO.pop(0)
        hist.append(f)
        if callback is not None:
            callback(it, f, g)
        if np.abs(g).max() < gtol:
            msg = "gtol"
            break
        if len(hist) > stall_window and (hist[-stall_window - 1] - hist[-1]) <= stall_tol * abs(hist[-1]):
            msg = "stall"
            break
    return x, f, g, {"nit": it, "nfev": nfev, "msg": msg, "nback": nback}
