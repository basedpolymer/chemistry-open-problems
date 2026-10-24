"""Equal-volume constrained problems (validation only):
  objective 'area' : minimise total area with V_i = V0 (Kelvin problem),
  objective 'dfm'  : minimise F with V_i = V0 ("Equal" rows of Reddy et al.).
Augmented Lagrangian, inner restricted L-BFGS, fixed box."""
import time
import numpy as np
from . import energy, relax
from .lbfgs import lbfgs

FOURPI = 4 * np.pi


def _obj_grad(M, X, objective, lam, mu, V0):
    Atot, V, Mo, S = energy.geometry(M, X)
    n = M.ncell
    c = V - V0
    if objective == "area":
        f = Atot
        gA = 1.0
        gV = lam + mu * c
        G = energy._grad(X, M.L, M.tri, M.tsh, M.tface, M.fi, M.fj, M.fsj, M.ref, n, gA, gV,
                         np.zeros((n, 3)), np.zeros(n))
    else:
        F, G0 = energy.dfm_F_and_grad(M, X)
        f = F
        G = G0 + energy.area_volume_grad(M, X, 0.0, lam + mu * c)
    val = f + lam @ c + 0.5 * mu * c @ c
    return val, G, c


def solve(M, objective="area", outer=30, mu=None, tol=1e-12, cycles=2, maxiter=20000):
    n = M.ncell
    Vbox = abs(np.linalg.det(M.L))
    V0 = Vbox / n
    lam = np.zeros(n)
    R0 = (3 * V0 / FOURPI) ** (1 / 3)
    if mu is None:
        mu = (10.0 / R0 ** 4) if objective == "area" else (10.0 / (V0 ** 2))
    t0 = time.time()
    for cyc in range(cycles):
        for it in range(outer):
            X0 = M.X.copy()
            import scipy.sparse as sp
            start, cnt, dirs = relax.restricted_basis(M)
            nd = int(cnt.sum())
            vidx = np.repeat(np.arange(len(X0)), cnt)
            kidx = np.concatenate([np.arange(c) for c in cnt])
            D = dirs[vidx, kidx]
            rows = (3 * vidx[:, None] + np.arange(3)[None, :]).ravel()
            cols = np.repeat(np.arange(nd), 3)
            B = sp.csr_matrix((D.ravel(), (rows, cols)), shape=(X0.size, nd))
            BT = B.T.tocsr()

            def fg(y):
                X = X0 + (B @ y).reshape(X0.shape)
                v, G, c = _obj_grad(M, X, objective, lam, mu, V0)
                return v, BT @ G.ravel()

            y, f, g, res = lbfgs(fg, np.zeros(nd), m=20, maxiter=maxiter, noise=1e-15)
            M.X = X0 + (B @ y).reshape(X0.shape)
            Atot, V, Mo, S = energy.geometry(M)
            c = V - V0
            lam = lam + mu * c
            if np.abs(c).max() < tol * V0:
                break
        if cyc < cycles - 1:
            relax.smooth_tangential(M, sweeps=3)
    Atot, V, Mo, S = energy.geometry(M)
    d = energy.dfm_values(M)
    R0 = (3 * V0 / FOURPI) ** (1 / 3)
    return {"A": d["A"], "I": d["I"], "F": d["F"], "maxdV": float(np.abs(V - V0).max() / V0),
            "time": time.time() - t0}
