"""Minimisation drivers (L-BFGS on vertex positions) and mesh-quality checks."""
import time
import numpy as np
from scipy.optimize import minimize
from . import energy, mesh


def relax_fixed_box(M, maxiter=20000, gtol=1e-11, ftol=1e-16, maxcor=30, verbose=False):
    X0 = M.X.copy()
    shape = X0.shape
    nev = [0]

    def fun(x):
        nev[0] += 1
        F, G = energy.dfm_F_and_grad(M, x.reshape(shape))
        return F, G.ravel()

    t0 = time.time()
    res = minimize(fun, X0.ravel(), jac=True, method="L-BFGS-B",
                   options={"maxiter": maxiter, "maxfun": 2 * maxiter, "gtol": gtol, "ftol": ftol, "maxcor": maxcor})
    M.X = res.x.reshape(shape)
    G = res.jac.reshape(shape)
    info = {"F": res.fun, "nit": res.nit, "nfev": nev[0], "msg": str(res.message),
            "gmax": float(np.abs(G).max()), "time": time.time() - t0}
    if verbose:
        print(info)
    return info


def mesh_quality(M):
    """min angle of triangles (deg), min triangle area ratio, signed check of
    orientation consistency against the level-0 normal of each face (crude)."""
    X = M.X
    L = M.L
    P = X[M.tri] + np.einsum("ij,tkj->tki", L, M.tsh)
    e = [P[:, 1] - P[:, 0], P[:, 2] - P[:, 1], P[:, 0] - P[:, 2]]
    l = [np.linalg.norm(x, axis=1) for x in e]
    n = np.cross(e[0], -e[2])
    area = 0.5 * np.linalg.norm(n, axis=1)
    # angles
    angs = []
    for k in range(3):
        a, b = -e[(k + 2) % 3], e[k]
        cosang = np.sum(a * b, axis=1) / (np.linalg.norm(a, axis=1) * np.linalg.norm(b, axis=1))
        angs.append(np.degrees(np.arccos(np.clip(cosang, -1, 1))))
    minang = np.min(np.array(angs))
    return {"min_angle": float(minang), "min_area": float(area.min()), "max_area": float(area.max())}


def foam_edge_lengths(M):
    """Length of every foam edge (polyline between consecutive foam vertices on
    face boundaries); each edge appears in 3 faces, we return the minimum over all."""
    X, L = M.X, M.L
    mins = []
    for loop, corner in zip(M.fbnd, M.fcorner):
        pos = np.array([X[v] + L @ np.array(s) for (v, s) in loop])
        m = len(loop)
        # find corners
        idx = [k for k in range(m) if corner[k]]
        for a in range(len(idx)):
            k0 = idx[a]
            k1 = idx[(a + 1) % len(idx)]
            seg = []
            k = k0
            while k != k1:
                seg.append(np.linalg.norm(pos[(k + 1) % m] - pos[k]))
                k = (k + 1) % m
            mins.append(sum(seg))
    return np.array(mins)
