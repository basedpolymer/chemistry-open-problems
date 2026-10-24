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


def restricted_basis(M):
    """Directions for restricted ('normal') motion: foam vertices 3 dof,
    foam-edge vertices 2 dof (plane normal to the edge), face-interior vertices
    1 dof (local surface normal).  Returns (start, cnt, dirs) with
    X_v = X0_v + sum_k y[start_v + k] * dirs[v, k]."""
    X, L = M.X, M.L
    nv = len(X)
    dirs = np.zeros((nv, 3, 3))
    cnt = np.zeros(nv, dtype=np.int64)
    # vertex normals from triangles (face-interior vertices only lie in one face)
    P = X[M.tri] + np.einsum("ij,tkj->tki", L, M.tsh)
    nrm = np.cross(P[:, 1] - P[:, 0], P[:, 2] - P[:, 0])
    vn = np.zeros((nv, 3))
    for k in range(3):
        np.add.at(vn, M.tri[:, k], nrm)
    # edge tangents from boundary loops
    tang = np.zeros((nv, 3))
    for loop in M.fbnd:
        m = len(loop)
        pos = np.array([X[v] + L @ np.array(s) for (v, s) in loop])
        for k in range(m):
            v = loop[k][0]
            if M.vtype[v] == 1:
                t = pos[(k + 1) % m] - pos[k - 1]
                # orientation of loops differs between faces; align sign
                if tang[v] @ t < 0:
                    t = -t
                tang[v] += t
    for v in range(nv):
        ty = M.vtype[v]
        if ty == 0:
            dirs[v] = np.eye(3); cnt[v] = 3
        elif ty == 1:
            t = tang[v] / np.linalg.norm(tang[v])
            tmp = np.array([1.0, 0, 0]) if abs(t[0]) < 0.9 else np.array([0, 1.0, 0])
            e1 = np.cross(t, tmp); e1 /= np.linalg.norm(e1)
            e2 = np.cross(t, e1)
            dirs[v, 0] = e1; dirs[v, 1] = e2; cnt[v] = 2
        else:
            n = vn[v] / np.linalg.norm(vn[v])
            dirs[v, 0] = n; cnt[v] = 1
    start = np.concatenate([[0], np.cumsum(cnt)[:-1]])
    return start, cnt, dirs


def relax(M, mode="full", maxiter=20000, gtol=1e-12, maxcor=30, stall_window=200, stall_tol=1e-13, verbose=False):
    """Minimise F over vertex positions.  mode 'full' (3 dof per vertex) or
    'restricted' (normal motion).  Stops on L-BFGS convergence or when the
    relative decrease of F over `stall_window` iterations is below stall_tol."""
    X0 = M.X.copy()
    if mode == "full":
        nd = X0.size
        def tox(y):
            return X0 + y.reshape(X0.shape)
        def tog(G):
            return G.ravel()
    else:
        start, cnt, dirs = restricted_basis(M)
        nd = int(cnt.sum())
        import scipy.sparse as sp
        vidx = np.repeat(np.arange(len(X0)), cnt)
        kidx = np.concatenate([np.arange(c) for c in cnt])
        D = dirs[vidx, kidx]               # (nd,3)
        rows = (3 * vidx[:, None] + np.arange(3)[None, :]).ravel()
        cols = np.repeat(np.arange(nd), 3)
        B = sp.csr_matrix((D.ravel(), (rows, cols)), shape=(X0.size, nd))
        BT = B.T.tocsr()
        def tox(y):
            return X0 + (B @ y).reshape(X0.shape)
        def tog(G):
            return BT @ G.ravel()
    hist = []

    def fun(y):
        F, G = energy.dfm_F_and_grad(M, tox(y))
        return F, tog(G)

    best = {"F": np.inf, "y": np.zeros(nd)}

    def fun_rec(y):
        F, g = fun(y)
        if F < best["F"]:
            best["F"] = F; best["y"] = y.copy()
        return F, g

    def callback(intermediate_result):
        hist.append(intermediate_result.fun)
        if len(hist) > stall_window and (hist[-stall_window - 1] - hist[-1]) < stall_tol * abs(hist[-1]):
            raise StopIteration

    t0 = time.time()
    res = minimize(fun_rec, np.zeros(nd), jac=True, method="L-BFGS-B", callback=callback,
                   options={"maxiter": maxiter, "maxfun": 3 * maxiter, "gtol": gtol, "ftol": 1e-17, "maxcor": maxcor})
    M.X = tox(best["y"])
    F, G = energy.dfm_F_and_grad(M, M.X)
    g = tog(G)
    info = {"F": F, "nit": res.nit, "msg": str(res.message), "gmax": float(np.abs(g).max()),
            "time": time.time() - t0, "ndof": nd}
    if verbose:
        print(info)
    return info
