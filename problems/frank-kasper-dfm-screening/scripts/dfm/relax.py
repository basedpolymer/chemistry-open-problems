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
        for c in range(3):
            vn[:, c] += np.bincount(M.tri[:, k], weights=nrm[:, c], minlength=nv)
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
    ty = M.vtype
    # foam vertices
    m0 = ty == 0
    dirs[m0] = np.eye(3)[None, :, :]
    cnt[m0] = 3
    # foam-edge vertices
    m1 = ty == 1
    t = tang[m1] / np.linalg.norm(tang[m1], axis=1)[:, None]
    tmp = np.where((np.abs(t[:, 0]) < 0.9)[:, None], np.array([1.0, 0, 0]), np.array([0, 1.0, 0]))
    e1 = np.cross(t, tmp); e1 /= np.linalg.norm(e1, axis=1)[:, None]
    e2 = np.cross(t, e1)
    dirs[m1, 0] = e1; dirs[m1, 1] = e2; cnt[m1] = 2
    # face-interior vertices
    m2 = ty == 2
    dirs[m2, 0] = vn[m2] / np.linalg.norm(vn[m2], axis=1)[:, None]
    cnt[m2] = 1
    start = np.concatenate([[0], np.cumsum(cnt)[:-1]])
    return start, cnt, dirs


def _sym_from6(e):
    return np.array([[e[0], e[3], e[4]], [e[3], e[1], e[5]], [e[4], e[5], e[2]]])


def relax(M, mode="restricted", box=False, maxiter=20000, gtol=1e-13, maxcor=20, stall_window=100,
          stall_tol=1e-14, verbose=False):
    """Minimise F over vertex positions ('full': 3 dof per vertex, 'restricted':
    normal motion) and, if box=True, over a symmetric homogeneous strain of the
    periodic box (shape and volume; F is scale invariant)."""
    from .lbfgs import lbfgs
    X0 = M.X.copy()
    L0 = M.L.copy()
    ref0 = M.ref.copy()
    if mode == "full":
        nd = X0.size
        def disp(y):
            return y.reshape(X0.shape)
        def tog(G):
            return G.ravel()
    else:
        import scipy.sparse as sp
        start_, cnt, dirs = restricted_basis(M)
        nd = int(cnt.sum())
        vidx = np.repeat(np.arange(len(X0)), cnt)
        kidx = np.concatenate([np.arange(c) for c in cnt])
        D = dirs[vidx, kidx]
        rows = (3 * vidx[:, None] + np.arange(3)[None, :]).ravel()
        cols = np.repeat(np.arange(nd), 3)
        B = sp.csr_matrix((D.ravel(), (rows, cols)), shape=(X0.size, nd))
        BT = B.T.tocsr()
        def disp(y):
            return (B @ y).reshape(X0.shape)
        def tog(G):
            return BT @ G.ravel()

    def apply(z):
        y = z[:nd]
        if box:
            A = np.eye(3) + _sym_from6(z[nd:])
        else:
            A = np.eye(3)
        X = (X0 + disp(y)) @ A.T
        return X, A

    def fun(z):
        X, A = apply(z)
        if box:
            M.L = A @ L0
            M.ref = ref0 @ A.T
            F, G, dE = energy.dfm_F_grad_virial(M, X)
            gy = tog(G @ A)
            dEt = dE @ np.linalg.inv(A).T
            ge = np.array([dEt[0, 0], dEt[1, 1], dEt[2, 2], dEt[0, 1] + dEt[1, 0], dEt[0, 2] + dEt[2, 0],
                           dEt[1, 2] + dEt[2, 1]])
            return F, np.concatenate([gy, ge])
        F, G = energy.dfm_F_and_grad(M, X)
        return F, tog(G)

    t0 = time.time()
    z0 = np.zeros(nd + (6 if box else 0))
    z, F, g, res = lbfgs(fun, z0, m=maxcor, maxiter=maxiter, gtol=gtol,
                         stall_window=stall_window, stall_tol=stall_tol)
    X, A = apply(z)
    M.X = X
    M.L = A @ L0
    M.ref = ref0 @ A.T
    F, G = energy.dfm_F_and_grad(M, M.X)
    info = {"F": F, "nit": res["nit"], "nfev": res["nfev"], "nback": res["nback"], "msg": res["msg"], "gmax": float(np.abs(g).max()),
            "time": time.time() - t0, "ndof": nd}
    if verbose:
        print(info)
    return info


def _neighbours(M):
    """neighbour lists for tangential smoothing: interior vertices -> all mesh
    neighbours (within their face); foam-edge vertices -> the two neighbours along
    the foam edge.  Neighbour references carry shifts relative to the vertex."""
    nv = len(M.X)
    nb = [dict() for _ in range(nv)]
    for t in range(len(M.tri)):
        for k in range(3):
            a = M.tri[t, k]
            if M.vtype[a] != 2:
                continue
            for kk in (1, 2):
                b = M.tri[t, (k + kk) % 3]
                d = tuple(M.tsh[t, (k + kk) % 3] - M.tsh[t, k])
                nb[a][(b, d)] = 1
    for loop in M.fbnd:
        m = len(loop)
        for k in range(m):
            a, sa = loop[k]
            if M.vtype[a] != 1:
                continue
            for kk in (-1, 1):
                b, sb = loop[(k + kk) % m]
                d = tuple(np.array(sb) - np.array(sa))
                nb[a][(b, d)] = 1
    idx, sh, ptr = [], [], [0]
    for v in range(nv):
        for (b, d) in nb[v]:
            idx.append(b); sh.append(d)
        ptr.append(len(idx))
    return np.array(idx, dtype=np.int64), np.array(sh, dtype=np.int64).reshape(-1, 3), np.array(ptr, dtype=np.int64)


def smooth_tangential(M, sweeps=3, omega=0.5):
    """Laplacian smoothing projected on the tangent plane (interior vertices) or
    on the edge tangent (foam-edge vertices).  Foam vertices are not moved."""
    if getattr(M, "_nbcache", None) is None or M._nbcache[0] != len(M.X):
        idx, sh, ptr = _neighbours(M)
        M._nbcache = (len(M.X), idx, sh, ptr)
    _, idx, sh, ptr = M._nbcache
    deg = np.diff(ptr)
    owner = np.repeat(np.arange(len(M.X)), deg)
    for _ in range(sweeps):
        start, cnt, dirs = restricted_basis(M)
        X = M.X
        nbpos = X[idx] + sh @ M.L.T
        avg = np.zeros_like(X)
        for c in range(3):
            avg[:, c] = np.bincount(owner, weights=nbpos[:, c], minlength=len(X))
        dg = np.maximum(deg, 1)[:, None]
        d = avg / dg - X
        ty = M.vtype
        n = dirs[:, 0]
        t = np.cross(dirs[:, 0], dirs[:, 1])
        dint = d - np.sum(d * n, axis=1)[:, None] * n
        dedge = np.sum(d * t, axis=1)[:, None] * t
        disp = np.where((ty == 2)[:, None], dint, np.where((ty == 1)[:, None], dedge, 0.0))
        M.X = X + omega * disp


def relax_cycles(M, cycles=3, sweeps=3, **kw):
    infos = []
    for c in range(cycles):
        infos.append(relax(M, mode="restricted", **kw))
        if c < cycles - 1:
            smooth_tangential(M, sweeps=sweeps)
    return infos


def cone_check(M):
    """Sufficient floating-point star-shapedness diagnostic at volume centroids.

    Positive boundary cones, with the manifold and degree-one hypotheses,
    certify a partition. A negative value does NOT prove self-intersection:
    another apex can lie in the cell's kernel when its centroid does not.
    See the exact Y04 counterexample in check_barrier_witnesses.py.
    Returns the minimum signed cone volume divided by mean cell volume.
    """
    from .energy import geometry
    Atot, V, Mo, S = geometry(M)
    cen = M.ref + Mo / V[:, None]
    L = M.L
    P = M.X[M.tri] + np.einsum("ij,tkj->tki", L, M.tsh)
    f = M.tface
    oi = cen[M.fi[f]]
    oj = cen[M.fj[f]] + M.fsj[f] @ L.T
    def vol(o, sg):
        a = P[:, 0] - o; b = P[:, 1] - o; c = P[:, 2] - o
        return sg * np.einsum("ij,ij->i", a, np.cross(b, c)) / 6.0
    vi = vol(oi, 1.0)
    vj = vol(oj, -1.0)
    vmean = abs(np.linalg.det(L)) / M.ncell
    return float(min(vi.min(), vj.min()) / vmean)
