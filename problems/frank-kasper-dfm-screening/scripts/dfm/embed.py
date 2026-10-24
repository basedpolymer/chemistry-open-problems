"""Re-embedding of a periodic FK contact graph (for catalogue entries whose
CGD coordinates do not reproduce the topology as a Voronoi diagram).

Minimise  E = sum_contacts (|x_j + L s - x_i| - (r_Zi + r_Zj))^2 - mu sum_tets log(vol_t / v_ref)
over the site positions and the lattice matrix (upper triangular), starting
from the given embedded triangulation.  Radii r_Z (units of the mean cell
volume^(1/3)) were fitted on the contact distances of 11 real FK crystals:
r12=0.519, r14=0.594, r15=0.630, r16=0.651.  The log barrier keeps every
tetrahedron positively oriented, so the result is again an embedded
triangulation of the same contact graph."""
import numpy as np
from scipy.optimize import minimize

RADII = {12: 0.519, 14: 0.594, 15: 0.630, 16: 0.651}


def reembed(L, P, edges, tets, coord, mu=1e-3, maxiter=5000):
    n = len(P)
    E = np.array([(i, j) for (i, j, s) in edges]); S = np.array([s for (i, j, s) in edges], dtype=float)
    ell = np.array([RADII[coord[i]] + RADII[coord[j]] for (i, j) in E])
    T = np.array([[k for (k, s) in key] for key in tets]); TS = np.array([[s for (k, s) in key] for key in tets], dtype=float)
    # orient every tetrahedron positively in the starting embedding
    Xt0 = (P[T] + TS) @ L.T
    v0 = np.einsum("ij,ij->i", Xt0[:, 1] - Xt0[:, 0], np.cross(Xt0[:, 2] - Xt0[:, 0], Xt0[:, 3] - Xt0[:, 0]))
    if np.any(v0 == 0):
        raise ValueError("degenerate starting tetrahedron")
    neg = v0 < 0
    T[neg, 2], T[neg, 3] = T[neg, 3].copy(), T[neg, 2].copy()
    TS[neg, 2], TS[neg, 3] = TS[neg, 3].copy(), TS[neg, 2].copy()
    # scale: unit mean cell volume
    sc = (abs(np.linalg.det(L)) / n) ** (1.0 / 3.0)
    L0 = L / sc
    iu = np.triu_indices(3)
    # make L upper triangular by QR-like rotation (LQ of L^T)
    q, rr = np.linalg.qr(L0)
    sg = np.sign(np.diag(rr)); rr = sg[:, None] * rr
    Lu = rr  # upper triangular, same metric as L0
    x0 = np.concatenate([P.ravel(), Lu[iu]])
    vref = 1.0 / 6.0 * 0.1

    def unpack(x):
        Pp = x[:3 * n].reshape(n, 3)
        Lm = np.zeros((3, 3)); Lm[iu] = x[3 * n:]
        return Pp, Lm

    def f(x):
        Pp, Lm = unpack(x)
        d = (Pp[E[:, 1]] + S - Pp[E[:, 0]]) @ Lm.T          # cartesian contact vectors
        ln = np.linalg.norm(d, axis=1)
        r = ln - ell
        e = np.sum(r * r)
        g_d = (2 * r / ln)[:, None] * d                          # dE/dd
        # tets
        Xt = Pp[T] + TS                                          # fractional (t,4,3)
        e1 = (Xt[:, 1] - Xt[:, 0]) @ Lm.T; e2 = (Xt[:, 2] - Xt[:, 0]) @ Lm.T; e3 = (Xt[:, 3] - Xt[:, 0]) @ Lm.T
        c23 = np.cross(e2, e3)
        vol = np.einsum("ij,ij->i", e1, c23) / 6.0
        if np.any(vol <= 0):
            return np.inf, np.zeros_like(x)
        e += -mu * np.sum(np.log(vol / vref))
        gv = -mu / vol / 6.0                                      # dE/d(6 vol)
        g_e1 = gv[:, None] * c23; g_e2 = gv[:, None] * np.cross(e3, e1); g_e3 = gv[:, None] * np.cross(e1, e2)
        # gradients w.r.t. fractional positions and L
        gP = np.zeros((n, 3)); gL = np.zeros((3, 3))
        gfr = g_d @ Lm                                           # d/d(fractional diff)
        np.add.at(gP, E[:, 1], gfr); np.add.at(gP, E[:, 0], -gfr)
        fd = Pp[E[:, 1]] + S - Pp[E[:, 0]]
        gL += g_d.T @ fd
        for (ge, a) in ((g_e1, 1), (g_e2, 2), (g_e3, 3)):
            gf = ge @ Lm
            np.add.at(gP, T[:, a], gf); np.add.at(gP, T[:, 0], -gf)
            fe = Xt[:, a] - Xt[:, 0]
            gL += ge.T @ fe
        return e, np.concatenate([gP.ravel(), gL[iu]])

    from .lbfgs import lbfgs
    x, fval, g, info = lbfgs(f, x0, m=20, maxiter=maxiter, gtol=1e-9, stall_window=200, stall_tol=1e-12, noise=1e-16)
    Pp, Lm = unpack(x)
    # positions are NOT wrapped into [0,1): the tetrahedra refer to images relative to them
    return Lm, Pp, fval, info["nit"]
