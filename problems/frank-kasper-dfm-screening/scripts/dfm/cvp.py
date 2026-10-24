"""Centroidal Voronoi partition by Lloyd iteration (flat cells), periodic box."""
import numpy as np
from . import topology, mesh, energy


def lloyd(L, P, maxit=2000, tol=1e-13):
    P = P.copy()
    Linv = np.linalg.inv(L)
    for it in range(maxit):
        tets, info = topology.periodic_delaunay(L, P, reach=1, check_degenerate=False)
        faces = topology.build_foam(L, P, tets)
        M = mesh.build_initial_mesh(L, P, tets, faces)
        Atot, V, Mo, S = energy.geometry(M)
        cen = M.ref + Mo / V[:, None]
        Pn = (cen @ Linv.T)
        Pn = Pn - np.floor(Pn)
        d = ((Pn - P + 0.5) % 1.0) - 0.5
        P = Pn
        if np.abs(L @ d.T).max() < tol:
            break
    d = energy.dfm_values(M)
    z = topology.coordination(len(P), faces)
    return P, d, z, it
