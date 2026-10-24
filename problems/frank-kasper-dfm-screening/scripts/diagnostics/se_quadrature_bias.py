#!/usr/bin/env python3
"""Diagnostic: evaluate, on our relaxed meshes, the second moment with the
one-point facet-centroid rule of Reddy et al. (Eqs. B2-B3, as used in their
Surface Evolver runs) and compare with the exact polyhedral value.
usage: se_quadrature_bias.py RUNDIR [RUNDIR...]   (reads state.pkl of each level run)
"""
import os, sys, pickle
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import dfm  # noqa
import numpy as np
from dfm import energy

FOURPI = 4 * np.pi


def se_style(M):
    L = M.L
    P = M.X[M.tri] + np.einsum("ij,tkj->tki", L, M.tsh)
    nrm = 0.5 * np.cross(P[:, 1] - P[:, 0], P[:, 2] - P[:, 0])  # area * unit normal (i -> j)
    Xc = P.mean(axis=1)
    f = M.tface
    n = M.ncell
    Atot, V, Mo, S = energy.geometry(M)
    Ise = np.zeros(n)
    for side, cells, org, sg in ((0, M.fi[f], M.ref[M.fi[f]], 1.0),
                                 (1, M.fj[f], M.ref[M.fj[f]] + M.fsj[f] @ L.T, -1.0)):
        r = Xc - org
        nd = sg * np.sum(nrm * r, axis=1)
        I0 = np.bincount(cells, weights=nd * np.sum(r * r, axis=1) / 5.0, minlength=n)
        d = np.zeros((n, 3))
        for c in range(3):
            d[:, c] = np.bincount(cells, weights=nd * r[:, c] / 4.0, minlength=n)
        if side == 0:
            I0s, ds = I0, d
        else:
            I0s = I0s + I0
            ds = ds + d
    xc = ds / V[:, None]
    Ise = I0s - V * np.sum(xc ** 2, axis=1)
    Iex = S - np.sum(Mo ** 2, axis=1) / V
    Vbox = abs(np.linalg.det(L))
    R0 = (3 * Vbox / (FOURPI * n)) ** (1 / 3)
    A = (2 * Atot / n) / (FOURPI * R0 ** 2)
    Ie = (Iex.sum() / n) / (FOURPI * R0 ** 5 / 5)
    Is = (Ise.sum() / n) / (FOURPI * R0 ** 5 / 5)
    return A, Ie, Is, (A * A * Ie) ** (1 / 3), (A * A * Is) ** (1 / 3)


if __name__ == "__main__":
    for d in sys.argv[1:]:
        M = pickle.load(open(os.path.join(d, "state.pkl"), "rb"))
        A, Ie, Is, Fe, Fs = se_style(M)
        print("%s L%d ntri=%d  A=%.6f  I_exact=%.6f  I_centroid_rule=%.6f  F_exact=%.6f  F_with_centroid_rule=%.6f"
              % (os.path.basename(d), M.level, len(M.tri), A, Ie, Is, Fe, Fs))
