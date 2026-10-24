"""Periodic Delaunay / Voronoi topology for a periodic point set.

Conventions
-----------
L : (3,3) array, columns are the lattice vectors a, b, c (Cartesian).
P : (n,3) fractional coordinates of the sites (cell centres), in [0,1).
A periodic object is referred to as (index, shift) where shift is an integer
3-vector: position = L @ (P[index] + shift).

The Delaunay complex is computed with scipy (qhull) on a block of periodic
images, then reduced to periodic classes.  The dual (Voronoi) complex gives the
foam topology used by the DFM: foam vertices = Delaunay tetrahedra, foam faces
= Delaunay edges (pairs of neighbouring cells).
"""
import itertools
import numpy as np
from scipy.spatial import Delaunay


def _canon_simplex(verts):
    """verts: list of (k, (s0,s1,s2)). Return canonical key (tuple) and the
    translation t such that canonical = verts translated by t."""
    best = None
    for (k0, s0) in verts:
        t = (-s0[0], -s0[1], -s0[2])
        tr = sorted((k, (s[0] + t[0], s[1] + t[1], s[2] + t[2])) for (k, s) in verts)
        key = tuple(tr)
        if best is None or key < best[0]:
            best = (key, t)
    return best


def circumcenter(pts):
    """pts (4,3) -> circumcentre (3,), radius"""
    a = pts[0]
    M = 2.0 * (pts[1:] - a)
    rhs = np.sum(pts[1:] ** 2 - a ** 2, axis=1)
    c = np.linalg.solve(M, rhs)
    return c, np.linalg.norm(c - a)


def periodic_delaunay(L, P, reach=1, check_degenerate=True, tol=1e-9):
    """Return list of canonical tetrahedra; each is a tuple of 4 (k, shift)."""
    n = len(P)
    shifts = list(itertools.product(range(-reach, reach + 1), repeat=3))
    pts, lab = [], []
    for s in shifts:
        for k in range(n):
            pts.append(L @ (P[k] + np.array(s)))
            lab.append((k, s))
    pts = np.array(pts)
    vol_scale = abs(np.linalg.det(L)) / n
    tri = Delaunay(pts, qhull_options="Qbb Qc Qz Q12")
    tets = {}
    for simp in tri.simplices:
        vs = [lab[i] for i in simp]
        if not any(s == (0, 0, 0) for (_, s) in vs):
            continue
        Xs = pts[simp]
        vol = abs(np.linalg.det(Xs[1:] - Xs[0])) / 6.0
        if vol < 1e-10 * vol_scale:
            continue
        key, t = _canon_simplex(vs)
        if key in tets:
            continue
        # keep only tetrahedra whose circumsphere is well inside the image block
        X = np.array([L @ (P[k] + np.array(s)) for (k, s) in key])
        c, R = circumcenter(X)
        tets[key] = (c, R)
    tets_list = list(tets.keys())
    info = {"n_tets": len(tets_list)}
    vsum = 0.0
    for key in tets_list:
        X = np.array([L @ (P[k] + np.array(s)) for (k, s) in key])
        vsum += abs(np.linalg.det(X[1:] - X[0])) / 6.0
    info["vol_ratio"] = vsum / abs(np.linalg.det(L))
    # sanity: each canonical tet must have been seen with all its vertices in range;
    # check degeneracy: count points on / inside circumsphere
    if check_degenerate:
        ndeg = 0
        worst = np.inf
        for key in tets_list:
            c, R = tets[key]
            d = np.linalg.norm(pts - c, axis=1)
            inside = np.sum(d < R * (1 - tol))
            on = np.sum(np.abs(d - R) <= R * tol)
            if inside > 0:
                raise RuntimeError("non-Delaunay tetrahedron found (image block too small?)")
            if on > 4:
                ndeg += 1
            # smallest relative gap to the next point outside
            dd = np.sort(d)
            worst = min(worst, (dd[4] - R) / R)
        info["n_degenerate"] = ndeg
        info["min_rel_gap"] = worst
    return tets_list, info


def build_foam(L, P, tets):
    """Build the periodic foam (dual of the periodic triangulation).

    Returns dict with
      tets      : list of canonical tets
      faces     : list of dicts {i, j, sj, ring: [(tet_index, shift)...]} -- the
                  Voronoi face between cell i (image 0) and cell j (image sj); the
                  ring lists foam vertices (tets) in cyclic order, with the shift that
                  places the canonical tet at that position in the frame of cell i.
                  Orientation: the ring is ordered counter-clockwise when seen from
                  cell j, i.e. the normal (right-hand rule) points from i to j.
    """
    n = len(P)
    tet_index = {key: idx for idx, key in enumerate(tets)}
    # map Delaunay edges -> list of (tet idx, translation)
    edge_tets = {}
    for ti, key in enumerate(tets):
        for a, b in itertools.combinations(range(4), 2):
            (ka, sa), (kb, sb) = key[a], key[b]
            # express in frame where the lower endpoint (by k, then by shift) is at image 0
            # canonical edge: (i, j, sj) with i <= j; if i == j choose sj lexicographically > 0
            for (k1, s1), (k2, s2) in (((ka, sa), (kb, sb)), ((kb, sb), (ka, sa))):
                t = tuple(-x for x in s1)
                sj = tuple(s2[c] + t[c] for c in range(3))
                if k1 < k2 or (k1 == k2 and sj > (0, 0, 0)):
                    ekey = (k1, k2, sj)
                    edge_tets.setdefault(ekey, []).append((ti, t))
    faces = []
    for (i, j, sj), lst in edge_tets.items():
        # order the ring of tets around the edge cyclically
        pi = L @ P[i]
        pj = L @ (P[j] + np.array(sj))
        axis = pj - pi
        axis /= np.linalg.norm(axis)
        # orthonormal frame
        tmp = np.array([1.0, 0, 0]) if abs(axis[0]) < 0.9 else np.array([0, 1.0, 0])
        e1 = np.cross(axis, tmp); e1 /= np.linalg.norm(e1)
        e2 = np.cross(axis, e1)
        mid = 0.5 * (pi + pj)
        angs = []
        for (ti, t) in lst:
            X = np.array([L @ (P[k] + np.array(s) + np.array(t)) for (k, s) in tets[ti]])
            c = X.mean(axis=0) - mid
            angs.append(np.arctan2(c @ e2, c @ e1))
        order = np.argsort(angs)
        ring = [lst[o] for o in order]
        # ring is counter-clockwise around axis (i->j) => normal points i->j
        faces.append({"i": i, "j": j, "sj": np.array(sj), "ring": ring})
    return faces


def coordination(n, faces):
    z = np.zeros(n, dtype=int)
    for f in faces:
        z[f["i"]] += 1
        z[f["j"]] += 1
    return z


def check_manifold(n, tets):
    """Each periodic triangle must lie in exactly two tetrahedra."""
    tri_count = {}
    for key in tets:
        for comb in itertools.combinations(range(4), 3):
            vs = [key[c] for c in comb]
            k, _ = _canon_simplex(vs)
            tri_count[k] = tri_count.get(k, 0) + 1
    bad = [k for k, v in tri_count.items() if v != 2]
    return len(bad) == 0, len(tri_count)
