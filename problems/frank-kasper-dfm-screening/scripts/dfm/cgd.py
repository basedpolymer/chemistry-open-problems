"""Reader for the Systre-style CGD files of the Frank-Kasper catalogue
(Dutour Sikiric, Delgado-Friedrichs, Deza, Acta Cryst. A 66 (2010) 602;
data: Zenodo 10.5281/zenodo.13691342, folder SpaceFullerene/ListDATApaper).

Each file contains one 'crystal' block: cell parameters, space group, atoms
(id, coordination number, fractional coordinates) and, for each atom, the
fractional coordinates of its neighbours ('edge' lines).  Atoms are the centres
of the Frank-Kasper cells; edges are the contacts (shared faces) between cells.
"""
import numpy as np
from fractions import Fraction
from .structures import lattice_from_params

# Systre group names -> gemmi names (Systre uses the first ITA setting, origin
# choice 2 for centrosymmetric groups with two origins; checked by edge matching)
_ALIAS = {}


def parse_cgd(path):
    txt = open(path).read().split("\n")
    blk = {"atoms": [], "edges": {}}
    for line in txt:
        w = line.split()
        if not w:
            continue
        k = w[0].lower()
        if k == "name":
            blk["name"] = w[1]
        elif k == "cell":
            blk["cell"] = tuple(float(x) for x in w[1:7])
        elif k == "group":
            blk["group"] = w[1]
        elif k == "atom":
            blk["atoms"].append((int(w[1]), int(w[2]), np.array([float(Fraction(x)) for x in w[3:6]])))
        elif k == "edge":
            blk["edges"].setdefault(int(w[1]), []).append(np.array([float(Fraction(x)) for x in w[2:5]]))
    return blk


def _gemmi_group(name):
    import gemmi
    cands = [name]
    if ":" not in name:
        cands = [name + ":2", name + ":H", name]
    for c in cands:
        try:
            sg = gemmi.find_spacegroup_by_name(c)
        except Exception:
            sg = None
        if sg is not None:
            return sg
    raise ValueError("unknown group " + name)


def expand_cgd(blk, tol=1e-3, group_override=None):
    """Return L, P (fractional sites), labels (coordination), site->(atom,op) map,
    and the periodic contact graph as a set of (k1, k2, shift) with k1<=k2."""
    sg = _gemmi_group(group_override or blk["group"])
    ops = list(sg.operations())
    L = lattice_from_params(*blk["cell"])
    pts, coord, orig = [], [], []
    for (aid, z, xyz) in blk["atoms"]:
        for oi, op in enumerate(ops):
            p = np.array(op.apply_to_xyz(list(xyz)))
            q = p - np.floor(p + 1e-9)
            found = False
            for k, r in enumerate(pts):
                d = (q - r + 0.5) % 1.0 - 0.5
                if np.linalg.norm(L @ d) < tol:
                    found = True
                    break
            if not found:
                pts.append(q)
                coord.append(z)
                orig.append((aid, oi, p - q))
    P = np.array(pts)
    n = len(P)

    def locate(x):
        for k in range(n):
            d = x - P[k]
            s = np.round(d)
            if np.linalg.norm(L @ (d - s)) < tol:
                return k, tuple(int(v) for v in s)
        return None

    edges = set()
    nmiss = 0
    for k in range(n):
        aid, oi, t = orig[k]
        op = ops[oi]
        for e in blk["edges"].get(aid, []):
            x = np.array(op.apply_to_xyz(list(e))) - t
            r = locate(x)
            if r is None:
                nmiss += 1
                continue
            j, s = r
            if (k < j) or (k == j and s > (0, 0, 0)):
                edges.add((k, j, s))
            else:
                edges.add((j, k, tuple(-v for v in s)))
    return L, P, coord, edges, nmiss, sg


def contact_degree(n, edges):
    deg = np.zeros(n, dtype=int)
    for (i, j, s) in edges:
        deg[i] += 1
        deg[j] += 1
    return deg
