"""Triangulated periodic foam mesh with exact (polyhedral) geometry.

A mesh vertex v has a stored Cartesian position X[v].  A reference to a vertex
inside a face is (v, shift) meaning position X[v] + L @ shift.  Each face of the
foam (interface between cell i at image 0 and cell j at image sj) carries a list
of triangles oriented with the normal pointing from i to j.
"""
import numpy as np


class Mesh:
    def __init__(self):
        self.X = None          # (Nv,3) float Cartesian positions
        self.vtype = None      # (Nv,) 0 foam vertex, 1 foam-edge vertex, 2 face-interior
        self.tri = None        # (T,3) int
        self.tsh = None        # (T,3,3) int shifts
        self.tface = None      # (T,) int
        self.fi = None         # (F,) int
        self.fj = None         # (F,) int
        self.fsj = None        # (F,3) int
        self.fbnd = None       # list of lists of (v, shift) : boundary loop of each face
        self.fcorner = None    # list of lists of bool : is the boundary vertex a foam vertex
        self.L = None
        self.ncell = 0
        self.ref = None        # (n,3) reference points of cells (Cartesian, image 0)
        self.level = 0

    def copy(self):
        m = Mesh()
        for k, v in self.__dict__.items():
            if isinstance(v, np.ndarray):
                setattr(m, k, v.copy())
            elif isinstance(v, list):
                setattr(m, k, [list(x) for x in v])
            else:
                setattr(m, k, v)
        return m


def _edge_key(a, ta, b, tb):
    """canonical key of the mesh edge between (a,ta) and (b,tb); returns key and the
    shift of the key's base vertex in the current frame."""
    d = (tb[0] - ta[0], tb[1] - ta[1], tb[2] - ta[2])
    if a < b or (a == b and d > (0, 0, 0)):
        return (a, b, d), ta
    md = (-d[0], -d[1], -d[2])
    return (b, a, md), tb


def build_initial_mesh(L, P, tets, faces, vertex_pos=None):
    """Level-0 mesh: foam vertices at given positions (default: circumcentres of
    the canonical tets), each face triangulated as a fan from its centroid."""
    from .topology import circumcenter
    n = len(P)
    nt = len(tets)
    X = []
    vtype = []
    for ti, key in enumerate(tets):
        if vertex_pos is not None:
            X.append(vertex_pos[ti])
        else:
            pts = np.array([L @ (P[k] + np.array(s)) for (k, s) in key])
            c, _ = circumcenter(pts)
            X.append(c)
        vtype.append(0)
    tri, tsh, tface = [], [], []
    fbnd, fcorner = [], []
    fi, fj, fsj = [], [], []
    for fidx, f in enumerate(faces):
        ring = f["ring"]
        pos = [X[ti] + L @ np.array(t) for (ti, t) in ring]
        c = np.mean(pos, axis=0)
        cid = len(X)
        X.append(c)
        vtype.append(2)
        m = len(ring)
        for k in range(m):
            a = ring[k]
            b = ring[(k + 1) % m]
            tri.append((cid, a[0], b[0]))
            tsh.append(((0, 0, 0), tuple(a[1]), tuple(b[1])))
            tface.append(fidx)
        fbnd.append([(ti, tuple(t)) for (ti, t) in ring])
        fcorner.append([True] * m)
        fi.append(f["i"]); fj.append(f["j"]); fsj.append(tuple(f["sj"]))
    M = Mesh()
    M.X = np.array(X, dtype=float)
    M.vtype = np.array(vtype, dtype=np.int64)
    M.tri = np.array(tri, dtype=np.int64)
    M.tsh = np.array(tsh, dtype=np.int64)
    M.tface = np.array(tface, dtype=np.int64)
    M.fi = np.array(fi, dtype=np.int64)
    M.fj = np.array(fj, dtype=np.int64)
    M.fsj = np.array(fsj, dtype=np.int64)
    M.fbnd = fbnd
    M.fcorner = fcorner
    M.L = np.array(L, dtype=float)
    M.ncell = n
    M.ref = np.array([L @ P[k] for k in range(n)])
    M.level = 0
    return M


def refine(M):
    """Split every triangle into 4 (edge midpoints).  Midpoints of edges lying on
    foam edges are shared between the three faces meeting there."""
    L = M.L
    X = list(M.X)
    vtype = list(M.vtype)
    mid = {}

    def midpoint(a, ta, b, tb, typ):
        key, tbase = _edge_key(a, ta, b, tb)
        if key not in mid:
            a0, b0, d = key
            pa = M.X[a0]
            pb = M.X[b0] + L @ np.array(d)
            X.append(0.5 * (pa + pb))
            vtype.append(typ)
            mid[key] = len(X) - 1
        return mid[key], tbase

    # boundary loops first (so that the midpoint type is "foam edge")
    nb, nc = [], []
    for fidx in range(len(M.fbnd)):
        loop = M.fbnd[fidx]
        corner = M.fcorner[fidx]
        new, newc = [], []
        m = len(loop)
        for k in range(m):
            a, ta = loop[k]
            b, tb = loop[(k + 1) % m]
            v, tv = midpoint(a, ta, b, tb, 1)
            new.append((a, ta)); newc.append(corner[k])
            new.append((v, tuple(tv))); newc.append(False)
        nb.append(new); nc.append(newc)
    tri, tsh, tface = [], [], []
    for t in range(len(M.tri)):
        a, b, c = M.tri[t]
        ta, tb, tc = (tuple(M.tsh[t, 0]), tuple(M.tsh[t, 1]), tuple(M.tsh[t, 2]))
        ab, tab = midpoint(a, ta, b, tb, 2)
        bc, tbc = midpoint(b, tb, c, tc, 2)
        ca, tca = midpoint(c, tc, a, ta, 2)
        f = M.tface[t]
        for (v0, s0), (v1, s1), (v2, s2) in (
            ((a, ta), (ab, tab), (ca, tca)),
            ((ab, tab), (b, tb), (bc, tbc)),
            ((ca, tca), (bc, tbc), (c, tc)),
            ((ab, tab), (bc, tbc), (ca, tca)),
        ):
            tri.append((v0, v1, v2))
            tsh.append((tuple(s0), tuple(s1), tuple(s2)))
            tface.append(f)
    N = M.copy()
    N.X = np.array(X)
    N.vtype = np.array(vtype, dtype=np.int64)
    N.tri = np.array(tri, dtype=np.int64)
    N.tsh = np.array(tsh, dtype=np.int64)
    N.tface = np.array(tface, dtype=np.int64)
    N.fbnd = nb
    N.fcorner = nc
    N.level = M.level + 1
    return N
