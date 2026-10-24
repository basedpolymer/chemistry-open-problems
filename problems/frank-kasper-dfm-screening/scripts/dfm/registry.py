"""Registry of structures.

Paper structures (Reddy et al., PNAS 2018): coordinates are crystallographic
(Dutour Sikiric et al. physical-phase files, ICSD-derived) and the periodic box
is the one stated in the captions of Figs. 7-18 (for Z: c/a = 0.9935, fitted to
the paper's flat-Voronoi values; caption says 0.993).  The relaxed DFM value
depends only on the cell topology and the box, not on the initial points.

Catalogue structures: 'cat:<ID>' with ID in E01..E11, N01..N60, Phys01..Phys27
(84 structures of Dutour Sikiric, Delgado-Friedrichs, Deza 2010) and
'cat2:<ID>' for the 21-cell structures of part II (folder ListDATApaperStr).
"""
import os
import numpy as np
from . import structures as S
from . import cgd, topology

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, "..", "..", "certificates", "sources", "dutour2010", "SpaceFullerene")


def _sg(name, sg, cell, sites):
    L = S.lattice_from_params(*cell)
    lab, P = S.expand(sg, sites)
    return L, P, lab


PAPER = {
    "BCC": lambda: S.bcc(),
    "A15": lambda: S.a15(),
    "C15": lambda: S.c15(),
    "sigma": lambda: _sg("sigma", "P 42/m n m", (1.9, 1.9, 1.0),
                         [("8j", (0.31733, 0.31733, 0.2481)), ("8i", (0.06611, 0.26065, 0)),
                          ("8i'", (0.13123, 0.53651, 0)), ("4g", (0.39875, 0.60125, 0)), ("2a", (0, 0, 0))]),
    "Z": lambda: _sg("Z", "P 6/m m m", (1, 1, 0.9935, 90, 90, 120),
                     [("3f", (0.5, 0, 0)), ("2e", (0, 0, 0.25)), ("2d", (1 / 3, 2 / 3, 0.5))]),
    "C14": lambda: _sg("C14", "P 63/m m c", (1, 1, 1.63, 90, 90, 120),
                       [("6h", (1 / 6, 1 / 3, 0.25)), ("4f", (1 / 3, 2 / 3, 0.5625)), ("2a", (0, 0, 0))]),
    "P": lambda: _sg("P", "P b n m", (1.91, 3.57, 1.0),
                     [("1", (0.0737, 0.1134, 0.25)), ("2", (0.1363, 0.2547, 0.25)), ("3", (0.3257, 0.1578, 0.25)),
                      ("4", (0.6058, 0.1819, 0.25)), ("5", (0.6650, 0.3253, 0.25)), ("6", (0.4746, 0.4536, 0.25)),
                      ("7", (0.1988, 0.4047, 0.25)), ("8", (0.8152, 0.0780, 0.25)), ("9", (0.9383, 0.3650, 0.25)),
                      ("10", (0.5202, 0.0355, 0.25)), ("11", (0.2504, 0.5375, 0.9986)), ("12", (0.3868, 0.2883, 0.0008))]),
    "delta": lambda: _sg("delta", "P 21 21 21", (1.03, 1.03, 1.0),
                         [("1", (0.4519, 0.1153, 0.5322)), ("2", (0.4424, 0.3662, 0.5971)), ("3", (0.3882, 0.0523, 0.2748)),
                          ("4", (0.1337, 0.0707, 0.2157)), ("5", (0.3768, 0.4358, 0.8567)), ("6", (0.0680, 0.1442, 0.9529)),
                          ("7", (0.2648, 0.1993, 0.7486)), ("8", (0.3136, 0.2464, 0.0740)), ("9", (0.0029, 0.1969, 0.6767)),
                          ("10", (0.1885, 0.0157, 0.4960)), ("11", (0.1031, 0.4192, 0.9133)), ("12", (0.1763, 0.4832, 0.6425)),
                          ("13", (0.0338, 0.3398, 0.1807)), ("14", (0.2289, 0.2865, 0.4098))]),
    "psigma": lambda: _sg("psigma", "P b a m", (1.95, 1.64, 1.0),
                          [("1", (0.0835, 0.1622, 0.5)), ("2", (0.1297, 0.5477, 0.5)), ("3", (0.0453, 0.3413, 0)),
                           ("4", (0.1830, 0.8416, 0.2566)), ("5", (0.2510, 0.0942, 0)), ("6", (0, 0, 0))]),
    "M": lambda: _sg("M", "P n m a", (1.89, 1.0, 3.3),
                     [("1", (0.1118, 0.0048, 0.2049)), ("2", (0.245, 0.0031, 0.455)), ("3", (0.0144, 0.25, 0.0518)),
                      ("4", (0.0242, 0.25, 0.541)), ("5", (0.0593, 0.25, 0.3506)), ("6", (0.1164, 0.25, 0.7068)),
                      ("7", (0.2996, 0.25, 0.1016)), ("8", (0.3168, 0.25, 0.5778)), ("9", (0.3255, 0.25, 0.8303)),
                      ("10", (0.3388, 0.25, 0.2938)), ("11", (0.5714, 0.25, 0.6225))]),
    "mu": lambda: _sg("mu", "R -3 m:H", (1, 1, 5.2, 90, 90, 120),
                      [("18h", (0.5, 0.5, 0.09)), ("6c", (0, 0, 0.052)), ("6c'", (0, 0, 0.154)),
                       ("6c''", (0, 0, 1 / 3)), ("3b", (0, 0, 0.5))]),
    "H": lambda: _sg("H", "C m m m", (1.0, 3.88, 1.0),
                     [("2a", (0, 0, 0)), ("4f", (0.25, 0.25, 0.5)), ("4i15", (0, 0.2892, 0)), ("4i14", (0, 0.4206, 0)),
                      ("4j", (0, 0.3677, 0.5)), ("4h", (0.2, 0, 0.5)), ("8p", (0, 0.125, 0.25))]),
}


def catalogue_ids(part=1):
    import glob
    d = os.path.join(SRC, "ListDATApaper" if part == 1 else "ListDATApaperStr")
    ids = sorted(os.path.basename(f)[:-4] for f in glob.glob(d + "/*.cgd") if "total" not in f)
    return ids


def get(name):
    """Return dict with L, P, coord/labels, tets, init ('circ' or 'centroid'), info."""
    if name in PAPER:
        L, P, lab = PAPER[name]()
        tets, info = topology.periodic_delaunay(L, P, reach=1)
        return {"name": name, "L": L, "P": P, "labels": lab, "tets": tets, "init": "circ", "info": info}
    if name.split(":")[0] in ("cat", "cat2", "catR", "cat2R"):
        pref = name.split(":")[0]
        part = 1 if pref in ("cat", "catR") else 2
        cid = name.split(":", 1)[1]
        d = os.path.join(SRC, "ListDATApaper" if part == 1 else "ListDATApaperStr")
        blk = cgd.parse_cgd(os.path.join(d, cid + ".cgd"))
        L, P, coord, edges, nmiss, sg = cgd.expand_cgd(blk)
        tets = topology.tets_from_graph(len(P), edges)
        ok, ntri = topology.check_manifold(len(P), tets)
        vols = topology.tet_orientations(L, P, tets)
        embedded = abs(np.abs(vols).sum() / abs(np.linalg.det(L)) - 1.0) < 1e-9
        # is the contact-graph triangulation the Delaunay triangulation of P?
        try:
            dt, dinfo = topology.periodic_delaunay(L, P, reach=1)
            delaunay = set(dt) == set(tets)
        except Exception:
            delaunay = False
        info = {"group": blk["group"], "nmiss": nmiss, "manifold": ok, "embedded": bool(embedded),
                "delaunay": bool(delaunay), "cell_cgd": blk["cell"]}
        deg = cgd.contact_degree(len(P), edges)
        if pref.endswith("R"):
            # spring re-embedding of the contact graph (see embed.py); tetrahedra stay positive
            from . import embed
            L, P, eS, nitS = embed.reembed(L, P, edges, tets, deg)
            info["reembedded"] = True
            info["spring_energy"] = float(eS)
            try:
                dt, dinfo = topology.periodic_delaunay(L, P, reach=1)
                delaunay = set(dt) == set(tets)
            except Exception:
                delaunay = False
            info["delaunay"] = bool(delaunay)
        info["coord_declared_doubled"] = bool(np.all(np.array(coord) == 2 * deg))
        return {"name": name, "L": L, "P": P, "labels": deg.tolist(), "tets": tets,
                "init": "circ" if delaunay else "centroid", "info": info}
    raise KeyError(name)
