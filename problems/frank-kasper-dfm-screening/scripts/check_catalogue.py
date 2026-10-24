#!/usr/bin/env python3
"""Checks on the 84 + 37 catalogue CGD files (Zenodo 10.5281/zenodo.13691342,
folders ListDATApaper and ListDATApaperStr):
  1. every edge endpoint lands on an expanded site (no miss);
  2. contact degrees equal the declared coordination (or exactly half of it:
     duplicated edge lines);
  3. the 4-cliques of the contact graph form a periodic triangulation of the
     3-torus: each triangle in exactly two tetrahedra, V - E + F - T = 0, and
     sum |vol| = det L with the file's coordinates (embedded);
  4. FK tiling: every foam face is a pentagon or hexagon and each cell has Z-12
     hexagons;
  5. Voronoi test (informative): is the contact triangulation the Delaunay
     triangulation of the given points?
  6. spring energy of the re-embedding (informative).
Prints one line per structure and PASS/FAIL for checks 1-4; writes
data/catalogue_checks.csv.  Runtime ~5-10 min, single thread."""
import os
import sys
import glob

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import dfm  # noqa: E402
import numpy as np  # noqa: E402
from dfm import cgd, topology, embed  # noqa: E402

ROOT = os.path.join(HERE, "..")
SRC = os.path.join(ROOT, "certificates", "sources", "dutour2010", "SpaceFullerene")
rows = ["id,group,n,Z12/Z14/Z15/Z16,miss,coord_ok,doubled,manifold,euler0,embedded,fk_faces,voronoi_realised,spring_energy_per_contact"]
nfail = 0
for part, sub in ((1, "ListDATApaper"), (2, "ListDATApaperStr")):
    files = [f for f in sorted(glob.glob(os.path.join(SRC, sub, "*.cgd"))) if "total" not in f]
    expected = 84 if part == 1 else 37
    if len(files) != expected:
        sys.exit("FAIL missing catalogue inputs: expected %d files in %s, found %d" % (expected, sub, len(files)))
    for f in files:
        cid = ("" if part == 1 else "II:") + os.path.basename(f)[:-4]
        blk = cgd.parse_cgd(f)
        L, P, coord, edges, nmiss, sg = cgd.expand_cgd(blk)
        n = len(P)
        deg = cgd.contact_degree(n, edges)
        coord = np.array(coord)
        coord_ok = bool(np.all(deg == coord) or np.all(2 * deg == coord))
        doubled = bool(np.all(2 * deg == coord))
        tets = topology.tets_from_graph(n, edges)
        manifold, ntri = topology.check_manifold(n, tets)
        euler0 = (n - len(edges) + ntri - len(tets)) == 0
        vols = np.abs(topology.tet_orientations(L, P, tets))
        embedded = abs(vols.sum() / abs(np.linalg.det(L)) - 1) < 1e-9 and vols.min() > 0
        faces = topology.build_foam(L, P, tets)
        hexn = np.zeros(n, int); z = np.zeros(n, int); sizes_ok = True
        for fc in faces:
            k = len(fc["ring"])
            sizes_ok &= k in (5, 6)
            for c in (fc["i"], fc["j"]):
                z[c] += 1; hexn[c] += (k == 6)
        fk = bool(sizes_ok and np.all(hexn == z - 12) and np.all(z == deg))
        try:
            dt, info = topology.periodic_delaunay(L, P, reach=1)
            vor = set(dt) == set(tets)
        except Exception:
            vor = False
        L2, P2, e, nit = embed.reembed(L, P, edges, tets, deg)
        comp = np.bincount(deg, minlength=17)[[12, 14, 15, 16]]
        ok = nmiss == 0 and coord_ok and manifold and euler0 and embedded and fk
        nfail += (not ok)
        print("%s %-8s %-10s n=%3d comp=%s doubled=%s voronoi=%s spring/contact=%.4f" % (
            "PASS" if ok else "FAIL", cid, blk["group"], n, "/".join(map(str, comp)), doubled, vor, e / len(edges)), flush=True)
        rows.append("%s,%s,%d,%s,%d,%s,%s,%s,%s,%s,%s,%s,%.5f" % (cid, blk["group"], n, "/".join(map(str, comp)), nmiss,
                    coord_ok, doubled, manifold, euler0, embedded, fk, vor, e / len(edges)))
open(os.path.join(ROOT, "data", "catalogue_checks.csv"), "w").write("\n".join(rows) + "\n")
print("%d FAIL" % nfail)
sys.exit(1 if nfail else 0)
