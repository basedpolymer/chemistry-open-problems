#!/usr/bin/env python3
"""For the 16 non-Voronoi catalogue entries: 𝓕 and embeddedness (cone test) of the
starting partitions (barycentric dual of the contact triangulation, flat faces) for
the catalogue coordinates and for the spring re-embedding.  These are valid upper
bounds for the topology (in the respective boxes) even when the relaxation
degenerates.  Writes data/frustrated_start_values.csv.  Runtime ~1 min."""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import dfm  # noqa
from dfm import registry, energy, relax
from run_dfm import initial_mesh
ROOT = os.path.join(HERE, "..")
ids = [l.strip() for l in open(os.path.join(ROOT, "data", "nonvoronoi.txt")) if l.strip()]
rows = ["id,start,F_start,cone_min_start"]
for nm in ids:
    for start, name in (("catalogue", nm), ("reembedded", nm.replace("cat:", "catR:").replace("cat2:", "cat2R:"))):
        st = registry.get(name)
        M = initial_mesh(st)
        d = energy.dfm_values(M)
        c = relax.cone_check(M)
        rows.append("%s,%s,%.8f,%.3e" % (nm, start, d["F"], c))
        print(rows[-1], flush=True)
open(os.path.join(ROOT, "data", "frustrated_start_values.csv"), "w").write("\n".join(rows) + "\n")
