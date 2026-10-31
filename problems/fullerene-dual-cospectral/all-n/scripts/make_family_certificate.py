#!/usr/bin/env python3
"""make_family_certificate.py -- writes ../certificates/family_pairs.json: for L = 20..70 (n = 60..160) the graph6 codes of
H_L, T_b(L), T_w(L), the hexagon histograms and the common integer characteristic polynomial (low to high)."""
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from family import *

out = os.path.join(os.path.dirname(HERE), 'certificates')
os.makedirs(out, exist_ok=True)
rows = []
for L in range(20, 71):
    faces, col = colour_faces(L)
    tb, _ = colour_graph(faces, col, 0)
    tw, _ = colour_graph(faces, col, 1)
    pb, pw = charpoly(tb), charpoly(tw)
    assert pb == pw
    rows.append(dict(L=L, n=2 * L + 20, H=to_graph6(adjacency_H(L)), T_black=to_graph6(tb), T_white=to_graph6(tw),
                     hexagons_with_j_pentagons_black=hex_histogram(tb), hexagons_with_j_pentagons_white=hex_histogram(tw),
                     charpoly_low_to_high=[str(pb[i]) for i in range(pb.degree() + 1)]))
json.dump(rows, open(os.path.join(out, 'family_pairs.json'), 'w'), indent=0)
print('PASS  wrote', len(rows), 'pairs')
