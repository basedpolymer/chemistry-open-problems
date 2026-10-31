#!/usr/bin/env python3
"""exploration: anatomy of the near-cospectral T^6 pairs found as hash collisions (n = 140)."""
import collections, json, os, subprocess, sys, tempfile
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from twin import *
N = int(sys.argv[1])
logdir = os.path.expanduser('~/fdc/t6')
pairs = []
for line in open(os.path.join(logdir, 'n%d.txt' % N)):
    if line.startswith('C '):
        _, ty, i, j = line.split()
        pairs.append((int(i), int(j)))
idx = sorted({x for p in pairs for x in p})
with tempfile.NamedTemporaryFile('w', suffix='.idx', delete=False) as f:
    f.write(' '.join(map(str, idx)))
    name = f.name
p1 = subprocess.Popen(['buckygen', '-q', '%dd' % N], stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
p2 = subprocess.run([os.path.expanduser('~/fdc/bin/pcselect'), name], stdin=p1.stdout, capture_output=True)
p1.stdout.close(); p1.kill(); p1.wait()
os.unlink(name)
graphs = dict(zip(idx, read_planar_code(p2.stdout)))
json.dump({str(k): to_graph6(v) for k, v in graphs.items()}, open(os.path.join(logdir, 'n%d.near.json' % N), 'w'))


def t5(rot):
    pent = [v for v in range(len(rot)) if len(rot[v]) == 5]
    G = nx.Graph(); G.add_nodes_from(pent)
    G.add_edges_from((u, v) for u in pent for v in rot[u] if v in pent)
    return sorted((len(c), G.subgraph(c).number_of_edges()) for c in nx.connected_components(G))


for i, j in pairs:
    r1, r2 = graphs[i], graphs[j]
    S1 = [v for v in range(len(r1)) if len(r1[v]) == 6]
    S2 = [v for v in range(len(r2)) if len(r2[v]) == 6]
    q1, q2 = charpoly(r1, S1), charpoly(r2, S2)
    d = q1 - q2
    cc, fac = d.factor()
    mp, core = best_common_patch(r1, r2)
    print(i, j, 'T5', t5(r1), t5(r2), 'core', len(core), 'of', len(r1), 'degree diff', d.degree(), 'of', q1.degree(),
          'content', cc, 'factors', [(f.degree(), e) for f, e in fac], 'charpoly T equal', charpoly(r1) == charpoly(r2), flush=True)
