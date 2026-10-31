#!/usr/bin/env python3
import sys, time
sys.path.insert(0, '..')
sys.path.insert(0, '../../../scripts')
from check_direct import *
# x12_b3 first line: 1772 b 3 a 396 m 187, gadgets (q,k,j): (5,5,2),(7,7,6),(11,11,7),(17,17,8),(29,29,18),(41,41,10), tree 32, 94
a, b = 396, 3
gad = [(5,2),(7,6),(11,7),(17,8),(29,18),(41,10)]
tree = (32, 94)
G, info = build(a, b, gad, tree)
print('n', G.number_of_nodes(), 'edges', G.number_of_edges(), flush=True)
t0 = time.time()
K0 = kf_exact(G)
print('Kf(G) time', time.time()-t0, 'Kf denominators digits', len(str(K0.q)), flush=True)
for name, v in info.items():
    Gv = G.copy(); Gv.remove_node(v)
    conn = nx.is_connected(Gv)
    t1 = time.time()
    if conn:
        k = kf_exact(Gv)
        print(name, 'good' if k == K0 else 'not good', 'Delta~', float(K0.p/K0.q - k.p/k.q) if False else '', round(time.time()-t1, 1), flush=True)
    else:
        print(name, 'cut vertex', flush=True)
