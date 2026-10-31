#!/usr/bin/env python3
import sys, time
sys.path.insert(0, '.')
from x6_instance import *
def check(a, b, m, gad, L, t, verts=None):
    G = build([b]*a)
    add_W(G, ('h', 0), gad, L, t)
    n = G.number_of_nodes()
    assert n == a*(b+1)+m, (n, a*(b+1)+m)
    t0 = time.time()
    K = kx.kf(G)
    gs = []
    tests = [('u', 1, 0)] if verts is None else verts
    for v in tests:
        Kv = kx.kf(kx.minus(G, v))
        gs.append(Kv == K)
    print('a', a, 'b', b, 'n', n, 'good-test', gs, 'time', round(time.time()-t0, 2), flush=True)
    return G, K
if __name__ == '__main__':
    G, K = check(5, 2, 7, [(3, 2)], 0, 0, verts=None)
    good = kx.good_vertices(G)
    print('n=22 all good vertices:', len(good), good)
    check(38, 2, 93, [(7, 5), (11, 4), (19, 5)], 8, 17)
    check(58, 2, 130, [(9, 3), (5, 0), (13, 2), (23, 21)], 8, 23)
