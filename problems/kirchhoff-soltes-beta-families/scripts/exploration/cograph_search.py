"""FLOATING-POINT EXPLORATION with exact Fraction re-check of hits: cographs given by cotrees with <= L twin classes of sizes <= N; prints those with at least THR good vertices. Usage: python3 cograph_search.py 4 14 0.5001. Runtime about 3 min.
Environment: Fedora 44 (WSL), Python 3.12, numpy, numba, networkx, nauty geng."""
# Search cographs (cotrees with <= L leaves, twin-class sizes <= N) for many good vertices.
# Floating-point screen, exact Fraction re-check of hits. Exploration.
import sys, itertools
from fractions import Fraction
def spec(node, sz, exact=False):
    # node: ('K',i) / ('I',i) leaf with size sz[i]; ('J',children) join; ('U',children) union
    t = node[0]
    if t == 'K':
        s = sz[node[1]]
        if s == 0: return 0, {}
        d = {0: 1}
        if s > 1: d[s] = d.get(s, 0) + s-1
        return s, d
    if t == 'I':
        s = sz[node[1]]
        if s == 0: return 0, {}
        return s, {0: s}
    parts = [spec(c, sz, exact) for c in node[1]]
    parts = [p for p in parts if p[0] > 0]
    if not parts: return 0, {}
    if t == 'U':
        n = 0; d = {}
        for m, e in parts:
            n += m
            for k, v in e.items(): d[k] = d.get(k, 0)+v
        return n, d
    # join: fold
    n, d = parts[0]
    for m, e in parts[1:]:
        nd = {0: 1}
        nd[n+m] = nd.get(n+m, 0)+1
        for k, v in d.items():
            vv = v-1 if k == 0 else v
            if vv: nd[k+m] = nd.get(k+m, 0)+vv
        for k, v in e.items():
            vv = v-1 if k == 0 else v
            if vv: nd[k+n] = nd.get(k+n, 0)+vv
        n, d = n+m, nd
    return n, d
def kf(node, sz, exact=False):
    n, d = spec(node, sz)
    if n == 0 or d.get(0, 0) != 1: return None
    if exact:
        return n*sum(Fraction(v, k) for k, v in d.items() if k)
    return n*sum(v/k for k, v in d.items() if k)
def leaves(node):
    if node[0] in 'KI': return [node[1]]
    out = []
    for c in node[1]: out += leaves(c)
    return out
def good_classes(node, sz, exact=False):
    K = kf(node, sz, exact)
    if K is None: return None
    res = []
    for i in leaves(node):
        if sz[i] == 0: continue
        s2 = list(sz); s2[i] -= 1
        K2 = kf(node, s2, exact)
        if K2 is None: res.append(False); continue
        res.append(K2 == K if exact else abs(K2-K) < 1e-9*(1+K))
    return K, res
# shapes: generate canonical-ish cotrees with up to L leaves, root join
def shapes(L):
    # trees with labelled leaf slots; node type alternates; leaves typed K or I
    from functools import lru_cache
    out = []
    def gen(nleaves, typ):
        # returns list of nodes (with placeholder leaf indices 0..nleaves-1 to be renumbered)
        res = []
        if nleaves == 1:
            return [('K', None), ('I', None)]
        # split nleaves into >=2 children (multiset, nonincreasing sizes)
        def parts(n, maxp):
            if n == 0: yield []; return
            for p in range(min(n, maxp), 0, -1):
                for rest in parts(n-p, p): yield [p]+rest
        child_typ = 'U' if typ == 'J' else 'J'
        for pr in parts(nleaves, nleaves-1):
            if len(pr) < 2: continue
            opts = [gen(p, child_typ) for p in pr]
            for combo in itertools.product(*opts):
                res.append((typ, list(combo)))
        return res
    for n in range(2, L+1):
        out += gen(n, 'J')
    # renumber leaves
    final = []
    for s in out:
        cnt = [0]
        def ren(nd):
            if nd[0] in 'KI':
                i = cnt[0]; cnt[0] += 1; return (nd[0], i)
            return (nd[0], [ren(c) for c in nd[1]])
        final.append(ren(s))
    return final
def size_of(node, sz):
    return sum(sz[i] for i in leaves(node))
if __name__ == '__main__':
    L = int(sys.argv[1]); N = int(sys.argv[2]); thr = float(sys.argv[3])
    shp = shapes(L)
    print('shapes', len(shp), file=sys.stderr)
    seen = set()
    for s in shp:
        nl = len(leaves(s))
        for sz in itertools.product(range(1, N+1), repeat=nl):
            r = good_classes(s, sz)
            if r is None: continue
            K, g = r
            ls = leaves(s)
            ng = sum(sz[i] for i, ok in zip(ls, g) if ok)
            n = sum(sz)
            if ng >= thr*n and ng > 0:
                re = good_classes(s, sz, exact=True)
                ng2 = sum(sz[i] for i, ok in zip(ls, re[1]) if ok)
                if ng2 >= thr*n:
                    print(ng2, n, '%.4f' % (ng2/n), s, sz, flush=True)
