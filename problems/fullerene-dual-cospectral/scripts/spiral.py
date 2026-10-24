"""spiral.py -- face spirals of fullerenes computed on the dual triangulation (Manolopoulos-Fowler spiral
algorithm), canonical spiral = lexicographically smallest 12-tuple of pentagon positions (1-based) over all
starting (vertex, neighbour, orientation) choices; the Atlas isomer number C_n:k is the rank of the canonical
spiral among all isomers of C_n (Fowler-Manolopoulos, An Atlas of Fullerenes, 1995; Bille et al. 2021 use
the same lexicographic numbering C_{n,k}).
"""


def spirals(adj):
    m = len(adj)
    pos = [{w: i for i, w in enumerate(adj[v])} for v in range(m)]
    deg = [len(a) for a in adj]
    for s in range(m):
        for t in adj[s]:
            for d in (1, -1):
                seq = [s, t]
                vis = [False] * m
                vis[s] = vis[t] = True
                nunvis = [deg[v] for v in range(m)]
                for v in (s, t):
                    for w in adj[v]:
                        nunvis[w] -= 1
                p = 0
                ok = True
                while len(seq) < m:
                    while p < len(seq) and nunvis[seq[p]] == 0:
                        p += 1
                    if p >= len(seq):
                        ok = False
                        break
                    a = seq[p]
                    last = seq[-1]
                    if last not in pos[a]:
                        ok = False
                        break
                    L = adj[a]
                    w = L[(pos[a][last] + d) % len(L)]
                    if vis[w]:
                        ok = False
                        break
                    vis[w] = True
                    seq.append(w)
                    for x in adj[w]:
                        nunvis[x] -= 1
                if ok:
                    yield tuple(k + 1 for k, v in enumerate(seq) if deg[v] == 5), (s, t, d)


def canonical_spiral(adj):
    best = None
    for sp, _ in spirals(adj):
        if best is None or sp < best:
            best = sp
    return best
