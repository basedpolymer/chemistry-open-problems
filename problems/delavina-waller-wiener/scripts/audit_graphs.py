#!/usr/bin/env python3
"""Graph-level audit of Sections 6-7 of M. Liu, Zenodo 10.5281/zenodo.22015517.

For every connected graph G produced by geng on n = 2d+1 vertices with diameter d, having exactly one
nontrivial block B (k = |B| >= 3), s = n-k >= 1, and some diametral path with its midpoint in B, check:
  L7.1  the height constraints (7.1), and 1 <= xbar <= min(s,r);
  (7.2) the exact decomposition of W(G);
  (7.3) W(C_n)-W(G) >= A/2 + (rho/2) sum t_i^2 + sum t_i q_i(q_i+eps) + sum_{i<j} t_i t_j e_ij - sum t_i c_{h_i};
  (7.6) W(C_n)-W(G) >= s(R(xbar)-cinterp(xbar)) + sum t_i q_i (q_i+eps).
Prints PASS or FAIL; exit status non-zero on FAIL.  Floating point with a tolerance of 1e-9 is used for
c_h (this is a test of the formulas, not part of a proof).
Usage: audit_graphs.py d [geng args]  Example: geng -c 9 | python audit_graphs.py 4
Environment: Python 3 with networkx. Deterministic.
"""
import math
import sys
import networkx as nx

d = int(sys.argv[1])
n = 2 * d + 1
r = d // 2
c = [1.0]
for _ in range(r + 2):
    c.append(n + c[-1] - 2 * math.sqrt(c[-1]))


def cint(x):
    if r == 1:
        return c[1]
    i = int(math.floor(x))
    if i >= r:
        return c[r]
    return c[i] + (x - i) * (c[i + 1] - c[i])


Wc = n * d * (d + 1) // 2
stats = dict(read=0, diam=0, used=0, skipped_no_midpoint=0)
bad = 0
minmargin = None


def fail(msg, G):
    global bad
    bad += 1
    if bad < 10:
        print("FAIL", msg, nx.to_graph6_bytes(G, header=False).decode().strip())


for line in sys.stdin:
    line = line.strip()
    if not line:
        continue
    stats["read"] += 1
    G = nx.from_graph6_bytes(line.encode())
    sp = dict(nx.all_pairs_shortest_path_length(G))
    D = max(max(v.values()) for v in sp.values())
    if D != d:
        continue
    stats["diam"] += 1
    blocks = [set(b) for b in nx.biconnected_components(G)]
    nontriv = [b for b in blocks if len(b) >= 3]
    if len(nontriv) != 1:
        continue
    B = nontriv[0]
    k = len(B)
    s = n - k
    if s < 1:
        continue
    BG = G.subgraph(B).copy()
    # midpoint test: does some diametral path have its midpoint in B?
    found = False
    for u in G:
        for v in G:
            if u < v and sp[u][v] == d:
                for path in nx.all_shortest_paths(G, u, v):
                    if d % 2 == 0:
                        if path[d // 2] in B:
                            found = True
                    else:
                        if path[r] in B and path[r + 1] in B:
                            found = True
                    if found:
                        break
            if found:
                break
        if found:
            break
    if not found:
        stats["skipped_no_midpoint"] += 1
        continue
    stats["used"] += 1
    # trees attached to B
    H = G.copy()
    H.remove_edges_from(BG.edges())
    comp = {}
    for cc in nx.connected_components(H):
        rootb = cc & B
        assert len(rootb) == 1
        x = next(iter(rootb))
        comp[x] = cc
    att = [x for x in B if len(comp[x]) > 1]
    t = {x: len(comp[x]) - 1 for x in att}
    h = {x: max(sp[x][w] for w in comp[x]) for x in att}
    rho = k // 2
    eps = k - 2 * rho
    a = d - rho
    assert s == 2 * a + 1 - eps
    A = 2 * a**3 + 2 * a**2 * rho + 3 * a**2 + 2 * a * rho**2 + 2 * a * rho + a + (1 - eps) * rho**2
    spB = dict(nx.all_pairs_shortest_path_length(BG))
    # Lemma 7.1
    for x in att:
        ecc = max(spB[x].values())
        if not (1 <= h[x] <= r and ecc <= min(rho, d - h[x])):
            fail("(7.1) eccentricity", G)
        for y in att:
            if x < y and spB[x][y] > min(rho, d - h[x] - h[y]):
                fail("(7.1) distance", G)
    xbar = sum(t[x] * h[x] for x in att) / s
    if not (1 - 1e-9 <= xbar <= min(s, r) + 1e-9):
        fail("xbar range", G)
    # (7.2)
    W = sum(sp[u][v] for u in G for v in G if u < v)
    tot = sum(spB[u][v] for u in B for v in B if u < v)
    for x in att:
        Tn = G.subgraph(comp[x])
        spT = dict(nx.all_pairs_shortest_path_length(Tn))
        WT = sum(spT[u][v] for u in Tn for v in Tn if u < v)
        sig = sum(spT[x].values())
        tot += WT + (n - len(comp[x])) * sig + t[x] * sum(spB[x].values())
    al = sorted(att)
    for i in range(len(al)):
        for j in range(i + 1, len(al)):
            tot += t[al[i]] * t[al[j]] * spB[al[i]][al[j]]
    if tot != W:
        fail("(7.2)", G)
    # (7.3)
    q = {x: max(h[x] - a, 0) for x in att}
    e = lambda x, y: max(h[x] + h[y] - a, 0)
    rhs = A / 2 + rho / 2 * sum(t[x] ** 2 for x in att) + sum(t[x] * q[x] * (q[x] + eps) for x in att)
    rhs += sum(t[al[i]] * t[al[j]] * e(al[i], al[j]) for i in range(len(al)) for j in range(i + 1, len(al)))
    rhs -= sum(t[x] * c[h[x]] for x in att)
    if Wc - W < rhs - 1e-9:
        fail("(7.3)", G)
    # (7.6)
    eta = lambda v: v * (rho - max(2 * v - a, 0))
    ell = (lambda x: eta(1)) if r == 1 else (lambda x: eta(1) + (eta(r) - eta(1)) * (x - 1) / (r - 1))
    R = A / (2 * s) + max(rho * s / (2 * k), rho * xbar / 2, s / 2 * max(2 * xbar - a, 0) + ell(xbar) / 2)
    rhs2 = s * (R - cint(xbar)) + sum(t[x] * q[x] * (q[x] + eps) for x in att)
    if Wc - W < rhs2 - 1e-9:
        fail("(7.6)", G)
    m = (Wc - W) - rhs2
    if minmargin is None or m < minmargin:
        minmargin = m

print("stats", stats, "smallest slack in (7.6):", minmargin)
print(("PASS" if bad == 0 else "FAIL"), f"d={d}: Lemma 7.1, (7.2), (7.3), (7.6) on {stats['used']} graphs, failures={bad}")
sys.exit(1 if bad else 0)
