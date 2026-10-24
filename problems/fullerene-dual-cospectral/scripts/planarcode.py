"""planarcode.py -- minimal reader/writer for planar_code (Brinkmann-McKay) and graph6, used by all scripts.

planar_code: header '>>planar_code<<' (optionally with 'le'/'be'), then for each graph:
  nv (1 byte; 0 followed by 2-byte little-endian count for nv > 255, which never happens here),
  then for each vertex v = 1..nv the neighbours of v in clockwise order (1-based), terminated by 0.
Buckygen (without -d) writes the triangulations = duals T_n of the fullerenes C_n (n = 2(m-2)).
"""
import sys

HEADER = b'>>planar_code<<'


def read_planar_code(data):
    """Yield adjacency lists (list of lists, 0-based, clockwise rotation kept) from planar_code bytes."""
    pos = 0
    if data.startswith(HEADER):
        pos = len(HEADER)
    elif data.startswith(b'>>planar_code le<<') or data.startswith(b'>>planar_code be<<'):
        raise ValueError('2-byte planar code not supported')
    L = len(data)
    while pos < L:
        nv = data[pos]
        pos += 1
        if nv == 0:
            raise ValueError('nv > 255 not supported')
        adj = []
        for _ in range(nv):
            nb = []
            while True:
                x = data[pos]
                pos += 1
                if x == 0:
                    break
                nb.append(x - 1)
            adj.append(nb)
        yield adj


def write_planar_code(adjs):
    out = bytearray(HEADER)
    for adj in adjs:
        out.append(len(adj))
        for nb in adj:
            out.extend(x + 1 for x in nb)
            out.append(0)
    return bytes(out)


def to_graph6(adj):
    n = len(adj)
    assert n <= 62
    bits = []
    for j in range(1, n):
        for i in range(j):
            bits.append(1 if j in adj[i] else 0)
    while len(bits) % 6:
        bits.append(0)
    s = chr(63 + n)
    for k in range(0, len(bits), 6):
        v = 0
        for b in bits[k:k + 6]:
            v = 2 * v + b
        s += chr(63 + v)
    return s


def from_graph6(s):
    s = s.strip()
    n = ord(s[0]) - 63
    assert 0 <= n <= 62
    bits = []
    for ch in s[1:]:
        v = ord(ch) - 63
        for k in range(5, -1, -1):
            bits.append((v >> k) & 1)
    adj = [[] for _ in range(n)]
    idx = 0
    for j in range(1, n):
        for i in range(j):
            if bits[idx]:
                adj[i].append(j)
                adj[j].append(i)
            idx += 1
    return adj


def check_triangulation_dual_of_fullerene(adj):
    """Return list of failed checks: symmetric simple graph, degrees in {5,6}, exactly 12 of degree 5,
    planar rotation system with all faces triangles and Euler characteristic 2 (so it is a plane triangulation)."""
    errs = []
    n = len(adj)
    for v in range(n):
        if len(set(adj[v])) != len(adj[v]) or v in adj[v]:
            errs.append('not simple at %d' % v)
        for w in adj[v]:
            if v not in adj[w]:
                errs.append('asymmetric %d-%d' % (v, w))
    deg = [len(a) for a in adj]
    if any(d not in (5, 6) for d in deg):
        errs.append('degree not in {5,6}')
    if deg.count(5) != 12:
        errs.append('number of degree-5 vertices != 12')
    # faces of the rotation system: next dart after (v->w) is (w -> successor of v in rotation of w)
    # (planar_code lists neighbours clockwise; face tracing works for either orientation)
    pos = [{w: i for i, w in enumerate(adj[v])} for v in range(n)]
    seen = set()
    nfaces = 0
    for v in range(n):
        for w in adj[v]:
            if (v, w) in seen:
                continue
            length = 0
            a, b = v, w
            while (a, b) not in seen:
                seen.add((a, b))
                length += 1
                i = pos[b][a]
                c = adj[b][(i + 1) % len(adj[b])]
                a, b = b, c
            nfaces += 1
            if length != 3:
                errs.append('non-triangular face')
    E = sum(deg) // 2
    if n - E + nfaces != 2:
        errs.append('Euler characteristic != 2')
    return errs
