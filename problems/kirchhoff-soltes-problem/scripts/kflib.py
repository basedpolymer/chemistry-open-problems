"""kflib.py -- exact Kirchhoff index utilities (python-flint), shared by the scripts of this folder.

Kf(G) = sum_{u<w} r(u,w) = n * sum_{i} 1/lambda_i  (nonzero Laplacian eigenvalues).
Two independent exact methods:
  kf_inv(G):      grounded Laplacian A (delete one vertex), Kf = n tr(A^{-1}) - 1^T A^{-1} 1  (fmpq_mat)
  kf_charpoly(G): Laplacian characteristic polynomial x^n + ... + c2 x^2 + c1 x, Kf = -n c2/c1 (fmpz_mat)
Graphs are dicts {vertex: set(neighbours)} or graph6 strings.
Environment: python 3.13, python-flint 0.9.
"""
from flint import fmpq_mat, fmpz_mat, fmpq, fmpz


def g6_decode(s):
    s = s.strip().split()[0]
    b = [ord(c) - 63 for c in s]
    if b[0] == 63:
        n = (b[1] << 12) | (b[2] << 6) | b[3]
        b = b[4:]
    else:
        n = b[0]
        b = b[1:]
    bits = []
    for x in b:
        for k in range(5, -1, -1):
            bits.append((x >> k) & 1)
    adj = {i: set() for i in range(n)}
    k = 0
    for j in range(1, n):
        for i in range(j):
            if bits[k]:
                adj[i].add(j)
                adj[j].add(i)
            k += 1
    return adj


def g6_encode(adj):
    vs = sorted(adj)
    idx = {v: i for i, v in enumerate(vs)}
    n = len(vs)
    bits = []
    for j in range(1, n):
        for i in range(j):
            bits.append(1 if vs[j] in adj[vs[i]] else 0)
    while len(bits) % 6:
        bits.append(0)
    out = [chr(n + 63)] if n < 63 else ['~', chr(((n >> 12) & 63) + 63), chr(((n >> 6) & 63) + 63), chr((n & 63) + 63)]
    for k in range(0, len(bits), 6):
        x = 0
        for t in range(6):
            x = (x << 1) | bits[k + t]
        out.append(chr(x + 63))
    return ''.join(out)


def induced(adj, S):
    S = set(S)
    return {u: adj[u] & S for u in S}


def delete(adj, v):
    return induced(adj, set(adj) - {v})


def is_connected(adj):
    vs = list(adj)
    if not vs:
        return True
    seen = {vs[0]}
    stack = [vs[0]]
    while stack:
        u = stack.pop()
        for w in adj[u]:
            if w not in seen:
                seen.add(w)
                stack.append(w)
    return len(seen) == len(vs)


def laplacian_rows(adj):
    vs = sorted(adj)
    idx = {v: i for i, v in enumerate(vs)}
    n = len(vs)
    L = [[0] * n for _ in range(n)]
    for u in vs:
        i = idx[u]
        L[i][i] = len(adj[u])
        for w in adj[u]:
            L[i][idx[w]] = -1
    return L


def kf_inv(adj):
    """Exact Kf via the inverse of the grounded Laplacian (requires connected graph)."""
    n = len(adj)
    if n <= 1:
        return fmpq(0)
    L = laplacian_rows(adj)
    A = fmpq_mat([row[:-1] for row in L[:-1]])
    Ai = A.inv()
    N = n - 1
    tr = sum((Ai[i, i] for i in range(N)), fmpq(0))
    tot = sum((Ai[i, j] for i in range(N) for j in range(N)), fmpq(0))
    return n * tr - tot


def kf_charpoly(adj):
    """Exact Kf via the Laplacian characteristic polynomial (requires connected graph)."""
    n = len(adj)
    if n <= 1:
        return fmpq(0)
    L = fmpz_mat(laplacian_rows(adj))
    cp = L.charpoly()           # monic, degree n
    c = [cp[k] for k in range(n + 1)]
    assert c[0] == 0
    c1, c2 = c[1], c[2]
    return fmpq(-n * c2, c1)


def spanning_trees(adj):
    n = len(adj)
    if n <= 1:
        return fmpz(1)
    L = laplacian_rows(adj)
    return fmpz_mat([row[:-1] for row in L[:-1]]).det()


def good_vertices(adj, method=kf_inv):
    """Return (Kf(G), {v: Kf(G-v) or None if G-v disconnected})."""
    k = method(adj)
    res = {}
    for v in sorted(adj):
        H = delete(adj, v)
        res[v] = method(H) if is_connected(H) else None
    return k, res


def is_kirchhoff_soltes(adj, method=kf_inv):
    if len(adj) < 3 or not is_connected(adj):
        return False
    k, res = good_vertices(adj, method)
    return all(x is not None and x == k for x in res.values())
