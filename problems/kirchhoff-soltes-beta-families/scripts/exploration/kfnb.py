"""kfnb.py -- numba helpers (grounded-Laplacian inverse, Kf, transmissions, gadget classes) for the floating-point
exploration scripts in this directory.  Copied from gadget_scan.py.  Exploration only."""
# batched numba scan of rooted gadgets (floating point, exploration only).
# For each graph B and root x: classes of internal non-cut vertices with equal (rho, delta) (tol 1e-7).
# Output classes with g/m > thr: g, m, x, rho, delta, a=R-m*rho, rho*=R/m.
import sys, numpy as np, numba as nb
@nb.njit(cache=True)
def inv_grounded(A, keep, nk, out):
    m = nk-1
    L = np.zeros((m, m))
    for i in range(m):
        d = 0.0
        for j in range(nk):
            d += A[keep[i], keep[j]]
        L[i, i] = d
        for j in range(m):
            if i != j:
                L[i, j] = -A[keep[i], keep[j]]
    # Gauss-Jordan
    I = np.eye(m)
    for c in range(m):
        p = L[c, c]
        for j in range(m):
            L[c, j] /= p; I[c, j] /= p
        for r in range(m):
            if r != c:
                f = L[r, c]
                if f != 0.0:
                    for j in range(m):
                        L[r, j] -= f*L[c, j]; I[r, j] -= f*I[c, j]
    for i in range(m):
        for j in range(m):
            out[i, j] = I[i, j]
@nb.njit(cache=True)
def kfrt(A, keep, nk, M, rt):
    inv_grounded(A, keep, nk, M)
    m = nk-1
    tr = 0.0; s = 0.0
    for i in range(m):
        tr += M[i, i]
        for j in range(m):
            s += M[i, j]
    for u in range(m):
        acc = M[u, u]
        for w in range(m):
            if w != u:
                acc += M[u, u]+M[w, w]-2*M[u, w]
        rt[u] = acc
    rt[m] = tr
    return nk*tr - s
@nb.njit(cache=True)
def conn(A, keep, nk):
    seen = np.zeros(nk, np.bool_); st = np.zeros(nk, np.int64); top = 1; st[0] = 0; seen[0] = True; cnt = 1
    while top > 0:
        top -= 1; u = st[top]
        for j in range(nk):
            if not seen[j] and A[keep[u], keep[j]] > 0:
                seen[j] = True; st[top] = j; top += 1; cnt += 1
    return cnt == nk
@nb.njit(cache=True)
def scan(As, n, thr, tol):
    G = As.shape[0]; m = n-1
    res = np.zeros((G*n*n, 8))  # gi, x, g, rho, delta, a, rhostar, flag (at most n classes per root)
    nr = 0
    M = np.zeros((n, n)); rt = np.zeros(n); rtv = np.zeros(n)
    allk = np.arange(n)
    rho = np.zeros((n, n)); dl = np.zeros(n); ok = np.zeros(n, np.bool_)
    keep = np.zeros(n-1, np.int64)
    for gi in range(G):
        A = As[gi]
        # ground at last vertex: rt from kfrt is indexed by keep order (keep = allk)
        kf = kfrt(A, allk, n, M, rt)
        for v in range(n):
            c = 0
            for u in range(n):
                if u != v:
                    keep[c] = u; c += 1
            if not conn(A, keep, n-1):
                ok[v] = False
                continue
            ok[v] = True
            kfv = kfrt(A, keep, n-1, M, rtv)
            dl[v] = kf - kfv
            for i in range(n-1):
                rho[v, keep[i]] = rtv[i] - rt[keep[i]]
        for x in range(n):
            used = np.zeros(n, np.bool_)
            for v in range(n):
                if v == x or not ok[v] or used[v]:
                    continue
                g = 0
                for w in range(n):
                    if w != x and ok[w] and not used[w]:
                        if abs(rho[w, x]-rho[v, x]) < tol and abs(dl[w]-dl[v]) < tol:
                            used[w] = True; g += 1
                if g > thr*m and nr < res.shape[0]:
                    r = rho[v, x]; R = rt[x]
                    res[nr, 0] = gi; res[nr, 1] = x; res[nr, 2] = g; res[nr, 3] = r
                    res[nr, 4] = dl[v]; res[nr, 5] = R - m*r; res[nr, 6] = R/m
                    nr += 1
    return res[:nr]
def parse(lines, n):
    G = len(lines)
    As = np.zeros((G, n, n))
    nb_ = n*(n-1)//2
    for gi, s in enumerate(lines):
        b = np.frombuffer(s[1:].encode(), dtype=np.uint8).astype(np.int64)-63
        bits = ((b[:, None] >> np.arange(5, -1, -1)) & 1).ravel()[:nb_]
        k = 0
        for j in range(1, n):
            for i in range(j):
                if bits[k]:
                    As[gi, i, j] = As[gi, j, i] = 1
                k += 1
    return As

