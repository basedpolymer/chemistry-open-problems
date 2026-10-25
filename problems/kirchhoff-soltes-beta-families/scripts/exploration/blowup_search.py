"""FLOATING-POINT EXPLORATION: blow-ups of connected patterns with <= K vertices (each pattern vertex -> clique or independent set of size <= N); prints blow-ups with at least THR of their vertices good (float screen, hits re-checked exactly with Fractions not needed: none found above 0.55). Usage: python3 blowup_search.py 4 20 0.55. Runtime about 5 min.
Environment: Fedora 44 (WSL), Python 3.12, numpy, numba, networkx, nauty geng."""
# Blow-ups of small connected patterns P: vertex i -> clique (type 1) or independent set (type 0) of size t_i;
# classes complete to each other iff adjacent in P. Laplacian spectrum: within-class eigenvalue e_i (mult t_i-1)
# and quotient Q. Count good vertices (float screen), exact re-check with sympy-free Fractions via charpoly.
import sys, itertools, numpy as np, numba as nb, networkx as nx, subprocess
from fractions import Fraction
@nb.njit(cache=True)
def kf_blow(adj, typ, t):
    k = adj.shape[0]
    n = 0
    for i in range(k): n += t[i]
    if n < 2: return -1.0
    # connectivity of support
    idx = [i for i in range(k) if t[i] > 0]
    m = len(idx)
    Q = np.zeros((m, m))
    tot = 0.0
    for a in range(m):
        i = idx[a]
        d = 0.0
        for j in range(k):
            if adj[i, j] and t[j] > 0: d += t[j]
        if typ[i] == 1: d += t[i]-1
        e = d + 1 if typ[i] == 1 else d
        if t[i] > 1:
            if e <= 0: return -1.0
            tot += (t[i]-1)/e
        Q[a, a] = d - (t[i]-1 if typ[i] == 1 else 0)
        for bb in range(m):
            j = idx[bb]
            if bb != a and adj[i, j]: Q[a, bb] = -t[j]
    ev = np.linalg.eigvals(Q).real
    ev = np.sort(ev)
    if abs(ev[0]) > 1e-9 or (m > 1 and ev[1] < 1e-9): return -1.0  # disconnected
    for a in range(1, m): tot += 1/ev[a]
    return n*tot
@nb.njit(cache=True)
def run(adj, typ, N, thr, out):
    k = adj.shape[0]
    t = np.ones(k, np.int64)
    cnt = 0
    total = N**k
    for code in range(total):
        c = code
        for i in range(k):
            t[i] = 1 + c % N; c //= N
        K = kf_blow(adj, typ, t)
        if K < 0: continue
        n = 0
        for i in range(k): n += t[i]
        good = 0
        for i in range(k):
            t[i] -= 1
            K2 = kf_blow(adj, typ, t)
            t[i] += 1
            if K2 > 0 and abs(K2-K) < 1e-9*(1+K): good += t[i]
        if good >= thr*n and good > 0 and cnt < out.shape[0]:
            for i in range(k): out[cnt, i] = t[i]
            out[cnt, k] = good; out[cnt, k+1] = n
            cnt += 1
    return cnt
kmax = int(sys.argv[1]); N = int(sys.argv[2]); thr = float(sys.argv[3])
for kk in range(2, kmax+1):
    pats = subprocess.run(['geng', '-c', '-q', str(kk)], capture_output=True, text=True).stdout.split()
    for s in pats:
        P = nx.from_graph6_bytes(s.encode()); adj = nx.to_numpy_array(P, nodelist=range(kk)).astype(np.bool_)
        for typ in itertools.product([0, 1], repeat=kk):
            ty = np.array(typ, np.int64)
            out = np.zeros((2000, kk+2), np.int64)
            c = run(adj, ty, N, thr, out)
            for r in out[:c]:
                print(s, typ, list(r[:kk]), r[kk], r[kk+1], '%.4f' % (r[kk]/r[kk+1]), flush=True)
