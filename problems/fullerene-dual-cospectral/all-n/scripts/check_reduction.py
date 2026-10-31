#!/usr/bin/env python3
"""check_reduction.py -- two small exact facts used in PROGRESS.md (PASS/FAIL).

  R1  literal 1(a): an isomer with n = 34 has T^6 = C6 + K1, cospectral with the subdivided claw and not isomorphic to it;
  R2  T^5(T_b(L)) and T^5(T_w(L)) are isomorphic for L = 20..LMAX (default 60): so sigma(T^5) is not injective for any
      even n >= 60 (for L >= 25 this is a local statement, the pentagons of the two caps being far apart);
  R3  tr A^5 of T^6(T_b(L)) and of T^6(T_w(L)) differ (by the constant -20 for L >= 25, a local count), so the
      hexagon graphs of the family are never cospectral.
  R4  A(T_b^6) - D_b and A(T_w^6) - D_w are cospectral, D_c(f) = number of pentagon pairs of which the hexagon f is a common
      neighbour (identity B_hh B_hh^T = 3I - D_b + A(T_b^6), B_hh^T B_hh = 3I - D_w + A(T_w^6), B_hh = black hexagons x white
      hexagons, edge adjacency in H); checked for L = 20..LMAX on the integer matrices and the integer polynomials.
Environment as in check_family.py.  Runtime < 1 min.
"""
import collections
import os
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from family import *

FAILED = []


def report(name, ok, detail=''):
    print('%s  %s%s' % ('PASS' if ok else 'FAIL', name, ('  ' + detail) if detail else ''), flush=True)
    if not ok:
        FAILED.append(name)


def induced(adj, S):
    idx = {v: i for i, v in enumerate(S)}
    out = [[] for _ in S]
    for v in S:
        for w in adj[v]:
            if w in idx:
                out[idx[v]].append(idx[w])
    return out


def main():
    lmax = int(sys.argv[1]) if len(sys.argv) > 1 else 60
    p = subprocess.run(['buckygen', '-q', '34d'], capture_output=True, check=True)
    hit = 0
    target = labelg_canon([[[1, 5], [0, 2], [1, 3], [2, 4], [3, 5], [4, 0], []]])[0]
    spider = [[1, 3, 5], [0, 2], [1], [0, 4], [3], [0, 6], [5]]  # subdivided claw
    ok = True
    for rot in read_planar_code(p.stdout):
        S = [v for v in range(len(rot)) if len(rot[v]) == 6]
        g = induced(rot, S)
        if len(g) == 7 and labelg_canon([g])[0] == target:
            hit += 1
            ok &= charpoly(g) == charpoly(spider) and labelg_canon([spider])[0] != target
    report('R1 n = 34: %d isomer(s) with T^6 = C6 + K1, cospectral with the (connected) subdivided claw' % hit, hit >= 1 and ok)
    ok = True
    for L in range(20, lmax + 1):
        faces, col = colour_faces(L)
        tb, _ = colour_graph(faces, col, 0)
        tw, _ = colour_graph(faces, col, 1)
        pb = [v for v in range(len(tb)) if len(tb[v]) == 5]
        pw = [v for v in range(len(tw)) if len(tw[v]) == 5]
        c = labelg_canon([induced(tb, pb), induced(tw, pw)])
        ok &= c[0] == c[1]
    report('R2 T^5(T_b(L)) and T^5(T_w(L)) isomorphic, L = 20..%d' % lmax, ok)
    ok = True
    vals = set()
    for L in range(20, lmax + 1):
        faces, col = colour_faces(L)
        tb, _ = colour_graph(faces, col, 0)
        tw, _ = colour_graph(faces, col, 1)
        tr = []
        for t in (tb, tw):
            S = [v for v in range(len(t)) if len(t[v]) == 6]
            g = induced(t, S)
            k = len(g)
            M = flint.fmpz_mat(k, k, [1 if j in g[i] else 0 for i in range(k) for j in range(k)])
            P = M * M * M * M * M
            tr.append(sum(int(P[i, i]) for i in range(k)))
        ok &= tr[0] != tr[1]
        if L >= 25:
            vals.add(tr[0] - tr[1])
    report('R3 tr A^5 of the hexagon graphs differ, L = 20..%d; difference for L >= 25: %s' % (lmax, sorted(vals)), ok and len(vals) == 1)
    ok = True
    for L in range(20, lmax + 1):
        faces, col = colour_faces(L)
        V = L + 14
        h = adjacency_H(L)
        four = {v for v in range(V) if len(h[v]) == 4}
        black = [f for f, c in zip(faces, col) if c == 0 and not set(f) & four]
        white = [f for f, c in zip(faces, col) if c == 1 and not set(f) & four]
        pwhite = [f for f, c in zip(faces, col) if c == 1 and set(f) & four]
        pblack = [f for f, c in zip(faces, col) if c == 0 and set(f) & four]
        k = len(black)
        ok &= k == len(white)
        Bhh = flint.fmpz_mat(k, k, [1 if len(set(f) & set(g)) == 2 else 0 for f in black for g in white])

        def potential(rows, others):
            return [sum(1 for g in others if len(set(f) & set(g)) == 2) for f in rows]
        Db, Dw = potential(black, pwhite), potential(white, pblack)

        def mat(rows, D):
            return flint.fmpz_mat(k, k, [(3 - D[i] if i == j else (1 if i != j and set(rows[i]) & set(rows[j]) else 0))
                                         for i in range(k) for j in range(k)])
        Mb, Mw = mat(black, Db), mat(white, Dw)
        ok &= (Bhh * Bhh.transpose()) == Mb and (Bhh.transpose() * Bhh) == Mw and Mb.charpoly() == Mw.charpoly()
    report('R4 A(T_b^6) - D_b and A(T_w^6) - D_w cospectral through B_hh, L = 20..%d' % lmax, ok)
    if FAILED:
        sys.exit(1)


if __name__ == '__main__':
    main()
