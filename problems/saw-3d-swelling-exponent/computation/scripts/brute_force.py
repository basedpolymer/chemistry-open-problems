"""brute_force.py -- independent brute-force enumerator of SAWs on Z^3 (pure Python, no symmetry reduction).

For n <= N it returns the same quantities as src/saw_enum.cpp: c_n, h_n (strict half-space walks), b_n (bridges,
0 < x1(w_i) <= x1(w_n) for 1 <= i <= n), the orbit-class endpoint distributions (A, H, B) and the intermediate-point
table M[(k,n)].  Used by check_enum.py.  Runtime: N = 8 about 10 s.
"""
from collections import defaultdict

STEPS = [(1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1)]


def orbit_size(a, b, c):
    perms = 1 if a == b == c else (3 if (a == b or b == c) else 6)
    return perms << ((a > 0) + (b > 0) + (c > 0))


def cls3(p):
    return tuple(sorted((abs(v) for v in p), reverse=True))


def cls1(p):
    return (p[0], max(abs(p[1]), abs(p[2])), min(abs(p[1]), abs(p[2])))


def stab_orbit(b, c):
    return (1 if b == c else 2) * (2 if b > 0 else 1) * (2 if c > 0 else 1)


def enumerate_all(N):
    c = defaultdict(int)
    h = defaultdict(int)
    br = defaultdict(int)
    A = defaultdict(lambda: defaultdict(int))
    H = defaultdict(lambda: defaultdict(int))
    B = defaultdict(lambda: defaultdict(int))
    M = defaultdict(lambda: defaultdict(int))
    walk = [(0, 0, 0)]
    seen = {(0, 0, 0)}

    def visit():
        m = len(walk) - 1
        e = walk[-1]
        c[m] += 1
        A[m][cls3(e)] += 1
        for k in range(1, m):
            M[(k, m)][cls3(walk[k])] += 1
        xs = [p[0] for p in walk[1:]]
        if all(x > 0 for x in xs):
            h[m] += 1
            H[m][cls1(e)] += 1
            if all(x <= e[0] for x in xs):
                br[m] += 1
                B[m][cls1(e)] += 1
        if m == N:
            return
        for s in STEPS:
            q = (e[0] + s[0], e[1] + s[1], e[2] + s[2])
            if q in seen:
                continue
            seen.add(q)
            walk.append(q)
            visit()
            walk.pop()
            seen.discard(q)

    visit()
    out = {"c": dict(c), "h": dict(h), "b": dict(br), "A": {}, "H": {}, "B": {}, "M": {}}
    for m in A:
        out["A"][m] = {k: (v, orbit_size(*k)) for k, v in A[m].items()}
    for m in H:
        out["H"][m] = {k: (v, stab_orbit(k[1], k[2])) for k, v in H[m].items()}
    for m in B:
        out["B"][m] = {k: (v, stab_orbit(k[1], k[2])) for k, v in B[m].items()}
    for km in M:
        out["M"][km] = {k: (v, orbit_size(*k)) for k, v in M[km].items()}
    return out


def gram_all(N):
    """G[(n,i,j)] = sum over n-step SAWs of w_i . w_j for 1 <= i <= j <= n (brute force, no symmetry)."""
    G = defaultdict(int)
    walk = [(0, 0, 0)]
    seen = {(0, 0, 0)}

    def visit():
        m = len(walk) - 1
        for i in range(1, m + 1):
            for j in range(i, m + 1):
                a, b = walk[i], walk[j]
                G[(m, i, j)] += a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
        if m == N:
            return
        e = walk[-1]
        for s in STEPS:
            q = (e[0] + s[0], e[1] + s[1], e[2] + s[2])
            if q in seen:
                continue
            seen.add(q)
            walk.append(q)
            visit()
            walk.pop()
            seen.discard(q)

    visit()
    return dict(G)


if __name__ == "__main__":
    import sys

    N = int(sys.argv[1]) if len(sys.argv) > 1 else 6
    r = enumerate_all(N)
    for n in range(N + 1):
        print(n, r["c"][n], r["h"].get(n, 0), r["b"].get(n, 0))
