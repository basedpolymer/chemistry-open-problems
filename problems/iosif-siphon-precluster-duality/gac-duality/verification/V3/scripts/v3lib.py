"""v3lib: small exact library written by verifier V3 (Claude Opus 5.5, effort xhigh), 27 October 2026.

Independent of the authors' code (crn.py, ecore.py, s1lib.py, d2lib.py are NOT imported here).
Exact arithmetic only (fractions.Fraction).  Standard library only.

A network is a pair (n, R) with R a list of reactions (src, prod), src/prod tuples of n nonnegative ints.
"""
from fractions import Fraction
from itertools import combinations, permutations
from math import gcd


# ---------------------------------------------------------------- linear algebra over Q
def rref(M):
    """Row-reduced echelon form of a list-of-lists matrix (Fractions); returns (R, pivots)."""
    A = [[Fraction(x) for x in row] for row in M]
    if not A:
        return A, []
    m, n = len(A), len(A[0])
    piv = []
    r = 0
    for c in range(n):
        p = next((i for i in range(r, m) if A[i][c] != 0), None)
        if p is None:
            continue
        A[r], A[p] = A[p], A[r]
        inv = 1 / A[r][c]
        A[r] = [x * inv for x in A[r]]
        for i in range(m):
            if i != r and A[i][c] != 0:
                f = A[i][c]
                A[i] = [a - f * b for a, b in zip(A[i], A[r])]
        piv.append(c)
        r += 1
        if r == m:
            break
    return A, piv


def rank(M):
    if not M or not M[0]:
        return 0
    return len(rref(M)[1])


def nullspace(M, ncols=None):
    """Basis of {x : M x = 0}."""
    if not M:
        return [[Fraction(int(i == j)) for i in range(ncols)] for j in range(ncols)]
    n = len(M[0])
    R, piv = rref(M)
    free = [c for c in range(n) if c not in piv]
    basis = []
    for f in free:
        v = [Fraction(0)] * n
        v[f] = Fraction(1)
        for i, p in enumerate(piv):
            v[p] = -R[i][f]
        basis.append(v)
    return basis


def transpose(M):
    return [list(r) for r in zip(*M)]


def matmul(A, B):
    return [[sum(Fraction(a) * b for a, b in zip(row, col)) for col in zip(*B)] for row in A]


# ---------------------------------------------------------------- exact LP feasibility (Phase I simplex, Bland)
def feasible(A, b):
    """Return x >= 0 (Fractions) with A x = b, or None.  A: m x n, b: length m."""
    m = len(A)
    n = len(A[0]) if m else 0
    A = [[Fraction(x) for x in row] for row in A]
    b = [Fraction(x) for x in b]
    for i in range(m):
        if b[i] < 0:
            A[i] = [-x for x in A[i]]
            b[i] = -b[i]
    # tableau with artificials n..n+m-1
    T = [A[i] + [Fraction(int(i == j)) for j in range(m)] + [b[i]] for i in range(m)]
    basis = [n + i for i in range(m)]
    N = n + m
    # objective: minimise sum of artificials -> reduced costs
    def obj_row():
        z = [Fraction(0)] * (N + 1)
        for i in range(m):
            if basis[i] >= n:
                z = [zz - t for zz, t in zip(z, T[i])]
        for j in range(n, N):
            z[j] += 1
        return z
    while True:
        z = obj_row()
        # entering: smallest index with negative reduced cost (Bland)
        e = next((j for j in range(N) if z[j] < 0 and j not in basis), None)
        if e is None:
            break
        rows = [i for i in range(m) if T[i][e] > 0]
        if not rows:
            break
        best = min(rows, key=lambda i: (T[i][-1] / T[i][e], basis[i]))
        piv = T[best][e]
        T[best] = [x / piv for x in T[best]]
        for i in range(m):
            if i != best and T[i][e] != 0:
                f = T[i][e]
                T[i] = [a - f * c for a, c in zip(T[i], T[best])]
        basis[best] = e
    x = [Fraction(0)] * N
    for i in range(m):
        x[basis[i]] = T[i][-1]
    if any(x[j] != 0 for j in range(n, N)):
        return None
    sol = x[:n]
    # verify exactly
    for i in range(m):
        assert sum(A[i][j] * sol[j] for j in range(n)) == b[i]
    return sol


def in_convex_hull(p, pts):
    """Exact test p in conv(pts)."""
    if not pts:
        return False
    d = len(p)
    A = [[Fraction(q[k]) for q in pts] for k in range(d)] + [[Fraction(1)] * len(pts)]
    return feasible(A, list(p) + [1]) is not None


def semipositive_in_kernel_rows(M, support):
    """Is there z >= 0, z != 0, supp z within `support`, with z^T M = 0 (M: rows indexed by species)?"""
    S = sorted(support)
    if not S:
        return False
    ncol = len(M[0])
    A = [[Fraction(M[i][j]) for i in S] for j in range(ncol)] + [[Fraction(1)] * len(S)]
    return feasible(A, [0] * ncol + [1]) is not None


def positive_in_kernel(M):
    """Is there v > 0 with M v = 0?  (substitute v = 1 + u, u >= 0: M u = -M 1; strictness by scaling)."""
    # v > 0 with Mv = 0  <=>  exists u >= 0 with M (1 + u) = 0 ... not exact in general (scaling); use:
    # v >= 1 componentwise (cone is homogeneous) <=> v = 1 + u, u >= 0.
    ncol = len(M[0])
    one = [Fraction(1)] * ncol
    rhs = [-sum(Fraction(M[i][j]) for j in range(ncol)) for i in range(len(M))]
    return feasible(M, rhs) is not None if M else True


# ---------------------------------------------------------------- networks
def Ye(n, R):
    return [[R[j][0][i] for j in range(len(R))] for i in range(n)]


def Yp(n, R):
    return [[R[j][1][i] for j in range(len(R))] for i in range(n)]


def Gamma(n, R):
    return [[R[j][1][i] - R[j][0][i] for j in range(len(R))] for i in range(n)]


def dual(n, R):
    """Iosif dual: r species, one reaction per species i: row i of Ye -> row i of Yp."""
    r = len(R)
    return r, [(tuple(R[j][0][i] for j in range(r)), tuple(R[j][1][i] for j in range(r))) for i in range(n)]


def complexes(R):
    seen = []
    for s, p in R:
        for c in (s, p):
            if c not in seen:
                seen.append(c)
    return seen


def scc(nodes, edges):
    """Tarjan-free simple SCC via reachability (small graphs)."""
    reach = {u: {u} for u in nodes}
    changed = True
    while changed:
        changed = False
        for u, v in edges:
            for w in list(reach):
                if u in reach[w] and v not in reach[w]:
                    reach[w].add(v)
                    changed = True
    comp = {}
    for u in nodes:
        comp[u] = frozenset(v for v in nodes if v in reach[u] and u in reach[v])
    return comp


def is_wr(R):
    C = complexes(R)
    comp = scc(C, [(s, p) for s, p in R])
    return all(p in comp[s] for s, p in R)


def linkage_classes(R):
    C = complexes(R)
    parent = {c: c for c in C}

    def f(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    for s, p in R:
        parent[f(s)] = f(p)
    cl = {}
    for c in C:
        cl.setdefault(f(c), []).append(c)
    return list(cl.values())


def counts(n, R):
    """(m, ell, s, delta)."""
    m = len(complexes(R))
    ell = len(linkage_classes(R))
    s = rank(Gamma(n, R))
    return m, ell, s, m - ell - s


def is_siphon(n, R, Z):
    Z = set(Z)
    for s, p in R:
        produces = any(p[i] > s[i] for i in Z)
        if produces and not any(s[i] > 0 for i in Z):
            return False
    return True


def minimal_siphons(n, R):
    sig = [frozenset(Z) for k in range(1, n + 1) for Z in combinations(range(n), k) if is_siphon(n, R, Z)]
    return [Z for Z in sig if not any(W < Z for W in sig)]


def conservative(n, R):
    """exists w > 0 with w^T Gamma = 0."""
    G = Gamma(n, R)
    return positive_in_kernel(transpose(G)) if R else True


def consistent(n, R):
    return positive_in_kernel(Gamma(n, R))


def H3(n, R):
    G = Gamma(n, R)
    return len(set(tuple(r) for r in G)) == n


def canonical(n, R):
    best = None
    for perm in permutations(range(n)):
        key = tuple(sorted((tuple(s[perm[i]] for i in range(n)), tuple(p[perm[i]] for i in range(n))) for s, p in R))
        if best is None or key < best:
            best = key
    return best


# ---------------------------------------------------------------- rank-one face classification (Theorem S1)
def rank_one_data(n, R):
    """For a rank-one network: primitive g (first nonzero entry of the first nonzero column made positive is NOT
    imposed; any primitive g spanning S) and c_j with Delta_j = c_j g."""
    G = Gamma(n, R)
    col = next(([G[i][j] for i in range(n)] for j in range(len(R)) if any(G[i][j] for i in range(n))), None)
    gg = 0
    for x in col:
        gg = gcd(gg, abs(x))
    g = [x // gg for x in col]
    c = []
    for j in range(len(R)):
        d = [G[i][j] for i in range(n)]
        k = next((i for i in range(n) if g[i] != 0), None)
        cj = Fraction(d[k], g[k])
        assert all(Fraction(d[i]) == cj * g[i] for i in range(n)), "not rank one"
        c.append(cj)
    return g, c


def face_types(n, R):
    """Return list of (sigma, Z, type) for every face of the rank-one conservative network (n, R).
    sigma=+1: upper end, Z subset of N- = {g_i < 0}; sigma=-1: lower end, Z subset of N+ = {g_i > 0}.
    Inward reaction j at sigma: sigma*c_j < 0."""
    g, c = rank_one_data(n, R)
    out = []
    for sigma, side in ((+1, [i for i in range(n) if g[i] < 0]), (-1, [i for i in range(n) if g[i] > 0])):
        for k in range(1, len(side) + 1):
            for Z in combinations(side, k):
                act = [j for j in range(len(R)) if c[j] != 0]
                order = {j: sum(R[j][0][i] for i in Z) for j in act}
                mu = min(order.values())
                J = [j for j in act if order[j] == mu]
                inward = [j for j in J if sigma * c[j] < 0]
                outward = [j for j in J if sigma * c[j] > 0]
                if not outward:
                    t = 'IN'
                elif not inward:
                    t = 'OUT'
                else:
                    off = [i for i in range(n) if i not in Z]
                    beta = lambda j: [R[j][0][i] for i in off]
                    ins = [beta(j) for j in inward]
                    if all(in_convex_hull(beta(j), ins) for j in outward):
                        t = 'DOM'
                    else:
                        t = 'ESC'
                out.append((sigma, Z, t))
    return out


def verdict_rank_one(n, R):
    ft = face_types(n, R)
    types = {t for _, _, t in ft}
    if types <= {'IN'}:
        return 'ALL', ft
    if 'OUT' in types or 'ESC' in types:
        return 'NONE', ft
    return 'DOM-only', ft


# ---------------------------------------------------------------- facet types in any rank (Theorem W)
def facet_types_any_rank(n, R):
    """For the dual of (n, R): every facet set Z = nonempty set of reactions of N with positively parallel
    nonzero reaction vectors; returns list of (Z, type) with types IN/OUT/ESC/DOM computed from the loads.
    Active dual reactions (species i of N): a_i != 0 where gamma_j = b_j a for j in Z."""
    G = Gamma(n, R)
    cols = [[G[i][j] for i in range(n)] for j in range(len(R))]
    r = len(R)

    def pos_parallel(u, v):
        if not any(u) or not any(v):
            return False
        k = next(i for i in range(n) if u[i] != 0)
        if v[k] == 0:
            return False
        t = Fraction(v[k], u[k])
        return t > 0 and all(Fraction(v[i]) == t * u[i] for i in range(n))
    classes = []
    for j in range(r):
        if not any(cols[j]):
            continue
        for cl in classes:
            if pos_parallel(cols[cl[0]], cols[j]):
                cl.append(j)
                break
        else:
            classes.append([j])
    out = []
    for cl in classes:
        a = cols[cl[0]]
        for k in range(1, len(cl) + 1):
            for Z in combinations(cl, k):
                act = [i for i in range(n) if a[i] != 0]
                load = {i: sum(R[j][0][i] for j in Z) for i in act}
                mu = min(load.values())
                J = [i for i in act if load[i] == mu]
                inward = [i for i in J if a[i] > 0]
                outward = [i for i in J if a[i] < 0]
                if not outward:
                    t = 'IN'
                elif not inward:
                    t = 'OUT'
                else:
                    off = [l for l in range(r) if l not in Z]
                    beta = lambda i: [R[l][0][i] for l in off]
                    ins = [beta(i) for i in inward]
                    t = 'DOM' if all(in_convex_hull(beta(o), ins) for o in outward) else 'ESC'
                out.append((Z, t))
    return out


def show(n, R, names=None):
    names = names or ['X%d' % (i + 1) for i in range(n)]

    def cx(y):
        terms = [(('%d' % y[i]) if y[i] > 1 else '') + names[i] for i in range(n) if y[i] > 0]
        return '+'.join(terms) if terms else '0'
    return ', '.join('%s -> %s' % (cx(s), cx(p)) for s, p in R)
