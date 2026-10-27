"""v2lib: independent exact helpers for verifier V2 (Claude Opus 5.5, effort xhigh, 27 Oct 2026).

Written from scratch for the verification; it does not import the authors' code (eclib.py, petri.py, e2lib.py).
Exact arithmetic with fractions.Fraction throughout.  A network is a list of reactions (src, prod), each a tuple of
n nonnegative integers.
"""
from fractions import Fraction as F
from itertools import combinations, product
import random


# ---------------------------------------------------------------------------------------------------- linear algebra
def rref(M):
    A = [[F(x) for x in row] for row in M]
    rows = len(A)
    cols = len(A[0]) if rows else 0
    piv = []
    r = 0
    for c in range(cols):
        p = next((i for i in range(r, rows) if A[i][c] != 0), None)
        if p is None:
            continue
        A[r], A[p] = A[p], A[r]
        pv = A[r][c]
        A[r] = [x / pv for x in A[r]]
        for i in range(rows):
            if i != r and A[i][c] != 0:
                f = A[i][c]
                A[i] = [a - f * b for a, b in zip(A[i], A[r])]
        piv.append(c)
        r += 1
        if r == rows:
            break
    return A, piv


def rank(M):
    if not M or not M[0]:
        return 0
    return len(rref(M)[1])


def nullspace(M, ncol=None):
    """Basis of {v : M v = 0}."""
    if not M:
        return [[F(int(i == j)) for i in range(ncol)] for j in range(ncol)]
    ncol = len(M[0])
    R, piv = rref(M)
    free = [c for c in range(ncol) if c not in piv]
    basis = []
    for f in free:
        v = [F(0)] * ncol
        v[f] = F(1)
        for i, c in enumerate(piv):
            v[c] = -R[i][f]
        basis.append(v)
    return basis


def transpose(M):
    return [list(r) for r in zip(*M)]


def matmul(A, B):
    Bt = transpose(B)
    return [[sum(a * b for a, b in zip(row, col)) for col in Bt] for row in A]


def det(M):
    A = [[F(x) for x in row] for row in M]
    n = len(A)
    d = F(1)
    for c in range(n):
        p = next((i for i in range(c, n) if A[i][c] != 0), None)
        if p is None:
            return F(0)
        if p != c:
            A[c], A[p] = A[p], A[c]
            d = -d
        d *= A[c][c]
        for i in range(c + 1, n):
            if A[i][c] != 0:
                f = A[i][c] / A[c][c]
                A[i] = [a - f * b for a, b in zip(A[i], A[c])]
    return d


def charpoly(M):
    """Coefficients [1, c1, ..., cn] of det(tI - M) (Faddeev-LeVerrier, exact)."""
    n = len(M)
    A = [[F(x) for x in row] for row in M]
    coeffs = [F(1)]
    Mk = [[F(0)] * n for _ in range(n)]
    I = [[F(int(i == j)) for j in range(n)] for i in range(n)]
    for k in range(1, n + 1):
        # Mk = A*M_{k-1} + c_{k-1} I ; c_k = -tr(A Mk)/k
        prev = Mk
        Mk = matmul(A, prev) if k > 1 else [row[:] for row in I]
        if k > 1:
            Mk = [[Mk[i][j] + coeffs[-1] * I[i][j] for j in range(n)] for i in range(n)]
        AM = matmul(A, Mk)
        ck = -sum(AM[i][i] for i in range(n)) / k
        coeffs.append(ck)
    return coeffs


def strip_zero_roots(coeffs, k):
    """Divide det(tI-M) by t^k; return None if t^k does not divide it."""
    if k == 0:
        return coeffs
    if any(c != 0 for c in coeffs[len(coeffs) - k:]):
        return None
    return coeffs[:len(coeffs) - k]


def is_hurwitz(coeffs):
    """coeffs = [1, a1, ..., ad] monic; True iff all roots have negative real part (Hurwitz determinants > 0)."""
    d = len(coeffs) - 1
    if d == 0:
        return True
    a = [F(c) for c in coeffs]
    if any(x <= 0 for x in a):
        return False
    H = [[F(0)] * d for _ in range(d)]
    for i in range(d):
        for j in range(d):
            idx = 2 * j - i + 1
            if 0 <= idx <= d:
                H[i][j] = a[idx]
    for k in range(1, d + 1):
        if det([row[:k] for row in H[:k]]) <= 0:
            return False
    return True


# ---------------------------------------------------------------------------------------------------- networks
def nsp(N):
    return len(N[0][0])


def Gamma(N):
    n = nsp(N)
    return [[N[t][1][i] - N[t][0][i] for t in range(len(N))] for i in range(n)]


def Ye(N):
    n = nsp(N)
    return [[N[t][0][i] for t in range(len(N))] for i in range(n)]


def Yp(N):
    n = nsp(N)
    return [[N[t][1][i] for t in range(len(N))] for i in range(n)]


def dual(N):
    E, P = Ye(N), Yp(N)
    return [(tuple(E[i]), tuple(P[i])) for i in range(nsp(N))]


def elementary_nonneg(M):
    """Extreme rays of {v >= 0 : M v = 0}: vectors whose support T has a one-dimensional kernel of M[:,T] spanned by
    a vector positive on T (minimal supports).  Brute force over supports (small networks only)."""
    ncol = len(M[0])
    rays = []
    supports = []
    for size in range(1, ncol + 1):
        for T in combinations(range(ncol), size):
            if any(set(S) <= set(T) for S in supports):
                continue
            sub = [[row[j] for j in T] for row in M]
            ns = nullspace(sub) if any(any(x != 0 for x in row) for row in sub) else \
                [[F(int(i == j)) for i in range(size)] for j in range(size)]
            if len(ns) != 1:
                continue
            v = ns[0]
            if all(x > 0 for x in v):
                pass
            elif all(x < 0 for x in v):
                v = [-x for x in v]
            else:
                continue
            full = [F(0)] * ncol
            for j, x in zip(T, v):
                full[j] = x
            rays.append(full)
            supports.append(T)
    return rays


def flux_rays(N):
    return elementary_nonneg(Gamma(N))


def cons_rays(N):
    return elementary_nonneg(transpose(Gamma(N)))


def conservative(N):
    R = cons_rays(N)
    return all(any(z[i] > 0 for z in R) for i in range(nsp(N)))


def consistent(N):
    R = flux_rays(N)
    return all(any(v[j] > 0 for v in R) for j in range(len(N)))


def rows_distinct(N):
    G = Gamma(N)
    return len(set(map(tuple, G))) == len(G)


def cols_distinct(N):
    G = transpose(Gamma(N))
    return len(set(map(tuple, G))) == len(G)


def H123(N):
    return conservative(N) and consistent(N) and rows_distinct(N)


def connected(N):
    n, r = nsp(N), len(N)
    parent = list(range(n + r))

    def f(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    for t, (a, b) in enumerate(N):
        for i in range(n):
            if a[i] or b[i]:
                parent[f(i)] = f(n + t)
    used = [i for i in range(n) if any(N[t][0][i] or N[t][1][i] for t in range(r))] + [n + t for t in range(r)]
    return len({f(x) for x in used}) == 1


def clusters(N):
    """Petri clusters: components of the bipartite species-reaction graph with reactant edges.
    Returns list of (species set, reaction set)."""
    n, r = nsp(N), len(N)
    parent = list(range(n + r))

    def f(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    for t, (a, b) in enumerate(N):
        for i in range(n):
            if a[i]:
                parent[f(i)] = f(n + t)
    comps = {}
    for x in range(n + r):
        comps.setdefault(f(x), []).append(x)
    out = []
    for c in comps.values():
        S = frozenset(x for x in c if x < n)
        T = frozenset(x - n for x in c if x >= n)
        if T:  # species never consumed form singleton components without reactions: ignore here
            out.append((S, T))
    return sorted(out, key=lambda st: min(st[1]))


def is_EC(N):
    for t, u in combinations(range(len(N)), 2):
        a, b = N[t][0], N[u][0]
        if any(x and y for x, y in zip(a, b)) and a != b:
            return False
    return True


def proportional(a, b):
    return rank([list(a), list(b)]) == 1 and all((x > 0) == (y > 0) for x, y in zip(a, b))


def is_PC(N):
    for t, u in combinations(range(len(N)), 2):
        a, b = N[t][0], N[u][0]
        if any(x and y for x, y in zip(a, b)) and not proportional(a, b):
            return False
    return True


def is_balanced(N):
    return all(len({x for x in a if x}) <= 1 for a, _ in N)


def is_siphon(N, S):
    S = set(S)
    for a, b in N:
        if any(b[i] for i in S) and not any(a[i] for i in S):
            return False
    return True


def minimal_siphons(N):
    n = nsp(N)
    sy = []
    for size in range(1, n + 1):
        for S in combinations(range(n), size):
            if any(set(T) <= set(S) for T in sy):
                continue
            if is_siphon(N, S):
                sy.append(S)
    return sy


def critical_siphons(N):
    """Minimal siphons containing no support of an extreme conservation law (enough: every semipositive law is a
    nonnegative combination of extreme ones, each with smaller support)."""
    R = cons_rays(N)
    sup = [frozenset(i for i, x in enumerate(z) if x > 0) for z in R]
    return [S for S in minimal_siphons(N) if not any(s <= set(S) for s in sup)]


def is_WFEC(N):
    if not is_EC(N) or not connected(N):
        return False
    if not (conservative(N) and consistent(N)):
        return False
    m = len({a for a, _ in N})
    return rank(Gamma(N)) == m - 1


def WF_cluster(N):
    """connected, (H1), (H2), s = #clusters - 1."""
    return connected(N) and conservative(N) and consistent(N) and rank(Gamma(N)) == len(clusters(N)) - 1


# ---------------------------------------------------------------------------------------------------- EC dynamics
def ec_structure(N):
    srcs = []
    for a, _ in N:
        if a not in srcs:
            srcs.append(a)
    ct = [srcs.index(a) for a, _ in N]
    return srcs, ct


def A_matrix(N, k):
    """n x m matrix with columns sum_{t in T_c} k_t (y'_t - y_c)."""
    srcs, ct = ec_structure(N)
    n, m = nsp(N), len(srcs)
    A = [[F(0)] * m for _ in range(n)]
    for t, (a, b) in enumerate(N):
        for i in range(n):
            A[i][ct[t]] += F(k[t]) * (b[i] - a[i])
    return A


def u_of_k(N, k):
    A = A_matrix(N, k)
    ns = nullspace(A)
    if len(ns) != 1:
        return None
    u = ns[0]
    if all(x < 0 for x in u):
        u = [-x for x in u]
    if not all(x > 0 for x in u):
        return None
    return u


def monomial(x, y):
    p = F(1)
    for xi, yi in zip(x, y):
        p *= F(xi) ** yi
    return p


def jacobian(N, k, x):
    """Exact Jacobian of mass action f(x) = Gamma diag(k) x^Ye at x."""
    n = nsp(N)
    J = [[F(0)] * n for _ in range(n)]
    for t, (a, b) in enumerate(N):
        w = F(k[t]) * monomial(x, a)
        for p in range(n):
            if a[p]:
                d = w * a[p] / F(x[p])
                for i in range(n):
                    J[i][p] += (b[i] - a[i]) * d
    return J


def field(N, k, x):
    n = nsp(N)
    f = [F(0)] * n
    for t, (a, b) in enumerate(N):
        w = F(k[t]) * monomial(x, a)
        for i in range(n):
            f[i] += (b[i] - a[i]) * w
    return f


def hurwitz_on_S(N, k, x):
    """True iff det(tI - J) = t^{n-s} q(t) with q Hurwitz."""
    n = nsp(N)
    s = rank(Gamma(N))
    cp = charpoly(jacobian(N, k, x))
    q = strip_zero_roots(cp, n - s)
    return q is not None and is_hurwitz(q)


def G_matrix(N, k, x):
    """G(x) = Dphi(x) A_k (m x m), phi_c = x^{y_c}."""
    srcs, ct = ec_structure(N)
    A = A_matrix(N, k)
    n, m = nsp(N), len(srcs)
    D = [[F(0)] * n for _ in range(m)]
    for c, y in enumerate(srcs):
        ph = monomial(x, y)
        for p in range(n):
            if y[p]:
                D[c][p] = y[p] * ph / F(x[p])
    return matmul(D, A)


def rand_pos(rng, size, lo=1, hi=9):
    return [F(rng.randint(lo, hi), rng.randint(lo, hi)) for _ in range(size)]


# ---------------------------------------------------------------------------------------------------- sign test
def maxminor_sign_condition(B1, B2):
    """B1, B2: s x n matrices of full row rank.  Returns True iff the products det(B1_I) det(B2_I) over all s-subsets I
    are all >= 0 or all <= 0 and not all 0 (Mueller et al. 2016: iff sigma(rowspace B1) cap sigma(rowspace(B2)^perp)
    = {0})."""
    s = len(B1)
    n = len(B1[0])
    signs = set()
    for I in combinations(range(n), s):
        d1 = det([[row[i] for i in I] for row in B1])
        d2 = det([[row[i] for i in I] for row in B2])
        p = d1 * d2
        if p > 0:
            signs.add(1)
        elif p < 0:
            signs.add(-1)
    return len(signs) == 1


def row_basis(M):
    R, piv = rref(M)
    return [row for row in R[:len(piv)]]


def parse(s, species):
    """'2A+B->3A, ...' -> network (list of (src, prod) tuples) over the species list."""
    N = []
    for rx in s.split(","):
        l, r = rx.split("->")

        def cpx(t):
            v = [0] * len(species)
            t = t.strip()
            if t == "0":
                return tuple(v)
            for term in t.split("+"):
                term = term.strip()
                i = 0
                while i < len(term) and term[i].isdigit():
                    i += 1
                c = int(term[:i]) if i else 1
                v[species.index(term[i:])] += c
            return tuple(v)
        N.append((cpx(l), cpx(r)))
    return N


def gen_wfec(rng, n, m, maxcoef=3, maxprod=3, maxper=3):
    """Random WF-EC network (same generator as gs_sim_v2.gen_wfec, with an explicit rng): m blocks, sources with
    coefficients in 1..maxcoef, products y' with Z y' = Z y_c for a random Z with n-m+1 rows (so rank <= m-1)."""
    from itertools import product as _prod
    while True:
        lab = [rng.randrange(m) for _ in range(n)]
        if len(set(lab)) == m:
            break
    B = [[p for p in range(n) if lab[p] == c] for c in range(m)]
    srcs = []
    for blk in B:
        y = [0] * n
        for p in blk:
            y[p] = rng.choice([1, 1] + list(range(2, maxcoef + 1)))
        srcs.append(tuple(y))
    rows = n - m + 1
    while True:
        Z = [[rng.choice([0, 0, 1, 1, 2]) for _ in range(n)] for _ in range(rows)]
        if all(any(Z[i][p] for i in range(rows)) for p in range(n)) and rank(Z) == rows:
            break
    cand = list(_prod(range(maxprod + 1), repeat=n))

    def zv(y):
        return tuple(sum(Z[i][p] * y[p] for p in range(n)) for i in range(rows))
    N = []
    for y in srcs:
        prods = [c for c in cand if c != y and zv(c) == zv(y)]
        if not prods:
            return None
        for _ in range(rng.randint(1, maxper)):
            N.append((y, rng.choice(prods)))
    N = sorted(set(N))
    return N if is_WFEC(N) else None
