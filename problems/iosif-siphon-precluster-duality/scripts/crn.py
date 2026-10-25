"""Exact combinatorics of mass-action networks for the Iosif siphon/precluster duality problem.

Environment: Python 3.12 (WSL Fedora 44), standard library only (fractions, itertools).
All linear algebra is exact over the rationals (fractions.Fraction); no floating point is used.

A network is given by a number of species n and a list of reactions (source, product), where source and
product are tuples of n nonnegative integers (stoichiometric coefficients).  Reactions are indexed 0..r-1,
species 0..n-1.

Objects implemented (definitions of Iosif, arXiv:2603.08767v2, and of Iosif's thesis, DOI 10.25673/13949):
  * stoichiometric matrix Gamma = Yp - Ye (n x r);
  * extreme rays of the flux cone ker(Gamma) cap R^r_{>=0} (internal cycles, Definition 20 / Prop. 1.19 of the
    thesis) and of the conservation cone left-ker(Gamma) cap R^n_{>=0};
  * the dual network N* (Section 3 of the paper): species of N* = reactions of N, reaction i of N* is
    k^{(Ye)_{i.}} -> k^{(Yp)_{i.}};
  * preclusters (Definition 19 of the paper = 'onestep'; Algorithm 4.4 of the thesis = 'thesis'; and the
    order-free closure 'closure');
  * maximal invariant polyhedral supports (MIPS): 'ray' reading (Definition 17 with each decorated type invariant,
    = Remark 18), 'semiflow' reading (supports of extreme conservation laws), and the 'set' reading of
    Definition 17 for conservation rank q <= 2 (the set of decorated types is invariant);
  * siphons (Definition 13) and minimal siphons (Definition 14).
"""
from fractions import Fraction
from itertools import combinations
from functools import lru_cache


# ----------------------------------------------------------------------------------------------
# exact linear algebra
# ----------------------------------------------------------------------------------------------
def rref(rows):
    """Reduced row echelon form of a list of rows (Fractions). Returns (R, pivots)."""
    M = [[Fraction(x) for x in row] for row in rows]
    if not M:
        return M, []
    m, ncol = len(M), len(M[0])
    piv = []
    r = 0
    for c in range(ncol):
        p = None
        for i in range(r, m):
            if M[i][c] != 0:
                p = i
                break
        if p is None:
            continue
        M[r], M[p] = M[p], M[r]
        inv = 1 / M[r][c]
        M[r] = [x * inv for x in M[r]]
        for i in range(m):
            if i != r and M[i][c] != 0:
                f = M[i][c]
                M[i] = [a - f * b for a, b in zip(M[i], M[r])]
        piv.append(c)
        r += 1
        if r == m:
            break
    return M[:r], piv


def rank(rows):
    if not rows or not rows[0]:
        return 0
    return len(rref(rows)[1])


def nullspace(rows, ncol):
    """Basis of {v : rows . v = 0} (list of Fraction vectors of length ncol)."""
    if not rows:
        return [[Fraction(int(i == j)) for i in range(ncol)] for j in range(ncol)]
    R, piv = rref(rows)
    free = [c for c in range(ncol) if c not in piv]
    basis = []
    for f in free:
        v = [Fraction(0)] * ncol
        v[f] = Fraction(1)
        for i, p in enumerate(piv):
            v[p] = -R[i][f]
        basis.append(v)
    return basis


def primitive_int(v):
    """Scale a rational vector to a primitive integer vector (same direction)."""
    from math import gcd
    den = 1
    for x in v:
        den = den * Fraction(x).denominator // gcd(den, Fraction(x).denominator)
    w = [int(Fraction(x) * den) for x in v]
    g = 0
    for x in w:
        g = gcd(g, abs(x))
    return tuple(x // g for x in w) if g else tuple(w)


def extreme_rays_nonneg_kernel(A, ncol):
    """Extreme rays of the pointed cone {v in R^ncol_{>=0} : A v = 0}, A given as list of rows.

    Uses the minimal-support characterisation: v is extreme iff its support S is minimal, iff the
    kernel of A restricted to the columns S is one-dimensional and spanned by a vector with full support S of
    constant sign.  Returns primitive integer vectors, sorted."""
    rays = []
    supports = []
    for size in range(1, ncol + 1):
        for S in combinations(range(ncol), size):
            # skip supersets of known extreme supports (their kernel cannot be 1-dim with full support)
            if any(set(T) <= set(S) for T in supports):
                continue
            sub = [[row[j] for j in S] for row in A]
            ker = nullspace(sub, len(S)) if A else [[Fraction(1)]] if len(S) == 1 else None
            if ker is None or len(ker) != 1:
                continue
            k = ker[0]
            if all(x > 0 for x in k) or all(x < 0 for x in k):
                v = [Fraction(0)] * ncol
                sgn = 1 if k[0] > 0 else -1
                for j, x in zip(S, k):
                    v[j] = sgn * x
                rays.append(primitive_int(v))
                supports.append(S)
    return sorted(rays)


# ----------------------------------------------------------------------------------------------
# networks
# ----------------------------------------------------------------------------------------------
class Network:
    def __init__(self, n, reactions, names=None):
        self.n = n
        self.reactions = [(tuple(a), tuple(b)) for a, b in reactions]
        self.r = len(self.reactions)
        self.names = names or [f"X{i+1}" for i in range(n)]
        for a, b in self.reactions:
            assert len(a) == n and len(b) == n

    # matrices -------------------------------------------------------------------------------
    @property
    def Ye(self):
        return [[self.reactions[j][0][i] for j in range(self.r)] for i in range(self.n)]

    @property
    def Yp(self):
        return [[self.reactions[j][1][i] for j in range(self.r)] for i in range(self.n)]

    @property
    def Gamma(self):
        return [[self.reactions[j][1][i] - self.reactions[j][0][i] for j in range(self.r)] for i in range(self.n)]

    def dual(self):
        """N*: species = reactions of N (k_1..k_r); reaction i of N* = k^{Ye row i} -> k^{Yp row i}."""
        Ye, Yp = self.Ye, self.Yp
        return Network(self.r, [(tuple(Ye[i]), tuple(Yp[i])) for i in range(self.n)],
                       names=[f"k{j+1}" for j in range(self.r)])

    # cones ----------------------------------------------------------------------------------
    def flux_rays(self):
        """Extreme rays of ker(Gamma) cap R^r_{>=0}; their supports are the internal cycles."""
        return extreme_rays_nonneg_kernel(self.Gamma, self.r)

    def conservation_rays(self):
        """Extreme rays of left-ker(Gamma) cap R^n_{>=0} (minimal semi-positive conservation laws)."""
        GT = [list(col) for col in zip(*self.Gamma)]
        return extreme_rays_nonneg_kernel(GT, self.n)

    def is_consistent(self):
        """Exists v > 0 with Gamma v = 0  <=>  some rate constants give a positive steady state."""
        rays = self.flux_rays()
        return all(any(v[j] > 0 for v in rays) for j in range(self.r))

    def is_conservative(self):
        rays = self.conservation_rays()
        return all(any(z[i] > 0 for z in rays) for i in range(self.n))

    def rank(self):
        return rank(self.Gamma)

    # hypotheses on rates --------------------------------------------------------------------
    def rows_distinct(self):
        """No two species with exactly the same rate: rows of Gamma pairwise distinct (xdot_i != xdot_j)."""
        G = self.Gamma
        return len(set(map(tuple, G))) == self.n

    def rows_pairwise_independent(self):
        """Stronger: no xdot_i is a constant multiple of xdot_j (rows pairwise linearly independent)."""
        G = self.Gamma
        return all(rank([G[i], G[j]]) == 2 for i, j in combinations(range(self.n), 2))

    def columns_distinct(self):
        G = self.Gamma
        cols = list(zip(*G))
        return len(set(cols)) == self.r

    def columns_pairwise_independent(self):
        cols = [list(c) for c in zip(*self.Gamma)]
        return all(rank([cols[i], cols[j]]) == 2 for i, j in combinations(range(self.r), 2))

    # complexes, weak reversibility, deficiency ---------------------------------------------
    def complexes(self):
        cs = []
        for a, b in self.reactions:
            for c in (a, b):
                if c not in cs:
                    cs.append(c)
        return cs

    def linkage_classes(self):
        cs = self.complexes()
        idx = {c: i for i, c in enumerate(cs)}
        parent = list(range(len(cs)))

        def find(x):
            while parent[x] != x:
                parent[x] = parent[parent[x]]
                x = parent[x]
            return x
        for a, b in self.reactions:
            parent[find(idx[a])] = find(idx[b])
        return len({find(i) for i in range(len(cs))})

    def deficiency(self):
        return len(self.complexes()) - self.linkage_classes() - self.rank()

    def weakly_reversible(self):
        cs = self.complexes()
        idx = {c: i for i, c in enumerate(cs)}
        m = len(cs)
        reach = [[i == j for j in range(m)] for i in range(m)]
        for a, b in self.reactions:
            reach[idx[a]][idx[b]] = True
        for k in range(m):
            for i in range(m):
                if reach[i][k]:
                    for j in range(m):
                        if reach[k][j]:
                            reach[i][j] = True
        return all(reach[idx[b]][idx[a]] for a, b in self.reactions)

    # siphons --------------------------------------------------------------------------------
    def is_siphon(self, S):
        S = set(S)
        if not S:
            return False
        for a, b in self.reactions:
            prod = {i for i in range(self.n) if b[i] > 0}
            src = {i for i in range(self.n) if a[i] > 0}
            if S & prod and not (S & src):
                return False
        return True

    def siphons(self):
        out = []
        for size in range(1, self.n + 1):
            for S in combinations(range(self.n), size):
                if self.is_siphon(S):
                    out.append(frozenset(S))
        return out

    def minimal_siphons(self):
        sy = self.siphons()
        return sorted([S for S in sy if not any(T < S for T in sy)], key=lambda s: (len(s), sorted(s)))

    # preclusters ----------------------------------------------------------------------------
    def doubling_pairs(self):
        """Pairs {i, j} of distinct reactions with the same source (edges of the graph G of Definition 19)."""
        return [(i, j) for i, j in combinations(range(self.r), 2) if self.reactions[i][0] == self.reactions[j][0]]

    def n_vectors(self):
        """n_1..n_r: rows of the matrix whose columns are the extreme rays of the flux cone."""
        rays = self.flux_rays()
        return [tuple(v[j] for v in rays) for j in range(self.r)]

    def preclusters(self, variant="onestep", return_graph=False):
        """Preclusters of N.

        variant 'onestep': Definition 19 of arXiv:2603.08767v2 -- adjoin {i,s},{j,s} whenever n_s is in
                           span{n_i, n_j} for an edge {i,j} of G (G = same-source graph), once;
        variant 'thesis' : Definition 4.3 / Algorithm 4.4 of Iosif's thesis -- starting from P = G, while there
                           is an edge {i,j} of P and s outside the component of i, j with n_s in span{n_i,n_j},
                           adjoin {i,s},{j,s} (edges and s scanned in increasing order);
        variant 'closure': like 'thesis' but without the component restriction, iterated to a fixed point
                           (the result does not depend on the order)."""
        nv = self.n_vectors()
        r = self.r
        G = self.doubling_pairs()

        @lru_cache(maxsize=None)
        def in_span(s, i, j):
            base = [list(nv[i]), list(nv[j])]
            if len(nv[0]) == 0:
                return True
            return rank(base + [list(nv[s])]) == rank(base)

        edges = set(G)
        if variant == "onestep":
            for (i, j) in G:
                for s in range(r):
                    if s not in (i, j) and in_span(s, i, j):
                        edges.add(tuple(sorted((i, s))))
                        edges.add(tuple(sorted((j, s))))
        elif variant == "closure":
            changed = True
            while changed:
                changed = False
                for (i, j) in sorted(edges):
                    for s in range(r):
                        if s in (i, j) or not in_span(s, i, j):
                            continue
                        for e in (tuple(sorted((i, s))), tuple(sorted((j, s)))):
                            if e not in edges:
                                edges.add(e)
                                changed = True
        elif variant == "thesis":
            while True:
                comp = components(r, edges)
                found = None
                for (i, j) in sorted(edges):
                    for s in range(r):
                        if comp[s] != comp[i] and in_span(s, i, j):
                            found = (i, j, s)
                            break
                    if found:
                        break
                if not found:
                    break
                i, j, s = found
                edges.add(tuple(sorted((i, s))))
                edges.add(tuple(sorted((j, s))))
        else:
            raise ValueError(variant)
        part = partition_from_components(r, edges)
        return (part, edges) if return_graph else part

    # MIPS -----------------------------------------------------------------------------------
    def conservation_columns(self):
        """For each species i, the vector (z^1_i, ..., z^m_i) of values of the extreme conservation laws."""
        rays = self.conservation_rays()
        return [tuple(z[i] for z in rays) for i in range(self.n)]

    def mips_ray(self):
        """MIPS, 'ray' reading: species i, j in the same block iff their columns in a matrix of conservation laws
        are positively proportional (Remark 18: the variables decorating one ray of the chamber decomposition).
        Requires N conservative (then the extreme conservation laws span the conservation space)."""
        cols = self.conservation_columns()
        key = {}
        for i, c in enumerate(cols):
            k = direction_key(c)
            key.setdefault(k, []).append(i)
        return normalize_partition(key.values())

    def mips_semiflow(self):
        """Alternative reading: supports of the extreme rays of the conservation cone."""
        return sorted({frozenset(i for i in range(self.n) if z[i] > 0) for z in self.conservation_rays()},
                      key=lambda s: (len(s), sorted(s)))

    def internal_cycles(self):
        return sorted({frozenset(j for j in range(self.r) if v[j] > 0) for v in self.flux_rays()},
                      key=lambda s: (len(s), sorted(s)))

    def describe(self):
        def cplx(c):
            terms = []
            for i, a in enumerate(c):
                if a == 1:
                    terms.append(self.names[i])
                elif a > 1:
                    terms.append(f"{a}{self.names[i]}")
            return " + ".join(terms) if terms else "0"
        return ", ".join(f"{cplx(a)} -> {cplx(b)}" for a, b in self.reactions)


# ----------------------------------------------------------------------------------------------
# helpers on partitions
# ----------------------------------------------------------------------------------------------
def direction_key(c):
    """Key of the open ray R_{>0} c (c a nonzero vector)."""
    c = [Fraction(x) for x in c]
    nz = [x for x in c if x != 0]
    if not nz:
        return ("zero",)
    f = abs(nz[0])
    return tuple(x / f for x in c)


def components(m, edges):
    parent = list(range(m))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    for i, j in edges:
        parent[find(i)] = find(j)
    return [find(i) for i in range(m)]


def partition_from_components(m, edges):
    comp = components(m, edges)
    blocks = {}
    for i in range(m):
        blocks.setdefault(comp[i], []).append(i)
    return normalize_partition(blocks.values())


def normalize_partition(blocks):
    return sorted((frozenset(b) for b in blocks), key=lambda s: (min(s), len(s)))


def refines(P, Q):
    """Every block of P is contained in a block of Q."""
    return all(any(p <= q for q in Q) for p in P)


def fmt(blocks, names):
    return "{" + ", ".join("{" + ",".join(names[i] for i in sorted(b)) + "}" for b in blocks) + "}"


# ----------------------------------------------------------------------------------------------
# 'set' reading of Definition 17 (only for conservation rank q <= 2)
# ----------------------------------------------------------------------------------------------
def vertex_supports(Zcols, c):
    """Vertex supports of P_c = {x >= 0 : Z x = c}: sets T with {Z e_t} linearly independent and
    c in the relative interior of cone(Z e_t : t in T) (all coefficients > 0).  c = 0 gives {emptyset}."""
    n = len(Zcols)
    q = len(c)
    if all(x == 0 for x in c):
        return frozenset([frozenset()])
    out = set()
    for size in range(1, q + 1):
        for T in combinations(range(n), size):
            M = [[Zcols[t][k] for t in T] for k in range(q)]  # q x |T|
            if rank(M) != size:
                continue
            # solve M a = c
            aug = [M[k] + [c[k]] for k in range(q)]
            R, piv = rref(aug)
            if size in piv:  # inconsistent
                continue
            a = [Fraction(0)] * size
            for row, p in zip(R, piv):
                a[p] = row[size]
            if all(x > 0 for x in a):
                out.add(frozenset(T))
    return frozenset(out)


def decorated_types_q2(N):
    """All decorated combinatorial types {vertex supports of P_c : c in the cone of conservation values},
    computed exactly when the conservation space has dimension q <= 2.  Returns None if q > 2."""
    rays = N.conservation_rays()
    # basis of the conservation space: rref of the extreme rays
    R, piv = rref([list(z) for z in rays])
    q = len(R)
    if q > 2:
        return None
    Zcols = [tuple(R[k][i] for k in range(q)) for i in range(N.n)]
    reps = [tuple(Fraction(0) for _ in range(q))]
    dirs = {}
    for c in Zcols:
        dirs[direction_key(c)] = c
    D = list(dirs.values())
    reps += D
    if q == 2:
        # sort directions by angle inside the pointed cone (all columns are in a pointed cone, since there is a
        # positive conservation law); consecutive pairs bound the 2-dimensional chambers
        def cross(u, v):
            return u[0] * v[1] - u[1] * v[0]
        # sort by angle exactly, comparing via cross products (valid since all directions lie in an open
        # half-plane: a positive conservation law w gives w-weighted coordinates > 0)
        from functools import cmp_to_key
        D_sorted = sorted(D, key=cmp_to_key(lambda u, v: -1 if cross(u, v) > 0 else (1 if cross(u, v) < 0 else 0)))
        for u, v in zip(D_sorted, D_sorted[1:]):
            reps.append(tuple(a + b for a, b in zip(u, v)))
    return {vertex_supports(Zcols, c) for c in reps}


def mips_set_q2(N):
    """'Set' reading of Definition 17: i ~ j iff the transposition (i j) maps the set of decorated types onto
    itself.  Only for q <= 2; returns None otherwise."""
    types = decorated_types_q2(N)
    if types is None:
        return None

    def swap(T, i, j):
        def s(x):
            return j if x == i else (i if x == j else x)
        return frozenset(frozenset(s(x) for x in U) for U in T)
    edges = []
    for i, j in combinations(range(N.n), 2):
        if {swap(T, i, j) for T in types} == types:
            edges.append((i, j))
    return partition_from_components(N.n, edges)
