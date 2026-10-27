"""V1 verifier library: an independent implementation (standard library only, exact Fractions).

Written by verifier V1 (Claude Opus 5.5, effort xhigh), 27 October 2026. It does NOT import crn.py, petri.py,
pclust.py or ecore.py; those are imported only by compare_*.py scripts, for comparison.

Conventions (manuscript / source): network N with species 0..n-1 and reactions y_j -> y'_j, Ye, Yp (n x r),
Gamma = Yp - Ye. Dual N*: species = reactions of N, one reaction per species i of N with source row i of Ye and
product row i of Yp, i.e. Ye(N*) = Ye^T, Yp(N*) = Yp^T.
"""
from fractions import Fraction as Fr
from itertools import combinations
import re

# ---------------------------------------------------------------- exact linear algebra

def rref(M):
    """Reduced row echelon form of a list-of-lists matrix (Fractions). Returns (R, pivots)."""
    A = [[Fr(x) for x in row] for row in M]
    if not A:
        return A, []
    rows, cols = len(A), len(A[0])
    piv = []
    r = 0
    for c in range(cols):
        p = None
        for i in range(r, rows):
            if A[i][c] != 0:
                p = i
                break
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


def nullspace(M, ncols=None):
    """Basis (list of vectors) of {x : M x = 0}."""
    if not M:
        n = ncols
        return [[Fr(int(i == j)) for i in range(n)] for j in range(n)]
    n = len(M[0])
    R, piv = rref(M)
    free = [c for c in range(n) if c not in piv]
    basis = []
    for f in free:
        v = [Fr(0)] * n
        v[f] = Fr(1)
        for i, p in enumerate(piv):
            v[p] = -R[i][f]
        basis.append(v)
    return basis


def transpose(M):
    return [list(x) for x in zip(*M)] if M else []


def in_span(v, vecs):
    if not vecs:
        return all(x == 0 for x in v)
    return rank([list(u) for u in vecs] + [list(v)]) == rank([list(u) for u in vecs])


def matmul(A, B):
    return [[sum(A[i][k] * B[k][j] for k in range(len(B))) for j in range(len(B[0]))] for i in range(len(A))]


# ---------------------------------------------------------------- networks

class Net:
    def __init__(self, Ye, Yp, species=None, name=""):
        self.Ye = [[int(x) for x in row] for row in Ye]
        self.Yp = [[int(x) for x in row] for row in Yp]
        self.n = len(Ye)
        self.r = len(Ye[0]) if Ye else 0
        self.species = species or ["X%d" % (i + 1) for i in range(self.n)]
        self.name = name
        self.G = [[self.Yp[i][j] - self.Ye[i][j] for j in range(self.r)] for i in range(self.n)]
        self._cache = {}

    def src(self, j):
        return tuple(self.Ye[i][j] for i in range(self.n))

    def prod(self, j):
        return tuple(self.Yp[i][j] for i in range(self.n))

    def dual(self):
        return Net(transpose(self.Ye), transpose(self.Yp), ["k%d" % (j + 1) for j in range(self.r)],
                   name=self.name + "*")

    def __repr__(self):
        return "; ".join(cplx_str(self.src(j), self.species) + " -> " + cplx_str(self.prod(j), self.species)
                         for j in range(self.r))


def cplx_str(y, names):
    t = []
    for i, c in enumerate(y):
        if c:
            t.append((str(c) if c > 1 else "") + names[i])
    return "+".join(t) if t else "0"


def parse_cplx(s, idx):
    s = s.strip()
    y = {}
    if s == "0":
        return y
    for term in s.split("+"):
        term = term.strip()
        m = re.match(r"^(\d*)\s*([A-Za-z_][A-Za-z_0-9']*)$", term)
        if not m:
            raise ValueError(term)
        c = int(m.group(1)) if m.group(1) else 1
        y[m.group(2)] = y.get(m.group(2), 0) + c
        if m.group(2) not in idx:
            idx[m.group(2)] = len(idx)
    return y


def net_from_str(text, species_order=None, name=""):
    """'A -> B, B -> 2A' ; '<->' or '<=>' gives both directions. Separators ',' or ';'."""
    idx = {}
    if species_order:
        for s in species_order:
            idx[s] = len(idx)
    reacs = []
    for part in re.split(r"[;,]", text):
        part = part.strip()
        if not part:
            continue
        if "<->" in part or "<=>" in part:
            a, b = re.split(r"<->|<=>", part)
            ya, yb = parse_cplx(a, idx), parse_cplx(b, idx)
            reacs += [(ya, yb), (yb, ya)]
        else:
            a, b = part.split("->")
            reacs.append((parse_cplx(a, idx), parse_cplx(b, idx)))
    names = [None] * len(idx)
    for s, i in idx.items():
        names[i] = s
    n, r = len(names), len(reacs)
    Ye = [[reacs[j][0].get(names[i], 0) for j in range(r)] for i in range(n)]
    Yp = [[reacs[j][1].get(names[i], 0) for j in range(r)] for i in range(n)]
    return Net(Ye, Yp, names, name)


# ---------------------------------------------------------------- cones

def extreme_rays(A, ncols):
    """Extreme rays of {x >= 0 : A x = 0} (A given as list of rows over ncols columns), by minimal supports.
    Returns list of (support tuple, vector with integer-normalised entries)."""
    rays = []
    found_supports = []
    cols = list(range(ncols))
    for size in range(1, ncols + 1):
        for S in combinations(cols, size):
            Sset = set(S)
            if any(set(F) <= Sset for F in found_supports):
                continue  # contains a smaller support: cannot be minimal
            sub = [[row[c] for c in S] for row in A] if A else []
            ns = nullspace(sub, ncols=len(S)) if sub else [[Fr(int(i == j)) for i in range(len(S))] for j in range(len(S))]
            if len(ns) != 1:
                continue
            v = ns[0]
            if all(x > 0 for x in v):
                pass
            elif all(x < 0 for x in v):
                v = [-x for x in v]
            else:
                continue
            full = [Fr(0)] * ncols
            for c, x in zip(S, v):
                full[c] = x
            m = min(x for x in full if x > 0)
            full = [x / m for x in full]
            rays.append((S, full))
            found_supports.append(S)
    return rays


def flux_rays(N):
    if "frays" not in N._cache:
        N._cache["frays"] = extreme_rays(N.G, N.r)
    return N._cache["frays"]


def cons_rays(N):
    if "crays" not in N._cache:
        N._cache["crays"] = extreme_rays(transpose(N.G), N.n)
    return N._cache["crays"]


def is_H1(N):
    """Conservative: exists z > 0 with z^T Gamma = 0 <=> union of supports of conservation rays = all species."""
    sup = set()
    for S, _ in cons_rays(N):
        sup |= set(S)
    return len(sup) == N.n


def is_H2(N):
    sup = set()
    for S, _ in flux_rays(N):
        sup |= set(S)
    return len(sup) == N.r


def is_H3(N):
    rows = [tuple(r) for r in N.G]
    return len(set(rows)) == len(rows)


def srank(N):
    return rank(N.G)


def nvec(N):
    """Rows n_t of the matrix whose columns are the extreme flux rays."""
    R = flux_rays(N)
    return [[R[k][1][t] for k in range(len(R))] for t in range(N.r)]


def wvec(N):
    """Columns w_i of the matrix whose rows are the extreme conservation laws."""
    R = cons_rays(N)
    return [[R[k][1][i] for k in range(len(R))] for i in range(N.n)]


# ---------------------------------------------------------------- partitions

def canon(blocks):
    return sorted(sorted(b) for b in blocks)


def components(nv, edges):
    par = list(range(nv))

    def f(x):
        while par[x] != x:
            par[x] = par[par[x]]
            x = par[x]
        return x
    for a, b in edges:
        ra, rb = f(a), f(b)
        if ra != rb:
            par[ra] = rb
    comp = {}
    for x in range(nv):
        comp.setdefault(f(x), []).append(x)
    return canon(comp.values())


def refines(A, B):
    """A refines B (every block of A inside a block of B)."""
    where = {}
    for k, b in enumerate(B):
        for x in b:
            where[x] = k
    return all(len({where[x] for x in a}) == 1 for a in A)


def join(A, B, nv):
    edges = []
    for P in (A, B):
        for b in P:
            edges += [(b[0], x) for x in b[1:]]
    return components(nv, edges)


def meet(A, B):
    wa = {x: k for k, b in enumerate(A) for x in b}
    wb = {x: k for k, b in enumerate(B) for x in b}
    d = {}
    for x in wa:
        d.setdefault((wa[x], wb[x]), []).append(x)
    return canon(d.values())


def prop_classes(vecs):
    """Classes of indices whose vectors are positively proportional (nonzero vectors). Zero vectors: singletons."""
    nv = len(vecs)
    edges = []
    for a in range(nv):
        for b in range(a + 1, nv):
            if any(vecs[a]) and any(vecs[b]) and positively_proportional(vecs[a], vecs[b]):
                edges.append((a, b))
    return components(nv, edges)


def positively_proportional(u, v):
    if rank([u, v]) != 1:
        return False
    # same direction?
    for x, y in zip(u, v):
        if x != 0:
            return (x > 0) == (y > 0)
    return False


# ---------------------------------------------------------------- Iosif's objects

def doubling_pairs(N):
    return [(i, j) for i, j in combinations(range(N.r), 2) if N.src(i) == N.src(j)]


def twin_pairs(N):
    """Pairs of species with identical rows of Ye (zero rows included, as in the manuscript's definition)."""
    return [(i, j) for i, j in combinations(range(N.n), 2) if N.Ye[i] == N.Ye[j]]


def preclusters_def19(N):
    """Definition 19 of arXiv:2603.08767v2, one step: G = doubling pairs; adjoin {i,s},{j,s} whenever n_s in
    span{n_i, n_j} for each edge {i,j} of G."""
    nv = nvec(N)
    E = list(doubling_pairs(N))
    G = list(E)
    for i, j in E:
        for s in range(N.r):
            if s in (i, j):
                continue
            if in_span(nv[s], [nv[i], nv[j]]):
                G += [(i, s), (j, s)]
    return components(N.r, G)


def preclusters_thesis(N):
    """Thesis Algorithm 4.4 (iterated; s outside the component of i and j)."""
    nv = nvec(N)
    edges = set(doubling_pairs(N))
    changed = True
    while changed:
        changed = False
        comps = components(N.r, list(edges))
        where = {x: k for k, b in enumerate(comps) for x in b}
        for (i, j) in sorted(edges):
            for s in range(N.r):
                if where[s] == where[i]:
                    continue
                if in_span(nv[s], [nv[i], nv[j]]):
                    edges |= {(min(i, s), max(i, s)), (min(j, s), max(j, s))}
                    changed = True
                    break
            if changed:
                break
    return components(N.r, list(edges))


def ray_mips(N):
    return prop_classes(wvec(N))


def full_coupling(N):
    return prop_classes(nvec(N))


# ---------------------------------------------------------------- siphons, traps, Petri notions

def siphons(N):
    """All siphons (nonempty sets S of species: every reaction with a product in S has a reactant in S)."""
    out = []
    for size in range(1, N.n + 1):
        for S in combinations(range(N.n), size):
            ok = True
            for j in range(N.r):
                if any(N.Yp[i][j] > 0 for i in S) and not any(N.Ye[i][j] > 0 for i in S):
                    ok = False
                    break
            if ok:
                out.append(frozenset(S))
    return out


def minimal(sets):
    sets = sorted(set(sets), key=len)
    out = []
    for s in sets:
        if not any(t <= s for t in out):
            out.append(s)
    return out


def min_siphons(N):
    if "msiph" not in N._cache:
        N._cache["msiph"] = minimal(siphons(N))
    return N._cache["msiph"]


def t_traps(N):
    """Nonempty reaction sets U with U-dot (species produced) subset of dot-U (species consumed)."""
    out = []
    for size in range(1, N.r + 1):
        for U in combinations(range(N.r), size):
            prod = {i for i in range(N.n) for j in U if N.Yp[i][j] > 0}
            cons = {i for i in range(N.n) for j in U if N.Ye[i][j] > 0}
            if prod <= cons:
                out.append(frozenset(U))
    return out


def t_siphons(N):
    out = []
    for size in range(1, N.r + 1):
        for U in combinations(range(N.r), size):
            prod = {i for i in range(N.n) for j in U if N.Yp[i][j] > 0}
            cons = {i for i in range(N.n) for j in U if N.Ye[i][j] > 0}
            if cons <= prod:
                out.append(frozenset(U))
    return out


def traps(N):
    out = []
    for size in range(1, N.n + 1):
        for S in combinations(range(N.n), size):
            ok = True
            for j in range(N.r):
                if any(N.Ye[i][j] > 0 for i in S) and not any(N.Yp[i][j] > 0 for i in S):
                    ok = False
                    break
            if ok:
                out.append(frozenset(S))
    return out


def critical_siphons(N):
    """Minimal siphons containing the support of no nonzero semi-positive conservation law."""
    sups = [frozenset(S) for S, _ in cons_rays(N)]
    return [R for R in min_siphons(N) if not any(Z <= R for Z in sups)]


def clusters(N):
    """Components of the bipartite graph species-reactions with an edge iff Ye[i][j] > 0.
    Returned as list of (species set, reaction set); vertices numbered species 0..n-1, reactions n..n+r-1."""
    edges = [(i, N.n + j) for i in range(N.n) for j in range(N.r) if N.Ye[i][j] > 0]
    comps = components(N.n + N.r, edges)
    return [(frozenset(x for x in c if x < N.n), frozenset(x - N.n for x in c if x >= N.n)) for c in comps]


def is_connected(N):
    edges = [(i, N.n + j) for i in range(N.n) for j in range(N.r) if N.Ye[i][j] > 0 or N.Yp[i][j] > 0]
    return len(components(N.n + N.r, edges)) == 1


def strongly_connected(N):
    """Petri graph: arcs i->j if Ye>0, j->i if Yp>0."""
    V = N.n + N.r
    adj = {v: set() for v in range(V)}
    radj = {v: set() for v in range(V)}
    for i in range(N.n):
        for j in range(N.r):
            if N.Ye[i][j] > 0:
                adj[i].add(N.n + j); radj[N.n + j].add(i)
            if N.Yp[i][j] > 0:
                adj[N.n + j].add(i); radj[i].add(N.n + j)

    def reach(start, g):
        seen = {start}
        st = [start]
        while st:
            x = st.pop()
            for y in g[x]:
                if y not in seen:
                    seen.add(y)
                    st.append(y)
        return seen
    return len(reach(0, adj)) == V and len(reach(0, radj)) == V


def is_ordinary(N):
    return all(x in (0, 1) for row in N.Ye + N.Yp for x in row)


def is_EC(N):
    for i, j in combinations(range(N.r), 2):
        a, b = N.src(i), N.src(j)
        if any(x > 0 and y > 0 for x, y in zip(a, b)) and a != b:
            return False
    return True


def balanced_joins(N):
    for j in range(N.r):
        nz = {x for x in N.src(j) if x}
        if len(nz) > 1:
            return False
    return True


def is_pure(N):
    return all(not (N.Ye[i][j] and N.Yp[i][j]) for i in range(N.n) for j in range(N.r))


def n_sources(N):
    return len({N.src(j) for j in range(N.r)})


def is_WFEC(N):
    return (is_EC(N) and is_connected(N) and is_H1(N) and is_H2(N) and srank(N) == n_sources(N) - 1)


def is_WFFC(N):
    return is_WFEC(N) and is_ordinary(N)


# ---------------------------------------------------------------- corrected preclusters (Q2)

def source_classes(N):
    d = {}
    for j in range(N.r):
        d.setdefault(N.src(j), []).append(j)
    keys = list(d.keys())
    return keys, [d[k] for k in keys]


def A_V(N, V):
    """n x m matrix with column c = sum over reactions t with source c of V_t (y'_t - y_t)."""
    keys, cls = source_classes(N)
    return [[sum(Fr(V[t]) * N.G[i][t] for t in cls[c]) for c in range(len(cls))] for i in range(N.n)]


def kappa_V(N, V):
    """kappa_V lifted to reactions: source classes c, c' equivalent iff e_c - e_c' in rowspace(A_V)
    (iff u_c = u_c' for every u in ker A_V)."""
    keys, cls = source_classes(N)
    m = len(cls)
    A = A_V(N, V)
    K = nullspace(A, ncols=m)
    edges = []
    for a in range(m):
        for b in range(a + 1, m):
            if all(u[a] == u[b] for u in K):
                edges.append((a, b))
    comps = components(m, edges)
    return canon([[t for c in comp for t in cls[c]] for comp in comps])


def flat_closure_partition(N):
    """pi(N): finest partition coarser than the source classes all of whose blocks are flats of the vector
    configuration (n_t)."""
    nv = nvec(N)
    keys, cls = source_classes(N)
    blocks = [list(c) for c in cls]
    changed = True
    while changed:
        changed = False
        for a in range(len(blocks)):
            B = blocks[a]
            vecs = [nv[t] for t in B]
            for b in range(len(blocks)):
                if a == b:
                    continue
                if any(in_span(nv[s], vecs) for s in blocks[b]):
                    blocks[a] = B + blocks[b]
                    del blocks[b]
                    changed = True
                    break
            if changed:
                break
    return canon(blocks)


def class_span_closure(N):
    """P_cs(N): merge a block C into a block B with |B| >= 2 when some s in C has n_s in span(B); iterate."""
    nv = nvec(N)
    keys, cls = source_classes(N)
    blocks = [list(c) for c in cls]
    changed = True
    while changed:
        changed = False
        for a in range(len(blocks)):
            B = blocks[a]
            if len(B) < 2:
                continue
            vecs = [nv[t] for t in B]
            for b in range(len(blocks)):
                if a == b:
                    continue
                if any(in_span(nv[s], vecs) for s in blocks[b]):
                    blocks[a] = B + blocks[b]
                    del blocks[b]
                    changed = True
                    break
            if changed:
                break
    return canon(blocks)


def interior_flux(N, weights):
    """V = sum_k weights[k] * ray_k (weights positive) : a point of C° = relative interior of the flux cone."""
    R = flux_rays(N)
    V = [Fr(0)] * N.r
    for w, (_, v) in zip(weights, R):
        V = [a + Fr(w) * b for a, b in zip(V, v)]
    return V


def linkage_classes(N):
    """Components of the reaction graph on complexes; returned as partition of reactions."""
    cpx = {}
    edges = []
    for j in range(N.r):
        a = cpx.setdefault(N.src(j), len(cpx))
        b = cpx.setdefault(N.prod(j), len(cpx))
        edges.append((a, b))
    comps = components(len(cpx), edges)
    where = {x: k for k, b in enumerate(comps) for x in b}
    d = {}
    for j in range(N.r):
        d.setdefault(where[cpx[N.src(j)]], []).append(j)
    return canon(d.values()), len(cpx), len(comps)


def is_WR(N):
    cpx = {}
    adj = {}
    for j in range(N.r):
        a = cpx.setdefault(N.src(j), len(cpx))
        b = cpx.setdefault(N.prod(j), len(cpx))
        adj.setdefault(a, set()).add(b)
        adj.setdefault(b, set())

    def reach(s):
        seen = {s}; st = [s]
        while st:
            x = st.pop()
            for y in adj[x]:
                if y not in seen:
                    seen.add(y); st.append(y)
        return seen
    for j in range(N.r):
        if cpx[N.src(j)] not in reach(cpx[N.prod(j)]):
            return False
    return True


def deficiency(N):
    lc, m, l = linkage_classes(N)
    return m - l - srank(N)


# ---------------------------------------------------------------- operations (Q2 section 3)

def duplicate(N, t):
    Ye = [row + [row[t]] for row in N.Ye]
    Yp = [row + [row[t]] for row in N.Yp]
    return Net(Ye, Yp, N.species, N.name + "dup%d" % t)


def split_choice(N, A):
    """A: set of reactions with a common source y. New species W; each t in A becomes W -> y'_t; add t_A: y -> W
    (appended as the last reaction). Returns the new network."""
    A = sorted(A)
    y = N.src(A[0])
    assert all(N.src(t) == y for t in A)
    Ye = [list(row) for row in N.Ye] + [[0] * N.r]
    Yp = [list(row) for row in N.Yp] + [[0] * N.r]
    for t in A:
        for i in range(N.n):
            Ye[i][t] = 0
        Ye[N.n][t] = 1
    for i in range(N.n):
        Ye[i].append(y[i]); Yp[i].append(0)
    Ye[N.n].append(0); Yp[N.n].append(1)
    return Net(Ye, Yp, N.species + ["W"], N.name + "split")
