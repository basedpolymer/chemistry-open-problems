#!/usr/bin/env python3
"""Certificates for the manuscript "Iosif's network duality is the Petri-net reverse dual" (result 17,
chemistry-open-problems, rewritten 27-28 October 2026).

Every explicit example and every finite claim stated in paper.tex is checked here in exact arithmetic:
networks, matrices, hypotheses (H1)-(H3), duals, conservation laws, fluxes, siphons, preclusters and the
corrected preclusters of Section 3, steady states, Jacobians, quadratic-form certificates, polynomial identities
used in the proofs of the examples, and the counterexamples (NC3, NC, NM, NB, NT, N0, N1, N2, U, the extinction
examples, the rank-one and rank-two examples, Theorems A, AF (u = 3..7, R = 1..7) and O).

Environment: Python >= 3.10, standard library only (fractions, itertools), plus the repository library
scripts/crn.py (imported, not modified).  No floating point is used in any claim.  Runtime: about 1 s
(Windows, Python 3.14; WSL Fedora 44, Python 3.12).
Usage, from the repository root:  py -3 problems/iosif-siphon-precluster-duality/scripts/certify_duality.py
The log is written to problems/iosif-siphon-precluster-duality/certificates/certify_duality.log and printed.
The environment, the runtime and the claim count are printed at the top of the log.  Exit status 1 on any FAIL.

Claim names start with the section of the paper (S2 ... S7) and the name of the example used there.
"""
import os
import platform
import sys
import time
from fractions import Fraction as Fr
from itertools import combinations, permutations, product

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.dont_write_bytecode = True
import crn  # noqa: E402  (exact library of result 17; not modified)

LINES = []
NPASS = 0
NFAIL = 0


def claim(name, ok, detail=""):
    global NPASS, NFAIL
    ok = bool(ok)
    if ok:
        NPASS += 1
    else:
        NFAIL += 1
    LINES.append(("PASS" if ok else "FAIL") + ": " + name + ((" -- " + detail) if detail else ""))


def section(title):
    LINES.append("")
    LINES.append("## " + title)


# =============================================================================================================
# 1. Exact multivariate polynomials with rational coefficients
# =============================================================================================================
class Poly:
    """Sparse polynomial: dict monomial -> Fraction, monomial = tuple of (variable name, exponent), sorted."""
    __slots__ = ("t",)

    def __init__(self, t=None):
        self.t = {}
        if t:
            for m, c in t.items():
                if c != 0:
                    self.t[m] = Fr(c)

    @staticmethod
    def const(c):
        return Poly({(): Fr(c)}) if c != 0 else Poly()

    @staticmethod
    def var(name):
        return Poly({((name, 1),): Fr(1)})

    @staticmethod
    def lift(x):
        return x if isinstance(x, Poly) else Poly.const(x)

    def __add__(self, o):
        o = Poly.lift(o)
        t = dict(self.t)
        for m, c in o.t.items():
            t[m] = t.get(m, 0) + c
        return Poly(t)

    __radd__ = __add__

    def __neg__(self):
        return Poly({m: -c for m, c in self.t.items()})

    def __sub__(self, o):
        return self + (-Poly.lift(o))

    def __rsub__(self, o):
        return Poly.lift(o) - self

    @staticmethod
    def _mulmono(a, b):
        d = dict(a)
        for v, e in b:
            d[v] = d.get(v, 0) + e
        return tuple(sorted((v, e) for v, e in d.items() if e))

    def __mul__(self, o):
        o = Poly.lift(o)
        t = {}
        for m1, c1 in self.t.items():
            for m2, c2 in o.t.items():
                m = Poly._mulmono(m1, m2)
                t[m] = t.get(m, 0) + c1 * c2
        return Poly(t)

    __rmul__ = __mul__

    def __pow__(self, k):
        r = Poly.const(1)
        for _ in range(k):
            r = r * self
        return r

    def is_zero(self):
        return not self.t

    def __eq__(self, o):
        return (self - Poly.lift(o)).is_zero()

    def __hash__(self):
        return hash(tuple(sorted(self.t.items())))

    def diff(self, v):
        t = {}
        for m, c in self.t.items():
            d = dict(m)
            e = d.get(v, 0)
            if e:
                d[v] = e - 1
                mm = tuple(sorted((x, y) for x, y in d.items() if y))
                t[mm] = t.get(mm, 0) + c * e
        return Poly(t)

    def subs(self, s):
        """s: dict variable -> Poly or number (or any ring element supporting + and *)."""
        acc = None
        for m, c in self.t.items():
            term = c
            for v, e in m:
                val = s.get(v)
                base = Poly.var(v) if val is None else val
                for _ in range(e):
                    term = term * base
            acc = term if acc is None else acc + term
        if acc is None:
            return Fr(0)
        return acc

    def value(self, s):
        """Evaluate with every variable given (returns a field element)."""
        r = self.subs(s)
        if isinstance(r, Poly):
            assert all(m == () for m in r.t), "not all variables substituted"
            return r.t.get((), Fr(0))
        return r

    def degree_in(self, v):
        return max((dict(m).get(v, 0) for m in self.t), default=-1)

    def low_order_in(self, v):
        """(lowest exponent of v, part of the polynomial with that exponent, divided by v^e)."""
        e = min(dict(m).get(v, 0) for m in self.t)
        part = {}
        for m, c in self.t.items():
            d = dict(m)
            if d.get(v, 0) == e:
                d.pop(v, None)
                part[tuple(sorted(d.items()))] = c
        return e, Poly(part)

    def coeffs(self):
        return list(self.t.values())

    def __repr__(self):
        if not self.t:
            return "0"
        out = []
        for m, c in sorted(self.t.items()):
            mono = "*".join(v if e == 1 else f"{v}^{e}" for v, e in m)
            out.append(f"{c}" + ("*" + mono if mono else ""))
        return " + ".join(out)


def V(*names):
    return [Poly.var(x) for x in names]


def det(M):
    """Determinant by Laplace expansion (small matrices, entries in any commutative ring)."""
    n = len(M)
    if n == 0:
        return Fr(1)
    if n == 1:
        return M[0][0]
    if n == 2:
        return M[0][0] * M[1][1] - M[0][1] * M[1][0]
    acc = None
    for j in range(n):
        minor = [row[:j] + row[j + 1:] for row in M[1:]]
        term = M[0][j] * det(minor)
        if j % 2:
            term = -term
        acc = term if acc is None else acc + term
    return acc


def matmul(A, B):
    return [[sum((A[i][k] * B[k][j] for k in range(len(B))), Fr(0)) for j in range(len(B[0]))]
            for i in range(len(A))]


def transpose(A):
    return [list(r) for r in zip(*A)]


# =============================================================================================================
# 2. The quadratic field Q(sqrt d), for the irrational steady states of NC3 and NB
# =============================================================================================================
class QS:
    __slots__ = ("a", "b", "d")

    def __init__(self, a, b=0, d=None):
        self.a, self.b, self.d = Fr(a), Fr(b), d

    def _c(self, o):
        if isinstance(o, QS):
            return o
        return QS(o, 0, self.d)

    def __add__(self, o):
        o = self._c(o)
        return QS(self.a + o.a, self.b + o.b, self.d or o.d)

    __radd__ = __add__

    def __neg__(self):
        return QS(-self.a, -self.b, self.d)

    def __sub__(self, o):
        return self + (-self._c(o))

    def __rsub__(self, o):
        return self._c(o) - self

    def __mul__(self, o):
        o = self._c(o)
        d = self.d or o.d
        return QS(self.a * o.a + (self.b * o.b * d if d else 0), self.a * o.b + self.b * o.a, d)

    __rmul__ = __mul__

    def inv(self):
        n = self.a * self.a - self.b * self.b * (self.d or 0)
        return QS(self.a / n, -self.b / n, self.d)

    def __truediv__(self, o):
        return self * self._c(o).inv()

    def __rtruediv__(self, o):
        return self._c(o) * self.inv()

    def sign(self):
        """Exact sign of a + b sqrt(d), d > 0 not a square."""
        a, b = self.a, self.b
        if b == 0:
            return (a > 0) - (a < 0)
        if a == 0:
            return (b > 0) - (b < 0)
        if (a > 0) == (b > 0):
            return 1 if a > 0 else -1
        # opposite signs: compare a^2 with b^2 d
        big_a = a * a > b * b * self.d
        return (1 if a > 0 else -1) if big_a else (1 if b > 0 else -1)

    def __repr__(self):
        return f"({self.a} + {self.b}*sqrt({self.d}))"


# =============================================================================================================
# 3. Networks: parser, mass-action fields, structural tests
# =============================================================================================================
def parse_complex(txt, species):
    txt = txt.strip()
    v = [0] * len(species)
    if txt == "0":
        return tuple(v)
    for term in txt.split("+"):
        term = term.strip()
        k = 0
        while k < len(term) and term[k].isdigit():
            k += 1
        coef = int(term[:k]) if k else 1
        name = term[k:].strip()
        v[species.index(name)] += coef
    return tuple(v)


def net(spec, species):
    """spec: 'A -> B; 2B -> A+C; ...' with species names from the list `species`."""
    reacs = []
    for part in spec.split(";"):
        part = part.strip()
        if not part:
            continue
        a, b = part.split("->")
        reacs.append((parse_complex(a, species), parse_complex(b, species)))
    return crn.Network(len(species), reacs, names=list(species))


def Xs(n):
    return [f"X{i + 1}" for i in range(n)]


def Ks(n):
    return [f"k{i + 1}" for i in range(n)]


def same_reactions(N, spec, species):
    """N equals the network given by spec, reactions in the same order."""
    M = net(spec, species)
    return N.n == M.n and N.reactions == M.reactions


def same_up_to_relabel(N, M):
    """Equal up to a permutation of species and of reactions."""
    if N.n != M.n or N.r != M.r:
        return False
    target = sorted(M.reactions)
    for p in permutations(range(N.n)):
        R = sorted((tuple(a[p[i]] for i in range(N.n)), tuple(b[p[i]] for i in range(N.n))) for a, b in N.reactions)
        if R == target:
            return True
    return False


def field(N, rates, x):
    """Mass-action field: list of n polynomials f_i = sum_j Gamma_ij kappa_j x^{y_j}."""
    G = N.Gamma
    f = [Poly() for _ in range(N.n)]
    for j, (a, b) in enumerate(N.reactions):
        mono = Poly.lift(rates[j])
        for i in range(N.n):
            if a[i]:
                mono = mono * (x[i] ** a[i])
        for i in range(N.n):
            if G[i][j]:
                f[i] = f[i] + G[i][j] * mono
    return f


def H123(N):
    return N.is_conservative(), N.is_consistent(), N.rows_distinct()


def H_all(N):
    return H123(N) == (True, True, True)


def counts(N):
    """(number of complexes, linkage classes, rank, deficiency)."""
    return len(N.complexes()), N.linkage_classes(), N.rank(), N.deficiency()


def connected(N):
    """Bipartite species-reaction graph (nonzero entries of Ye or Yp) connected."""
    m = N.n + N.r
    edges = []
    for j, (a, b) in enumerate(N.reactions):
        for i in range(N.n):
            if a[i] or b[i]:
                edges.append((i, N.n + j))
    comp = crn.components(m, edges)
    return len(set(comp)) == 1


def clusters(N):
    """Petri clusters: components of the graph with an edge {X_i, R_j} iff (Ye)_ij > 0.
    Returns list of (species set, reaction set)."""
    m = N.n + N.r
    edges = [(i, N.n + j) for j, (a, b) in enumerate(N.reactions) for i in range(N.n) if a[i]]
    comp = crn.components(m, edges)
    out = {}
    for v in range(m):
        out.setdefault(comp[v], (set(), set()))
        if v < N.n:
            out[comp[v]][0].add(v)
        else:
            out[comp[v]][1].add(v - N.n)
    return [(frozenset(s), frozenset(t)) for s, t in out.values()]


def sources(N):
    out = []
    for a, b in N.reactions:
        if a not in out:
            out.append(a)
    return out


def is_EC(N):
    """Two reactions sharing a reactant species have the same source."""
    for (a, _), (c, _) in combinations(N.reactions, 2):
        if a != c and any(a[i] and c[i] for i in range(N.n)):
            return False
    return True


def is_PC(N):
    """Two reactions sharing a reactant species have positively proportional sources."""
    for (a, _), (c, _) in combinations(N.reactions, 2):
        if any(a[i] and c[i] for i in range(N.n)) and crn.rank([list(a), list(c)]) != 1:
            return False
    return True


def is_BH(N):
    """Block-homogeneous: every source is a positive multiple of the indicator of a block of a partition."""
    blocks = []
    for a in sources(N):
        B = frozenset(i for i in range(N.n) if a[i])
        if len({a[i] for i in B}) != 1:
            return False
        blocks.append(B)
    for B, C in combinations(set(blocks), 2):
        if B & C:
            return False
    return True


def balanced_joins(N):
    return all(len({a[i] for i in range(N.n) if a[i]}) <= 1 for a in sources(N))


def ordinary(N):
    return all(x in (0, 1) for a, b in N.reactions for x in a + b)


def pure(N):
    return all(not (a[i] and b[i]) for a, b in N.reactions for i in range(N.n))


def rank_eq_clusters(N):
    return N.rank() == len(clusters(N)) - 1


def WF_EC(N):
    return (is_EC(N) and connected(N) and N.is_conservative() and N.is_consistent()
            and N.rank() == len(sources(N)) - 1)


def cluster_rank_class(N):
    return connected(N) and N.is_conservative() and N.is_consistent() and rank_eq_clusters(N)


def is_T_trap(N, U):
    """Set U of reactions: every species produced by U is consumed by U."""
    prod = {i for j in U for i in range(N.n) if N.reactions[j][1][i]}
    cons = {i for j in U for i in range(N.n) if N.reactions[j][0][i]}
    return bool(U) and prod <= cons


def minimal_T_traps(N):
    tr = [frozenset(U) for k in range(1, N.r + 1) for U in combinations(range(N.r), k) if is_T_trap(N, U)]
    return sorted([U for U in tr if not any(V_ < U for V_ in tr)], key=lambda s: (len(s), sorted(s)))


def law_in(N, S):
    """Some nonzero z >= 0 with z^T Gamma = 0 has support inside S (via extreme rays)."""
    return any(all(z[i] == 0 for i in range(N.n) if i not in S) for z in N.conservation_rays())


def critical_siphons(N):
    return [S for S in N.siphons() if not law_in(N, S)]


def conservation_basis(N):
    """Basis of the left kernel of Gamma (rows)."""
    GT = transpose(N.Gamma)
    return crn.nullspace(GT, N.n)


def twins(N):
    """Pairs of distinct species with equal rows of Ye (zero rows included)."""
    Ye = N.Ye
    return [(i, j) for i, j in combinations(range(N.n), 2) if Ye[i] == Ye[j]]


def conservation_coupled(N, i, j):
    cols = N.conservation_columns()
    return crn.direction_key(cols[i]) == crn.direction_key(cols[j]) and any(cols[i])


def fully_coupled(N, s, t):
    nv = N.n_vectors()
    return crn.direction_key(nv[s]) == crn.direction_key(nv[t])


def full_coupling_partition(N):
    """MIPS(N*) = classes of fully coupled reactions (under (H2))."""
    return N.dual().mips_ray()


def source_classes(N):
    cls = {}
    for j, (a, b) in enumerate(N.reactions):
        cls.setdefault(a, []).append(j)
    return crn.normalize_partition(cls.values())


def in_span(nv, s, B):
    base = [list(nv[t]) for t in B]
    return crn.rank(base + [list(nv[s])]) == crn.rank(base)


def P_cs(N):
    """Class-span closure: merge a block C into a block B with |B| >= 2 if some s in C has n_s in span(B)."""
    nv = N.n_vectors()
    blocks = [set(b) for b in source_classes(N)]
    changed = True
    while changed:
        changed = False
        for B in blocks:
            if len(B) < 2:
                continue
            for C in blocks:
                if C is B:
                    continue
                if any(in_span(nv, s, sorted(B)) for s in C):
                    B |= C
                    blocks.remove(C)
                    changed = True
                    break
            if changed:
                break
    return crn.normalize_partition(blocks)


def flat_closure(N):
    """pi(N): finest partition coarser than the source classes whose blocks are flats of the n-vectors."""
    nv = N.n_vectors()
    blocks = [set(b) for b in source_classes(N)]
    changed = True
    while changed:
        changed = False
        for B in blocks:
            for C in blocks:
                if C is B:
                    continue
                if any(in_span(nv, s, sorted(B)) for s in C):
                    B |= C
                    blocks.remove(C)
                    changed = True
                    break
            if changed:
                break
    return crn.normalize_partition(blocks)


def A_matrix(N, V_):
    """A_V = Gamma diag(V) Q (n x m), m = number of source classes, columns in the order of source_classes."""
    G = N.Gamma
    cls = source_classes(N)
    return [[sum((G[i][t] * V_[t] for t in c), Fr(0)) for c in cls] for i in range(N.n)], cls


def kappa_V(N, V_):
    """Partition of reactions: source classes c, c' equivalent iff u_c = u_c' on ker A_V (1 in ker A_V)."""
    A, cls = A_matrix(N, V_)
    K = crn.nullspace(A, len(cls))
    m = len(cls)
    edges = [(c, d) for c, d in combinations(range(m), 2) if all(u[c] == u[d] for u in K)]
    comp = crn.components(m, edges)
    blocks = {}
    for c in range(m):
        blocks.setdefault(comp[c], set()).update(cls[c])
    return crn.normalize_partition(blocks.values())


def partition_names(P, names):
    return crn.fmt(P, names)


def Rnames(r):
    return [f"R{j + 1}" for j in range(r)]


def is_wr(N):
    return N.weakly_reversible()


# =============================================================================================================
# 4. Rank-one face classification (Theorem S1) and drain certificates
# =============================================================================================================
def in_convex_hull(p, pts):
    """Exact test p in conv(pts) by Caratheodory: solve on affinely independent subsets."""
    if not pts:
        return False
    d = len(p)
    for k in range(1, min(len(pts), d + 1) + 1):
        for sub in combinations(pts, k):
            rows = [[Fr(q[i]) for q in sub] + [Fr(p[i])] for i in range(d)] + [[Fr(1)] * k + [Fr(1)]]
            R, piv = crn.rref(rows)
            if any(all(x == 0 for x in row[:-1]) and row[-1] != 0 for row in R):
                continue  # inconsistent
            if len(piv) < k:
                continue  # affinely dependent subset: covered by a smaller one
            lam = [Fr(0)] * k
            for row, c in zip(R, piv):
                lam[c] = row[-1]
            if all(x >= 0 for x in lam):
                return True
    return False


def rank_one_faces(M):
    """For a conservative rank-one network M: list of (sigma, Z, type) over all faces; type in IN/OUT/ESC/DOM."""
    G = M.Gamma
    cols = transpose(G)
    g = None
    for c in cols:
        if any(c):
            g = c
            break
    from math import gcd
    gg = 0
    for x in g:
        gg = gcd(gg, abs(x))
    g = [x // gg for x in g]
    c = []
    for col in cols:
        k = next((i for i in range(M.n) if g[i] != 0), None)
        cj = Fr(col[k], g[k])
        assert all(Fr(col[i]) == cj * g[i] for i in range(M.n))
        c.append(cj)
    Nm = [i for i in range(M.n) if g[i] < 0]
    Np = [i for i in range(M.n) if g[i] > 0]
    out = []
    for sigma, side in ((1, Nm), (-1, Np)):
        for k in range(1, len(side) + 1):
            for Z in combinations(side, k):
                act = [j for j in range(M.r) if c[j] != 0]
                order = {j: sum(M.reactions[j][0][i] for i in Z) for j in act}
                mu = min(order.values())
                J = [j for j in act if order[j] == mu]
                inward = [j for j in J if sigma * c[j] < 0]
                outward = [j for j in J if sigma * c[j] > 0]
                comp = [i for i in range(M.n) if i not in Z]
                beta = {j: tuple(M.reactions[j][0][i] for i in comp) for j in J}
                if not outward:
                    t = "IN"
                elif not inward:
                    t = "OUT"
                else:
                    pts = [beta[q] for q in inward]
                    t = "DOM" if all(in_convex_hull(beta[o], pts) for o in outward) else "ESC"
                out.append((sigma, Z, t))
    return out


def rank_one_verdict(M):
    ft = rank_one_faces(M)
    types = [t for _, _, t in ft]
    if all(t == "IN" for t in types):
        return "ALL", ft
    if any(t in ("OUT", "ESC") for t in types):
        return "NONE", ft
    return "DOM", ft


def drain_species(N):
    """Prop. (drain species): species i and d0 with a reaction of source exactly d0*X_i consuming X_i, and every
    reaction producing X_i has X_i-coefficient > d0 in its source.  Returns list of (i, d0)."""
    out = []
    for i in range(N.n):
        for a, b in N.reactions:
            if a[i] >= 1 and all(a[l] == 0 for l in range(N.n) if l != i) and b[i] < a[i]:
                d0 = a[i]
                if all(a2[i] > d0 for a2, b2 in N.reactions if b2[i] > a2[i]):
                    out.append((i, d0))
    return sorted(set(out))


def is_drain_set(N, Z, w, d):
    """(D1) each i in Z has a reaction with source exactly d_i X_i that strictly decreases V = <w, x_Z>;
    (D2) every reaction strictly increasing V has a source y with sum_{i in Z} y_i/d_i > 1."""
    def dV(a, b):
        return sum(Fr(w[i]) * (b[i] - a[i]) for i in Z)
    for i in Z:
        if not any(a[i] == d[i] and all(a[l] == 0 for l in range(N.n) if l != i) and dV(a, b) < 0
                   for a, b in N.reactions):
            return False
    for a, b in N.reactions:
        if dV(a, b) > 0 and not sum(Fr(a[i], d[i]) for i in Z) > 1:
            return False
    return True


def restricted_jacobian(f, xs, point, basis_cols, rows_idx):
    """Jacobian J of f at `point`, restricted to S = span(basis_cols): J B = B J_S, solved on rows rows_idx where
    the corresponding square submatrix of B is invertible.  Entries may lie in any field."""
    n = len(xs)
    J = [[f[i].diff(xs_name(xs[j])).value(point) for j in range(n)] for i in range(n)]
    B = basis_cols  # n x s
    JB = matmul(J, B)
    s = len(B[0])
    Bs = [[Fr(B[i][j]) for j in range(s)] for i in rows_idx]
    JBs = [JB[i] for i in rows_idx]
    if s == 2:
        dB = Bs[0][0] * Bs[1][1] - Bs[0][1] * Bs[1][0]
        inv = [[Bs[1][1] / dB, -Bs[0][1] / dB], [-Bs[1][0] / dB, Bs[0][0] / dB]]
        JS = matmul(inv, JBs)
    else:
        raise NotImplementedError
    # check J B = B J_S on all rows
    ok = all(JB[i][j] - sum((B[i][k] * JS[k][j] for k in range(s)), Fr(0)) == 0 or
             (hasattr(JB[i][j] - sum((B[i][k] * JS[k][j] for k in range(s)), Fr(0)), "sign") and
              (JB[i][j] - sum((B[i][k] * JS[k][j] for k in range(s)), Fr(0))).sign() == 0)
             for i in range(n) for j in range(s))
    return JS, ok


def xs_name(p):
    """Name of a polynomial variable."""
    (m, c), = p.t.items()
    (v, e), = m
    return v


def sgn(x):
    if isinstance(x, QS):
        return x.sign()
    return (x > 0) - (x < 0)


# =============================================================================================================
# 5. The networks of the paper
# =============================================================================================================
X2, X3, X4, X5 = Xs(2), Xs(3), Xs(4), Xs(5)
NETS = {
    "R": net("X1 -> X2; X2 -> X1", X2),
    "T": net("X1 -> X2; X2 -> X3; X3 -> X1", X3),
    "W": net("X1 -> X2; X1 -> X3; X2 -> X1; X3 -> X1", X3),
    "U": net("X4 -> X2+X3; X3+X4 -> X1; X1+X2 -> 2X4", X4),
    "D": net("X1 -> X2; X1 -> X3; X2+X3 -> 2X1", X3),
    "F": net("S0+E -> ES0; ES0 -> S0+E; ES0 -> S1+E; S1+F -> FS1; FS1 -> S1+F; FS1 -> S0+F",
             ["S0", "S1", "E", "F", "ES0", "FS1"]),
    "N3": net("2X1 -> X1+X2; X1+X2 -> 2X1; X1+X2 -> 2X2", X2),
    "N412": net("2X1 -> X1+X2; X1+X2 -> 2X1", X2),
    "M": net("3A -> 2A+B; 2A+B -> 3A", ["A", "B"]),
}
NETS["Wstar"] = NETS["W"].dual()


def B_k(k, dup=False):
    sp_ = ["X"] + [f"Y{a}" for a in range(1, k + 1)] + ["Z"]
    spec = "; ".join(f"X -> Y{a}" for a in range(1, k + 1)) + "; " + \
        "; ".join(f"Y{a} -> Z" for a in range(1, k + 1)) + "; Z -> X" + ("; Z -> X" if dup else "")
    return net(spec, sp_)


# =============================================================================================================
# 6. Section 2: the dual and the Petri dictionary
# =============================================================================================================
def sec2():
    section("S2 The dual network and the Petri dictionary")
    ex = ["R", "T", "W", "U", "D", "F", "N3", "N412", "M"]
    ok_t = ok_inv = ok_h = ok_h3 = ok_pre = ok_sip = ok_trap = ok_cl = ok_ic = ok_mips = ok_rank = True
    for name in ex:
        N = NETS[name]
        Nd = N.dual()
        ok_t &= Nd.Gamma == transpose(N.Gamma)
        ok_inv &= Nd.dual().reactions == N.reactions
        ok_rank &= Nd.rank() == N.rank() and Nd.n == N.r and Nd.r == N.n
        ok_h &= (N.is_conservative() == Nd.is_consistent()) and (N.is_consistent() == Nd.is_conservative())
        ok_h3 &= Nd.rows_distinct() == N.columns_distinct()
        # reverse-dual with weights: Pre(N*) = Pre(N)^T, Post(N*) = Post(N)^T, i.e. every arc reversed
        ok_pre &= Nd.Ye == transpose(N.Ye) and Nd.Yp == transpose(N.Yp)
        # siphons of N = T-traps of N*, traps of N = T-siphons of N*
        for k in range(1, N.n + 1):
            for S in combinations(range(N.n), k):
                ok_sip &= N.is_siphon(S) == is_T_trap(Nd, set(S))
                prod_d = {j for i in S for j in range(Nd.n) if Nd.reactions[i][1][j]}
                cons_d = {j for i in S for j in range(Nd.n) if Nd.reactions[i][0][j]}
                is_trap = all(any(N.reactions[j][1][i] for i in S) for j in range(N.r)
                              if any(N.reactions[j][0][i] for i in S))
                ok_trap &= is_trap == (cons_d <= prod_d)
        # clusters are exchanged
        cN = {(s, t) for s, t in clusters(N)}
        cD = {(t, s) for s, t in clusters(Nd)}
        ok_cl &= cN == cD
        # internal cycles of N* = supports of extreme conservation laws of N
        if N.is_conservative():
            ok_ic &= set(Nd.internal_cycles()) == set(N.mips_semiflow())
        # MIPS(N*) = full-coupling classes of N (under (H2))
        if N.is_consistent():
            nv = N.n_vectors()
            key = {}
            for j in range(N.r):
                key.setdefault(crn.direction_key(nv[j]), []).append(j)
            ok_mips &= Nd.mips_ray() == crn.normalize_partition(key.values())
    tag = "on R, T, W, U, D, F, N3, N4.12, M"
    claim("S2 dual: Gamma(N*) = Gamma^T " + tag, ok_t)
    claim("S2 dual: N** = N " + tag, ok_inv)
    claim("S2 dual: N* has r species, n reactions and the same rank s " + tag, ok_rank)
    claim("S2 dual: (H1) for N <=> (H2) for N*, and (H2) for N <=> (H1) for N* " + tag, ok_h)
    claim("S2 dual: (H3) for N* <=> columns of Gamma pairwise distinct " + tag, ok_h3)
    claim("S2 Prop reverse-dual: Pre(N*) = Pre(N)^T and Post(N*) = Post(N)^T (every arc reversed, weight kept) "
          + tag, ok_pre)
    claim("S2 dictionary: siphons of N = T-traps of N* (all subsets) " + tag, ok_sip)
    claim("S2 dictionary: traps of N = T-siphons of N* (all subsets) " + tag, ok_trap)
    claim("S2 dictionary: clusters of N* = clusters of N with species and reactions exchanged " + tag, ok_cl)
    claim("S2 dictionary: internal cycles of N* = supports of extreme conservation laws of N " + tag, ok_ic)
    claim("S2 dictionary: MIPS(N*) = classes of fully coupled reactions of N " + tag, ok_mips)
    # EC is not self-dual; EC with balanced joins is (on the examples)
    E = net("2X1+X2 -> X3; X3 -> 2X1+X2", X3)
    claim("S2 EC not self-dual: 2X1+X2 <-> X3 is EC, its dual 2K1 -> 2K2, K1 -> K2, K2 -> K1 is not EC",
          is_EC(E) and not is_EC(E.dual()) and
          same_reactions(E.dual(), "2k1 -> 2k2; k1 -> k2; k2 -> k1", Ks(2)))
    okb = all((not (is_EC(NETS[x]) and balanced_joins(NETS[x]))) or
              (is_EC(NETS[x].dual()) and balanced_joins(NETS[x].dual())) for x in ex)
    claim("S2 EC with balanced joins is preserved by the dual " + tag, okb)
    # the figure example: 3A <-> 2A+B and its dual
    M = NETS["M"]
    claim("S2 figure 1: M = 3A <-> 2A+B has Ye = [[3,2],[0,1]], Yp = [[2,3],[1,0]]; M* = 3U+2V -> 2U+3V, V -> U",
          M.Ye == [[3, 2], [0, 1]] and M.Yp == [[2, 3], [1, 0]] and
          same_reactions(M.dual(), "3U+2V -> 2U+3V; V -> U", ["U", "V"]))


# =============================================================================================================
# 7. Section 3: preclusters, MIPS and siphons
# =============================================================================================================
def sec3_archive():
    section("S3.1 Precise forms (condensed from the archived report)")
    R, T, W, U, D, F, N3, N412 = (NETS[x] for x in ("R", "T", "W", "U", "D", "F", "N3", "N412"))
    Wd = NETS["Wstar"]
    allH = all(H_all(N) for N in (R, T, W, U, D, F, N3, N412, Wd))
    claim("S3 table: R, T, W, W*, U, D, F, N3, N4.12 satisfy (H1)-(H3)", allH)
    strong = all(N.rows_pairwise_independent() and N.columns_pairwise_independent() for N in (T, U, D)) and \
        W.rows_pairwise_independent()
    claim("S3 table: T, U, D have rows and columns of Gamma pairwise linearly independent; so do the rows of W",
          strong)
    one = lambda r: [frozenset(range(r))]  # noqa: E731
    sing = lambda r: [frozenset([j]) for j in range(r)]  # noqa: E731
    Td = T.dual()
    claim("S3 T: P(T) singletons, MIPS(T*) one block, P(T*) singletons, MIPS(T) one block",
          T.preclusters() == sing(3) and Td.mips_ray() == one(3) and Td.preclusters() == sing(3)
          and T.mips_ray() == one(3))
    claim("S3 T: the only siphon of T and of T* is the whole species set; T* is again a 3-cycle",
          T.siphons() == [frozenset(range(3))] and Td.siphons() == [frozenset(range(3))] and
          same_up_to_relabel(Td, T))
    claim("S3 T: the internal cycle {R1,R2,R3} lies in no precluster (Remark 23 fails)",
          T.internal_cycles() == [frozenset(range(3))] and
          not any(frozenset(range(3)) <= b for b in T.preclusters()))
    claim("S3 T: weakly reversible, deficiency 0", is_wr(T) and T.deficiency() == 0)
    claim("S3 R: P(R) singletons, MIPS(R*) one block, R* isomorphic to R",
          R.preclusters() == sing(2) and R.dual().mips_ray() == one(2) and same_up_to_relabel(R.dual(), R))
    claim("S3 W: P(W) one block, MIPS(W*) = {R1,R3},{R2,R4}; doubling pair {R1,R2} not fully coupled; "
          "W weakly reversible with deficiency 0",
          W.preclusters() == one(4) and W.dual().mips_ray() == [frozenset({0, 2}), frozenset({1, 3})]
          and not fully_coupled(W, 0, 1) and is_wr(W) and W.deficiency() == 0)
    claim("S3 W*: = K1+K2 -> K3+K4, K3 -> K1, K4 -> K2, not weakly reversible",
          same_reactions(Wd, "k1+k2 -> k3+k4; k3 -> k1; k4 -> k2", Ks(4)) and not is_wr(Wd))
    claim("S3 U: MIPS(U) singletons, P(U*) one block, twins X1,X2 not conservation-coupled",
          U.mips_ray() == sing(4) and U.dual().preclusters() == one(4) and (0, 1) in twins(U)
          and not conservation_coupled(U, 0, 1))
    claim("S3 D: P(D) = MIPS(D*) = one block, P(D*) = MIPS(D) = one block, D* isomorphic to D",
          D.preclusters() == one(3) and D.dual().mips_ray() == one(3) and D.dual().preclusters() == one(3)
          and D.mips_ray() == one(3) and same_up_to_relabel(D.dual(), D))
    claim("S3 D: minimal siphons {X1,X2}, {X1,X3}; the single precluster of D* lies in none (Q3 fails)",
          D.minimal_siphons() == [frozenset({0, 1}), frozenset({0, 2})])
    claim("S3 F: P(F) one block, MIPS(F*) = {R3,R6} + singletons, P(F*) one block, MIPS(F) = {S0,S1} + singletons",
          F.preclusters() == one(6) and
          set(F.dual().mips_ray()) == {frozenset({2, 5}), frozenset({0}), frozenset({1}), frozenset({3}),
                                       frozenset({4})}
          and F.dual().preclusters() == one(6) and
          set(F.mips_ray()) == {frozenset({0, 1}), frozenset({2}), frozenset({3}), frozenset({4}), frozenset({5})})
    claim("S3 F: minimal siphons {E,ES0}, {F,FS1}, {S0,S1,ES0,FS1}",
          set(F.minimal_siphons()) == {frozenset({2, 4}), frozenset({3, 5}), frozenset({0, 1, 4, 5})})
    claim("S3 N3: P(N3) one block, MIPS(N3*) singletons (no two reactions fully coupled)",
          N3.preclusters() == one(3) and N3.dual().mips_ray() == sing(3))
    claim("S3 N4.12: P = {R1},{R2}, MIPS(N4.12*) one block; internal cycle {R1,R2} in no precluster",
          N412.preclusters() == sing(2) and N412.dual().mips_ray() == one(2) and
          N412.internal_cycles() == [frozenset({0, 1})])
    agree = all(N.preclusters("onestep") == N.preclusters("thesis") == N.preclusters("closure")
                for N in (R, T, W, U, D, F, N3, N412, Wd) + tuple(N.dual() for N in (R, T, W, U, D, F, N3, N412)))
    claim("S3 the three precluster constructions agree on all these networks and their duals", agree)
    wffc = all(ordinary(N) and WF_EC(N) and ordinary(N.dual()) and WF_EC(N.dual()) for N in (R, T, W, Wd, F))
    claim("S3 R, T, W, W*, F are well-formed free-choice networks with well-formed free-choice duals", wffc)
    claim("S3 W: P(W*) singletons (no twins), as for every WR network with (H1)-(H3)",
          W.dual().preclusters() == sing(3) and twins(W) == [])


def one_block_certificate(N):
    """CF(N) is one block if some (m-1)-minor of A_V, V = E lambda, is a nonzero polynomial in lambda with all
    coefficients of one sign (then rank A_V = m - 1 for every lambda > 0, so ker A_V = R 1 for every V in C°)."""
    rays = N.flux_rays()
    lam = V(*[f"l{i + 1}" for i in range(len(rays))])
    Vsym = [sum((lam[i] * rays[i][j] for i in range(len(rays))), Poly()) for j in range(N.r)]
    G = N.Gamma
    cls = source_classes(N)
    m = len(cls)
    A = [[sum((Vsym[t] * G[i][t] for t in c), Poly()) for c in cls] for i in range(N.n)]
    if m == 1:
        return True
    for rows in combinations(range(N.n), m - 1):
        for cols in combinations(range(m), m - 1):
            d = Poly.lift(det([[A[i][j] for j in cols] for i in rows]))
            cs = d.coeffs()
            if cs and (all(c > 0 for c in cs) or all(c < 0 for c in cs)):
                return True
    return False


def sec3_phat():
    section("S3.2 Definition 19 is not invariant; corrected preclusters")
    B2, B3, B4 = B_k(2), B_k(3), B_k(4)
    one = lambda r: [frozenset(range(r))]  # noqa: E731
    claim("S3 B2, B3, B4 satisfy (H1)-(H3) and are well-formed free-choice",
          all(H_all(B) and WF_EC(B) and ordinary(B) for B in (B2, B3, B4)))
    claim("S3 P(B2) is one block", B2.preclusters() == one(5))
    claim("S3 P(B3) = {all but Z->X}, {Z->X}; P(B4) has two blocks",
          B3.preclusters() == [frozenset(range(6)), frozenset({6})] and len(B4.preclusters()) == 2)
    claim("S3 class-span closure P_cs is one block for B2, B3 and B4",
          P_cs(B2) == one(5) and P_cs(B3) == one(7) and P_cs(B4) == one(9))
    B3s = net("X -> Y1; X -> W; W -> Y2; W -> Y3; Y1 -> Z; Y2 -> Z; Y3 -> Z; Z -> X",
              ["X", "Y1", "Y2", "Y3", "Z", "W"])
    claim("S3 splitting the choice of B3 through a new species W gives one precluster (Definition 19 is not "
          "invariant under splitting)", H_all(B3s) and B3s.preclusters() == one(8))
    N412 = NETS["N412"]
    N412d = net("2X1 -> X1+X2; 2X1 -> X1+X2; X1+X2 -> 2X1", X2)
    claim("S3 duplicating R1 of N4.12 merges its two preclusters (Definition 19 is not invariant under "
          "duplication)", N412.preclusters() == [frozenset({0}), frozenset({1})] and N412d.preclusters() == one(3))
    B2d = B_k(2, dup=True)
    claim("S3 B2dup (Z -> X duplicated): P_cs = pi = {R1..R4}, {R5,R6}, while P_cs(B2) = pi(B2) = one block",
          P_cs(B2d) == [frozenset({0, 1, 2, 3}), frozenset({4, 5})] == flat_closure(B2d) and
          P_cs(B2) == flat_closure(B2) == one(5))
    B2ds = net("X -> W; W -> Y1; W -> Y2; Y1 -> Z; Y2 -> Z; Z -> X; Z -> X", ["X", "Y1", "Y2", "Z", "W"])
    claim("S3 splitting the choice of B2dup through W: P_cs and pi become one block",
          P_cs(B2ds) == one(7) == flat_closure(B2ds))
    claim("S3 CF(B2dup) is one block: certificate minor sign-definite on the open orthant of ray weights",
          one_block_certificate(B2d))
    R = NETS["R"]
    claim("S3 hierarchy: B3 has P (two blocks) strictly finer than P_cs (one block)",
          len(B3.preclusters()) == 2 and P_cs(B3) == one(7))
    claim("S3 hierarchy: R has P_cs = singletons strictly finer than pi = MIPS(R*) = one block",
          P_cs(R) == [frozenset({0}), frozenset({1})] and flat_closure(R) == one(2) == R.dual().mips_ray())
    K3 = net("X -> Y; X -> Z; Y -> X; Y -> Z; Z -> X; Z -> Y", ["X", "Y", "Z"])
    claim("S3 hierarchy: K3 (complete digraph on X, Y, Z) is EC, connected, (H1), (H2), weakly reversible, "
          "deficiency 0; pi(K3) = the three source classes, while CF(K3) is one block",
          is_EC(K3) and connected(K3) and K3.is_conservative() and K3.is_consistent() and is_wr(K3)
          and K3.deficiency() == 0 and flat_closure(K3) == source_classes(K3) and len(source_classes(K3)) == 3
          and one_block_certificate(K3))
    XY4 = net("X+Y -> 2X; X+Y -> 2Y; X+2Y -> 2X+Y; X+2Y -> 3Y", ["X", "Y"])
    rays = XY4.flux_rays()
    Vgen = [sum((Fr(w) * v[j] for w, v in zip((1, 2, 3, 5), rays)), Fr(0)) for j in range(4)]
    claim("S3 hierarchy: XY4 satisfies (H1)-(H3); kappa_V is one block at a generic V, two blocks at "
          "V = (1,1,2,2), pi has two blocks, so CF = pi = {R1,R2},{R3,R4} is strictly finer than kappa_V(generic)",
          H_all(XY4) and kappa_V(XY4, Vgen) == [frozenset(range(4))] and
          kappa_V(XY4, [1, 1, 2, 2]) == [frozenset({0, 1}), frozenset({2, 3})] == flat_closure(XY4) and
          all(x > 0 for x in Vgen) and all(x == 0 for row in matmul(XY4.Gamma, [[v] for v in Vgen]) for x in row),
          "V_generic = %s" % [str(x) for x in Vgen])
    Z4 = net("2X2 -> 2X1; 2X1 -> 2X2; X1+X2 -> 2X1; X1+X2 -> 2X2", X2)
    rays = Z4.flux_rays()
    Vg = [sum((Fr(w) * v[j] for w, v in zip((1, 2, 3, 5), rays)), Fr(0)) for j in range(4)]
    claim("S3 Z4: kappa_V = {R1},{R2},{R3,R4} at a generic V but {R1,R2},{R3,R4} at V = (1,1,2,2); "
          "pi = {R1},{R2},{R3,R4}, hence CF = {R1},{R2},{R3,R4}",
          H_all(Z4) and kappa_V(Z4, Vg) == [frozenset({0}), frozenset({1}), frozenset({2, 3})] ==
          flat_closure(Z4) and kappa_V(Z4, [1, 1, 2, 2]) == [frozenset({0, 1}), frozenset({2, 3})])
    claim("S3 Remark 23 for P-hat fails for XY4 and Z4 (an internal cycle meets two blocks of CF), and holds "
          "for T (CF(T) is one block)",
          not all(any(c <= b for b in flat_closure(XY4)) for c in XY4.internal_cycles()) and
          not all(any(c <= b for b in flat_closure(Z4)) for c in Z4.internal_cycles()) and
          one_block_certificate(NETS["T"]))
    okEC = True
    for N in (NETS["Wstar"], NETS["F"], B3, K3):
        rays = N.flux_rays()
        V_ = [sum((Fr(i + 1) * v[j] for i, v in enumerate(rays)), Fr(0)) for j in range(N.r)]
        A, cls = A_matrix(N, V_)
        okEC &= crn.rank(A) == len(cls) - 1
    claim("S3 Theorem (EC): ker A_V = R 1 at a sample flux V for W*, F, B3, K3", okEC)
    claim("S3 WR deficiency 0: kappa_V is the linkage-class partition (one block) for W and T at the flux 1",
          kappa_V(NETS["W"], [1, 1, 1, 1]) == [frozenset(range(4))] and
          kappa_V(NETS["T"], [1, 1, 1]) == [frozenset(range(3))])
    U, F, W, N3, Wd = NETS["U"], NETS["F"], NETS["W"], NETS["N3"], NETS["Wstar"]
    tw_ok = all(any(not conservation_coupled(N, i, j) for i, j in twins(N)) for N in (U, F, Wd))
    dp_ok = all(any(not fully_coupled(N, i, j) for i, j in N.doubling_pairs()) for N in (W, N3, F))
    claim("S3 P-hat, species side: U, F, W* have a twin pair that is not conservation-coupled", tw_ok)
    claim("S3 P-hat, reaction side: W, N3, F have a doubling pair that is not fully coupled", dp_ok)


def sec3_siphons():
    section("S3.3 Siphons: transversality, Theorem B, MIPS and siphons, K8")
    D, N412 = NETS["D"], NETS["N412"]

    def thmB(N):
        Nd = N.dual()
        lhs = all(any(b <= S for S in N.minimal_siphons()) for b in Nd.preclusters())
        rhs = twins(N) == [] and all(any(i in S for S in N.minimal_siphons()) for i in range(N.n))
        return lhs == rhs
    ex = [NETS[x] for x in ("R", "T", "W", "U", "D", "F", "N3", "N412", "M", "Wstar")]
    claim("S3 Theorem B (twins with zero rows included) holds on both sides of R, T, W, W*, U, D, F, N3, N4.12, M",
          all(thmB(N) and thmB(N.dual()) for N in ex))
    Dg = net("X1 -> X2; X2 -> X1", X4)
    Dgd = Dg.dual()
    claim("S3 degenerate case: X1 <-> X2 with X3, X4 in no reaction: P(N*) = {X1},{X2},{X3,X4}, minimal siphons "
          "{X3},{X4},{X1,X2}; {X3,X4} is a twin pair of zero rows of Ye",
          Dgd.preclusters() == [frozenset({0}), frozenset({1}), frozenset({2, 3})] and
          set(Dg.minimal_siphons()) == {frozenset({2}), frozenset({3}), frozenset({0, 1})} and
          (2, 3) in twins(Dg) and thmB(Dg))

    def transversal(N):
        cons = {i: frozenset(j for j in range(N.r) if N.reactions[j][0][i]) for i in range(N.n)}
        return all(not (i in S and j in S) for S in N.minimal_siphons() for i, j in combinations(range(N.n), 2)
                   if cons[i] == cons[j])
    claim("S3 Theorem A (transversality) on both sides of the same networks",
          all(transversal(N) and transversal(N.dual()) for N in ex))
    claim("S3 Prop G: M(D) is one block while D has the two minimal siphons {X1,X2}, {X1,X3} (readings (i)-(iii) "
          "fail); N4.12 has the law X1+X2 but the only minimal siphon {X1} (reading (iv) fails)",
          D.mips_ray() == [frozenset(range(3))] and len(D.minimal_siphons()) == 2 and
          N412.conservation_rays() == [(1, 1)] and N412.minimal_siphons() == [frozenset({0})])

    def propH(N):
        if not N.is_conservative():
            return True
        M_ = N.mips_ray()
        return all(all(m <= b or not (m & b) for m in M_) for b in N.dual().preclusters() if len(b) >= 2)
    claim("S3 Prop H: every precluster of N* with >= 2 elements is a union of MIPS of N, both sides of the examples",
          all(propH(N) and propH(N.dual()) for N in ex))
    M3 = net("X4 -> 2X3; X4 -> X2+X3; 2X3 -> 2X2; X2+X3 -> 2X1; 2X2 -> X4; 2X1 -> 2X3", X4)
    M4 = net("X4 -> X2+X3; 2X3 -> X4; 2X3 -> 2X2; X2+X3 -> 2X1; 2X2 -> X4; 2X1 -> 2X3", X4)
    hyp = all(H_all(N) and H_all(N.dual()) and is_wr(N) and pure(N) and N.rows_pairwise_independent() and
              N.columns_pairwise_independent() for N in (M3, M4))
    claim("S3 K8: M3, M4 weakly reversible, pure, (H1)-(H3) on both sides, rows and columns of Gamma pairwise "
          "linearly independent", hyp)

    def key(N):
        dcl = [c for c in source_classes(N) if len(c) > 1]
        return dcl, minimal_T_traps(N)
    k3, k4 = key(M3), key(M4)
    claim("S3 K8: doubling classes {R1,R2} (M3), {R2,R3} (M4); minimal T-traps {R1,R3,R5},{R2,R3,R5},{R4,R6} (M3), "
          "{R1,R2,R5},{R1,R3,R5},{R4,R6} (M4); P(M3) = {R1,R2,R5},{R3},{R4},{R6}; P(M4) = {R2,R3} + singletons",
          k3[0] == [frozenset({0, 1})] and k4[0] == [frozenset({1, 2})] and
          set(k3[1]) == {frozenset({0, 2, 4}), frozenset({1, 2, 4}), frozenset({3, 5})} and
          set(k4[1]) == {frozenset({0, 1, 4}), frozenset({0, 2, 4}), frozenset({3, 5})} and
          set(M3.preclusters()) == {frozenset({0, 1, 4}), frozenset({2}), frozenset({3}), frozenset({5})} and
          set(M4.preclusters()) == {frozenset({1, 2}), frozenset({0}), frozenset({3}), frozenset({4}),
                                    frozenset({5})})
    match = differ = 0
    P3, P4 = set(M3.preclusters()), set(M4.preclusters())
    for p in permutations(range(6)):
        img = lambda S: frozenset(p[j] for j in S)  # noqa: E731
        if {img(c) for c in k3[0]} == set(k4[0]) and {img(U_) for U_ in k3[1]} == set(k4[1]):
            match += 1
            if {img(b) for b in P3} != P4:
                differ += 1
    claim("S3 K8: exactly 8 relabellings of the reactions match the doubling classes and minimal T-traps of M3 "
          "with those of M4, and the preclusters differ under all 8", match == 8 and differ == 8,
          f"matching {match}, differing {differ}")

    def thmD(N):
        laws = {frozenset(i for i in range(N.n) if z[i]) for z in N.conservation_rays()}
        ind = all(all(x in (0, 1) for x in z) for z in N.conservation_rays())
        fl = {frozenset(j for j in range(N.r) if v[j]) for v in N.flux_rays()}
        ind2 = all(all(x in (0, 1) for x in v) for v in N.flux_rays())
        return ind and ind2 and laws == set(N.minimal_siphons()) and fl == set(minimal_T_traps(N)) and \
            critical_siphons(N) == [] and critical_siphons(N.dual()) == []
    claim("S3 Theorem D on B3 and F: extreme conservation laws = 0/1 indicators of the minimal siphons, extreme "
          "fluxes = indicators of the minimal T-traps, no critical siphon on either side",
          thmD(B_k(3)) and thmD(NETS["F"]))


# =============================================================================================================
# 8. Section 4: what the duality transfers
# =============================================================================================================
def mcounts(N):
    """(number of distinct sources, number of Petri clusters, rank)."""
    return len(sources(N)), len(clusters(N)), N.rank()


def no_crit_both(N):
    return critical_siphons(N) == [] and critical_siphons(N.dual()) == []


def sym_field(N, kname="k", xnames=None):
    ks = V(*[f"{kname}{j + 1}" for j in range(N.r)])
    xn = xnames or [f"x{i + 1}" for i in range(N.n)]
    xs = V(*xn)
    return field(N, ks, xs), ks, xs, xn


def basis_of_S(N):
    """Columns of Gamma forming a basis of S = im Gamma (greedy)."""
    G = N.Gamma
    cols = transpose(G)
    chosen = []
    for c in cols:
        if crn.rank(chosen + [c]) > len(chosen):
            chosen.append(c)
    return transpose(chosen)


def J_restricted(N, rates, point):
    """Exact Jacobian of the mass-action field restricted to S at `point` (dict name -> value)."""
    xn = [f"x{i + 1}" for i in range(N.n)]
    f = field(N, rates, V(*xn))
    B = basis_of_S(N)
    s = len(B[0])
    J = [[f[i].diff(xn[j]).value(point) for j in range(N.n)] for i in range(N.n)]
    JB = matmul(J, B)
    for rows in combinations(range(N.n), s):
        sub = [[Fr(B[i][j]) for j in range(s)] for i in rows]
        if det(sub) != 0:
            break
    if s == 2:
        a, b_, c, d = sub[0][0], sub[0][1], sub[1][0], sub[1][1]
        dt = a * d - b_ * c
        inv = [[d / dt, -b_ / dt], [-c / dt, a / dt]]
    elif s == 1:
        inv = [[1 / sub[0][0]]]
    else:
        raise NotImplementedError
    JS = matmul(inv, [JB[i] for i in rows])
    resid = [[JB[i][j] - sum((B[i][k] * JS[k][j] for k in range(s)), Fr(0)) for j in range(s)]
             for i in range(N.n)]
    ok = all(sgn(x) == 0 for row in resid for x in row)
    return JS, ok, f


def is_steady(f, point):
    return all(sgn(fi.value(point)) == 0 for fi in f)


def sec4():
    section("S4 Transfers: cluster rank equation, well-formed equal conflict, proportional conflict")
    Wd, F, B3 = NETS["Wstar"], NETS["F"], B_k(3)
    claim("S4 W*, F, B3 are WF-EC (EC, connected, (H1), (H2), s = m - 1); W* and F are not weakly reversible "
          "(hence complex balanced for no rate vector)",
          WF_EC(Wd) and WF_EC(F) and WF_EC(B3) and not is_wr(Wd) and not is_wr(F))
    claim("S4 B2, B3, B4 are weakly reversible with deficiency 0 (hence complex balanced for every rate vector)",
          all(is_wr(B_k(k)) and B_k(k).deficiency() == 0 for k in (2, 3, 4)))
    Wj = net("A -> C; 2C -> 2A", ["A", "C"])
    claim("S4 W' = A -> C, 2C -> 2A is WF-EC with (H1)-(H3) and ker Gamma = R(2,1) (a circuit with unequal entries: "
          "no translation of deficiency zero)", WF_EC(Wj) and H_all(Wj) and Wj.flux_rays() == [(2, 1)])
    # N0 and its dual
    N0 = net("2A+B -> 3A; 2A+B -> A+2B", ["A", "B"])
    N0d = N0.dual()
    f, ks, xs, xn = sym_field(N0, xnames=["a", "b"])
    a, b = xs
    claim("S4 N0 = 2A+B -> 3A, 2A+B -> A+2B: EC, unbalanced, connected and (H1)-(H3) on both sides, s = m = 1; "
          "a' = (k1 - k2) a^2 b",
          is_EC(N0) and not balanced_joins(N0) and connected(N0) and connected(N0d) and H_all(N0) and H_all(N0d)
          and N0.rank() == 1 and len(sources(N0)) == 1 and f[0] == (ks[0] - ks[1]) * a ** 2 * b)
    fd, kd, xd, _ = sym_field(N0d, "c", ["K1", "K2"])
    K1, K2 = xd
    claim("S4 N0* = 2K1+2K2 -> 3K1+K2, K1+K2 -> 2K2; K1' = cA (K1 K2)^2 - cB K1 K2, so (1, cB/cA) is a positive "
          "steady state for every rate vector",
          same_reactions(N0d, "2k1+2k2 -> 3k1+k2; k1+k2 -> 2k2", Ks(2)) and
          fd[0] == kd[0] * (K1 * K2) ** 2 - kd[1] * K1 * K2 and fd[1] == -fd[0])
    N1 = net("2A+B -> C; C -> 3A; C -> 3B", ["A", "B", "C"])
    N1d = N1.dual()
    claim("S4 N1 = 2A+B -> C, C -> 3A, C -> 3B: EC, (H1)-(H3), s = m = 2, ker Gamma = R(3,2,1) (so a positive "
          "steady state needs k2 = 2 k3); N1* = 2K1 -> 3K2, K1 -> 3K3, K2+K3 -> K1",
          is_EC(N1) and H_all(N1) and N1.rank() == 2 and len(sources(N1)) == 2 and N1.flux_rays() == [(3, 2, 1)]
          and same_reactions(N1d, "2k1 -> 3k2; k1 -> 3k3; k2+k3 -> k1", Ks(3)))
    fd, kd, xd, xdn = sym_field(N1d, "c", ["K1", "K2", "K3"])
    cA, cB, cC = kd
    # steady state (cB/cA, 3 cA K1^2/cC, 1): check after clearing denominators via substitution with values
    ok1 = True
    for vals in ((1, 1, 1), (2, 3, 5), (Fr(1, 3), 7, Fr(2, 5))):
        A_, B_, C_ = (Fr(v) for v in vals)
        K1v = B_ / A_
        pt = {"c1": A_, "c2": B_, "c3": C_, "K1": K1v, "K2": 3 * A_ * K1v ** 2 / C_, "K3": Fr(1)}
        ok1 &= all(fi.value(pt) == 0 for fi in fd)
    # symbolic: with K1 = cB/cA, K2 = 3 cA K1^2/cC, K3 = 1, the field vanishes identically: multiply out
    claim("S4 N1*: (cB/cA, 3 cA K1^2/cC, 1) is a positive steady state (checked at three rational rate vectors)",
          ok1)
    # N2 and tor*
    N2 = net("2A+B -> C; C -> 2A+B", ["A", "B", "C"])
    N2d = N2.dual()
    fd, kd, xd, _ = sym_field(N2d, "c", ["K1", "K2"])
    K1, K2 = xd
    c1, c2, c3 = kd
    claim("S4 N2 = 2A+B <-> C is WF-EC with unbalanced joins; N2* = 2K1 -> 2K2, K1 -> K2, K2 -> K1; "
          "K1' = -2 c1 K1^2 - c2 K1 + c3 K2 (steady states K2 = (2 c1 K1^2 + c2 K1)/c3)",
          WF_EC(N2) and not balanced_joins(N2) and same_reactions(N2d, "2k1 -> 2k2; k1 -> k2; k2 -> k1", Ks(2)) and
          fd[0] == -2 * c1 * K1 ** 2 - c2 * K1 + c3 * K2)
    claim("S4 tor*: with h = log(2 c1 K1^2 + c2 K1), dh/dlogK1 = (4c1K1 + c2)/(2c1K1 + c2) and "
          "K1 d/dK1 of it = 2 c1 c2 K1/(2 c1 K1 + c2)^2 (nonzero: the steady-state curve is not log-affine)",
          K1 * (4 * c1 * (2 * c1 * K1 + c2) - 2 * c1 * (4 * c1 * K1 + c2)) == 2 * c1 * c2 * K1)
    claim("S4 N2* is block-homogeneous (blocks {K1}, {K2}), not EC, connected, (H1), (H2), s = #blocks - 1",
          is_BH(N2d) and not is_EC(N2d) and connected(N2d) and N2d.is_conservative() and N2d.is_consistent()
          and N2d.rank() == 1 and len(clusters(N2d)) == 2)
    NT = net("X1+X2 -> 2X2; X1+X2 -> 2X3+X4; 2X3+X4 -> 2X1", X4)
    NTd = NT.dual()
    laws = conservation_basis(NT)
    claim("S4 NT = X1+X2 -> 2X2, X1+X2 -> 2X3+X4, 2X3+X4 -> 2X1: EC, unbalanced, connected and (H1)-(H3) on both "
          "sides, s = m = 2; every conservation law has z1 = z2; the reactions of NT* dual to X1 and X2 have the "
          "same source K1+K2 (so NT* has no positive steady state when their rate constants differ)",
          is_EC(NT) and not balanced_joins(NT) and connected(NT) and connected(NTd) and H_all(NT) and H_all(NTd)
          and NT.rank() == 2 and len(sources(NT)) == 2 and all(z[0] == z[1] for z in laws) and
          NTd.reactions[0][0] == NTd.reactions[1][0] == (1, 1, 0))
    NTc = net("2A+B+C -> 2A+2B; 2A+B+C -> 2A+2C", ["A", "B", "C"])
    NTcd = NTc.dual()
    claim("S4 NTc = 2A+B+C -> 2A+2B, 2A+B+C -> 2A+2C: EC with one cluster, (H1)-(H3) on both sides; every law has "
          "z_B = z_C; the dual reactions of B and C share the source K1+K2",
          is_EC(NTc) and H_all(NTc) and H_all(NTcd) and all(z[1] == z[2] for z in conservation_basis(NTc)) and
          NTcd.reactions[1][0] == NTcd.reactions[2][0] == (1, 1))
    # U and U*
    U = NETS["U"]
    Ud = U.dual()
    f, ks, xs, xn = sym_field(U)
    k1, k2, k3 = ks
    okU = True
    for vals in ((1, 1, 1), (2, 3, 5), (Fr(1, 7), Fr(9, 2), 4)):
        pt = {"k1": Fr(vals[0]), "k2": Fr(vals[1]), "k3": Fr(vals[2])}
        pt.update({"x1": pt["k1"] / pt["k3"], "x2": Fr(1), "x3": pt["k1"] / pt["k2"], "x4": Fr(1)})
        okU &= is_steady(f, pt)
    claim("S4 U: (k1/k3, 1, k1/k2, 1) is a positive steady state (three rational rate vectors; all fluxes equal k1)",
          okU)
    fd, kd, xd, _ = sym_field(Ud, "c", ["K1", "K2", "K3"])
    claim("S4 U* = K3 -> K2, K3 -> K1, K2 -> K1, K1+K2 -> 2K3; f_K1 - f_K2 = (c2 - c1) K3 + 2 c3 K2 (> 0 on the "
          "orthant when c2 >= c1); U and U* satisfy (H1)-(H3) and are not EC",
          same_reactions(Ud, "k3 -> k2; k3 -> k1; k2 -> k1; k1+k2 -> 2k3", Ks(3)) and
          fd[0] - fd[1] == (kd[1] - kd[0]) * xd[2] + 2 * kd[2] * xd[1] and H_all(U) and H_all(Ud)
          and not is_EC(U) and not is_EC(Ud))
    # NP: proportional conflict
    NP = net("A+2B -> C; 2A+4B -> C+D; C -> A+2B; D -> A+2B", ["A", "B", "C", "D"])
    NPd = NP.dual()
    claim("S4 NP = A+2B -> C, 2A+4B -> C+D, C -> A+2B, D -> A+2B: PC, neither EC nor BH, connected, (H1), (H2), "
          "s = q - 1 = 2; NP* is PC and satisfies the same; no critical siphon on either side",
          is_PC(NP) and not is_EC(NP) and not is_BH(NP) and cluster_rank_class(NP) and NP.rank() == 2 and
          is_PC(NPd) and cluster_rank_class(NPd) and no_crit_both(NP))
    # NC3
    sp4 = ["A", "B", "C", "D"]
    NC3 = net("D -> B+C; C -> A; A+B -> B+C; 2A+B -> C+D", sp4)
    claim("S4 NC3 = D -> B+C, C -> A, A+B -> B+C, 2A+B -> C+D: connected, (H1)-(H3) on NC3 and NC3*, laws B+D and "
          "A+C+D, flux rays (0,1,1,0), (1,2,0,1), s = 2, clusters {D}, {C}, {A,B} (q = 3 = s + 1)",
          connected(NC3) and H_all(NC3) and H_all(NC3.dual()) and
          set(NC3.conservation_rays()) == {(0, 1, 0, 1), (1, 0, 1, 1)} and
          set(NC3.flux_rays()) == {(0, 1, 1, 0), (1, 2, 0, 1)} and NC3.rank() == 2 and len(clusters(NC3)) == 3)
    tefc = all(frozenset(i for i in range(4) if a_[i]) == frozenset(i for i in range(4) if c_[i])
               for (a_, _), (c_, _) in combinations(NC3.reactions, 2) if any(a_[i] and c_[i] for i in range(4)))
    claim("S4 NC3: TEFC (sources sharing a species have the same support), not PC, not EC; no critical siphon on "
          "either side (as the cluster-rank theorem says)",
          tefc and not is_PC(NC3) and not is_EC(NC3) and no_crit_both(NC3))
    kNC3 = [Fr(1, 6), Fr(1, 3), Fr(6), Fr(6)]
    pt = {"x1": Fr(1, 2), "x2": Fr(1, 3), "x3": Fr(6), "x4": Fr(3)}
    JS, okJ, f = J_restricted(NC3, kNC3, pt)
    dJ = JS[0][0] * JS[1][1] - JS[0][1] * JS[1][0]
    claim("S4 NC3, k = (1/6, 1/3, 6, 6): x = (1/2, 1/3, 6, 3) is a positive steady state with det(J|_S) = -7/9 "
          "(a saddle in its class B+D = 10/3, A+C+D = 19/2)",
          is_steady(f, pt) and okJ and dJ == Fr(-7, 9) and pt["x2"] + pt["x4"] == Fr(10, 3) and
          pt["x1"] + pt["x3"] + pt["x4"] == Fr(19, 2), f"J_S = {[[str(x) for x in r] for r in JS]}")
    # exact count: parametrisation and the cubic 72A^3 - 204A^2 + 122A - 19 = (2A - 1)(36A^2 - 84A + 19)
    A_, B_ = V("x1", "x2")
    D_ = 36 * A_ ** 2 * B_
    C_ = 18 * A_ * B_ + 36 * A_ ** 2 * B_
    sub = {"x3": C_, "x4": D_}
    fA, fB, fC, fD = f
    claim("S4 NC3: f_B = -f_D and f_A + f_C + f_D = 0, so steady states are f_A = f_D = 0; on D = 36 A^2 B, "
          "C = 18 A B + 36 A^2 B both vanish identically",
          (fB + fD).is_zero() and (fA + fC + fD).is_zero() and Poly.lift(fA.subs(sub)).is_zero() and
          Poly.lift(fD.subs(sub)).is_zero())
    cub = 72 * A_ ** 3 - 204 * A_ ** 2 + 122 * A_ - 19
    lhs = 2 * (A_ * (1 + 36 * A_ ** 2) + Fr(10, 3) * (18 * A_ + 72 * A_ ** 2) - Fr(19, 2) * (1 + 36 * A_ ** 2))
    claim("S4 NC3: with B = (10/3)/(1 + 36 A^2) the class equation A + C + D = 19/2 becomes "
          "72A^3 - 204A^2 + 122A - 19 = (2A - 1)(36A^2 - 84A + 19) = 0",
          lhs == cub and cub == (2 * A_ - 1) * (36 * A_ ** 2 - 84 * A_ + 19))
    okst = True
    kinds = []
    for Aval in (QS(Fr(1, 2), 0, 30), QS(Fr(7, 6), Fr(1, 6), 30), QS(Fr(7, 6), Fr(-1, 6), 30)):
        okst &= Aval.sign() > 0
        Bval = QS(Fr(10, 3), 0, 30) / (1 + 36 * Aval * Aval)
        Cval = 18 * Aval * Bval + 36 * Aval * Aval * Bval
        Dval = 36 * Aval * Aval * Bval
        okst &= Bval.sign() > 0 and Cval.sign() > 0 and Dval.sign() > 0
        p = {"x1": Aval, "x2": Bval, "x3": Cval, "x4": Dval}
        okst &= all(sgn(fi.value(p)) == 0 for fi in f)
        okst &= (Bval + Dval - Fr(10, 3)).sign() == 0 and (Aval + Cval + Dval - Fr(19, 2)).sign() == 0
        JS, okJ, _ = J_restricted(NC3, kNC3, p)
        okst &= okJ
        tr = JS[0][0] + JS[1][1]
        dt = JS[0][0] * JS[1][1] - JS[0][1] * JS[1][0]
        kinds.append("stable" if (sgn(tr) < 0 and sgn(dt) > 0) else ("saddle" if sgn(dt) < 0 else "other"))
    claim("S4 NC3: exactly three positive steady states in the class, A = 1/2 (saddle) and A = (7 +- sqrt 30)/6 "
          "(trace < 0 < det on S: exponentially stable)", okst and kinds == ["saddle", "stable", "stable"],
          str(kinds))
    # NC
    NC = net("C -> B; 2A -> D; 2B+D -> 2A+2B; 2B+2D -> 2A+2C+D", sp4)
    tefc = all(frozenset(i for i in range(4) if a_[i]) == frozenset(i for i in range(4) if c_[i])
               for (a_, _), (c_, _) in combinations(NC.reactions, 2) if any(a_[i] and c_[i] for i in range(4)))
    claim("S4 NC = C -> B, 2A -> D, 2B+D -> 2A+2B, 2B+2D -> 2A+2C+D: connected, (H1)-(H3), s = 2, three clusters, "
          "TEFC, not PC; no critical siphon on either side",
          connected(NC) and H_all(NC) and NC.rank() == 2 and len(clusters(NC)) == 3 and tefc and not is_PC(NC)
          and no_crit_both(NC))
    kNC = [Fr(2, 3), Fr(5, 18), Fr(24), Fr(64)]
    pt = {"x1": Fr(3), "x2": Fr(1, 2), "x3": Fr(3), "x4": Fr(1, 4)}
    JS, okJ, f = J_restricted(NC, kNC, pt)
    dJ = JS[0][0] * JS[1][1] - JS[0][1] * JS[1][0]
    claim("S4 NC, k = (2/3, 5/18, 24, 64): x = (3, 1/2, 3, 1/4) is a positive steady state with det(J|_S) = -88/9 "
          "in the class B + C = 7/2, A + 2D = 7/2", is_steady(f, pt) and okJ and dJ == Fr(-88, 9) and
          pt["x2"] + pt["x3"] == Fr(7, 2) and pt["x1"] + 2 * pt["x4"] == Fr(7, 2))
    # NM
    NM = net("2X3 -> X1+X2; 3X3 -> X2+2X3; X1+2X2 -> 3X3", X3)
    f, ks, xs, xn = sym_field(NM)
    x1, x2, x3 = xs
    r1, r2, r3 = ks[0] * x3 ** 2, ks[1] * x3 ** 3, ks[2] * x1 * x2 ** 2
    claim("S4 NM = 2X3 -> X1+X2, 3X3 -> X2+2X3, X1+2X2 -> 3X3: connected, (H1)-(H3) on both sides, PC, not EC, "
          "m_src = 3, q = 2, s = 2; critical siphons {X1,X3} and {X2,X3}",
          connected(NM) and H_all(NM) and H_all(NM.dual()) and is_PC(NM) and not is_EC(NM) and
          mcounts(NM) == (3, 2, 2) and set(S for S in critical_siphons(NM) if len(S) == 2) ==
          {frozenset({0, 2}), frozenset({1, 2})})
    p = Poly.var("t")
    pc = p ** 3 - 2 * p ** 2 + Fr(1, 4)
    vals = [pc.value({"t": Fr(v)}) for v in (-1, 0, 1, 2)]
    claim("S4 NM, k = (1,2,1): steady states satisfy x3 = 1/2 and x1 x2^2 = 1/4 (f1 = r1 - r3, f2 = r1 + r2 - 2r3); "
          "in the class x1+x2+x3 = 5/2, x2 is a root of t^3 - 2t^2 + 1/4 in (0,2), which has exactly two "
          "(sign changes at -1, 0, 1, 2: three real roots, one negative)",
          f[0] == r1 - r3 and f[1] == r1 + r2 - 2 * r3 and (f[0] + f[1] + f[2]).is_zero() and
          vals[0] < 0 < vals[1] and vals[2] < 0 < vals[3])
    NB = net("X3 -> X2; 2X3 -> X1; X1+X2 -> 3X3", X3)
    f, ks, xs, xn = sym_field(NB)
    x1, x2, x3 = xs
    okNB = True
    for sgn_ in (1, -1):
        x1v = QS(Fr(5, 4), Fr(sgn_, 4), 17)
        x2v = 5 - 2 * x1v
        pt = {"x1": x1v, "x2": x2v, "x3": QS(1, 0, 17), "k1": 1, "k2": 1, "k3": 1}
        okNB &= x1v.sign() > 0 and x2v.sign() > 0 and all(sgn(fi.value(pt)) == 0 for fi in f) and \
            (2 * x1v + x2v + 1 - 6).sign() == 0
    claim("S4 NB = X3 -> X2, 2X3 -> X1, X1+X2 -> 3X3: BH (blocks {X3}, {X1,X2}), q = 2 = s; for k = 1 the class "
          "2x1+x2+x3 = 6 contains the two positive steady states x3 = 1, x1 = (5 +- sqrt 17)/4, x2 = 5 - 2 x1",
          is_BH(NB) and H_all(NB) and H_all(NB.dual()) and connected(NB) and mcounts(NB) == (3, 2, 2) and okNB
          and (2, 1, 1) in NB.conservation_rays())
    # boundary table: the cluster rank equation fails for each counterexample
    M = NETS["M"]
    E4 = net("X4 -> 3X3; 3X3 -> X1+2X2; X1+2X2 -> 3X2; 3X2 -> X4", X4)
    rows = {"U*": (Ud, (3, 2, 2)), "U": (U, (3, 2, 2)), "3A<->2A+B": (M, (2, 1, 1)), "4-cycle": (E4, (4, 3, 3)),
            "NM": (NM, (3, 2, 2)), "NB": (NB, (3, 2, 2)), "NT": (NT, (2, 2, 2))}
    okT = all(mcounts(N) == c and N.rank() != len(clusters(N)) - 1 and connected(N) and H_all(N) and
              H_all(N.dual()) for N, c in rows.values())
    claim("S4 table of failures: (m_src, m_cl, s) = U* (3,2,2), U (3,2,2), 3A<->2A+B (2,1,1), 4-cycle (4,3,3), "
          "NM (3,2,2), NB (3,2,2), NT (2,2,2); each is connected with (H1)-(H3) on both sides and violates "
          "s = m_cl - 1", okT)
    # Theorem E / E2.5': WR and (H1) give im Gamma in im Ye
    NA = net("3X3 -> X1+X2; X1+X2 -> X2+2X3+3X4; X2+2X3+3X4 -> 3X3", X4)
    claim("S4 im Gamma is contained in im Ye for the network of Theorem A and for W, T, M (weakly reversible)",
          all(crn.rank(transpose(N.Ye) + transpose(N.Gamma)) == crn.rank(transpose(N.Ye)) for N in
              (NA, NETS["W"], NETS["T"], M)))


# =============================================================================================================
# 9. Section 5: what the duality does not transfer
# =============================================================================================================
def sec5():
    section("S5 Non-transfers: weak reversibility, ADLS, persistence, convergence")
    S2 = net("2X2 -> X1; X1 -> 2X2", X2)
    claim("S5 2X2 <-> X1 is conservative and WR; its dual k2 -> k1, 2k1 -> 2k2 is not WR",
          S2.is_conservative() and is_wr(S2) and same_reactions(S2.dual(), "k2 -> k1; 2k1 -> 2k2", Ks(2))
          and not is_wr(S2.dual()))
    H = net("2A -> 2B; 2B -> 2A; 2A -> A+B; A+B -> 2A", ["A", "B"])
    Hd = H.dual()
    claim("S5 H = 2A <-> 2B, 2A <-> A+B: WR, (H1)-(H3), only minimal siphon {A,B} = support of a+b (no critical "
          "siphon); H* = 2U1+2U3+U4 -> 2U2+U3+2U4, 2U2+U4 -> 2U1+U3 has (H1)-(H3) and the critical siphon {U4}",
          is_wr(H) and H_all(H) and H.minimal_siphons() == [frozenset({0, 1})] and critical_siphons(H) == [] and
          same_reactions(Hd, "2U1+2U3+U4 -> 2U2+U3+2U4; 2U2+U4 -> 2U1+U3", ["U1", "U2", "U3", "U4"]) and
          H_all(Hd) and frozenset({3}) in critical_siphons(Hd))
    E = net("X4 -> X2+X3; X2+X3 -> X4; 2X4 -> X1+X3; X1+X3 -> 2X4", X4)
    claim("S5 converse (agent E): X4 <-> X2+X3, 2X4 <-> X1+X3 is WR with (H1)-(H3) and has the critical siphons "
          "{X3,X4}, {X2,X3,X4}; its dual satisfies (H1)-(H3) and has no critical siphon",
          is_wr(E) and H_all(E) and H_all(E.dual()) and
          set(critical_siphons(E)) == {frozenset({2, 3}), frozenset({1, 2, 3})} and critical_siphons(E.dual()) == [])
    M = NETS["M"]
    D = M.dual()
    claim("S5 M = 3A <-> 2A+B: (m, l, s, delta) = (2,1,1,0), reversible; D = M* = 3U+2V -> 2U+3V, V -> U: "
          "(4,2,1,1), not WR; both satisfy (H1)-(H3) with Gamma = [[-1,1],[1,-1]]",
          counts(M) == (2, 1, 1, 0) and counts(D) == (4, 2, 1, 1) and is_wr(M) and not is_wr(D) and
          H_all(M) and H_all(D) and M.Gamma == D.Gamma == [[-1, 1], [1, -1]])
    fD, kD, xD, _ = sym_field(D, "c", ["u", "v"])
    u, v = xD
    c1, c2 = kD
    claim("S5 D, all rates: v' = v (c1 u^3 v - c2); with u <= T and v <= min(T, c2/(2 c1 T^3)): v' <= -(c2/2) v",
          fD[1] == v * (c1 * u ** 3 * v - c2))
    fM, kM, xM, _ = sym_field(M, "k", ["a", "b"])
    a, b = xM
    claim("S5 M, all rates: a' = a^2 (k2 b - k1 a) (one globally attracting equilibrium on each segment)",
          fM[0] == a ** 2 * (kM[1] * b - kM[0] * a))
    t = Poly.var("t")
    h = (2 - t) ** 3 * t
    cub = t ** 3 - 5 * t ** 2 + 7 * t - 1
    disc = 18 * (-5) * 7 * (-1) - 4 * (-5) ** 3 * (-1) + (-5) ** 2 * 7 ** 2 - 4 * 7 ** 3 - 27 * (-1) ** 2
    claim("S5 D, unit rates, class u+v = 2: v' = v (h(v) - 1), h = (2-v)^3 v, h' = (2-v)^2 (2-4v), "
          "h(1/8) = 3375/4096, so v' <= -(721/4096) v on (0,1/8]; h - 1 = -(v-1)(v^3 - 5v^2 + 7v - 1) and the cubic "
          "has negative discriminant (-44) and one root in (1/8, 1/2): exactly two positive equilibria",
          h.diff("t") == (2 - t) ** 2 * (2 - 4 * t) and h.value({"t": Fr(1, 8)}) == Fr(3375, 4096) and
          1 - Fr(3375, 4096) == Fr(721, 4096) and h - 1 == -(t - 1) * cub and disc == -44 and
          cub.value({"t": Fr(1, 8)}) < 0 < cub.value({"t": Fr(1, 2)}) and
          fD[1].subs({"c1": 1, "c2": 1, "u": 2 - v}) == v * (((2 - v) ** 3) * v - 1))
    # drain species
    claim("S5 drain species: V (the dual species of 2A+B -> 3A) is a drain species of D with d0 = 1; M has none",
          drain_species(D) == [(1, 1)] and drain_species(M) == [])
    E4 = net("X4 -> 3X3; 3X3 -> X1+2X2; X1+2X2 -> 3X2; 3X2 -> X4", X4)
    E4d = E4.dual()
    claim("S5 4-cycle X4 -> 3X3 -> X1+2X2 -> 3X2 -> X4: WR, (H1)-(H3), (m,l,s,delta) = (4,1,3,0); the dual has "
          "(H1)-(H3) and exactly one drain species, k3 (dual to X1+2X2 -> 3X2), d0 = 1; N has none",
          is_wr(E4) and H_all(E4) and counts(E4) == (4, 1, 3, 0) and H_all(E4d) and drain_species(E4d) == [(2, 1)]
          and drain_species(E4) == [])
    E5 = net("X5 -> 3X4; 3X4 -> X2+X4; X2+X4 -> X1+2X3; X1+2X3 -> 3X3; 3X3 -> X5", X5)
    E5d = E5.dual()
    ds5 = drain_species(E5d)
    claim("S5 5-cycle X5 -> 3X4 -> X2+X4 -> X1+2X3 -> 3X3 -> X5: WR, (H1)-(H3), (m,l,s,delta) = (5,1,4,0); the "
          "dual has (H1)-(H3) and the drain species k4 (dual to X1+2X3 -> 3X3)",
          is_wr(E5) and H_all(E5) and counts(E5) == (5, 1, 4, 0) and H_all(E5d) and ds5 == [(3, 1)], str(ds5))
    # drain set example C and the pure WR example
    C = net("2X3 -> X2+X3; X2+X3 -> X1+X3; X1+X3 -> 2X3", X3)
    Cd = C.dual()
    fC, kC, xC, _ = sym_field(Cd, "c", ["k1", "k2", "k3"])
    k1, k2, k3 = xC
    claim("S5 C = 2X3 -> X2+X3 -> X1+X3 -> 2X3: WR, (H1)-(H3), (3,1,2,0); C* = k3 -> k2, k2 -> k1, "
          "2k1+k2+k3 -> k1+k2+2k3; V = k2 + 2k3 has V' = -c1 k3 - c2 k2 + 2 c3 k1^2 k2 k3; {k2,k3} is a drain set "
          "(w = (1,2), d = (1,1)); C* has no drain species",
          is_wr(C) and H_all(C) and counts(C) == (3, 1, 2, 0) and
          same_reactions(Cd, "k3 -> k2; k2 -> k1; 2k1+k2+k3 -> k1+k2+2k3", Ks(3)) and
          fC[1] + 2 * fC[2] == -kC[0] * k3 - kC[1] * k2 + 2 * kC[2] * k1 ** 2 * k2 * k3 and
          is_drain_set(Cd, [1, 2], {1: 1, 2: 2}, {1: 1, 2: 1}) and drain_species(Cd) == [])
    Q6 = net("X4 -> X2+X3; X2+X3 -> X4; 3X3 -> X1+2X4; X1+2X4 -> 3X3", X4)
    Q6d = Q6.dual()
    claim("S5 X4 <-> X2+X3, 3X3 <-> X1+2X4: pure, WR, (H1)-(H3), (4,2,2,0); its dual has the drain set "
          "of two species with d = (1,1) (pure weakly reversible networks can have non-persistent duals)",
          pure(Q6) and is_wr(Q6) and H_all(Q6) and counts(Q6) == (4, 2, 2, 0) and
          any(is_drain_set(Q6d, [i, j], {i: wi, j: wj}, {i: 1, j: 1})
              for i, j in combinations(range(4), 2) for wi in range(1, 5) for wj in range(1, 5)),
          str([((i + 1, j + 1), (wi, wj)) for i, j in combinations(range(4), 2) for wi in range(1, 5)
               for wj in range(1, 5) if is_drain_set(Q6d, [i, j], {i: wi, j: wj}, {i: 1, j: 1})][:2]))


# =============================================================================================================
# 10. Section 6: persistence of duals in low rank
# =============================================================================================================
def faces_str(ft, names):
    return ", ".join("{" + ",".join(names[i] for i in Z) + "}:" + t for _, Z, t in ft)


def single_facet_type(N, z):
    """Type of the facet {k_z = 0} from the field of species z: lowest power of x_z among the terms of f_z;
    IN (all such terms positive), OUT (all negative), ESC/DOM (mixed: is every negative exponent off z in the
    convex hull of the positive ones?).  Rates symbolic, so terms with different reactions do not cancel."""
    f, ks, xs, xn = sym_field(N)
    e, part = f[z].low_order_in(xn[z])
    pos, neg = [], []
    for m, c in part.t.items():
        d = dict(m)
        expo = tuple(d.get(xn[i], 0) for i in range(N.n) if i != z)
        (pos if c > 0 else neg).append(expo)
    if not neg:
        return "IN", e
    if not pos:
        return "OUT", e
    return ("DOM" if all(in_convex_hull(q, pos) for q in neg) else "ESC"), e


def sec6():
    section("S6 Persistence of duals in low rank")
    X = net("2X2 -> X1+X2; X1+X2 -> 2X2", X2)
    Xd = X.dual()
    f, ks, xs, xn = sym_field(Xd, "c", ["k1", "k2"])
    v_, ft = rank_one_verdict(Xd)
    claim("S6 X = 2X2 <-> X1+X2: WR, (H1)-(H3), deficiency 0; X* = k2 -> k1, 2k1+k2 -> k1+2k2; "
          "k2' = k2 (c2 k1^2 - c1); faces of X*: " + faces_str(ft, ["k1", "k2"]) + "; verdict NONE",
          is_wr(X) and H_all(X) and X.deficiency() == 0 and
          same_reactions(Xd, "k2 -> k1; 2k1+k2 -> k1+2k2", Ks(2)) and
          f[1] == xs[1] * (ks[1] * xs[0] ** 2 - ks[0]) and v_ == "NONE" and
          any(Z == (1,) and t == "ESC" for _, Z, t in ft))
    claim("S6 X*: no drain species and {k2} with d = 1 is not a drain set", drain_species(Xd) == [] and
          not is_drain_set(Xd, [1], {1: 1}, {1: 1}))
    Md = NETS["M"].dual()
    v_, ft = rank_one_verdict(Md)
    claim("S6 (3A <-> 2A+B)*: faces " + faces_str(ft, ["U", "V"]) + "; verdict NONE (an OUT face)",
          v_ == "NONE" and any(t == "OUT" for _, _, t in ft))
    P = net("2X1+X3 -> 3X1+2X2; 3X1+2X2 -> 2X1+X3; 2X2+X3 -> X1+4X2; X1+4X2 -> 2X2+X3", X3)
    Pd = P.dual()
    v_, ft = rank_one_verdict(Pd)
    claim("S6 P = 2X1+X3 <-> 3X1+2X2, 2X2+X3 <-> X1+4X2: WR, s = 1, (H1)-(H3); every singleton face of P* is IN, "
          "the face {k1,k3} of the two forward reactions is ESC: verdict NONE",
          is_wr(P) and H_all(P) and P.rank() == 1 and v_ == "NONE" and
          all(t == "IN" for _, Z, t in ft if len(Z) == 1) and
          [Z for _, Z, t in ft if t == "ESC"] == [(0, 2)], faces_str(ft, Ks(4)))
    S = net("X1+X2 -> 2X2+2X3; X1+X2+2X3 -> 2X1", X3)
    Sd = S.dual()
    f, ks, xs, xn = sym_field(Sd, "c", ["k1", "k2"])
    v_, ft = rank_one_verdict(Sd)
    claim("S6 S = X1+X2 -> 2X2+2X3, X1+X2+2X3 -> 2X1: (H1)-(H3), s = 1, not WR; S* = k1+k2 -> 2k2, "
          "k1+k2 -> 2k1, 2k2 -> 2k1; k1' = k2((c2 - c1) k1 + 2 c3 k2); DOM faces at one end only, others IN",
          H_all(S) and S.rank() == 1 and not is_wr(S) and
          same_reactions(Sd, "k1+k2 -> 2k2; k1+k2 -> 2k1; 2k2 -> 2k1", Ks(2)) and
          f[0] == xs[1] * ((ks[1] - ks[0]) * xs[0] + 2 * ks[2] * xs[1]) and v_ == "DOM" and
          len({sg for sg, Z, t in ft if t == "DOM"}) == 1, faces_str(ft, Ks(2)))
    Bn = net("X1+X2 -> 2X2; X1+X2 -> 2X1", X2)
    Bd = Bn.dual()
    f, ks, xs, xn = sym_field(Bd, "c", ["k1", "k2"])
    v_, ft = rank_one_verdict(Bd)
    claim("S6 B = X1+X2 -> 2X2, X1+X2 -> 2X1: (H1)-(H3); B* = k1+k2 -> 2k2, k1+k2 -> 2k1, k1' = (c2 - c1) k1 k2; "
          "DOM faces at both ends (persistent only on the hypersurface c1 = c2)",
          H_all(Bn) and f[0] == (ks[1] - ks[0]) * xs[0] * xs[1] and v_ == "DOM" and
          len({sg for sg, Z, t in ft if t == "DOM"}) == 2)
    Z0 = net("0 -> A; A -> 0", ["A"])
    Z0d = Z0.dual()
    v_, ft = rank_one_verdict(Z0d)
    claim("S6 0 <-> A: pure, WR, consistent, s = 1, not conservative (unbounded classes); its dual k2 -> k1 is "
          "conservative with an OUT face (NONE): the corollary for pure networks needs (H1)",
          pure(Z0) and is_wr(Z0) and Z0.is_consistent() and not Z0.is_conservative() and Z0.rank() == 1 and
          same_reactions(Z0d, "k2 -> k1", Ks(2)) and Z0d.is_conservative() and v_ == "NONE" and
          any(t == "OUT" for _, _, t in ft))
    Zx = net("A+B -> 2A; 2A -> A+B; 3A+B -> 4A", ["A", "B"])
    f, ks, xs, xn = sym_field(Zx, "k", ["a", "b"])
    v_, ft = rank_one_verdict(Zx)
    # w_I-endotacticity for the indicator of the siphon {A}: the w-maximal source 3A+B points outward
    w = (1, 0)
    hts = [sum(w[i] * a[i] for i in range(2)) for a, b in Zx.reactions]
    dl = [sum(w[i] * (b[i] - a[i]) for i in range(2)) for a, b in Zx.reactions]
    hm = max(h for h, d_ in zip(hts, dl) if d_ != 0)
    claim("S6 A+B -> 2A, 2A -> A+B, 3A+B -> 4A: conservative, rank one, every face IN (persistent for all rates), "
          "a' = a (k1 b - k2 a + k3 a^2 b); {A} is a siphon and the w_{A}-maximal source 3A+B points outward, so "
          "the network is not W_I-endotactic (Zhang-Gao-Dochain's sufficient condition does not apply)",
          Zx.is_conservative() and Zx.rank() == 1 and v_ == "ALL" and
          f[0] == xs[0] * (ks[0] * xs[1] - ks[1] * xs[0] + ks[2] * xs[0] ** 2 * xs[1]) and
          Zx.is_siphon((0,)) and any(d_ > 0 and h == hm for h, d_ in zip(hts, dl)))
    Qn = net("X1 -> X2; X1 -> X3; X2+X3 -> 2X1", X3)
    Qd = Qn.dual()
    f, ks, xs, xn = sym_field(Qd, "c", ["k1", "k2", "k3"])
    claim("S6 Q = X1 -> X2, X1 -> X3, X2+X3 -> 2X1 (pure, (H1)-(H3), s = 2, not WR): Q* = k1+k2 -> 2k3, k3 -> k1, "
          "k3 -> k2 and (k2 - k1)' = (c3 - c2) k3",
          pure(Qn) and H_all(Qn) and Qn.rank() == 2 and not is_wr(Qn) and
          same_reactions(Qd, "k1+k2 -> 2k3; k3 -> k1; k3 -> k2", Ks(3)) and
          f[1] - f[0] == (ks[2] - ks[1]) * xs[2])
    # rank two examples (D2)
    AL = net("2X3 -> X2; X2 -> 2X1+X3; 2X1+X3 -> 2X3", X3)
    ALd = AL.dual()
    f, ks, xs, xn = sym_field(ALd, "c", ["k1", "k2", "k3"])
    sips1 = [Z for Z in ((0,), (1,), (2,)) if ALd.is_siphon(Z)]
    sips2 = [Z for Z in ((0, 1), (0, 2), (1, 2)) if ALd.is_siphon(Z)]
    claim("S6 rank two, ALL example 2X3 -> X2 -> 2X1+X3 -> 2X3: WR, (H1)-(H3), s = 2; dual 2k3 -> 2k2, k2 -> k1, "
          "2k1+k3 -> k2+2k3 with k3' = k3 (c3 k1^2 - 2 c1 k3); siphon edges {k3} only, siphon vertex {k2,k3} only",
          is_wr(AL) and H_all(AL) and AL.rank() == 2 and
          same_reactions(ALd, "2k3 -> 2k2; k2 -> k1; 2k1+k3 -> k2+2k3", Ks(3)) and
          f[2] == xs[2] * (ks[2] * xs[0] ** 2 - 2 * ks[0] * xs[2]) and sips1 == [(2,)] and sips2 == [(1, 2)])
    EE = net("2X3 -> 2X2+X3; 2X2+X3 -> X1+X2+X3; X1+X2+X3 -> 2X3", X3)
    EEd = EE.dual()
    t_, e_ = single_facet_type(EEd, 2)
    claim("S6 rank two, NONE example (edge equilibrium) 2X3 -> 2X2+X3 -> X1+X2+X3 -> 2X3: WR, (H1)-(H3), s = 2; the "
          "facet {k3} of the dual is a siphon of type ESC with minimal order 1 and zero edge flow (every reaction "
          "of the dual that changes the state has k3 in its source)",
          is_wr(EE) and H_all(EE) and EE.rank() == 2 and EEd.is_siphon((2,)) and t_ == "ESC" and e_ == 1 and
          all(a[2] >= 1 for a, b in EEd.reactions if a != b))
    CM = net("2X3 -> X2+X3; X2+X3 -> X1+X2; X1+X2 -> 2X3", X3)
    CMd = CM.dual()
    f, ks, xs, xn = sym_field(CMd, "c", ["k1", "k2", "k3"])
    Tt, aa, bb, h2 = V("T", "a", "b", "h2")
    sub = {"k1": aa, "k3": bb, "k2": Tt - aa - bb}
    fa, fb = Poly.lift(f[0].subs(sub)), Poly.lift(f[2].subs(sub))

    def coef(p, ea, eb):
        out = {}
        for m, c in p.t.items():
            d = dict(m)
            if d.get("a", 0) == ea and d.get("b", 0) == eb:
                d.pop("a", None)
                d.pop("b", None)
                out[tuple(sorted(d.items()))] = c
        return Poly(out)
    c1, c2, c3 = ks
    lin_ok = (coef(fa, 0, 1) == c2 * Tt and coef(fa, 2, 0) == -c3 * Tt and coef(fa, 1, 0).is_zero() and
              coef(fb, 0, 1) == -(c1 + c2 * Tt) and coef(fb, 2, 0) == 2 * c3 * Tt and coef(fb, 1, 0).is_zero())
    # centre manifold b = h2 a^2 with (c1 + c2 T) h2 = 2 c3 T; reduced a' = (c2 T h2 - c3 T) a^2 + ...
    red = (c2 * Tt * (2 * c3 * Tt) - c3 * Tt * (c1 + c2 * Tt))
    claim("S6 rank two, centre-manifold example 2X3 -> X2+X3 -> X1+X2 -> 2X3: dual k3 -> k2, k2+k3 -> k1+k2, "
          "2k1+k2 -> k1+2k3; at (0,T,0) with a = k1, b = k3: a' = c2 T b - c3 T a^2 + ..., "
          "b' = -(c1 + c2 T) b + 2 c3 T a^2 + ...; reduced coefficient c3 T (c2 T - c1)/(c1 + c2 T)",
          same_reactions(CMd, "k3 -> k2; k2+k3 -> k1+k2; 2k1+k2 -> k1+2k3", Ks(3)) and lin_ok and
          red == c3 * Tt * (c2 * Tt - c1))
    LL = net("2X3 -> 2X2+X3; 2X2+X3 -> 2X1+X3; 2X1+X3 -> 2X3", X3)
    LLd = LL.dual()
    f, ks, xs, xn = sym_field(LLd, "c", ["k1", "k2", "k3"])
    claim("S6 rank two, Lyapunov example 2X3 -> 2X2+X3 -> 2X1+X3 -> 2X3: V = k2 + 2k3 has "
          "V' = -2c1 k3^2 - 2c2 k2^2 + 2c3 k1^2 k2 k3; the form is negative definite iff c3^2 T^4 < 4 c1 c2",
          f[1] + 2 * f[2] == -2 * ks[0] * xs[2] ** 2 - 2 * ks[1] * xs[1] ** 2 + 2 * ks[2] * xs[0] ** 2 * xs[1] * xs[2]
          and det([[-2 * ks[1], ks[2] * Tt ** 2], [ks[2] * Tt ** 2, -2 * ks[0]]]) ==
          4 * ks[0] * ks[1] - ks[2] ** 2 * Tt ** 4)
    E27 = net("2X3 -> X2; X2 -> X1+X3; X1+X3 -> 2X3", X3)
    E27d = E27.dual()
    f, ks, xs, xn = sym_field(E27d, "c", ["k1", "k2", "k3"])
    sub = {"k1": Tt - Poly.var("k2") - Poly.var("k3")}
    g2, g3 = Poly.lift(f[1].subs(sub)), Poly.lift(f[2].subs(sub))
    Jv = [[g.diff(v_).subs({"k2": 0, "k3": 0}) for v_ in ("k2", "k3")] for g in (g2, g3)]
    claim("S6 vertex scaling example 2X3 -> X2 -> X1+X3 -> 2X3: dual k3 -> k2, k2 -> k1, 2k1+k3 -> k2+2k3; "
          "linearisation at (T,0,0) in (k2,k3) = [[-c2, c1 + c3 T^2], [0, c3 T^2 - c1]]",
          same_reactions(E27d, "k3 -> k2; k2 -> k1; 2k1+k3 -> k2+2k3", Ks(3)) and H_all(E27) and is_wr(E27) and
          Poly.lift(Jv[0][0]) == -ks[1] and Poly.lift(Jv[0][1]) == ks[0] + ks[2] * Tt ** 2 and
          Poly.lift(Jv[1][0]).is_zero() and Poly.lift(Jv[1][1]) == ks[2] * Tt ** 2 - ks[0])
    UU = net("3X3 -> 2X1+X2; X2+2X3 -> 3X3; 2X1+X2 -> X2+2X3", X3)
    UUd = UU.dual()
    f, ks, xs, xn = sym_field(UUd, "c", ["k1", "k2", "k3"])
    k1v, a_, b_ = xs  # a = k2, b = k3
    c1, c2, c3 = ks
    K = c3 * k1v ** 3
    # z = a/b: z' b^2 = a' b - a b'  ;  claim z' = a[(2c1 - c2) + K z(1 - 2z)]  => times b^3
    lhs = f[1] * b_ - a_ * f[2]
    rhs = a_ * ((2 * c1 - c2) * b_ ** 2 + K * a_ * (b_ - 2 * a_))
    lhs_u = f[2] * a_ - b_ * f[1]
    rhs_u = a_ * ((c2 - 2 * c1) * b_ ** 2 - K * b_ * a_ + 2 * K * a_ ** 2)
    claim("S6 example U 3X3 -> 2X1+X2, X2+2X3 -> 3X3, 2X1+X2 -> X2+2X3 (WR, (H1)-(H3), s = 2): dual "
          "2k3 -> 2k1, k2+k3 -> k1+k3, 3k1+2k2 -> 3k2+2k3; with a = k2, b = k3, K = c3 k1^3: "
          "z = a/b has z' = a[(2c1 - c2) + K z(1 - 2z)], u = b/a has u' = a[(c2 - 2c1)u^2 - K u + 2K], "
          "a' = a^2 (K - c2 u), b' = b^2(-2c1 + 2K z^2)",
          is_wr(UU) and H_all(UU) and UU.rank() == 2 and
          same_reactions(UUd, "2k3 -> 2k1; k2+k3 -> k1+k3; 3k1+2k2 -> 3k2+2k3", Ks(3)) and
          lhs == rhs and lhs_u == rhs_u and f[1] == a_ * (K * a_ - c2 * b_) and
          f[2] == -2 * c1 * b_ ** 2 + 2 * K * a_ ** 2)


def thmA_network(u, R):
    return net(f"{u}X3 -> X1+X2; X1+X2 -> X2+2X3+{R}X4; X2+2X3+{R}X4 -> {u}X3", X4)


def sec6_A():
    section("S6.2 Theorems A, AF (refutation of the OUT-facet conjecture) and O (periodic orbit)")
    okfam = True
    detail = []
    a, b, c = V("a", "b", "c")
    al, be, ga, de = V("al", "be", "ga", "de")
    for u in range(3, 8):
        for R in range(1, 8):
            N = thmA_network(u, R)
            z = (u * R + 2, u - 2, R + 1, u - 2)
            G = N.Gamma
            law = all(sum(z[i] * G[i][j] for i in range(4)) == 0 for j in range(3))
            val = all(sum(z[i] * cpx[i] for i in range(4)) == u * (R + 1) for cpx in N.complexes())
            flux = all(sum(G[i][j] for j in range(3)) == 0 for i in range(4))
            hyp = law and val and flux and N.rank() == 2 and N.rows_distinct() and is_wr(N) and \
                N.deficiency() == 0 and N.is_conservative() and N.is_consistent()
            Nd = N.dual()
            dual_ok = same_reactions(Nd, f"k2 -> k1; k2+k3 -> k1+k2; {u}k1+2k3 -> 2k2+{u}k3; {R}k3 -> {R}k2", Ks(3))
            f = field(Nd, [al, be, ga, de], [a, b, c])
            eq_ok = (f[1] == -al * b + 2 * ga * a ** u * c ** 2 + R * de * c ** R and
                     f[2] == -be * b * c + (u - 2) * ga * a ** u * c ** 2 - R * de * c ** R)
            # (AF.1): al * (c' - c (be/al) b') = c^2 [ ... ] times al
            ident = (al * f[2] - be * c * f[1] ==
                     c ** 2 * (al * (u - 2) * ga * a ** u - al * R * de * c ** (R - 2) - 2 * be * ga * a ** u * c
                               - R * be * de * c ** (R - 1))) if R >= 2 else \
                (al * f[2] - be * c * f[1] ==
                 c * (al * (u - 2) * ga * a ** u * c - al * de - 2 * be * ga * a ** u * c ** 2 - be * de * c))
            e_, part = f[2].low_order_in("c")
            out_ok = e_ == 1 and all(cc < 0 for cc in part.coeffs())
            okfam &= hyp and dual_ok and eq_ok and ident and out_ok
            if not (hyp and dual_ok and eq_ok and ident and out_ok):
                detail.append((u, R, hyp, dual_ok, eq_ok, ident, out_ok))
    claim("S6 Theorem AF, u = 3..7, R = 1..7: N(u,R) = uX3 -> X1+X2 -> X2+2X3+RX4 -> uX3 is WR, deficiency 0, rank 2, "
          "(H1)-(H3) (law (uR+2, u-2, R+1, u-2) with value u(R+1) on every complex, flux (1,1,1)); its dual is "
          "b -> a, b+c -> a+b, ua+2c -> 2b+uc, Rc -> Rb; b', c' as stated; identity (AF.1) for V = c exp(-be b/al); "
          "the facet c = 0 is OUT (all terms of c' of order 1 in c are negative)", okfam, str(detail[:3]))
    # Theorem A specifics (u = R = 3)
    N = thmA_network(3, 3)
    Nd = N.dual()
    f = field(Nd, [al, be, ga, de], [a, b, c])
    claim("S6 Theorem A: a' = al b + be b c - 3 ga a^3 c^2, b' = -al b + 2 ga a^3 c^2 + 3 de c^3, "
          "c' = -be b c + ga a^3 c^2 - 3 de c^3; on a = 0: a' = al b + be b c; on b = 0: b' = 2 ga a^3 c^2 + 3 de c^3; "
          "on c = 0: b' = -al b",
          f[0] == al * b + be * b * c - 3 * ga * a ** 3 * c ** 2 and
          f[1] == -al * b + 2 * ga * a ** 3 * c ** 2 + 3 * de * c ** 3 and
          f[2] == -be * b * c + ga * a ** 3 * c ** 2 - 3 * de * c ** 3 and
          f[0].subs({"a": 0}) == al * b + be * b * c and f[1].subs({"b": 0}) == 2 * ga * a ** 3 * c ** 2 + 3 * de * c ** 3
          and f[1].subs({"c": 0}) == -al * b and (f[0] + f[1] + f[2]).is_zero())
    claim("S6 Theorem A, identity (A.1): al c' - be c b' = c^2 [al ga a^3 - (3 al de + 2 be ga a^3) c - 3 be de c^2]",
          al * f[2] - be * c * f[1] == c ** 2 * (al * ga * a ** 3 - (3 * al * de + 2 * be * ga * a ** 3) * c
                                                - 3 * be * de * c ** 2))
    # Theorem O: uniqueness of the equilibrium per class
    claim("S6 Theorem O: b' - 2c' = (2 be c - al) b + 9 de c^3",
          f[1] - 2 * f[2] == (2 * be * c - al) * b + 9 * de * c ** 3)
    # parametrisation: with b (al - 2 be c) = 9 de c^3 and ga a^3 (al - 2 be c) = 3 de c (al + be c), all three
    # equations hold: check by clearing the denominator (al - 2 be c)
    Dn = al - 2 * be * c
    bb_ = 9 * de * c ** 3        # = b * Dn
    a3D = 3 * de * c * (al + be * c)  # = ga a^3 * Dn
    # c' * Dn = -be (b Dn) c + (ga a^3 Dn) c^2 - 3 de c^3 Dn
    cprime = -be * bb_ * c + a3D * c ** 2 - 3 * de * c ** 3 * Dn
    bprime = -al * bb_ + 2 * a3D * c ** 2 + 3 * de * c ** 3 * Dn
    claim("S6 Theorem O: on 0 < c < al/(2 be), b = 9 de c^3/(al - 2 be c) and ga a^3 = 3 de c (al + be c)/(al - 2 be c) "
          "solve b' = c' = 0 (identities after clearing the denominator)", cprime.is_zero() and bprime.is_zero())
    # monotonicity: d/dc [c^3/(al - 2be c)] numerator = c^2 (3 al - 4 be c); d/dc [c (al + be c)/(al - 2 be c)]
    # numerator = al^2 + 2 al be c - 2 be^2 c^2
    n1 = (3 * c ** 2) * Dn - c ** 3 * (-2 * be)
    n2 = (al + 2 * be * c) * Dn - c * (al + be * c) * (-2 * be)
    claim("S6 Theorem O: d/dc of c^3/(al - 2be c) has numerator c^2 (3al - 4be c) > 0 and d/dc of "
          "c(al + be c)/(al - 2be c) has numerator al^2 + 2 al be c - 2 be^2 c^2 > 0 on (0, al/(2be))",
          n1 == c ** 2 * (3 * al - 4 * be * c) and n2 == al ** 2 + 2 * al * be * c - 2 * be ** 2 * c ** 2)
    # exact source certificate
    rates = {"al": Fr(17), "be": Fr(400), "ga": Fr(70000), "de": Fr(1000000)}
    pt = {"a": Fr(1), "b": Fr(1), "c": Fr(1, 100)}
    fl = [rates["al"] * pt["b"], rates["be"] * pt["b"] * pt["c"], rates["ga"] * pt["a"] ** 3 * pt["c"] ** 2,
          rates["de"] * pt["c"] ** 3]
    eq0 = all(fi.value({**rates, **pt}) == 0 for fi in f)
    T_ = pt["a"] + pt["b"] + pt["c"]
    sub = {"a": Poly.const(T_) - Poly.var("b") - Poly.var("c")}
    gb, gc = Poly.lift(f[1].subs(sub)), Poly.lift(f[2].subs(sub))
    J = [[g.diff(v_).value({**rates, "b": pt["b"], "c": pt["c"]}) for v_ in ("b", "c")] for g in (gb, gc)]
    tr = J[0][0] + J[1][1]
    dt = J[0][0] * J[1][1] - J[0][1] * J[1][0]
    Hm = [[Fr(93655), Fr(-290457)], [Fr(-290457), Fr(13471234)]]
    JtH = matmul(transpose(J), Hm)
    HJ = matmul(Hm, J)
    S_ = [[JtH[i][j] + HJ[i][j] for j in range(2)] for i in range(2)]
    claim("S6 Theorem O: at (al,be,ga,de) = (17,400,70000,10^6), (1,1,1/100) is an equilibrium (fluxes 17,4,7,1) in "
          "the class T = 201/100; Jacobian in (b,c) = [[-59,3658],[-25,79]], trace 20, det 86789 "
          "(discriminant 20^2 - 4*86789 < 0: complex eigenvalues 10 +- i sqrt(86689))",
          eq0 and fl == [17, 4, 7, 1] and T_ == Fr(201, 100) and J == [[-59, 3658], [-25, 79]] and tr == 20 and
          dt == 86789 and tr ** 2 - 4 * dt == -4 * 86689)
    claim("S6 Theorem O: H = [[93655,-290457],[-290457,13471234]] is positive definite (det 1177283151421) and "
          "J^T H + H J = 3471560 I", Hm[0][0] > 0 and det(Hm) == 1177283151421 and
          S_ == [[3471560, 0], [0, 3471560]])


# =============================================================================================================
# 11. Section 7: the GAC and small codimension
# =============================================================================================================
def beta_delta(N):
    """(beta, delta) of a network: beta = dim ker of the complex incidence map = r - (m - l)."""
    m, l, s, d = counts(N)
    return N.r - (m - l), d


def sec7():
    section("S7 The GAC: dimension, codimension, boundary dictionary")
    okc = True
    for n in range(2, 7):
        sp_ = Xs(n)
        N = net("; ".join(f"X{i + 1} -> X{(i + 1) % n + 1}" for i in range(n)), sp_)
        Nd = N.dual()
        okc &= (N.n - N.rank() == 1 and N.rank() == n - 1 and H_all(N) and is_wr(N) and N.deficiency() == 0 and
                ordinary(Nd) and all(sum(a) == 1 and sum(b) == 1 for a, b in Nd.reactions) and is_wr(Nd) and
                Nd.linkage_classes() == 1 and len(Nd.complexes()) == n)
    claim("S7 the monomolecular n-cycle (n = 2..6) has codimension q = n - s = 1, s = n - 1, (H1)-(H3), WR, "
          "deficiency 0, and its dual is again a monomolecular n-cycle", okc)
    okq = True
    for name in ("R", "T", "W", "U", "D", "F", "N3", "N412", "M"):
        N = NETS[name]
        Nd = N.dual()
        bs, ds = beta_delta(Nd)
        okq &= N.n - N.rank() == bs + ds
    claim("S7 q = beta* + delta* (q = n - s of N; beta*, delta* of the complex graph of N*) on R, T, W, U, D, F, "
          "N3, N4.12, M", okq)
    M = NETS["M"]
    AB = net("A -> B; B -> A", ["A", "B"])
    claim("S7 M and A <-> B have the same Gamma, m, l, delta, but (m*, l*, s*, delta*) = (4,2,1,1) for M* and "
          "(2,1,1,0) for (A <-> B)*", M.Gamma == AB.Gamma and counts(M) == counts(AB) and
          counts(M.dual()) == (4, 2, 1, 1) and counts(AB.dual()) == (2, 1, 1, 0))
    # five-species example
    R5 = []
    for mult in (1, 2):
        for i in range(5):
            a = [0] * 5
            b = [0] * 5
            a[i] += mult
            a[(i + 2) % 5] += mult
            b[(i + 1) % 5] += mult
            b[(i + 3) % 5] += mult
            R5.append((tuple(a), tuple(b)))
    N5 = crn.Network(5, R5)
    N5d = N5.dual()
    G5 = N5.Gamma
    Sig = (0, 1, 2)
    claim("S7 five-species example (X_i+X_{i+2} -> X_{i+1}+X_{i+3} and twice it, indices mod 5): 10 reactions, "
          "(m,l,s,delta) = (10,2,4,4), q = 1, WR, (H1)-(H3); the dual is one 5-cycle with deficiency 0; "
          "{X1,X2,X3} is a critical siphon, realised at (0,0,0,2,3) in the class of (1,1,1,1,1), with face "
          "codimension rank Gamma[Sigma,:] = 3; no blind reaction and dim S0 = 1",
          len(R5) == 10 and counts(N5) == (10, 2, 4, 4) and is_wr(N5) and H_all(N5) and is_wr(N5d) and
          N5d.linkage_classes() == 1 and N5d.deficiency() == 0 and len(N5d.complexes()) == 5 and
          N5.is_siphon(Sig) and not law_in(N5, set(Sig)) and sum((0, 0, 0, 2, 3)) == 5 and
          crn.rank([G5[i] for i in Sig]) == 3 and
          [j for j, (a, b) in enumerate(R5) if all(a[i] == 0 for i in Sig)] == [] and
          4 - crn.rank([G5[i] for i in Sig]) == 1)
    R4 = net("X1 -> X2; X2 -> X3; X3 -> X4; X4 -> X1; X1+X3 -> X2+X4; X2+X4 -> X1+X3", X4)
    w = (3, 2, 1, 2)
    hts = [sum(w[i] * a[i] for i in range(4)) for a, b in R4.reactions]
    dl = [sum(w[i] * (b[i] - a[i]) for i in range(4)) for a, b in R4.reactions]
    claim("S7 four-species example A -> B -> C -> D -> A, A+C <-> B+D: q = 1, (H3), WR with WR dual; for "
          "w = (3,2,1,2) the maximal source height 4 carries no strictly inward reaction (not strongly endotactic)",
          R4.n - R4.rank() == 1 and R4.rows_distinct() and is_wr(R4) and is_wr(R4.dual()) and max(hts) == 4 and
          not any(d_ < 0 and h == 4 for h, d_ in zip(hts, dl)), f"heights {hts}, displacements {dl}")
    # Example 8.1
    E8 = net("3A+2C -> 2A+3B; 2A+3B -> 2B+3C; 2B+3C -> 3A+2C", ["A", "B", "C"])
    G8 = E8.Gamma
    alt = net("10A+2C -> 9A+3B; 2A+B -> 3C; 10B+C -> 3A+8B", ["A", "B", "C"])
    xa, xb, xc = V("xa", "xb", "xc")
    g = [xa ** 10 * xc ** 2, xa ** 2 * xb, xb ** 10 * xc]
    Gd = [[gi.diff(v_).value({"xa": 1, "xb": 1, "xc": 1}) for v_ in ("xa", "xb", "xc")] for gi in g]
    J = matmul(G8, Gd)
    B = [[1, 0], [0, 1], [-1, -1]]
    JB = matmul(J, B)
    JS = [JB[0], JB[1]]  # rows 1, 2 of B form the identity
    okS = all(JB[i][j] == sum(B[i][k] * JS[k][j] for k in range(2)) for i in range(3) for j in range(2))
    claim("S7 Example (structural stability): Gamma = [[-1,-2,3],[3,-1,-2],[-2,3,-1]], rank 2, deficiency 0, "
          "WR with WR dual of deficiency 0; (Ye)_ij > 0 iff Gamma_ij < 0; the realisation 10A+2C -> 9A+3B, "
          "2A+B -> 3C, 10B+C -> 3A+8B has the same Gamma; G = [[10,0,2],[2,1,0],[0,10,1]]; "
          "J_S = [[-15,27],[24,-25]], det -273",
          G8 == [[-1, -2, 3], [3, -1, -2], [-2, 3, -1]] and counts(E8) == (3, 1, 2, 0) and is_wr(E8) and
          is_wr(E8.dual()) and E8.dual().deficiency() == 0 and
          all((E8.Ye[i][j] > 0) == (G8[i][j] < 0) for i in range(3) for j in range(3)) and alt.Gamma == G8 and
          Gd == [[10, 0, 2], [2, 1, 0], [0, 10, 1]] and okS and JS == [[-15, 27], [24, -25]] and det(JS) == -273)


# =============================================================================================================
# 12. Main
# =============================================================================================================
def main():
    t0 = time.time()
    sec2()
    sec3_archive()
    sec3_phat()
    sec3_siphons()
    sec4()
    sec5()
    sec6()
    sec6_A()
    sec7()
    dt = time.time() - t0
    head = [
        "# certify_duality.py -- certificates for the manuscript of result 17 (Iosif's duality)",
        f"# Python {platform.python_version()} ({platform.python_implementation()}) on {platform.system()} "
        f"{platform.machine()}; standard library only (fractions) + scripts/crn.py; exact arithmetic",
        f"# runtime {dt:.1f} s; claims {NPASS + NFAIL}: {NPASS} PASS, {NFAIL} FAIL",
    ]
    out = "\n".join(head + LINES + ["", "ALL CLAIMS PASS" if NFAIL == 0 else f"{NFAIL} CLAIM(S) FAILED"]) + "\n"
    logdir = os.path.join(os.path.dirname(HERE), "certificates")
    os.makedirs(logdir, exist_ok=True)
    with open(os.path.join(logdir, "certify_duality.log"), "w", encoding="utf-8", newline="\n") as fh:
        fh.write(out)
    sys.stdout.write(out)
    sys.exit(1 if NFAIL else 0)


if __name__ == "__main__":
    main()
