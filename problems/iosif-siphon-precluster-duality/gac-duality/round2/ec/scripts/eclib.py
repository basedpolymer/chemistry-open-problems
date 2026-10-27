"""Library for agent Q3 (round 2): exact objects of equal-conflict (EC) reaction networks.

Agent Q3, Claude Opus 5.5 (effort xhigh), 27 October 2026.
Environment: Python 3.12 (WSL Fedora 44 or Windows), standard library only, exact rationals (fractions.Fraction).
Imports (never edits) `problems/iosif-siphon-precluster-duality/scripts/crn.py` and agent P's `petri.py`.

Notation (round2/ec/notes.md): N EC means two reactions sharing a reactant species have the same source complex,
i.e. the distinct source complexes y_1..y_m have pairwise disjoint supports S_1..S_m.  T_c = reactions with source
y_c.  A_k = n x m matrix with columns a_c(k) = sum_{t in T_c} k_t (y'_t - y_c).  WF-EC = EC + connected + (H1) +
(H2) + s = m - 1 (every species consumed, which (H2) forces).
"""
import os
import sys
from fractions import Fraction
from itertools import combinations

HERE = os.path.dirname(os.path.abspath(__file__))
PETRI = os.path.normpath(os.path.join(HERE, "..", "..", "..", "petri", "scripts"))
CRN_DIR = os.path.normpath(os.path.join(HERE, "..", "..", "..", "..", "scripts"))
for p in (CRN_DIR, PETRI):
    if p not in sys.path:
        sys.path.insert(0, p)

from crn import Network, rank, nullspace, rref  # noqa: E402,F401
import petri  # noqa: E402


# ------------------------------------------------------------------------------------------------
# linear algebra helpers (exact)
# ------------------------------------------------------------------------------------------------
def matmul(A, B):
    return [[sum((A[i][k] * B[k][j] for k in range(len(B))), Fraction(0)) for j in range(len(B[0]))]
            for i in range(len(A))]


def transpose(A):
    return [list(r) for r in zip(*A)]


def det(M):
    M = [[Fraction(x) for x in row] for row in M]
    n = len(M)
    d = Fraction(1)
    for c in range(n):
        p = next((r for r in range(c, n) if M[r][c] != 0), None)
        if p is None:
            return Fraction(0)
        if p != c:
            M[c], M[p] = M[p], M[c]
            d = -d
        d *= M[c][c]
        for r in range(c + 1, n):
            f = M[r][c] / M[c][c]
            if f:
                for k in range(c, n):
                    M[r][k] -= f * M[c][k]
    return d


def charpoly(A):
    """Coefficients c_0..c_n of det(tI - A) (Faddeev-LeVerrier, exact)."""
    n = len(A)
    c = [Fraction(0)] * (n + 1)
    c[n] = Fraction(1)
    M = [[Fraction(0)] * n for _ in range(n)]
    for k in range(1, n + 1):
        AM = matmul(A, M) if k > 1 else [[Fraction(0)] * n for _ in range(n)]
        M = [[AM[i][j] + (c[n - k + 1] if i == j else 0) for j in range(n)] for i in range(n)]
        AM = matmul(A, M)
        c[n - k] = -sum(AM[i][i] for i in range(n)) / k
    return c


def is_hurwitz(coeffs):
    """coeffs = b_0..b_d (b_d != 0).  True iff all roots have negative real part (Routh array, exact; any zero in
    the first column counts as 'not strictly Hurwitz')."""
    b = [Fraction(x) for x in coeffs]
    if b[-1] < 0:
        b = [-x for x in b]
    d = len(b) - 1
    if d == 0:
        return True
    if any(x <= 0 for x in b):
        return False
    hi = b[::-1]  # b_d, b_{d-1}, ..., b_0
    r0 = hi[0::2]
    r1 = hi[1::2]
    rows = [r0, r1]
    for _ in range(d - 1):
        a, bb = rows[-2], rows[-1]
        if not bb or bb[0] == 0:
            return False
        new = []
        for i in range(len(a) - 1):
            x = a[i + 1]
            y = bb[i + 1] if i + 1 < len(bb) else Fraction(0)
            new.append((bb[0] * x - a[0] * y) / bb[0])
        while new and new[-1] == 0 and len(new) > 1:
            new.pop()
        rows.append(new if new else [Fraction(0)])
    first = [r[0] if r else Fraction(0) for r in rows[:d + 1]]
    return all(x > 0 for x in first)


# ------------------------------------------------------------------------------------------------
# EC structure
# ------------------------------------------------------------------------------------------------
class ECStructure:
    def __init__(self, N):
        self.N = N
        srcs = []
        for a, _ in N.reactions:
            if a not in srcs:
                srcs.append(a)
        self.sources = srcs                                   # y_c
        self.m = len(srcs)
        self.T = [[t for t in range(N.r) if N.reactions[t][0] == y] for y in srcs]
        self.S = [frozenset(i for i in range(N.n) if y[i] > 0) for y in srcs]
        self.cl_of_reaction = [srcs.index(N.reactions[t][0]) for t in range(N.r)]
        self.cl_of_species = [None] * N.n
        for c, Sc in enumerate(self.S):
            for i in Sc:
                self.cl_of_species[i] = c

    def is_ec(self):
        if any(not any(y) for y in self.sources):
            return False
        for c, d in combinations(range(self.m), 2):
            if self.S[c] & self.S[d]:
                return False
        return True

    def every_species_consumed(self):
        return all(c is not None for c in self.cl_of_species)

    def A(self, k):
        N = self.N
        A = [[Fraction(0)] * self.m for _ in range(N.n)]
        for t, (a, b) in enumerate(N.reactions):
            c = self.cl_of_reaction[t]
            for i in range(N.n):
                A[i][c] += Fraction(k[t]) * (b[i] - a[i])
        return A

    def cluster_successor_graph_strong(self):
        """c -> c' if some reaction of T_c produces a species of S_c'."""
        N = self.N
        succ = {c: set() for c in range(self.m)}
        for t, (a, b) in enumerate(N.reactions):
            c = self.cl_of_reaction[t]
            for i in range(N.n):
                if b[i] > 0 and self.cl_of_species[i] is not None:
                    succ[c].add(self.cl_of_species[i])
        for start in range(self.m):
            seen, stack = {start}, [start]
            while stack:
                x = stack.pop()
                for y in succ[x]:
                    if y not in seen:
                        seen.add(y)
                        stack.append(y)
            if len(seen) != self.m:
                return False
        return True


def positive_flux(N):
    rays = N.flux_rays()
    v = [Fraction(0)] * N.r
    for ray in rays:
        for j in range(N.r):
            v[j] += ray[j]
    return v


def positive_conservation(N):
    rays = N.conservation_rays()
    z = [Fraction(0)] * N.n
    for ray in rays:
        for i in range(N.n):
            z[i] += ray[i]
    return z


def is_wf_ec(N):
    E = ECStructure(N)
    if not E.is_ec() or not E.every_species_consumed():
        return False
    return (N.rank() == E.m - 1 and petri.is_connected(N) and N.is_consistent() and N.is_conservative())


def is_pure(N):
    return all(not (a[i] > 0 and b[i] > 0) for a, b in N.reactions for i in range(N.n))


def is_balanced(N):
    return petri.is_balanced_joins(N)


# ------------------------------------------------------------------------------------------------
# Theorem C (no purity): constructive certificate
# ------------------------------------------------------------------------------------------------
def theorem_c_certificate(N, E=None):
    """For each minimal siphon R: (i) R meets each S_c at most once; (ii) B = A_v[R, C_R] is Metzler with zero row
    sums; (iii) the left kernel of B contains z >= 0, z != 0 (computed); (iv) z extended by 0 is a conservation
    law of N with support exactly R.  Returns (ok, message)."""
    E = E or ECStructure(N)
    v = positive_flux(N)
    A = E.A(v)
    G = N.Gamma
    for R in N.minimal_siphons():
        R = sorted(R)
        cls = [E.cl_of_species[i] for i in R]
        if len(set(cls)) != len(cls):
            return False, f"siphon {R} meets a cluster twice"
        B = [[A[i][c] for c in cls] for i in R]
        k = len(R)
        if any(B[a][b] < 0 for a in range(k) for b in range(k) if a != b):
            return False, f"B not Metzler for {R}"
        if any(sum(B[a]) != 0 for a in range(k)):
            return False, f"B row sums nonzero for {R}"
        # columns of A outside C_R vanish on the rows R
        others = [c for c in range(E.m) if c not in cls]
        if any(A[i][c] != 0 for i in R for c in others):
            return False, f"A[R, outside] != 0 for {R}"
        ker = nullspace(transpose(B), k)
        cand = None
        for w in ker:
            if all(x >= 0 for x in w) or all(x <= 0 for x in w):
                cand = [abs(x) for x in w]
                break
        if cand is None or not any(cand):
            return False, f"no nonnegative left kernel vector for {R} (dim {len(ker)})"
        z = [Fraction(0)] * N.n
        for i, x in zip(R, cand):
            z[i] = x
        for j in range(N.r):
            if sum(z[i] * G[i][j] for i in range(N.n)) != 0:
                return False, f"z not a conservation law for {R}"
        if {i for i in range(N.n) if z[i] != 0} != set(R):
            return False, f"support of z != R for {R}"
    return True, ""


# ------------------------------------------------------------------------------------------------
# sign condition sigma(S) cap sigma(D^perp) = {0}
# ------------------------------------------------------------------------------------------------
def sign_condition_minors(N, E=None):
    """Mueller-Feliu-Regensburger-Conradi-Shiu-Dickenstein (2016) criterion: for bases M of S = im Gamma and Mt of
    D = span{y_c - y_c0}, sigma(S) cap sigma(D^perp) = {0} iff the products det(M_I) det(Mt_I) over s-subsets I of
    rows are all >= 0 or all <= 0, and not all zero.  Returns (holds, products)."""
    E = E or ECStructure(N)
    s = N.rank()
    G = N.Gamma
    Mcols = independent_columns(G)
    y0 = E.sources[0]
    Dcols = [[Fraction(E.sources[c][i] - y0[i]) for i in range(N.n)] for c in range(1, E.m)]
    assert len(Mcols) == s and len(Dcols) == s
    prods = []
    for I in combinations(range(N.n), s):
        a = det([[Mcols[j][i] for j in range(s)] for i in I])
        b = det([[Dcols[j][i] for j in range(s)] for i in I])
        prods.append(a * b)
    pos = any(p > 0 for p in prods)
    neg = any(p < 0 for p in prods)
    return (pos != neg), prods


def independent_columns(M):
    cols = [list(map(Fraction, c)) for c in zip(*M)]
    chosen = []
    for c in cols:
        if rank(chosen + [c]) > len(chosen):
            chosen.append(c)
    return chosen


def _fm_feasible(cons, d):
    """cons: list of (a, b) meaning a.t >= b.  Exact Fourier-Motzkin feasibility (copied from agent E's endo.py,
    computation/scripts/endo.py, with no change of logic)."""
    cons = [([Fraction(x) for x in a], Fraction(b)) for a, b in cons]
    for v in range(d):
        pos = [c for c in cons if c[0][v] > 0]
        neg = [c for c in cons if c[0][v] < 0]
        new = [c for c in cons if c[0][v] == 0]
        for ap, bp in pos:
            for an, bn in neg:
                fp, fn = ap[v], -an[v]
                a = [x / fp + y / fn for x, y in zip(ap, an)]
                a[v] = Fraction(0)
                new.append((a, bp / fp + bn / fn))
        seen = set()
        cons = []
        for a, b in new:
            key = (tuple(a), b)
            if key not in seen:
                seen.add(key)
                cons.append((a, b))
    return all(b <= 0 for a, b in cons)


def sign_vector_in_span(basis_cols, tau):
    """tau in {-1,0,1}^n is the sign vector of some vector of span(basis_cols) (columns, list of R^n vectors)."""
    n = len(tau)
    d = len(basis_cols)
    if d == 0:
        return all(t == 0 for t in tau)
    cons = []
    for i in range(n):
        row = [basis_cols[j][i] for j in range(d)]
        if tau[i] > 0:
            cons.append((row, 1))
        elif tau[i] < 0:
            cons.append(([-x for x in row], 1))
        else:
            cons.append((row, 0))
            cons.append(([-x for x in row], 0))
    return _fm_feasible(cons, d)


def sign_condition_direct(N, E=None):
    """Direct check: no nonzero tau in sigma(S) cap sigma(D^perp) (enumerates tau; small n only)."""
    from itertools import product
    E = E or ECStructure(N)
    Mcols = independent_columns(N.Gamma)
    y0 = E.sources[0]
    Drows = [[Fraction(E.sources[c][i] - y0[i]) for i in range(N.n)] for c in range(1, E.m)]
    Kperp = nullspace(Drows, N.n) if Drows else [[Fraction(int(i == j)) for i in range(N.n)] for j in range(N.n)]
    for tau in product((-1, 0, 1), repeat=N.n):
        if not any(tau):
            continue
        # symmetry: tau and -tau together
        first = next(t for t in tau if t)
        if first < 0:
            continue
        if sign_vector_in_span(Mcols, tau) and sign_vector_in_span(Kperp, tau):
            return False, tau
    return True, None


# ------------------------------------------------------------------------------------------------
# Jacobian at a constructed steady state
# ------------------------------------------------------------------------------------------------
def monomial(x, y):
    p = Fraction(1)
    for xi, yi in zip(x, y):
        p *= Fraction(xi) ** yi
    return p


def rates_making_steady(N, x, v=None):
    """k_t = v_t / x^{y_t} with v > 0 in ker Gamma: then x is a steady state for k."""
    v = v or positive_flux(N)
    return [v[t] / monomial(x, N.reactions[t][0]) for t in range(N.r)]


def jacobian(N, k, x):
    n = N.n
    J = [[Fraction(0)] * n for _ in range(n)]
    for t, (a, b) in enumerate(N.reactions):
        rate = k[t] * monomial(x, a)
        for p in range(n):
            if a[p] == 0:
                continue
            d = rate * a[p] / Fraction(x[p])      # d rate / d x_p
            for i in range(n):
                J[i][p] += (b[i] - a[i]) * d
    return J


def stability_certificate(N, x, E=None):
    """x > 0 rational.  Builds k making x a steady state, checks: G = Y_s^T X^{-1} A_k diag(phi) is Metzler with zero
    row sums and irreducible; det(tI - J) = t^{n-s} q(t) with q(0) != 0 and q Hurwitz.  Returns (ok, msg)."""
    E = E or ECStructure(N)
    k = rates_making_steady(N, x)
    A = E.A(k)
    phi = [monomial(x, y) for y in E.sources]
    f = [sum(A[i][c] * phi[c] for c in range(E.m)) for i in range(N.n)]
    if any(f):
        return False, "x not a steady state"
    G = [[sum(Fraction(E.sources[c][p]) / Fraction(x[p]) * A[p][d] * phi[d] for p in range(N.n))
          for d in range(E.m)] for c in range(E.m)]
    if any(G[c][d] < 0 for c in range(E.m) for d in range(E.m) if c != d):
        return False, "G not Metzler"
    if any(sum(G[c]) != 0 for c in range(E.m)):
        return False, "G row sums nonzero"
    J = jacobian(N, k, x)
    cp = charpoly(J)
    s = N.rank()
    low = N.n - s
    if any(cp[i] != 0 for i in range(low)):
        return False, f"charpoly has t^j, j < n - s, coefficient nonzero: {cp}"
    q = cp[low:]
    if q[0] == 0:
        return False, "degenerate steady state (q(0) = 0)"
    if not is_hurwitz(q):
        return False, f"q not Hurwitz: {q}"
    return True, ""
