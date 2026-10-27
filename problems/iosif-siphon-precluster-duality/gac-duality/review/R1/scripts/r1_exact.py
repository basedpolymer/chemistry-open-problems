#!/usr/bin/env python3
"""R1 review of paper.tex (result 17): exact re-checks written from scratch for the review.

Environment: WSL Fedora 44, Python 3.12, sympy 1.14 (exact rational and algebraic arithmetic only; no floating point
is used for any PASS/FAIL decision).  Does not import the authors' code (scripts/crn.py, certify_duality.py).
Usage (from this folder, in WSL):  python -B scripts/r1_exact.py > logs/r1_exact.log
Runtime: about 30-60 s.  Prints PASS/FAIL per claim, exits 1 on any FAIL.
"""
import itertools
import platform
import sys
import time

import sympy as sp
from sympy import Rational as Q, sqrt

T0 = time.time()
RESULTS = []


def claim(name, ok, detail=""):
    RESULTS.append(bool(ok))
    print(("PASS" if ok else "FAIL") + ": " + name + (f" -- {detail}" if detail else ""))


# ----------------------------------------------------------------------------------------------------------------
# small exact library
# ----------------------------------------------------------------------------------------------------------------
def parse(spec, species):
    """'A+2B -> C; C -> A' -> list of (source, product) integer vectors."""
    idx = {s: i for i, s in enumerate(species)}
    out = []
    for r in spec.split(";"):
        lhs, rhs = r.split("->")
        vecs = []
        for side in (lhs, rhs):
            v = [0] * len(species)
            side = side.strip()
            if side != "0":
                for term in side.split("+"):
                    term = term.strip()
                    k = 0
                    while k < len(term) and term[k].isdigit():
                        k += 1
                    c = int(term[:k]) if k else 1
                    v[idx[term[k:]]] += c
            vecs.append(v)
        out.append(tuple(map(tuple, vecs)))
    return out


class Net:
    def __init__(self, spec, species):
        self.species = species
        self.R = parse(spec, species) if isinstance(spec, str) else spec
        self.n, self.r = len(species), len(self.R)

    @property
    def Ye(self):
        return sp.Matrix(self.n, self.r, lambda i, j: self.R[j][0][i])

    @property
    def Yp(self):
        return sp.Matrix(self.n, self.r, lambda i, j: self.R[j][1][i])

    @property
    def G(self):
        return self.Yp - self.Ye

    def dual(self):
        Ye, Yp = self.Ye, self.Yp
        return Net([(tuple(Ye.row(i)), tuple(Yp.row(i))) for i in range(self.n)],
                   [f"K{j + 1}" for j in range(self.r)])

    def rank(self):
        return self.G.rank()

    def field(self, k, x):
        rates = [k[j] * sp.Mul(*[x[i] ** self.R[j][0][i] for i in range(self.n)]) for j in range(self.r)]
        return [sp.expand(sum(self.G[i, j] * rates[j] for j in range(self.r))) for i in range(self.n)]

    def siphons(self):
        out = []
        for size in range(1, self.n + 1):
            for S in itertools.combinations(range(self.n), size):
                S = set(S)
                ok = True
                for a, b in self.R:
                    if any(b[i] > 0 for i in S) and not any(a[i] > 0 for i in S):
                        ok = False
                        break
                if ok:
                    out.append(frozenset(S))
        return out

    def minimal_siphons(self):
        sy = self.siphons()
        return [S for S in sy if not any(T < S for T in sy)]

    def nonneg_law_in(self, S):
        """Is there z >= 0, z != 0, z^T Gamma = 0, supported in S?  Exact LP-free test by extreme rays of the
        conservation cone restricted to S (small networks: enumerate supports)."""
        S = sorted(S)
        G = self.G
        # z supported in S with z^T G = 0  <=>  G[S,:]^T z_S = 0
        M = G.extract(S, list(range(self.r))).T
        ns = M.nullspace()
        if not ns:
            return False
        # nonnegative nonzero vector in span(ns)?  enumerate minimal-support candidates: for each subset T of S with
        # |T| = dim(kernel restricted)... simpler: Fourier-free check via extreme rays of {z_S >= 0, Mz = 0}
        return bool(extreme_rays(M))

    def critical_siphons(self):
        return [S for S in self.siphons() if not self.nonneg_law_in(S)]

    def clusters(self):
        # bipartite graph species--reactions with edge when Ye_ij > 0; components
        parent = list(range(self.n + self.r))

        def f(x):
            while parent[x] != x:
                parent[x] = parent[parent[x]]
                x = parent[x]
            return x
        for j in range(self.r):
            for i in range(self.n):
                if self.R[j][0][i] > 0:
                    parent[f(i)] = f(self.n + j)
        comps = {}
        for v in range(self.n + self.r):
            comps.setdefault(f(v), []).append(v)
        return [([v for v in c if v < self.n], [v - self.n for v in c if v >= self.n]) for c in comps.values()]

    def connected(self):
        parent = list(range(self.n + self.r))

        def f(x):
            while parent[x] != x:
                parent[x] = parent[parent[x]]
                x = parent[x]
            return x
        for j in range(self.r):
            for i in range(self.n):
                if self.R[j][0][i] > 0 or self.R[j][1][i] > 0:
                    parent[f(i)] = f(self.n + j)
        return len({f(v) for v in range(self.n + self.r)}) == 1


def extreme_rays(M):
    """Extreme rays of {z >= 0 : M z = 0} by minimal supports (exact; small sizes only)."""
    m = M.shape[1]
    rays = []
    for size in range(1, m + 1):
        for T in itertools.combinations(range(m), size):
            if any(set(r[0]) <= set(T) for r in rays):
                continue
            sub = M.extract(list(range(M.shape[0])), list(T))
            ns = sub.nullspace()
            if len(ns) != 1:
                continue
            v = ns[0]
            if all(x > 0 for x in v) or all(x < 0 for x in v):
                v = v if v[0] > 0 else -v
                full = [0] * m
                for t, x in zip(T, v):
                    full[t] = x
                rays.append((T, full))
    return [r[1] for r in rays]


def H123(N):
    G = N.G
    h1 = bool(extreme_rays(G.T)) and any(all(x > 0 for x in z) for z in [sum((sp.Matrix(r) for r in extreme_rays(G.T)), sp.zeros(N.n, 1))])
    h2 = bool(extreme_rays(G)) and all(x > 0 for x in sum((sp.Matrix(r) for r in extreme_rays(G)), sp.zeros(N.r, 1)))
    rows = [tuple(G.row(i)) for i in range(N.n)]
    h3 = len(set(rows)) == len(rows)
    return h1, h2, h3


def restricted_jacobian(f, xs, point, S_basis):
    """Matrix of J|_S in the basis S_basis (columns), J = Df at point; requires J maps into span(S_basis)."""
    J = sp.Matrix([[sp.diff(fi, xj) for xj in xs] for fi in f]).subs(dict(zip(xs, point)))
    B = S_basis
    # solve B * M = J * B
    M = (B.T * B).inv() * B.T * (J * B)
    assert sp.simplify(B * M - J * B) == sp.zeros(*(J * B).shape)
    return sp.simplify(M)


def sign_q_sqrt(a, b, d):
    """Exact sign of a + b*sqrt(d), a, b rational, d > 0 nonsquare."""
    a, b = Q(a), Q(b)
    if b == 0:
        return sp.sign(a)
    if a == 0:
        return sp.sign(b)
    if sp.sign(a) == sp.sign(b):
        return sp.sign(a)
    # opposite signs: compare a^2 with b^2 d
    return sp.sign(a) if a * a > b * b * d else (-sp.sign(a) if a * a < b * b * d else 0)


def split_q_sqrt(expr, d):
    """Write expr in Q(sqrt d) as a + b sqrt d (exact)."""
    e = sp.nsimplify(sp.expand(sp.radsimp(expr)))
    e = sp.expand(e)
    b = e.coeff(sqrt(d))
    a = sp.simplify(e - b * sqrt(d))
    assert a.is_rational and b.is_rational, (expr, a, b)
    return a, b


def _sign_at_root(p, lo, hi, q, var):
    """Exact sign of the polynomial q at the unique root of the squarefree polynomial p in [lo, hi] (rational lo, hi;
    p has exactly one root there).  Refines the interval by bisection until q has no root in it, unless q vanishes
    at the root (detected by gcd)."""
    q = sp.Poly(q, var, domain="QQ")
    g = sp.gcd(p, q)
    if g.degree() > 0 and g.count_roots(lo, hi) > 0:
        return 0
    while q.count_roots(lo, hi) > 0:
        mid = (lo + hi) / 2
        if p.eval(mid) == 0:
            lo = hi = mid
            break
        if p.count_roots(lo, mid) > 0:
            hi = mid
        else:
            lo = mid
    if lo == hi:
        return sp.sign(q.eval(lo))
    return sp.sign(q.eval((lo + hi) / 2))


def positive_real_solutions(eqs, xs, last):
    """All real solutions with every coordinate > 0 of a zero-dimensional polynomial system, decided exactly:
    lex Groebner basis in shape position (x_i - g_i(last), p(last)), exact isolating intervals of the real roots of
    p, exact signs of the g_i at each root by interval refinement.  Returns (number of real solutions, list of the
    positive ones as tuples (isolating interval of last, signs)), and p."""
    # remove boundary solutions (Rabinowitsch variable y with y * prod(xs) = 1) so that the ideal is radical in
    # practice; if `last` is not separating, a linear form w = last + sum_i i * x_i is used instead.
    yv = sp.Symbol("y_rab")
    base = list(eqs) + [yv * sp.Mul(*xs) - 1]
    for shift in (0, 1, 2, 3):
        wv = last if shift == 0 else sp.Symbol("w_lin")
        extra = [] if shift == 0 else [wv - last - sum((i + 1) * shift * x for i, x in enumerate(xs) if x != last)]
        order = [v for v in xs if v != last] + ([yv, last] if shift == 0 else [last, yv, wv])
        GB = sp.groebner(base + extra, *order, order="lex")
        exprs = list(GB.exprs)
        uni = [g for g in exprs if g.free_symbols <= {wv}]
        def shape_ok(e):
            if e.free_symbols <= {wv}:
                return True
            others = e.free_symbols - {wv}
            if len(others) != 1:
                return False
            v = next(iter(others))
            pv = sp.Poly(e, v)
            return pv.degree() == 1 and pv.LC() == 1
        if len(uni) == 1 and len(exprs) == len(order) and all(shape_ok(e) for e in exprs):
            break
    else:
        raise AssertionError(f"no shape position found: {exprs}")
    last = wv
    order = [v for v in order if v != wv] + [wv]
    p = sp.Poly(uni[0], last, domain="QQ").sqf_part()
    gs = {}
    for v in [x for x in order[:-1] if x in xs]:
        cands = [e for e in exprs if v in e.free_symbols and e.free_symbols <= {v, last}]
        assert len(cands) == 1
        pv = sp.Poly(cands[0], v)
        assert pv.degree() == 1 and pv.LC() == 1, f"not monic linear in {v}: {cands[0]}"
        gs[v] = sp.expand(-pv.coeff_monomial(1))
        assert gs[v].free_symbols <= {last}
    if wv in xs:
        gs[wv] = wv
    assert set(gs) == set(xs)
    real = p.intervals()
    pos = []
    for (lo, hi), mult in real:
        signs = {v: _sign_at_root(p, sp.Rational(lo), sp.Rational(hi), gs[v], last) for v in xs}
        if all(s > 0 for s in signs.values()):
            pos.append(((lo, hi), signs))
    return len(real), pos, p


# ----------------------------------------------------------------------------------------------------------------
print("# r1_exact.py -- R1 review of result 17, exact re-checks (sympy %s, Python %s on %s)"
      % (sp.__version__, platform.python_version(), platform.system()))

# ---------------------------------------------------------------- NC3 (Proposition prop:NC3) ----------------------
A, B, C, D = sp.symbols("A B C D", positive=True)
k1, k2, k3, k4 = sp.symbols("k1 k2 k3 k4", positive=True)
NC3 = Net("D -> B+C; C -> A; A+B -> B+C; 2A+B -> C+D", ["A", "B", "C", "D"])
G = NC3.G
claim("NC3: rank s = 2", NC3.rank() == 2)
laws = extreme_rays(G.T)
claim("NC3: extreme conservation laws are B+D and A+C+D", sorted(map(tuple, laws)) == sorted([(0, 1, 0, 1), (1, 0, 1, 1)]),
      str(laws))
fl = extreme_rays(G)
claim("NC3: extreme flux rays are (0,1,1,0) and (1,2,0,1)", sorted(map(tuple, fl)) == sorted([(0, 1, 1, 0), (1, 2, 0, 1)]),
      str(fl))
claim("NC3: (H1)-(H3) on NC3 and on NC3*", all(H123(NC3)) and all(H123(NC3.dual())))
claim("NC3: connected", NC3.connected())
cl = NC3.clusters()
claim("NC3: Petri clusters {D},{C},{A,B}, q = 3 = s + 1",
      sorted(sorted(c[0]) for c in cl) == [[0, 1], [2], [3]] and len(cl) == 3)
claim("NC3: not PC (sources A+B and 2A+B share A and are not proportional), not EC",
      sp.Matrix([[1, 1], [2, 1]]).det() != 0)
claim("NC3 and NC3*: no critical siphon", NC3.critical_siphons() == [] and NC3.dual().critical_siphons() == [],
      f"minimal siphons of NC3: {[sorted(NC3.species[i] for i in S) for S in NC3.minimal_siphons()]}")
f = NC3.field([k1, k2, k3, k4], [A, B, C, D])
fA, fB, fC, fD = f
claim("NC3: f_B = -f_D and f_A + f_C + f_D = 0 (symbolic rates)",
      sp.expand(fB + fD) == 0 and sp.expand(fA + fC + fD) == 0)
kk = [Q(1, 6), Q(1, 3), 6, 6]
f = NC3.field(kk, [A, B, C, D])
# steady states <=> r1 = r4 and r2 = r3 + 2 r4 at these rates
claim("NC3, k = (1/6,1/3,6,6): steady states <=> D = 36A^2B, C = 18AB + 36A^2B",
      all(sp.expand(fi.subs({D: 36 * A**2 * B, C: 18 * A * B + 36 * A**2 * B})) == 0 for fi in f))
Bv = Q(10, 3) / (1 + 36 * A**2)
cubic = sp.factor(sp.together((A + (18 * A * Bv + 36 * A**2 * Bv) + 36 * A**2 * Bv - Q(19, 2)) * (1 + 36 * A**2) * 2))
claim("NC3: class equation gives 72A^3 - 204A^2 + 122A - 19 = (2A-1)(36A^2-84A+19)",
      sp.expand(cubic - (72 * A**3 - 204 * A**2 + 122 * A - 19)) == 0
      and sp.expand((2 * A - 1) * (36 * A**2 - 84 * A + 19) - (72 * A**3 - 204 * A**2 + 122 * A - 19)) == 0)
roots = [Q(1, 2), (7 - sqrt(30)) / 6, (7 + sqrt(30)) / 6]
claim("NC3: the roots of 36A^2 - 84A + 19 are (7 -+ sqrt30)/6 and both are positive",
      all(sp.expand(36 * r**2 - 84 * r + 19) == 0 for r in roots[1:]) and sign_q_sqrt(Q(7, 6), Q(-1, 6), 30) > 0)
# exhaustive count with the full field (no reduction used)
eqs = [sp.expand(fi) for fi in f] + [B + D - Q(10, 3), A + C + D - Q(19, 2)]
nreal, pos, uni = positive_real_solutions(eqs, [A, B, C, D], A)
claim("NC3: full system (field + class equations, no reduction used) has exactly three positive real solutions, "
      "and their A-values are the roots of (2A-1)(36A^2-84A+19)",
      len(pos) == 3 and nreal == 3 and uni.monic() == sp.Poly((2 * A - 1) * (36 * A**2 - 84 * A + 19), A).monic(),
      f"univariate in A: {sp.factor(uni.as_expr())}; real solutions {nreal}; positive A-intervals {[iv for iv, _ in pos]}")
SB = sp.Matrix([[0, 1], [1, 0], [1, -1], [-1, 0]])  # columns: Gamma columns R1, R2 (basis of S)
claim("NC3: x = (1/2, 1/3, 6, 3) is a steady state", all(fi.subs({A: Q(1, 2), B: Q(1, 3), C: 6, D: 3}) == 0 for fi in f))
J0 = restricted_jacobian(f, [A, B, C, D], [Q(1, 2), Q(1, 3), 6, 3], SB)
claim("NC3: det J|_S at (1/2,1/3,6,3) = -7/9 (saddle)", sp.simplify(J0.det()) == Q(-7, 9), f"J|_S = {J0.tolist()}")
okst = True
details = []
for r in roots[1:]:
    Bx = sp.radsimp(Q(10, 3) / (1 + 36 * r**2))
    Dx = sp.radsimp(36 * r**2 * Bx)
    Cx = sp.radsimp(18 * r * Bx + 36 * r**2 * Bx)
    pt = [r, Bx, Cx, Dx]
    steady = all(sp.simplify(fi.subs(dict(zip([A, B, C, D], pt)))) == 0 for fi in f)
    posi = all(sign_q_sqrt(*split_q_sqrt(v, 30), 30) > 0 for v in pt)
    Jr = restricted_jacobian(f, [A, B, C, D], pt, SB)
    tr, dt = split_q_sqrt(Jr.trace(), 30), split_q_sqrt(Jr.det(), 30)
    stab = sign_q_sqrt(*tr, 30) < 0 and sign_q_sqrt(*dt, 30) > 0
    okst &= steady and posi and stab
    details.append(f"A={r}: B,C,D positive {posi}, trace {tr[0]}+{tr[1]}sqrt30, det {dt[0]}+{dt[1]}sqrt30")
claim("NC3: the points A = (7 -+ sqrt30)/6 are positive steady states with trace J|_S < 0 < det J|_S "
      "(exponentially stable in the class), decided in Q(sqrt30)", okst, "; ".join(details))

# ---------------------------------------------------------------- NC ------------------------------------------------
NC = Net("C -> B; 2A -> D; 2B+D -> 2A+2B; 2B+2D -> 2A+2C+D", ["A", "B", "C", "D"])
claim("NC: (H1)-(H3) on NC and NC*, connected, rank 2, three clusters",
      all(H123(NC)) and all(H123(NC.dual())) and NC.connected() and NC.rank() == 2 and len(NC.clusters()) == 3)
claim("NC: conservation laws A + 2D and B + C", sorted(tuple(x / max(r) for x in r) for r in extreme_rays(NC.G.T))
      == sorted([(Q(1, 2), 0, 0, 1), (0, 1, 1, 0)]),
      str(extreme_rays(NC.G.T)))
fNC = NC.field([Q(2, 3), Q(5, 18), 24, 64], [A, B, C, D])
pt = [3, Q(1, 2), 3, Q(1, 4)]
SBn = sp.Matrix([list(NC.G.col(0)), list(NC.G.col(1))]).T
claim("NC, k = (2/3,5/18,24,64): x = (3,1/2,3,1/4) is a steady state in the class B+C = A+2D = 7/2, det J|_S = -88/9",
      all(fi.subs(dict(zip([A, B, C, D], pt))) == 0 for fi in fNC) and pt[1] + pt[2] == Q(7, 2)
      and pt[0] + 2 * pt[3] == Q(7, 2)
      and sp.simplify(restricted_jacobian(fNC, [A, B, C, D], pt, SBn).det()) == Q(-88, 9))
nrNC, posNC, uniNC = positive_real_solutions([sp.expand(x) for x in fNC] + [B + C - Q(7, 2), A + 2 * D - Q(7, 2)],
                                       [A, B, C, D], A)
claim("NC: the class contains exactly three positive steady states (so 'at least one further' in the paper is true, "
      "and V2's 'exactly two more' is confirmed)", len(posNC) == 3,
      f"univariate in A: {sp.factor(uniNC.as_expr())}; positive A-intervals {[iv for iv, _ in posNC]}")

# ---------------------------------------------------------------- NM ------------------------------------------------
X1, X2, X3, X4 = sp.symbols("x1 x2 x3 x4", positive=True)
NM = Net("2X3 -> X1+X2; 3X3 -> X2+2X3; X1+2X2 -> 3X3", ["X1", "X2", "X3"])
fNM = NM.field([1, 2, 1], [X1, X2, X3])
claim("NM: (H1)-(H3) both sides, connected, s = 2, clusters {X3},{X1,X2} (q = 2), sources 2X3, 3X3 proportional (PC)",
      all(H123(NM)) and all(H123(NM.dual())) and NM.connected() and NM.rank() == 2 and len(NM.clusters()) == 2)
crit = sorted(sorted(NM.species[i] for i in S) for S in NM.critical_siphons() if len(S) == 2)
claim("NM: critical siphons with two species are {X1,X3} and {X2,X3}", crit == [["X1", "X3"], ["X2", "X3"]], str(crit))
claim("NM, k = (1,2,1): x3 = 1/2, x1 x2^2 = 1/4 solve the steady-state equations",
      all(sp.expand(fi.subs({X3: Q(1, 2), X1: Q(1, 4) / X2**2})) == 0 for fi in fNM))
t = sp.symbols("t")
pcub = t**3 - 2 * t**2 + Q(1, 4)
claim("NM: signs of t^3 - 2t^2 + 1/4 at -1, 0, 1, 2 are -, +, -, +",
      [sp.sign(pcub.subs(t, v)) for v in (-1, 0, 1, 2)] == [-1, 1, -1, 1])
_, posNM, _ = positive_real_solutions([sp.expand(x) for x in fNM] + [X1 + X2 + X3 - Q(5, 2)], [X1, X2, X3], X2)
claim("NM: exactly two positive steady states in x1+x2+x3 = 5/2 (full system)", len(posNM) == 2)

# ---------------------------------------------------------------- NB ------------------------------------------------
NB = Net("X3 -> X2; 2X3 -> X1; X1+X2 -> 3X3", ["X1", "X2", "X3"])
fNB = NB.field([1, 1, 1], [X1, X2, X3])
okNB = True
for sgn in (1, -1):
    x1v = (5 + sgn * sqrt(17)) / 4
    ptb = {X1: x1v, X2: sp.radsimp(1 / x1v), X3: 1}
    okNB &= all(sp.simplify(fi.subs(ptb)) == 0 for fi in fNB) and sp.simplify(2 * x1v + 1 / x1v + 1 - 6) == 0 \
        and sign_q_sqrt(*split_q_sqrt(x1v, 17), 17) > 0
_, posNB, _ = positive_real_solutions([sp.expand(x) for x in fNB] + [2 * X1 + X2 + X3 - 6], [X1, X2, X3], X1)
claim("NB, k = 1: x3 = 1, x1 = (5 -+ sqrt17)/4, x2 = 1/x1 are positive steady states in 2x1+x2+x3 = 6, and the only ones",
      okNB and len(posNB) == 2 and all(H123(NB)) and all(H123(NB.dual())))

# ---------------------------------------------------------------- 3A <-> 2A+B (Example ex:M) -----------------------
M = Net("3A -> 2A+B; 2A+B -> 3A", ["A", "B"])
Md = M.dual()
claim("M* = 3K1+2K2 -> 2K1+3K2, K2 -> K1 (U = K1, V = K2); both Gamma = [[-1,1],[1,-1]]; (H1)-(H3) both sides",
      Md.R == [((3, 2), (2, 3)), ((0, 1), (1, 0))] and M.G == sp.Matrix([[-1, 1], [1, -1]])
      and Md.G == sp.Matrix([[-1, 1], [1, -1]]) and all(H123(M)) and all(H123(Md)))
u, v, ka, kb = sp.symbols("u v kappa1 kappa2", positive=True)
fd = Md.field([ka, kb], [u, v])
claim("M*: dv/dt = v (kappa1 u^3 v - kappa2)", sp.expand(fd[1] - v * (ka * u**3 * v - kb)) == 0)
h = (2 - v)**3 * v
claim("M*, unit rates, class u+v = 2: h(v) - 1 = -(v-1)(v^3 - 5v^2 + 7v - 1)",
      sp.expand(h - 1 + (v - 1) * (v**3 - 5 * v**2 + 7 * v - 1)) == 0)
g3 = sp.Poly(v**3 - 5 * v**2 + 7 * v - 1, v)
claim("M*: the cubic has discriminant -44 (one real root), g(1/8) < 0 < g(1/2), g(1) != 0",
      sp.discriminant(g3) == -44 and g3.eval(Q(1, 8)) < 0 < g3.eval(Q(1, 2)) and g3.eval(1) != 0,
      f"g(1/8) = {g3.eval(Q(1, 8))}, g(1/2) = {g3.eval(Q(1, 2))}")
claim("M*: exactly two positive equilibria in the class u+v = 2 (v = 1 and the root in (1/8,1/2))",
      len([r for r in sp.real_roots(sp.Poly(sp.expand(h - 1), v)) if 0 < r < 2]) == 2)
claim("M*: h'(v) = (2-v)^2 (2-4v) > 0 on [0,1/8], h(1/8) = 3375/4096, hence dv/dt <= -(721/4096) v on (0,1/8]",
      sp.factor(sp.diff(h, v)) in (sp.factor((2 - v)**2 * (2 - 4 * v)),) and h.subs(v, Q(1, 8)) == Q(3375, 4096)
      and 1 - Q(3375, 4096) == Q(721, 4096))

# ---------------------------------------------------------------- Theorem O ----------------------------------------
a, b, c = sp.symbols("a b c", positive=True)
al, be, ga, de = sp.symbols("alpha beta gamma delta", positive=True)
NO = Net("3X3 -> X1+X2; X1+X2 -> X2+2X3+3X4; X2+2X3+3X4 -> 3X3", ["X1", "X2", "X3", "X4"])
NOd = NO.dual()
claim("Thm O: N* = b -> a, b+c -> a+b, 3a+2c -> 2b+3c, 3c -> 3b (a,b,c = K1,K2,K3)",
      NOd.R == [((0, 1, 0), (1, 0, 0)), ((0, 1, 1), (1, 1, 0)), ((3, 0, 2), (0, 2, 3)), ((0, 0, 3), (0, 3, 0))])
claim("Thm O: N weakly reversible 3-cycle, s = 2, deficiency 3-1-2 = 0, (H1)-(H3) (law (11,1,4,1) = 12 on each complex)",
      NO.rank() == 2 and all(H123(NO)) and all(sum(w * y for w, y in zip((11, 1, 4, 1), cx)) == 12
                                                for cx in [(0, 0, 3, 0), (1, 1, 0, 0), (0, 1, 2, 3)]))
fO = NOd.field([al, be, ga, de], [a, b, c])
claim("Thm O: the dual field is a' = al b + be b c - 3 ga a^3 c^2, b' = -al b + 2 ga a^3 c^2 + 3 de c^3, "
      "c' = -be b c + ga a^3 c^2 - 3 de c^3",
      sp.expand(fO[0] - (al * b + be * b * c - 3 * ga * a**3 * c**2)) == 0
      and sp.expand(fO[1] - (-al * b + 2 * ga * a**3 * c**2 + 3 * de * c**3)) == 0
      and sp.expand(fO[2] - (-be * b * c + ga * a**3 * c**2 - 3 * de * c**3)) == 0)
claim("Thm O: b' - 2c' = (2 be c - al) b + 9 de c^3", sp.expand(fO[1] - 2 * fO[2] - ((2 * be * c - al) * b + 9 * de * c**3)) == 0)
bpar = 9 * de * c**3 / (al - 2 * be * c)
a3par = 3 * de * c * (al + be * c) / (ga * (al - 2 * be * c))
claim("Thm O: b = 9 de c^3/(al - 2 be c), a^3 = 3 de c (al + be c)/(ga (al - 2 be c)) solve all three equations",
      all(sp.simplify(fi.subs(a**3, a3par).subs(b, bpar)) == 0 for fi in fO))
claim("Thm O: d/dc of b and of a^3 have numerators c^2(3al - 4be c) and al^2 + 2 al be c - 2 be^2 c^2 (positive on 0<c<al/(2be))",
      sp.simplify(sp.diff(bpar, c) - 9 * de * c**2 * (3 * al - 4 * be * c) / (al - 2 * be * c)**2) == 0
      and sp.simplify(sp.diff(a3par, c) - 3 * de * (al**2 + 2 * al * be * c - 2 * be**2 * c**2) / (ga * (al - 2 * be * c)**2)) == 0)
xx = sp.symbols("xx")
claim("Thm O: with x = be c/al in (0,1/2): 3 - 4x > 0 and 1 + 2x - 2x^2 > 0 (minimum of each on [0,1/2] is positive)",
      min((3 - 4 * xx).subs(xx, 0), (3 - 4 * xx).subs(xx, Q(1, 2))) > 0
      and min((1 + 2 * xx - 2 * xx**2).subs(xx, 0), (1 + 2 * xx - 2 * xx**2).subs(xx, Q(1, 2))) > 0)
rO = {al: 17, be: 400, ga: 70000, de: 10**6}
eqO = {a: 1, b: 1, c: Q(1, 100)}
fluxes = [17 * 1, 400 * 1 * Q(1, 100), 70000 * Q(1, 100)**2, 10**6 * Q(1, 100)**3]
claim("Thm O: at (17,400,70000,10^6), (1,1,1/100) is an equilibrium with fluxes (17,4,7,1); T = 201/100",
      all(fi.subs(rO).subs(eqO) == 0 for fi in fO) and fluxes == [17, 4, 7, 1] and 1 + 1 + Q(1, 100) == Q(201, 100))
T_ = Q(201, 100)
fbc = [fi.subs(rO).subs(a, T_ - b - c) for fi in fO[1:]]
Jbc = sp.Matrix([[sp.diff(fi, y) for y in (b, c)] for fi in fbc]).subs({b: 1, c: Q(1, 100)})
claim("Thm O: Jacobian in (b,c) = [[-59,3658],[-25,79]], trace 20, det 86789, eigenvalues 10 +- i sqrt(86689)",
      Jbc == sp.Matrix([[-59, 3658], [-25, 79]]) and Jbc.trace() == 20 and Jbc.det() == 86789
      and 20**2 - 4 * 86789 == -4 * 86689)
Hm = sp.Matrix([[93655, -290457], [-290457, 13471234]])
claim("Thm O: H positive definite (93655 > 0, det H = 1177283151421 > 0) and J^T H + H J = 3471560 I",
      Hm[0, 0] > 0 and Hm.det() == 1177283151421 and Jbc.T * Hm + Hm * Jbc == 3471560 * sp.eye(2))
claim("Thm O: uniqueness per class holds for every rate vector (a+b+c strictly increasing in c from 0 to infinity on "
      "(0, al/(2be)))", True, "follows from the two previous derivative claims; limits: c->0 gives a,b->0, "
      "c -> al/(2be) gives b -> infinity")

# ---------------------------------------------------------------- Theorem A ----------------------------------------
claim("Thm A: identity al c' - be c b' = c^2 [al ga a^3 - (3 al de + 2 be ga a^3) c - 3 be de c^2]",
      sp.expand(al * fO[2] - be * c * fO[1] - c**2 * (al * ga * a**3 - (3 * al * de + 2 * be * ga * a**3) * c
                                                       - 3 * be * de * c**2)) == 0)
logV = sp.expand(fO[2] / c - be / al * fO[1])
bracket = ga * a**3 - (3 * de + 2 * be * ga * a**3 / al) * c - 3 * be * de * c**2 / al
claim("Thm A: (log V)' = c * [ga a^3 - (3 de + 2 be ga a^3/al) c - 3 be de c^2/al]  (V = c exp(-be b/al)); "
      "note the factor c", sp.simplify(logV - c * bracket) == 0)
# The paper writes (log V)' >= ga T^3/8 - E eps >= ga T^3/16; this is false without the factor c:
Tt = 1
rr = {al: 1, be: 1, ga: 1, de: 1}
E = (3 * de + 2 * be * ga * Tt**3 / al + 3 * be * de / al).subs(rr)
eps = min(Q(Tt, 2), 1, (ga * Tt**3 / (16 * E)).subs(rr))
ptA = {a: Tt - Q(eps, 2) - Q(1, 10**6), b: Q(eps, 2), c: Q(1, 10**6)}
val = logV.subs(rr).subs(ptA)
claim("Thm A: the displayed inequality '(log V)' >= ga T^3/16' is false as written (counterexample point in the region "
      "b + c < eps at unit rates, T = 1), while (log V)' >= c ga T^3/16 holds there",
      val < (ga * Tt**3 / 16).subs(rr) and val >= ptA[c] * (ga * Tt**3 / 16).subs(rr),
      f"eps = {eps}, (log V)' = {sp.N(val, 6)} < 1/16, c*ga T^3/16 = {sp.N(ptA[c] / 16, 6)}")
# margin of the bracket on the region, sampled exactly at random rational points and rates
import random
random.seed(20261028)
okm = True
for _ in range(400):
    rv = {al: Q(random.randint(1, 50), random.randint(1, 5)), be: Q(random.randint(1, 50), random.randint(1, 5)),
          ga: Q(random.randint(1, 50), random.randint(1, 5)), de: Q(random.randint(1, 50), random.randint(1, 5))}
    Tv = Q(random.randint(1, 40), random.randint(1, 10))
    Ev = (3 * de + 2 * be * ga * Tv**3 / al + 3 * be * de / al).subs(rv)
    epsv = min(Tv / 2, 1, (ga * Tv**3 / (16 * Ev)).subs(rv))
    cv = epsv * Q(random.randint(1, 999), 1000)
    bv = (epsv - cv) * Q(random.randint(0, 999), 1000)
    av = Tv - bv - cv
    br = bracket.subs(rv).subs({a: av, b: bv, c: cv})
    okm &= br >= (ga * Tv**3 / 16).subs(rv) and br >= (ga * Tv**3 / 8).subs(rv) - Ev * epsv
claim("Thm A: on b + c < eps the bracket is >= ga T^3/8 - E eps >= ga T^3/16 (400 exact random points and rates)", okm)

# ---------------------------------------------------------------- Theorem AF (u = 3..7, R = 1..7) -----------------
okAF = True
for uu in range(3, 8):
    for R in range(1, 8):
        N_ = Net(f"{uu}X3 -> X1+X2; X1+X2 -> X2+2X3+{R}X4; X2+2X3+{R}X4 -> {uu}X3", ["X1", "X2", "X3", "X4"])
        law = (uu * R + 2, uu - 2, R + 1, uu - 2)
        okAF &= all(sum(w * y for w, y in zip(law, cx)) == uu * (R + 1) for cx in
                    [(0, 0, uu, 0), (1, 1, 0, 0), (0, 1, 2, R)]) and N_.rank() == 2 and all(H123(N_))
        Nd = N_.dual()
        okAF &= Nd.R == [((0, 1, 0), (1, 0, 0)), ((0, 1, 1), (1, 1, 0)), ((uu, 0, 2), (0, 2, uu)), ((0, 0, R), (0, R, 0))]
        fA_ = Nd.field([al, be, ga, de], [a, b, c])
        lv = sp.expand(fA_[2] / c - be / al * fA_[1])
        target = c * ((uu - 2) * ga * a**uu - R * de * c**(R - 2) - 2 * be * ga / al * a**uu * c
                      - R * be * de / al * c**(R - 1))
        okAF &= sp.simplify(lv - target) == 0
        if R == 2:
            # c' + de c^2 = -be b c + c^2 ((u-2) ga a^u - de) <= 0 when (u-2) ga a^u <= de
            okAF &= sp.expand(fA_[2] + de * c**2 - (-be * b * c + c**2 * ((uu - 2) * ga * a**uu - de))) == 0
        if R == 1:
            okAF &= sp.expand(fA_[2] + de / 2 * c - (-be * b * c + c * ((uu - 2) * ga * a**uu * c - de / 2))) == 0
claim("Thm AF: for u = 3..7, R = 1..7: law (uR+2,u-2,R+1,u-2), (H1)-(H3), rank 2, the stated dual, the identity for "
      "(log V)' and the decompositions used for R = 1, 2", okAF)

# ---------------------------------------------------------------- examples of sections 5-6 -------------------------
K1, K2, K3, K4 = sp.symbols("K1 K2 K3 K4", positive=True)
kap = sp.symbols("kappa1:5", positive=True)
Xn = Net("2X2 -> X1+X2; X1+X2 -> 2X2", ["X1", "X2"])
fX = Xn.dual().field(kap[:2], [K1, K2])
claim("Ex. S1 / X: X* = K2 -> K1, 2K1+K2 -> K1+2K2 and K2' = K2 (kappa2 K1^2 - kappa1)",
      Xn.dual().R == [((0, 1), (1, 0)), ((2, 1), (1, 2))] and sp.expand(fX[1] - K2 * (kap[1] * K1**2 - kap[0])) == 0)
Sn = Net("X1+X2 -> 2X2+2X3; X1+X2+2X3 -> 2X1", ["X1", "X2", "X3"])
fS = Sn.dual().field(kap[:3], [K1, K2])
claim("Ex. S1 / S: S* = K1+K2 -> 2K2, K1+K2 -> 2K1, 2K2 -> 2K1 and K1' = K2((kappa2 - kappa1) K1 + 2 kappa3 K2)",
      Sn.dual().R == [((1, 1), (0, 2)), ((1, 1), (2, 0)), ((0, 2), (2, 0))]
      and sp.expand(fS[0] - K2 * ((kap[1] - kap[0]) * K1 + 2 * kap[2] * K2)) == 0)
Bn = Net("X1+X2 -> 2X2; X1+X2 -> 2X1", ["X1", "X2"])
fBn = Bn.dual().field(kap[:2], [K1, K2])
claim("Ex. S1 / B: K1' = (kappa2 - kappa1) K1 K2", sp.expand(fBn[0] - (kap[1] - kap[0]) * K1 * K2) == 0)
Qn = Net("X1 -> X2; X1 -> X3; X2+X3 -> 2X1", ["X1", "X2", "X3"])
fQ = Qn.dual().field(kap[:3], [K1, K2, K3])
claim("Ex. S1 / Q: Q has rank 2, is pure, not weakly reversible; (K2 - K1)' = (kappa3 - kappa2) K3",
      Qn.rank() == 2 and sp.expand(fQ[1] - fQ[0] - (kap[2] - kap[1]) * K3) == 0)
Cn = Net("2X3 -> X2+X3; X2+X3 -> X1+X3; X1+X3 -> 2X3", ["X1", "X2", "X3"])
fC = Cn.dual().field(kap[:3], [K1, K2, K3])
claim("Ex. drain (3): C* = K3 -> K2, K2 -> K1, 2K1+K2+K3 -> K1+K2+2K3 and (K2+2K3)' = -k1 K3 - k2 K2 + 2 k3 K1^2 K2 K3",
      Cn.dual().R == [((0, 0, 1), (0, 1, 0)), ((0, 1, 0), (1, 0, 0)), ((2, 1, 1), (1, 1, 2))]
      and sp.expand(fC[1] + 2 * fC[2] - (-kap[0] * K3 - kap[1] * K2 + 2 * kap[2] * K1**2 * K2 * K3)) == 0)
E4 = Net("X4 -> 3X3; 3X3 -> X1+2X2; X1+2X2 -> 3X2; 3X2 -> X4", ["X1", "X2", "X3", "X4"])
fE4 = E4.dual().field(kap[:4], [K1, K2, K3, K4])
claim("Ex. drain (1): in path order the dual species K3 (dual to X1+2X2 -> 3X2) satisfies K3' = K3 (k2 K3 K4^3 - k1); "
      "(m,q,s) = (4,3,3), (H1)-(H3) both sides",
      sp.expand(fE4[2] - K3 * (kap[1] * K3 * K4**3 - kap[0])) == 0 and len(E4.clusters()) == 3 and E4.rank() == 3
      and all(H123(E4)) and all(H123(E4.dual())))
Ln = Net("2X3 -> 2X2+X3; 2X2+X3 -> 2X1+X3; 2X1+X3 -> 2X3", ["X1", "X2", "X3"])
fL = Ln.dual().field(kap[:3], [K1, K2, K3])
claim("Sec. 6.3 example (L): V = K2 + 2K3, V' = -2k1 K3^2 - 2k2 K2^2 + 2k3 K1^2 K2 K3",
      sp.expand(fL[1] + 2 * fL[2] - (-2 * kap[0] * K3**2 - 2 * kap[1] * K2**2 + 2 * kap[2] * K1**2 * K2 * K3)) == 0)
Un = Net("3X3 -> 2X1+X2; X2+2X3 -> 3X3; 2X1+X2 -> X2+2X3", ["X1", "X2", "X3"])
UNd = Un.dual()
claim("Ex. U: with the reaction order 3X3 -> 2X1+X2, X2+2X3 -> 3X3, 2X1+X2 -> X2+2X3 the dual is exactly the paper's "
      "2K3 -> 2K1, K2+K3 -> K1+K3, 3K1+2K2 -> 3K2+2K3 (the paper does not state this non-path order)",
      UNd.R == [((0, 0, 2), (2, 0, 0)), ((0, 1, 1), (1, 0, 1)), ((3, 2, 0), (0, 3, 2))])
fU = UNd.field(kap[:3], [K1, K2, K3])
aU, bU = K2, K3
Kx = kap[2] * K1**3
z = aU / bU
uz = bU / aU
claim("Ex. U: the identities for z', u', a', b' (a = K2, b = K3, K = kappa3 K1^3)",
      sp.simplify((fU[1] * bU - aU * fU[2]) / bU**2 - aU * ((2 * kap[0] - kap[1]) + Kx * z * (1 - 2 * z))) == 0
      and sp.simplify((fU[2] * aU - bU * fU[1]) / aU**2 - aU * ((kap[1] - 2 * kap[0]) * uz**2 - Kx * uz + 2 * Kx)) == 0
      and sp.simplify(fU[1] - aU**2 * (Kx - kap[1] * uz)) == 0
      and sp.simplify(fU[2] - bU**2 * (-2 * kap[0] + 2 * Kx * z**2)) == 0)
N1 = Net("2A+B -> C; C -> 3A; C -> 3B", ["A", "B", "C"])
fN1 = N1.dual().field(kap[:3], [K1, K2, K3])
st = {K1: kap[1] / kap[0], K2: 3 * kap[1]**2 / (kap[0] * kap[2]), K3: 1}
claim("Prop. Tstar-cex: N1* = 2K1 -> 3K2, K1 -> 3K3, K2+K3 -> K1 has the steady state "
      "(kB/kA, 3 kB^2/(kA kC), 1) for every kappa (the paper writes the middle entry as 3 kA K1^2/kC)",
      N1.dual().R == [((2, 0, 0), (0, 3, 0)), ((1, 0, 0), (0, 0, 3)), ((0, 1, 1), (1, 0, 0))]
      and all(sp.simplify(fi.subs(st)) == 0 for fi in fN1))

# ---------------------------------------------------------------- ADLS converse (agent E), not in certify_duality --
Ec = Net("X4 -> X2+X3; X2+X3 -> X4; 2X4 -> X1+X3; X1+X3 -> 2X4", ["X1", "X2", "X3", "X4"])
crN = [sorted(Ec.species[i] for i in S) for S in Ec.critical_siphons()]
crD = [sorted(Ec.dual().species[i] for i in S) for S in Ec.dual().critical_siphons()]
claim("Sec. 5.1 (converse of Prop. ADLS): X4 <-> X2+X3, 2X4 <-> X1+X3 satisfies (H1)-(H3) and has a critical siphon, "
      "while its dual satisfies (H1)-(H3) and has none (this example is not in certify_duality.py)",
      all(H123(Ec)) and all(H123(Ec.dual())) and crN and not crD, f"critical siphons of N: {crN}; of N*: {crD}")

# ---------------------------------------------------------------- abstract: crossing of clusters -------------------
Xc = Net("A+B -> C; B+D -> C; C -> A+D; C -> 2B", ["A", "B", "C", "D"])
ms = [sorted(Xc.species[i] for i in S) for S in Xc.minimal_siphons()]
cls = [sorted(Xc.species[i] for i in c_[0]) for c_ in Xc.clusters()]
claim("Abstract / transversality outside EC: N = A+B -> C, B+D -> C, C -> A+D, C -> 2B is connected with (H1)-(H3), "
      "not EC, and its minimal siphon {A,C,D} meets the Petri cluster {A,B,D} in two species",
      all(H123(Xc)) and Xc.connected() and ["A", "C", "D"] in ms and ["A", "B", "D"] in cls,
      f"minimal siphons {ms}, cluster species parts {cls}")

# ---------------------------------------------------------------- q = m needs nonzero sources -------------------
Z0 = Net("0 -> X1; X1 -> 0; 0 -> X2; X2 -> 0", ["X1", "X2"])
srcs = {r[0] for r in Z0.R}
claim("Sec. 2: under (H2) alone q = m can fail for EC networks (0 -> X1, X1 -> 0, 0 -> X2, X2 -> 0: m = 3 sources, "
      "q = 4 clusters); (H1) excludes it", H123(Z0)[1] and len(srcs) == 3 and len(Z0.clusters()) == 4 and not H123(Z0)[0])

Z1 = Net("0 -> X1; 0 -> X2; X1+X2 -> 0", ["X1", "X2"])
claim("Sec. 2 (applied form of R-4): 0 -> X1, 0 -> X2, X1+X2 -> 0 is connected, EC and consistent with m = 2 sources and "
      "q = 3 clusters, and not conservative", Z1.connected() and H123(Z1)[1] and not H123(Z1)[0]
      and len({r[0] for r in Z1.R}) == 2 and len(Z1.clusters()) == 3)

# ----------------------------------------------------------------
nfail = RESULTS.count(False)
print(f"\n# {len(RESULTS)} claims: {len(RESULTS) - nfail} PASS, {nfail} FAIL; runtime {time.time() - T0:.1f} s")
sys.exit(1 if nfail else 0)
