"""Exact face analysis of conservative mass-action networks of stoichiometric rank one (agent Q1, round 2).

Environment: Python 3.12, standard library only; exact rational arithmetic (fractions.Fraction).

A network M is given as (p, rx): p species, rx a list of reactions (source, product), each an integer p-tuple.
If rank Gamma(M) = 1, every reaction vector is c_j * g with g a primitive integer vector (sign fixed by the first
nonzero entry being positive is NOT imposed: g is the primitive vector of the first nonzero reaction vector).
M is conservative iff g has entries of both signs.  Then every positive compatibility class is a segment
x(xi) = x0 + xi g, and xi' = F(xi) = sum_j c_j k_j x^{y_j}.

Faces.  An *upper face* is a nonempty set Z of species with g_i < 0 for all i in Z; a *lower face* is a nonempty
Z with g_i > 0 for all i in Z.  The upper endpoint of a class is a point x+ with zero set Z (an upper face) and
x+_i > 0 off Z; every such point is the upper endpoint of a class.  Near it, x = x+ - eps g, and
F = eps^mu * L(k, x+) + O(eps^(mu+1)), with mu = min over active reactions (c_j != 0) of o_Z(j) = sum_{i in Z}
y_j(i), and L = sum_{j : o_Z(j) = mu} c_j k_j prod_{i in Z} |g_i|^{y_j(i)} prod_{i not in Z} (x+_i)^{y_j(i)}.
The endpoint attracts nearby interior points iff sigma * L > 0 (sigma = +1 upper, -1 lower; if L = 0 look at
higher order).  A reaction j of minimal order is *inward* if sigma * c_j < 0 and *outward* if sigma * c_j > 0.

Face types:
  IN   all minimal-order active reactions are inward      -> repelling for all rates and all classes
  OUT  all are outward ("drain face")                     -> attracting for all rates and all classes
  ESC  mixed, and some outward exponent beta_j = y_j|_{Z^c} is not in the convex hull of the inward ones
       -> for every rate vector some class reaching the face has an attracting endpoint
  DOM  mixed, and every outward exponent lies in the convex hull of the inward ones
       -> repelling for all classes when the outward rate constants are small enough; attracting somewhere when
          one outward rate constant is large enough.
"""
from fractions import Fraction
from itertools import combinations
from math import gcd


def deltas(rx):
    return [tuple(b - a for a, b in zip(s, t)) for s, t in rx]


def rank1(p, rx):
    """Return (g, c) with delta_j = c_j g, g primitive, or None if the rank is not 1."""
    ds = deltas(rx)
    g = None
    for d in ds:
        if any(d):
            G = 0
            for x in d:
                G = gcd(G, abs(x))
            g = tuple(x // G for x in d)
            break
    if g is None:
        return None
    k0 = next(i for i in range(p) if g[i] != 0)
    c = []
    for d in ds:
        if d[k0] % g[k0] != 0:
            return None
        cj = d[k0] // g[k0]
        if tuple(cj * x for x in g) != d:
            return None
        c.append(cj)
    return g, c


def solve_exact(A, b):
    """Solve A x = b (A: list of rows) exactly; return a solution if the system is consistent and the solution is
    unique, else None."""
    m = len(A)
    ncol = len(A[0]) if m else 0
    M = [[Fraction(x) for x in row] + [Fraction(bb)] for row, bb in zip(A, b)]
    piv = []
    r = 0
    for col in range(ncol):
        pr = next((i for i in range(r, m) if M[i][col] != 0), None)
        if pr is None:
            return None  # not unique (free variable)
        M[r], M[pr] = M[pr], M[r]
        pv = M[r][col]
        M[r] = [x / pv for x in M[r]]
        for i in range(m):
            if i != r and M[i][col] != 0:
                f = M[i][col]
                M[i] = [x - f * y for x, y in zip(M[i], M[r])]
        piv.append(col)
        r += 1
    for i in range(r, m):
        if M[i][ncol] != 0:
            return None
    return [M[i][ncol] for i in range(ncol)]


def in_conv(beta, Q):
    """Exact test beta in conv(Q) (Q a list of integer vectors).  Returns barycentric weights (dict index->Fraction)
    on an affinely independent subset, or None.  Caratheodory: it suffices to try affinely independent subsets."""
    if not Q:
        return None
    d = len(beta)
    for k in range(1, min(len(Q), d + 1) + 1):
        for sub in combinations(range(len(Q)), k):
            A = [[Q[q][i] for q in sub] for i in range(d)] + [[1] * k]
            b = list(beta) + [1]
            lam = solve_exact(A, b)
            if lam is not None and all(x >= 0 for x in lam):
                return {q: x for q, x in zip(sub, lam)}
    return None


def faces(p, g):
    neg = [i for i in range(p) if g[i] < 0]
    pos = [i for i in range(p) if g[i] > 0]
    out = []
    for sigma, S in ((1, neg), (-1, pos)):
        for k in range(1, len(S) + 1):
            for Z in combinations(S, k):
                out.append((sigma, Z))
    return out


def face_analysis(p, rx, sigma, Z, g, c):
    """Classify the face (sigma, Z).  Returns a dict with type, minimal order mu, the minimal reactions split
    into inward / outward, and the exponents beta_j = y_j restricted to the complement of Z."""
    Zs = set(Z)
    act = [j for j in range(len(rx)) if c[j] != 0]
    order = {j: sum(rx[j][0][i] for i in Z) for j in act}
    mu = min(order.values())
    J = [j for j in act if order[j] == mu]
    inw = [j for j in J if sigma * c[j] < 0]
    outw = [j for j in J if sigma * c[j] > 0]
    comp = [i for i in range(p) if i not in Zs]
    beta = {j: tuple(rx[j][0][i] for i in comp) for j in J}
    res = {"sigma": sigma, "Z": tuple(Z), "mu": mu, "in": inw, "out": outw, "beta": beta, "comp": comp}
    if not outw:
        res["type"] = "IN"
    elif not inw:
        res["type"] = "OUT"
    else:
        Q = [beta[j] for j in inw]
        esc = [j for j in outw if in_conv(beta[j], Q) is None]
        res["type"] = "ESC" if esc else "DOM"
        res["esc"] = esc
    return res


def classify(p, rx):
    """Classify a rank-one network.  Returns (verdict, data) with verdict in
    'not-rank-1', 'not-conservative', 'ALL', 'NONE', 'SOME', 'TWO-SIDED' (DOM faces on both sides, decided by
    tropical_some below), and data = (g, c, list of face analyses)."""
    rc = rank1(p, rx)
    if rc is None:
        return "not-rank-1", None
    g, c = rc
    if not (any(x > 0 for x in g) and any(x < 0 for x in g)):
        return "not-conservative", (g, c, [])
    fa = [face_analysis(p, rx, s, Z, g, c) for s, Z in faces(p, g)]
    types = [f["type"] for f in fa]
    if all(t == "IN" for t in types):
        return "ALL", (g, c, fa)
    if any(t in ("OUT", "ESC") for t in types):
        return "NONE", (g, c, fa)
    sides = {f["sigma"] for f in fa if f["type"] == "DOM"}
    if len(sides) == 1:
        return "SOME", (g, c, fa)
    return "TWO-SIDED", (g, c, fa)


# ------------------------------------------------------------------------------------------------------------
# exact endpoint expansion (independent check of the face formula)
# ------------------------------------------------------------------------------------------------------------
def poly_mul(a, b):
    out = [Fraction(0)] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                out[i + j] += x * y
    return out


def poly_pow(a, e):
    out = [Fraction(1)]
    for _ in range(e):
        out = poly_mul(out, a)
    return out


def endpoint_expansion(p, rx, g, c, kappa, xplus, sigma):
    """F along x = xplus - sigma*eps*g as a polynomial in eps (list of Fractions, constant term first)."""
    total = [Fraction(0)]
    for j, (src, _) in enumerate(rx):
        if c[j] == 0:
            continue
        term = [Fraction(c[j]) * Fraction(kappa[j])]
        for i in range(p):
            if src[i]:
                lin = [Fraction(xplus[i]), Fraction(-sigma * g[i])]
                term = poly_mul(term, poly_pow(lin, src[i]))
        if len(term) > len(total):
            total += [Fraction(0)] * (len(term) - len(total))
        for k, x in enumerate(term):
            total[k] += x
    return total


def lowest(poly):
    for k, x in enumerate(poly):
        if x != 0:
            return k, x
    return None, Fraction(0)


def leading_L(p, rx, fa, g, c, kappa, xplus):
    """Leading coefficient L predicted by the face formula at the point xplus of the face fa."""
    L = Fraction(0)
    for j in fa["in"] + fa["out"]:
        t = Fraction(c[j]) * Fraction(kappa[j])
        for i in range(p):
            e = rx[j][0][i]
            if e:
                t *= (Fraction(abs(g[i])) ** e) if i in fa["Z"] else Fraction(xplus[i]) ** e
        L += t
    return L


# ------------------------------------------------------------------------------------------------------------
# duals and pretty printing
# ------------------------------------------------------------------------------------------------------------
def dual(n, rx):
    """Iosif dual: species = reactions of N, one reaction per species of N (source row i of Ye, product row i of
    Yp)."""
    r = len(rx)
    return r, [(tuple(rx[j][0][i] for j in range(r)), tuple(rx[j][1][i] for j in range(r))) for i in range(n)]


def cplx(y, names):
    parts = []
    for i, e in enumerate(y):
        if e:
            parts.append((str(e) if e > 1 else "") + names[i])
    return "+".join(parts) if parts else "0"


def show(n, rx, names=None):
    if names is None:
        names = [f"X{i+1}" for i in range(n)]
    return ", ".join(f"{cplx(a, names)} -> {cplx(b, names)}" for a, b in rx)
