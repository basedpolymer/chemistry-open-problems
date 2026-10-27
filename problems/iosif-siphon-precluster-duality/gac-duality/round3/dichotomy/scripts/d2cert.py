"""Agent D2 (round 3): rate-independent certificates for the persistence of the rank-2 dual N* (s = 2).

Environment: Python 3.12, standard library only; exact rational arithmetic (fractions.Fraction); LPs by exact
Fourier-Motzkin elimination (E2's e2lib.fm_point / lp_point, imported unchanged).  Statements and proofs: notes.md.

Boundary pieces of the classes of N* (all realizable ones, ties included):
  facets   Z, a nonempty subset of a direction class D of N (reactions with positively parallel reaction vectors);
  vertices Zv, a set of reactions of N of rank 2 whose reduced vectors g_j lie in an open half-plane.
Certificates:
  ALL  (Theorem P2 + Proposition R): every siphon facet is IN with a one-way, nonzero edge flow; every vertex
       equilibrium has a repelling function a^lam b^mu or lam a + mu b passing the domination test; the directions
       that carry a siphon facet do not positively span the plane (no boundary cycle).
  NONE (Proposition H, tropical Hurwitz vertex; Proposition E, tropical edge equilibrium on an OUT/ESC facet).
"""
import os
import sys
from fractions import Fraction
from itertools import combinations
from math import comb

HERE = os.path.dirname(os.path.abspath(__file__))
E2DIR = os.path.normpath(os.path.join(HERE, "..", "..", "..", "round2", "computation", "scripts"))
if E2DIR not in sys.path:
    sys.path.insert(0, E2DIR)
if HERE not in sys.path:
    sys.path.insert(0, HERE)
import e2lib  # noqa: E402
import d2lib  # noqa: E402

LAMS = [(1, 0), (0, 1), (1, 1), (1, 2), (2, 1), (1, 3), (3, 1), (2, 3), (3, 2)]


def cross(u, v):
    return u[0] * v[1] - u[1] * v[0]


def direction_classes(N):
    cl = []
    for j in range(N.r):
        for c in cl:
            if N.parallel(c[0], j):
                c.append(j)
                break
        else:
            cl.append([j])
    return cl


# ---------------------------------------------------------------------------------------------- facets
def facets(N):
    """All realizable facets Z (nonempty subsets of a direction class).  For siphon facets: type, minimal order mu,
    minimal-order species J (inward = net produced by the reactions of Z), exponents beta_i off Z, edge flow."""
    out = []
    for D in direction_classes(N):
        j0 = D[0]
        ch = [i for i in range(N.n) if N.G[i][j0] != 0]
        sgn = {i: (1 if N.G[i][j0] > 0 else -1) for i in ch}
        d = (-N.gh[j0][1], N.gh[j0][0])
        for size in range(1, len(D) + 1):
            for Z in combinations(D, size):
                load = [sum(N.Y[i][j] for j in Z) for i in range(N.n)]
                siphon = all(load[i] > 0 for i in ch)
                rec = {"Z": Z, "D": tuple(D), "siphon": siphon}
                if siphon:
                    mu = min(load[i] for i in ch)
                    J = [i for i in ch if load[i] == mu]
                    off = [l for l in range(N.r) if l not in Z]
                    beta = {i: tuple(N.Y[i][l] for l in off) for i in range(N.n)}
                    inw = [i for i in J if sgn[i] > 0]
                    outw = [i for i in J if sgn[i] < 0]
                    if not outw:
                        typ = "IN"
                    elif not inw:
                        typ = "OUT"
                    else:
                        hull = [list(beta[i]) for i in inw]
                        esc = any(e2lib.in_convex_hull(list(beta[o]), hull) is None for o in outw)
                        typ = "ESC" if esc else "DOM"
                    flow = []
                    for i in range(N.n):
                        if load[i] == 0 and any(N.pi[i]):
                            p = N.pi[i]
                            c = p[0] / d[0] if d[0] != 0 else p[1] / d[1]
                            assert (c * d[0], c * d[1]) == p
                            flow.append((i, c))
                    signs = {1 if c > 0 else -1 for i, c in flow}
                    fl = "zero" if not flow else ("uni" if len(signs) == 1 else "bi")
                    rec.update(mu=mu, J=J, inw=inw, outw=outw, type=typ, beta=beta, off=off, flow=flow, fl=fl)
                out.append(rec)
    return out


def edge_equilibrium_cert(N, f):
    """Proposition E: for an OUT/ESC siphon facet, a direction theta (coordinates off Z) along which the
    minimal-order transverse terms are dominated by outward ones and the edge flow has a two-coloured leading face.
    Returns a witness or None.  Zero edge flow: requires mu = 1."""
    if f["type"] not in ("OUT", "ESC"):
        return None
    dim = len(f["off"])
    beta = f["beta"]
    outs = f["outw"] if f["inw"] else [None]
    for o in outs:
        base = []
        if o is not None:
            for i in f["inw"]:
                base.append(([-(x - y) for x, y in zip(beta[o], beta[i])], -1))
        if f["fl"] == "zero":
            if f["mu"] != 1:
                return None
            th = e2lib.fm_point(base, dim) if base else [Fraction(0)] * dim
            if th is not None:
                return ("zero-flow", o, th)
            continue
        pos = [i for i, c in f["flow"] if c > 0]
        neg = [i for i, c in f["flow"] if c < 0]
        allf = [i for i, c in f["flow"]]
        for ip in pos:
            for im in neg:
                if beta[ip] == beta[im]:
                    continue
                eqs = [([x - y for x, y in zip(beta[ip], beta[im])], 0)]
                ineqs = list(base)
                for i in allf:
                    if i in (ip, im):
                        continue
                    ineqs.append(([-(x - y) for x, y in zip(beta[ip], beta[i])], -1))
                th = e2lib.lp_point(eqs, ineqs, dim)
                if th is not None:
                    return ("bicoloured", o, ip, im, th)
    return None


# ---------------------------------------------------------------------------------------------- vertices
def polymul(P, Q):
    R = {}
    for (a, b), c in P.items():
        for (x, y), e in Q.items():
            k = (a + x, b + y)
            R[k] = R.get(k, 0) + c * e
    return R


def binpow(al, be, m):
    P = {}
    for k in range(m + 1):
        c = comb(m, k) * al ** k * be ** (m - k)
        if c != 0:
            P[(k, m - k)] = c
    return P


def vertices(N):
    """All realizable vertex zero sets Zv with the corner data: extreme reactions e1, e2 (a = k_e1, b = k_e2),
    k_l = al_l a + be_l b (l in Zv), equilibrium flag, local monomial support of each changed species."""
    out = []
    for size in range(2, N.r + 1):
        for Zv in combinations(range(N.r), size):
            g = [N.gh[j] for j in Zv]
            if all(cross(g[0], x) == 0 for x in g):
                continue
            X0 = e2lib.fm_point([([x[0], x[1]], -1) for x in g], 2)
            if X0 is None:
                continue
            m = (-X0[0], -X0[1])  # <g_j, m> > 0 for all j in Zv
            ks = {j: cross(m, N.gh[j]) / (m[0] * N.gh[j][0] + m[1] * N.gh[j][1]) for j in Zv}
            e1 = min(Zv, key=lambda j: ks[j])
            e2 = max(Zv, key=lambda j: ks[j])
            G1, G2 = N.gh[e1], N.gh[e2]
            det = cross(G1, G2)
            ab = {}
            for l in Zv:
                gl = N.gh[l]
                al = cross(gl, G2) / det
                be = cross(G1, gl) / det
                assert al >= 0 and be >= 0
                ab[l] = (al, be)
            eq = all(sum(N.Y[i][l] for l in Zv) >= 1 for i in N.changed)
            rec = {"Zv": Zv, "e1": e1, "e2": e2, "ab": ab, "eq": eq}
            if eq:
                loc = {}
                for i in N.changed:
                    P = {(0, 0): Fraction(1)}
                    for l in Zv:
                        if N.Y[i][l]:
                            P = polymul(P, binpow(ab[l][0], ab[l][1], N.Y[i][l]))
                    loc[i] = {k: v for k, v in P.items() if v != 0}
                rec["loc"] = loc
            out.append(rec)
    return out


def newton_dom(e, base):
    """Is the exponent e in conv(base) + R^2_{>=0}?  Returns 'strict' if some point p of conv(base) satisfies
    p <= e, p != e (then x^e = o(sum of the base monomials) near 0 in the open quadrant, by weighted AM-GM), a list
    of supports [(e1, e2)] if only boundary points p = e exist, or None.  (Pairs suffice in the plane.)"""
    best = None
    for e1 in base:
        for e2 in base:
            # p(l) = e2 + l (e1 - e2) <= e, l in [0, 1]
            lo, hi = Fraction(0), Fraction(1)
            ok = True
            for c in (0, 1):
                d = e1[c] - e2[c]
                rhs = e[c] - e2[c]
                if d == 0:
                    if rhs < 0:
                        ok = False
                elif d > 0:
                    hi = min(hi, Fraction(rhs, d))
                else:
                    lo = max(lo, Fraction(rhs, d))
            if not ok or lo > hi:
                continue
            for l in {lo, hi, (lo + hi) / 2}:
                p = (e2[0] + l * (e1[0] - e2[0]), e2[1] + l * (e1[1] - e2[1]))
                if p != (e[0], e[1]):
                    return "strict"
            best = (best or []) + [(e1, e2)]
    return best


def dominated(terms):
    """terms: list of (exponent (x, y), sign).  True if there is a positive term and every negative term is
    strictly Newton-dominated by the positive terms (then the positive terms win near 0 for all coefficients)."""
    pos = sorted({e for e, s in terms if s > 0})
    neg = {e for e, s in terms if s < 0}
    if not pos:
        return False
    for e in neg:
        if newton_dom(e, pos) != "strict":
            return False
    return True


def vertex_repelling(N, v):
    """Proposition R: a function V = a^lam b^mu or lam a + mu b with V' > 0 near the vertex on the interior of the
    class, for every kappa and every class (domination test).  Returns the function or None."""
    e1, e2 = v["e1"], v["e2"]
    for lam, mu in LAMS:
        tm, tl = [], []
        for i, P in v["loc"].items():
            g1, g2 = N.G[i][e1], N.G[i][e2]
            for (x, y) in P:
                if lam and g1:
                    tm.append(((x - 1, y), g1))
                if mu and g2:
                    tm.append(((x, y - 1), g2))
                c = lam * g1 + mu * g2
                if c:
                    tl.append(((x, y), c))
        if dominated(tm):
            return ("monomial", lam, mu)
        if dominated(tl):
            return ("linear", lam, mu)
    return None


def vertex_hurwitz_cert(N, v):
    """Proposition H: a direction theta for the coordinates off Zv such that, for v_off = exp(t theta), t large, the
    linearisation at the vertex has trace < 0 and determinant > 0 for every kappa.  Returns a witness or None."""
    Zv = v["Zv"]
    off = [l for l in range(N.r) if l not in Zv]
    dim = len(off)
    e1, e2 = v["e1"], v["e2"]
    o1 = []
    for i in N.changed:
        tot = sum(N.Y[i][l] for l in Zv)
        if tot != 1:
            continue
        l = next(l for l in Zv if N.Y[i][l] == 1)
        u = (Fraction(N.G[i][e1]), Fraction(N.G[i][e2]))
        vv = v["ab"][l]
        w = tuple(N.Y[i][x] for x in off)
        o1.append((i, u, vv, w))
    if not o1:
        return None
    tr = {i: u[0] * vv[0] + u[1] * vv[1] for i, u, vv, w in o1}
    W = {i: w for i, u, vv, w in o1}
    pairs = []
    for (i, u, vv, w), (i2, u2, vv2, w2) in combinations(o1, 2):
        m = cross(u, u2) * cross(vv, vv2)
        if m != 0:
            pairs.append(((i, i2), m, tuple(a + b for a, b in zip(w, w2))))
    negt = [i for i in tr if tr[i] < 0]
    post = [i for i in tr if tr[i] > 0]
    posp = [p for p in pairs if p[1] > 0]
    negp = [p for p in pairs if p[1] < 0]
    for it in negt:
        for P in posp:
            ineqs = []
            for i in post:
                ineqs.append(([-(a - b) for a, b in zip(W[it], W[i])], -1))
            for Q in negp:
                ineqs.append(([-(a - b) for a, b in zip(P[2], Q[2])], -1))
            if not ineqs:
                return ("theta", it, P[0], [Fraction(0)] * dim)
            if dim == 0:
                continue
            th = e2lib.fm_point(ineqs, dim)
            if th is not None:
                return ("theta", it, P[0], th)
    return None


def _lead_neg(terms, dim):
    """terms: list of (w, sign) with w in Z^dim.  Inequalities on theta saying that some negative term strictly
    dominates (by >= 1 in <theta, w>) every positive term.  Returns a list of alternatives (lists of ineqs)."""
    negs = [w for w, s in terms if s < 0]
    poss = [w for w, s in terms if s > 0]
    alts = []
    for wn in negs:
        alts.append([([-(a - b) for a, b in zip(wn, wp)], -1) for wp in poss])
    return alts


def vertex_cm_cert(N, v):
    """Proposition CM: one column of the linearisation at the vertex vanishes identically (no local monomial a, say),
    the other diagonal entry delta is < 0, and the reduced equation on the center manifold, a' = c a^m + ..., has
    c < 0, for v_off = exp(t theta), t large, every kappa.  Returns a witness or None."""
    Zv = v["Zv"]
    off = [l for l in range(N.r) if l not in Zv]
    dim = len(off)
    W = {i: tuple(N.Y[i][x] for x in off) for i in N.changed}
    for side in (0, 1):
        ea, eb = (v["e1"], v["e2"]) if side == 0 else (v["e2"], v["e1"])

        def mono(x, y):
            return (x, y) if side == 0 else (y, x)
        loc = v["loc"]
        # first column zero: no species with the monomial "a" (and nonzero effect)
        if any(mono(1, 0) in loc[i] and (N.G[i][ea] or N.G[i][eb]) for i in N.changed):
            continue
        Sd = [(W[i], 1 if N.G[i][eb] > 0 else -1) for i in N.changed if mono(0, 1) in loc[i] and N.G[i][eb]]
        if not any(s < 0 for w, s in Sd):
            continue
        ms = [x for i in N.changed for (x, y) in [((m[0], m[1]) if side == 0 else (m[1], m[0])) for m in loc[i]]
              if y == 0 and N.G[i][ea]]
        ks = [x for i in N.changed for (x, y) in [((m[0], m[1]) if side == 0 else (m[1], m[0])) for m in loc[i]]
              if y == 0 and N.G[i][eb]]
        if not ms:
            continue
        m = min(ms)
        k = min(ks) if ks else None
        alts_d = _lead_neg(Sd, dim)
        if k is None or m < k:
            Sf = [(W[i], 1 if N.G[i][ea] > 0 else -1) for i in N.changed if mono(m, 0) in loc[i] and N.G[i][ea]]
            alts_c = _lead_neg(Sf, dim)
        elif m == k:
            A = [i for i in N.changed if mono(k, 0) in loc[i] and (N.G[i][ea] or N.G[i][eb])]
            B = [i for i in N.changed if mono(0, 1) in loc[i] and (N.G[i][ea] or N.G[i][eb])]
            prs = []
            for i in A:
                for i2 in B:
                    cr = N.G[i][ea] * N.G[i2][eb] - N.G[i2][ea] * N.G[i][eb]
                    if cr:
                        prs.append((tuple(a + b for a, b in zip(W[i], W[i2])), -1 if cr > 0 else 1))
            alts_c = _lead_neg(prs, dim)  # need sum C C' cross > 0: "negative" here means cross > 0
        else:
            continue
        for ad in alts_d:
            for ac in alts_c:
                ineqs = ad + ac
                if not ineqs:
                    return ("cm", side, m, k, [Fraction(0)] * dim)
                if dim == 0:
                    continue
                th = e2lib.fm_point(ineqs, dim)
                if th is not None:
                    return ("cm", side, m, k, th)
    return None


LAMS_L = [(1, 1), (1, 2), (2, 1), (1, 3), (3, 1), (2, 3), (3, 2), (1, 4), (4, 1)]


def vertex_lyap_cert(N, v, cap=400):
    """Proposition L: V = lam a + mu b (lam, mu > 0) with V' < 0 on the open quadrant near the vertex, for
    v_off = exp(t theta), t large, every kappa: every positive term of V' is Newton-dominated by the negative terms,
    strictly, or on the boundary with the negative coefficients tropically larger.  Returns a witness or None."""
    Zv = v["Zv"]
    off = [l for l in range(N.r) if l not in Zv]
    dim = len(off)
    e1, e2 = v["e1"], v["e2"]
    W = {i: tuple(N.Y[i][x] for x in off) for i in N.changed}
    for lam, mu in LAMS_L:
        terms = []
        for i, P in v["loc"].items():
            c = lam * N.G[i][e1] + mu * N.G[i][e2]
            if c:
                for e in P:
                    terms.append((e, 1 if c > 0 else -1, i))
        negs = [(e, i) for e, s, i in terms if s < 0]
        if not negs:
            continue
        nexp = sorted({e for e, i in negs})
        options = []
        ok = True
        for e, s, i in terms:
            if s < 0:
                continue
            r = newton_dom(e, nexp)
            if r is None:
                ok = False
                break
            if r == "strict":
                continue
            opts = []
            for (x1, x2) in r:
                for n1 in [j for ee, j in negs if ee == x1]:
                    for n2 in [j for ee, j in negs if ee == x2]:
                        cons = [([-(a - b) for a, b in zip(W[n1], W[i])], -1)]
                        if n2 != n1:
                            cons.append(([-(a - b) for a, b in zip(W[n2], W[i])], -1))
                        opts.append(cons)
            options.append(opts)
        if not ok:
            continue
        combos = [[]]
        for opts in options:
            combos = [c + o for c in combos for o in opts][:cap]
        for ineqs in combos:
            if not ineqs:
                return ("lyap", lam, mu, [Fraction(0)] * dim)
            if dim == 0:
                continue
            th = e2lib.fm_point(ineqs, dim)
            if th is not None:
                return ("lyap", lam, mu, th)
    return None


def vertex_sector_cert(N, v, cap=3000):
    """Proposition S: an attracting sector.  With b = a u (or a = b u) and, for t large, v_off = exp(t theta),
    u in [sig1, sig2] exp(-t omega) (or u in [0, sig2] exp(-t omega)): a' < 0 in the sector, u' > 0 on the lower
    ray and u' < 0 on the upper ray, near the vertex, for every kappa.  Exact LP in (theta, omega)."""
    Zv = v["Zv"]
    off = [l for l in range(N.r) if l not in Zv]
    dim = len(off) + 1  # theta, omega
    for side in (0, 1):
        ea, eb = (v["e1"], v["e2"]) if side == 0 else (v["e2"], v["e1"])
        ta, tu = [], []
        for i in N.changed:
            g1, g2 = N.G[i][ea], N.G[i][eb]
            w = [N.Y[i][x] for x in off]
            for m in v["loc"][i]:
                x, y = (m[0], m[1]) if side == 0 else (m[1], m[0])
                if g1:
                    ta.append((x + y, y, w, 1 if g1 > 0 else -1))
                if g2:
                    tu.append((x + y - 1, y, w, 1 if g2 > 0 else -1))
                if g1:
                    tu.append((x + y - 1, y + 1, w, -1 if g1 > 0 else 1))
        if not ta or not tu:
            continue
        da = min(t[0] for t in ta)
        Ga = [t for t in ta if t[0] == da]
        du = min(t[0] for t in tu)
        Gu = [t for t in tu if t[0] == du]

        def tau(t):  # coefficient vector of <theta, w> - omega * pow
            return list(t[2]) + [-t[1]]

        def ge1(p, q):  # tau_p - tau_q >= 1  as  a.x <= b
            return ([-(a - b) for a, b in zip(tau(p), tau(q))], -1)
        negA = [t for t in Ga if t[3] < 0]
        posA = [t for t in Ga if t[3] > 0]
        for axis in (False, True):
            altsA = []
            for nA in negA:
                if axis and any(p[1] < nA[1] for p in posA):
                    continue
                altsA.append([ge1(nA, p) for p in posA])
            if not altsA:
                continue
            pos_u = [t for t in Gu if t[3] > 0]
            neg_u = [t for t in Gu if t[3] < 0]
            pairs = []
            for Q in neg_u:
                if axis:
                    pairs.append((None, Q))
                for P in pos_u:
                    if P[1] < Q[1]:
                        pairs.append((P, Q))
            tried = 0
            for P, Q in pairs:
                ref = Q
                eqs = []
                if P is not None:
                    eqs.append(([a - b for a, b in zip(tau(P), tau(Q))], 0))
                ineqs = []
                bad = False
                for R in Gu:
                    if R is P or R is Q:
                        continue
                    if tau(R) == tau(ref):
                        lo = P[1] if P is not None else -1
                        if R[1] <= lo and R[3] < 0:
                            bad = True
                        if R[1] >= Q[1] and R[3] > 0:
                            bad = True
                        if P is None and R[1] < Q[1] and R[3] > 0:
                            pass
                        continue
                    ineqs.append(ge1(ref, R))
                if bad:
                    continue
                for aA in altsA:
                    tried += 1
                    if tried > cap:
                        break
                    pt = e2lib.lp_point(eqs, ineqs + aA, dim)
                    if pt is not None:
                        return ("sector", side, axis, pt)
    return None


# ---------------------------------------------------------------------------------------------- verdict
def positively_spanning(vecs):
    if not vecs:
        return False
    for g in vecs:
        for X in ((-g[1], g[0]), (g[1], -g[0])):
            if all(h[0] * X[0] + h[1] * X[1] >= 0 for h in vecs):
                return False
    return True


def flow_nondegenerate(N, f):
    """For a siphon facet with a two-signed edge flow F: True if F cannot vanish identically on an edge (the terms of
    highest degree in the edge parameter have one sign of c_i * prod_l d_l^{s_il}, for every kappa and class)."""
    D0 = f["D"][0]
    d = (-N.gh[D0][1], N.gh[D0][0])
    dl = {l: N.gh[l][0] * d[0] + N.gh[l][1] * d[1] for l in f["off"]}
    best, signs = None, set()
    for i, c in f["flow"]:
        deg = sum(N.Y[i][l] for l in f["off"] if dl[l] != 0)
        sg = 1 if c > 0 else -1
        for l in f["off"]:
            if dl[l] < 0 and N.Y[i][l] % 2 == 1:
                sg = -sg
        if best is None or deg > best:
            best, signs = deg, {sg}
        elif deg == best:
            signs.add(sg)
    return len(signs) == 1


def facet_ok_for_all(N, f):
    """Condition (A1) of Theorem P2: no equilibrium in the open edge (one-way flow), or IN with finitely many edge
    equilibria (two-way flow that cannot vanish identically)."""
    if f["fl"] == "uni":
        return True
    if f["fl"] == "bi" and f["type"] == "IN" and flow_nondegenerate(N, f):
        return True
    return False


def blocks(N):
    """Connected components of the dual N* (species k_j linked when they occur in a common dual reaction)."""
    parent = list(range(N.r))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    for i in N.changed:
        sup = [j for j in range(N.r) if N.Y[i][j] or N.Yp[i][j]]
        for j in sup[1:]:
            parent[find(j)] = find(sup[0])
    comp = {}
    for j in range(N.r):
        comp.setdefault(find(j), []).append(j)
    return list(comp.values())


def product_verdict(N):
    """Proposition B: if N* splits into two blocks of rank one, its classes are products of segments and N* is
    persistent iff both blocks are; Q1's Theorem S1 (s1lib.classify, imported unchanged) decides each block."""
    B = blocks(N)
    if len(B) < 2:
        return None
    QDIR = os.path.normpath(os.path.join(HERE, "..", "..", "..", "round2", "persistence", "scripts"))
    if QDIR not in sys.path:
        sys.path.insert(0, QDIR)
    import s1lib
    verdicts = []
    for blk in B:
        rxb = []
        for i in N.changed:
            if any(N.Y[i][j] or N.Yp[i][j] for j in blk):
                rxb.append((tuple(N.Y[i][j] for j in blk), tuple(N.Yp[i][j] for j in blk)))
        verdicts.append(s1lib.classify(len(blk), rxb)[0])
    if all(v == "ALL" for v in verdicts):
        return "ALL", verdicts
    if any(v == "NONE" for v in verdicts):
        return "NONE", verdicts
    return None


def analyse(n, rx):
    N = d2lib.Net(n, rx)
    pv = product_verdict(N)
    if pv is not None:
        return {"verdict": pv[0], "none": [("B", pv[1])] if pv[0] == "NONE" else [], "uf": [], "uv": [], "rv": [],
                "cycle": None, "F": [], "V": [], "product": pv[1]}
    F = facets(N)
    V = vertices(N)
    res = {"none": [], "uf": [], "uv": [], "rv": [], "cycle": None, "F": F, "V": V}
    for f in F:
        if f["siphon"] and f["type"] in ("OUT", "ESC"):
            c = edge_equilibrium_cert(N, f)
            if c is not None:
                res["none"].append(("E", f["Z"], f["type"], c))
    for v in V:
        if v["eq"]:
            c = vertex_hurwitz_cert(N, v)
            if c is not None:
                res["none"].append(("H", v["Zv"], c))
                continue
            c = vertex_cm_cert(N, v)
            if c is not None:
                res["none"].append(("CM", v["Zv"], c))
                continue
            c = vertex_lyap_cert(N, v)
            if c is not None:
                res["none"].append(("L", v["Zv"], c))
                continue
            c = vertex_sector_cert(N, v)
            if c is not None:
                res["none"].append(("S", v["Zv"], c))
    for f in F:
        if f["siphon"] and not facet_ok_for_all(N, f):
            res["uf"].append((f["Z"], f["type"], f["fl"], f["mu"]))
    for v in V:
        if v["eq"]:
            rp = vertex_repelling(N, v)
            if rp is None:
                res["uv"].append(v["Zv"])
            else:
                res["rv"].append((v["Zv"], rp))
    dsip = []
    for D in direction_classes(N):
        if any(f["siphon"] and f["D"] == tuple(D) for f in F):
            dsip.append(N.gh[D[0]])
    res["cycle"] = positively_spanning(dsip)
    all_ok = not res["uf"] and not res["uv"] and not res["cycle"]
    if res["none"] and all_ok:
        res["verdict"] = "CONFLICT"
    elif res["none"]:
        res["verdict"] = "NONE"
    elif all_ok:
        res["verdict"] = "ALL"
    else:
        res["verdict"] = "UNDECIDED"
    return res
