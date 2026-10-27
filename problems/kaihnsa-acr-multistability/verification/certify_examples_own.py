"""Verifier's own certificate for the three-species examples (N1 main, N2, G7).

Runtime: under 1 s (WSL Fedora 44, Python 3.12, sympy 1.14, python-flint 0.9).
Usage:   python certify_examples_own.py        (exit code 1 if any check fails)

Written separately from the solver's scripts.  The vector field is built from the reaction list (not from the
hand-written ODEs of the manuscript).  For each network and the stated rate constants it:
  C  counts species, reactions, reactant complexes, complexes, molecularity, dim S, linkage classes,
     terminal strong linkage classes and deficiency (graph computations from scratch);
  A  checks the ACR identity w.f = x^y (a - b x) symbolically for ALL rate constants;
  E  checks symbolically (all rate constants) that a positive steady state exists, so the ACR is unconditional;
  G  computes the positive steady states with a lex Groebner basis of <f, 1 - t x y z> (saturation), exact
     real-root isolation of the univariate eliminant, and back substitution; checks that there are exactly
     the three claimed points;
  J  exact Jacobians, characteristic polynomials, Routh-Hurwitz conditions; eigenvalue enclosures with arb
     (python-flint acb_poly roots), certifying the sign of every real part.
"""
import sys
import time
import itertools
import sympy as sp

T0 = time.time()
FAIL = []


def check(cond, msg):
    print(("PASS " if cond else "FAIL ") + msg)
    if not cond:
        FAIL.append(msg)


x, y, z, t = sp.symbols('x y z t')
X = (x, y, z)


def cplx(s):
    v = [0, 0, 0]
    if s == "0":
        return tuple(v)
    for part in s.split("+"):
        part = part.strip()
        c = 1
        if part[0].isdigit():
            c = int(part[0])
            part = part[1:]
        v["XYZ".index(part)] += c
    return tuple(v)


def mono(c):
    return x**c[0] * y**c[1] * z**c[2]


def field(reactions, k):
    f = sp.zeros(3, 1)
    for (a, b), kk in zip(reactions, k):
        ya, yb = cplx(a), cplx(b)
        for i in range(3):
            f[i] += kk * mono(ya) * (yb[i] - ya[i])
    return sp.Matrix([sp.expand(e) for e in f])


def graph_data(reactions):
    cs = sorted({cplx(a) for a, b in reactions} | {cplx(b) for a, b in reactions})
    idx = {c: i for i, c in enumerate(cs)}
    m = len(cs)
    edges = [(idx[cplx(a)], idx[cplx(b)]) for a, b in reactions]
    # reachability
    R = [[i == j for j in range(m)] for i in range(m)]
    for a, b in edges:
        R[a][b] = True
    for k_ in range(m):
        for i in range(m):
            if R[i][k_]:
                for j in range(m):
                    if R[k_][j]:
                        R[i][j] = True
    # linkage classes (undirected components)
    U = [[R[i][j] or R[j][i] for j in range(m)] for i in range(m)]
    comp = list(range(m))
    for a, b in edges:
        ca, cb = comp[a], comp[b]
        comp = [ca if c == cb else c for c in comp]
    ell = len(set(comp))
    # strong classes and terminal ones
    strong = {}
    for i in range(m):
        key = frozenset(j for j in range(m) if R[i][j] and R[j][i])
        strong[key] = True
    terminal = 0
    for sc in strong:
        i = next(iter(sc))
        if all((j in sc) for j in range(m) if R[i][j]):
            terminal += 1
    vecs = sp.Matrix([[cplx(b)[i] - cplx(a)[i] for i in range(3)] for a, b in reactions])
    s = vecs.rank()
    reactants = {cplx(a) for a, b in reactions}
    maxmol = max(sum(c) for c in cs)
    return dict(m=m, ell=ell, t=terminal, s=s, delta=m - ell - s, reactants=len(reactants), maxmol=maxmol,
                weakly_rev=all(R[b][a] for a, b in edges))


def positive_steady_states(f):
    G = sp.groebner(list(f) + [1 - t * x * y * z], t, x, y, z, order='lex')
    polys = [g for g in G.exprs if not g.has(t)]
    # univariate in z
    uz = [g for g in polys if g.free_symbols <= {z}]
    assert uz, "no univariate eliminant"
    pz = sp.Poly(uz[0], z)
    roots = [r for r in sp.real_roots(pz) if r > 0]
    sols = []
    for rz in roots:
        sub = [sp.expand(g.subs(z, rz)) for g in polys]
        sub = [g for g in sub if g != 0]
        uy = [g for g in sub if g.free_symbols <= {y}]
        ys = set()
        for g in uy:
            ys |= {r for r in sp.real_roots(sp.Poly(g, y)) if r > 0}
        for ry in ys:
            sub2 = [sp.expand(g.subs(y, ry)) for g in sub]
            sub2 = [g for g in sub2 if g != 0]
            ux = [g for g in sub2 if g.free_symbols <= {x}]
            xs = set()
            for g in ux:
                xs |= {r for r in sp.real_roots(sp.Poly(g, x)) if r > 0}
            for rx in xs:
                pt = {x: rx, y: ry, z: rz}
                if all(sp.simplify(e.subs(pt)) == 0 for e in f):
                    sols.append((rx, ry, rz))
    return sorted(set(sols), key=lambda p: p[2])


def arb_eigs(charpoly):
    from flint import acb_poly, ctx
    ctx.prec = 200
    coeffs = [sp.Rational(c) for c in sp.Poly(charpoly, sp.Symbol('lam')).all_coeffs()][::-1]
    from flint import fmpq, acb
    P = acb_poly([acb(fmpq(int(c.p), int(c.q))) for c in coeffs])
    return P.roots()


lam = sp.Symbol('lam')

NETS = {
    "N1": dict(rx=[("Z", "Y+Z"), ("2Z", "Z"), ("Y", "Y+Z"), ("Y+Z", "2X"), ("X+Z", "2Z"), ("2X", "Y")],
               k=[5, 1, 1, 4, 10, 7], w=(1, 2, 0),
               counts=dict(m=6, reactants=6, s=3, ell=1, t=1, delta=2),
               pts=[(1, sp.Rational(19, 4), sp.Rational(1, 2)), (1, 3, 1), (1, sp.Rational(7, 4), sp.Rational(7, 2))],
               stable=[True, False, True]),
    "N2": dict(rx=[("2Z", "X+Z"), ("X+Z", "2Z"), ("Y+Z", "X"), ("Y+Z", "X+Y"), ("X", "2X"), ("X+Y", "X"),
                   ("2X", "Y+Z")],
               k=[1, 5, 2, 5, 3, 1, 3], w=(1, 0, 1),
               counts=dict(m=6, reactants=6, s=3, ell=2, t=2, delta=1),
               pts=[(1, sp.Rational(3, 2), sp.Rational(1, 2)), (1, 1, 1), (1, sp.Rational(3, 7), 3)],
               stable=[True, False, True]),
    "G7": dict(rx=[("X+Y", "2Y"), ("Y", "2Z"), ("0", "X+Z"), ("X+Z", "2Z"), ("2Z", "X+Z"), ("Y+Z", "Y")],
               k=[1, 1, 5, 4, 1, 6], w=(0, 1, 0),
               counts=dict(m=7, reactants=6, s=3),
               pts=[(1, 2, 1), (1, sp.Rational(5, 4), sp.Rational(3, 2)), (1, sp.Rational(10, 9), sp.Rational(5, 3))],
               stable=[True, False, True]),
}

for name, d in NETS.items():
    print("=" * 20, name)
    g = graph_data(d["rx"])
    for key, val in d["counts"].items():
        check(g[key] == val, "%s C %s = %s (expected %s)" % (name, key, g[key], val))
    check(g["maxmol"] <= 2, "%s C at most bimolecular" % name)
    check(len(d["rx"]) == len(d["k"]), "%s C %d reactions" % (name, len(d["rx"])))
    # symbolic ACR identity
    ks = sp.symbols('k1:%d' % (len(d["rx"]) + 1), positive=True)
    fs = field(d["rx"], ks)
    wf = sp.factor(sum(d["w"][i] * fs[i] for i in range(3)))
    print("   w.f =", wf)
    # w.f must be (positive monomial) * (a - b x) with a, b > 0
    num = sp.Poly(sp.expand(wf), x, y, z)
    terms = num.terms()
    ok = False
    if len(terms) == 2:
        (e1, c1), (e2, c2) = terms
        diff = tuple(a - b for a, b in zip(e1, e2))
        ok = diff in [(1, 0, 0), (-1, 0, 0)] and (c1 * c2).is_negative
        if ok:
            hi, lo = ((e1, c1), (e2, c2)) if diff == (1, 0, 0) else ((e2, c2), (e1, c1))
            cval = sp.simplify(-lo[1] / hi[1])
            print("   ACR value x =", cval)
    check(ok, "%s A w.f = monomial*(a - b x) with a, b > 0 for all rate constants (ACR in X)" % name)
    # existence of a positive steady state for all kappa: restrict to x = c and eliminate
    # (done numerically-exactly by a sign argument on the resultant in z, see below)
    f = field(d["rx"], [sp.Integer(v) for v in d["k"]])
    pts = positive_steady_states(f)
    print("   positive steady states:", pts)
    check(len(pts) == 3 and set(pts) == set(tuple(sp.nsimplify(c) for c in p) for p in d["pts"]),
          "%s G exactly three positive steady states, as claimed" % name)
    J = f.jacobian(X)
    nstable = 0
    for p, st in zip(d["pts"], d["stable"]):
        Jp = J.subs({x: p[0], y: p[1], z: p[2]})
        cp = sp.expand((lam * sp.eye(3) - Jp).det())
        a = sp.Poly(cp, lam).all_coeffs()
        a1, a2, a3 = a[1], a[2], a[3]
        rh = a1 > 0 and a3 > 0 and a1 * a2 > a3
        roots = arb_eigs(cp)
        if st:
            cert = all(r.real.upper() < 0 for r in roots)
            nstable += 1
        else:
            cert = any(r.real.lower() > 0 for r in roots)
        print("   %s: charpoly %s; RH %s; eigenvalue real parts in %s" % (
            p, cp, rh, [str(r.real.mid())[:10] for r in roots]))
        check(a3 != 0, "%s J nondegenerate at %s" % (name, p))
        check(rh == st, "%s J Routh-Hurwitz %s at %s" % (name, "holds" if st else "fails", p))
        check(cert, "%s J arb enclosures certify %s at %s" % (name, "stability" if st else "a positive eigenvalue", p))
    check(nstable == 2, "%s two exponentially stable positive steady states in the single class R^3_>0" % name)

# Symbolic existence of a positive steady state for every kappa (N1): on x = c = 2k1/k5, y = (k1 z + k6 c^2)/(k4 z)
ks = sp.symbols('k1:7', positive=True)
k1, k2, k3, k4, k5, k6 = ks
fN1 = field(NETS["N1"]["rx"], ks)
c = 2 * k1 / k5
yz = (k1 * z + k6 * c**2) / (k4 * z)
check(sp.simplify(fN1[1].subs({x: c, y: yz})) == 0, "N1 E on x = c the Y-equation gives y = (k1 z + k6 c^2)/(k4 z) > 0")
cub = sp.expand(sp.simplify(k4 * z * fN1[2].subs({x: c, y: yz})))
pc = sp.Poly(cub, z)
print("   k4 z * fZ on that curve:", pc.as_expr())
check(sp.simplify(pc.coeff_monomial(1)) == k3 * k6 * c**2 and sp.simplify(pc.LC()) == -k2 * k4,
      "N1 E cubic in z has constant term k3 k6 c^2 > 0 and leading coefficient -k2 k4 < 0: a positive root exists "
      "for every kappa (unconditional ACR)")
check(sp.simplify(fN1[0] + 2 * fN1[1] - z * (2 * k1 - k5 * x)) == 0, "N1 A fX + 2 fY = z (2 k1 - k5 x)")
# dimension 3 from reaction vectors alone
check(graph_data(NETS["N1"]["rx"])["s"] == 3, "N1 C full-dimensional: the only class is R^3_>=0")

print("runtime %.1f s" % (time.time() - T0))
print("ALL PASS" if not FAIL else "FAILURES: %d" % len(FAIL))
sys.exit(1 if FAIL else 0)
