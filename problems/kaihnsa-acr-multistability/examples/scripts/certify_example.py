#!/usr/bin/env python3
"""
certify_example.py -- exact/interval certificate for a 3-species at-most-bimolecular network with ACR and
multistability (task E1).

Runtime: a few seconds (WSL, Python 3.12, sympy 1.14, python-flint 0.9).
Usage:   python certify_example.py [network-name]      (default: the network chosen in report.md)

What is certified (each line prints PASS/FAIL; exit code 1 on any FAIL):
  1. the network data: bimolecular, number of species / reactions / reactant complexes / complexes, linkage
     classes, rank of the stoichiometric matrix (= 3: full-dimensional, so the only compatibility class meeting
     the positive orthant is R^3_{>0} itself), deficiency;
  2. ACR (structural): dy/dt = y*(k1*x - k2) identically, so x = k2/k1 at every positive steady state, for every
     vector of rate constants;
  3. the positive steady states for the given rational kappa: exact lex Groebner basis of f_X, f_Z restricted to
     x = k2/k1 (shape {y - q(z), p(z)}), exact real-root isolation of p, certified sign of y = q(z) on each
     isolating interval (arb ball arithmetic) -> the complete list of positive steady states with enclosures;
  4. stability: Routh-Hurwitz coefficients a1, a2, a3 of det(lambda I - J) and a1*a2 - a3 enclosed with arb
     at each steady state; exponentially stable iff a1 > 0, a3 > 0, a1*a2 - a3 > 0; a3 < 0 gives a real positive
     eigenvalue (unstable); also eigenvalue enclosures (acb) are printed;
  5. robustness: the steady states are nondegenerate (a3 != 0) and simple roots of p, so by the implicit function
     theorem they persist, with the same stability type, for all kappa in an open neighbourhood; ACR holds for
     every kappa.
"""
import sys, time
import sympy as sp
import flint

flint.ctx.prec = 256
T0 = time.time()
x, y, z, lam = sp.symbols('x y z lam')
SPECIES = (x, y, z)
NAMES = ('X', 'Y', 'Z')

# network: list of (reactant, product) exponent vectors over (X, Y, Z); rate constants as rationals.
NETWORKS = {
    # chosen example (see report.md): 6 reactions, 6 reactant complexes, 7 complexes
    'G7': dict(rx=[((1, 1, 0), (0, 2, 0)),   # X+Y -> 2Y   k1
                   ((0, 1, 0), (0, 0, 2)),   # Y   -> 2Z   k2
                   ((0, 0, 0), (1, 0, 1)),   # 0   -> X+Z  k3
                   ((1, 0, 1), (0, 0, 2)),   # X+Z -> 2Z   k4
                   ((0, 0, 2), (1, 0, 1)),   # 2Z  -> X+Z  k5
                   ((0, 1, 1), (0, 1, 0))],  # Y+Z -> Y    k6
               k=["1", "1", "5", "4", "1", "6"]),
}
FAIL = []


def check(cond, msg):
    print(("PASS " if cond else "FAIL ") + msg)
    if not cond:
        FAIL.append(msg)


def cname(c):
    parts = []
    for i, a in enumerate(c):
        if a == 1:
            parts.append(NAMES[i])
        elif a == 2:
            parts.append('2' + NAMES[i])
    return '+'.join(parts) if parts else '0'


def to_arb(q):
    q = sp.Rational(q)
    return flint.arb(flint.fmpq(int(q.p), int(q.q)))


def arb_interval(lo, hi):
    return to_arb(lo).union(to_arb(hi))


def eval_arb(expr, env):
    """evaluate a polynomial expression with rational coefficients at arb values env = {sym: arb}."""
    P = sp.Poly(sp.expand(expr), *env.keys())
    tot = flint.arb(0)
    syms = list(env.keys())
    for mon, coef in P.terms():
        t = to_arb(coef)
        for s, e in zip(syms, mon):
            if e:
                t = t * env[s] ** e
        tot += t
    return tot


def refine(p, lo, hi, bits=200):
    """bisect the isolating interval [lo, hi] of the square-free polynomial p (exact rationals)."""
    lo, hi = sp.Rational(lo), sp.Rational(hi)
    slo = sp.sign(p.eval(lo))
    if slo == 0:
        return lo, lo
    while hi - lo > sp.Rational(1, 2 ** bits):
        mid = (lo + hi) / 2
        sm = sp.sign(p.eval(mid))
        if sm == 0:
            return mid, mid
        if sm == slo:
            lo = mid
        else:
            hi = mid
    return lo, hi


def main(name, kappa):
    net = NETWORKS[name]
    rx = net['rx']
    r = len(rx)
    ks = sp.symbols('k1:%d' % (r + 1), positive=True)
    print("network %s:" % name)
    for i, (a, b) in enumerate(rx):
        print("   R%d: %s -> %s   k%d = %s" % (i + 1, cname(a), cname(b), i + 1, kappa[i]))
    # 1. network data
    bimol = all(sum(a) <= 2 and sum(b) <= 2 for a, b in rx)
    reac = sorted(set(a for a, _ in rx))
    comps = sorted(set([a for a, _ in rx] + [b for _, b in rx]))
    G = sp.Matrix([[b[i] - a[i] for (a, b) in rx] for i in range(3)])
    rank = G.rank()
    # linkage classes
    import networkx as nx
    gr = nx.Graph()
    gr.add_nodes_from(comps)
    gr.add_edges_from(rx)
    ell = nx.number_connected_components(gr)
    dg = nx.DiGraph()
    dg.add_nodes_from(comps)
    dg.add_edges_from(rx)
    cond = nx.condensation(dg)
    term = sum(1 for v in cond.nodes if cond.out_degree(v) == 0)
    defi = len(comps) - ell - rank
    check(bimol, "at most bimolecular")
    print("     species 3, reactions %d, reactant complexes %d, complexes %d, linkage classes %d, terminal strong "
          "linkage classes %d, rank %d, deficiency %d" % (r, len(reac), len(comps), ell, term, rank, defi))
    check(rank == 3, "stoichiometric subspace = R^3 (full-dimensional; one compatibility class R^3_{>0})")
    # ODEs
    f = [0, 0, 0]
    for k, (a, b) in zip(ks, rx):
        mon = k * x ** a[0] * y ** a[1] * z ** a[2]
        for i in range(3):
            f[i] += (b[i] - a[i]) * mon
    f = [sp.expand(e) for e in f]
    for i in range(3):
        print("     d%s/dt = %s" % (NAMES[i].lower(), f[i]))
    # 2. ACR (structural)
    check(sp.expand(f[1] - y * (ks[0] * x - ks[1])) == 0,
          "structural ACR: dy/dt = y*(k1*x - k2) for all kappa, so x = k2/k1 at every positive steady state")
    # 3. steady states for the given kappa
    sub = {ks[i]: sp.Rational(kappa[i]) for i in range(r)}
    fk = [e.subs(sub) for e in f]
    xs = sp.Rational(kappa[1]) / sp.Rational(kappa[0])
    print("     ACR value x* = k2/k1 = %s" % xs)
    FX = sp.expand(fk[0].subs(x, xs))
    FZ = sp.expand(fk[2].subs(x, xs))
    GB = sp.groebner([FX, FZ], y, z, order='lex')
    polys = list(GB.exprs)
    pz = [g for g in polys if not g.has(y)]
    py = [g for g in polys if g.has(y)]
    shape = (len(pz) == 1 and len(py) == 1 and sp.Poly(py[0], y).degree() == 1)
    check(shape, "lex Groebner basis of {f_X, f_Z} at x = x* has shape {y - q(z), p(z)}")
    p = sp.Poly(pz[0], z)
    lc = sp.Poly(py[0], y).coeffs()[0]
    check(not lc.has(z) and lc != 0, "leading coefficient of the y-element is a nonzero constant")
    q = sp.expand(-(py[0] - lc * y) / lc)
    print("     p(z) = %s" % p.as_expr())
    print("     y = q(z) = %s" % q)
    psf = sp.Poly(sp.quo(p, sp.gcd(p, p.diff(z))), z)
    check(psf.degree() == p.degree(), "p is square-free (all roots simple)")
    ivs = psf.intervals(eps=sp.Rational(1, 2 ** 200))
    check(all(hi1 < lo2 for ((lo1, hi1), _), ((lo2, hi2), _) in zip(ivs, ivs[1:])),
          "isolating intervals (width < 2^-200) are pairwise disjoint")
    check(p.eval(0) != 0, "z = 0 is not a root of p")
    pos = []
    for (lo, hi), mult in ivs:
        if hi <= 0:
            continue
        check(lo > 0, "root interval does not straddle 0")
        if lo != hi:
            check(sp.sign(psf.eval(lo)) * sp.sign(psf.eval(hi)) < 0, "sign change of p on [%s, %s]" % (sp.N(lo, 8), sp.N(hi, 8)))
        else:
            check(psf.eval(lo) == 0, "exact rational root z = %s" % lo)
        zb = arb_interval(lo, hi)
        yb = eval_arb(q, {z: zb})
        if yb > 0:
            pos.append((lo, hi, zb, yb))
        elif yb < 0:
            print("     root z in [%s, %s] gives y < 0 (not positive)" % (sp.N(lo, 12), sp.N(hi, 12)))
        else:
            check(False, "sign of y undetermined")
    print("     number of real roots of p: %d; positive steady states: %d" % (len(ivs), len(pos)))
    # 4. stability
    J = sp.Matrix(fk).jacobian([x, y, z])
    cp = sp.Poly((lam * sp.eye(3) - J).det(), lam)
    co = cp.all_coeffs()
    assert co[0] == 1
    a1, a2, a3 = [sp.expand(c.subs(x, xs)) for c in co[1:]]
    nstable = 0
    for idx, (lo, hi, zb, yb) in enumerate(pos):
        env = {y: yb, z: zb}
        A1, A2, A3 = eval_arb(a1, env), eval_arb(a2, env), eval_arb(a3, env)
        H = A1 * A2 - A3
        stable = (A1 > 0) and (A3 > 0) and (H > 0)
        unstable = (A3 < 0) or (H < 0 and A1 > 0 and A3 > 0) or (A1 < 0)
        # eigenvalue enclosures
        cpol = flint.acb_poly([flint.acb(A3), flint.acb(A2), flint.acb(A1), flint.acb(1)])
        try:
            ev = cpol.roots()
        except Exception:
            ev = []
        print("   steady state %d: x = %s, y in %s, z in %s" % (idx + 1, xs, yb.str(12, radius=False),
                                                                  zb.str(12, radius=False)))
        print("       a1 = %s, a2 = %s, a3 = %s, a1*a2-a3 = %s" % (A1.str(15), A2.str(15), A3.str(15), H.str(15)))
        if ev:
            print("       eigenvalues: " + ", ".join(e.str(15) for e in ev))
        check(A3 != 0 and not A3.contains(0), "steady state %d nondegenerate (det J != 0, certified)" % (idx + 1))
        if stable:
            nstable += 1
            check(True, "steady state %d exponentially stable (Routh-Hurwitz a1>0, a3>0, a1*a2>a3 certified)" % (idx + 1))
        elif unstable:
            print("     steady state %d unstable (certified: %s)" % (idx + 1, "a3<0, real positive eigenvalue" if A3 < 0 else "RH violated"))
        else:
            check(False, "stability of steady state %d undecided" % (idx + 1))
    check(nstable >= 2, "two exponentially stable positive steady states in the single class R^3_{>0} "
                        "(multistability), together with ACR in X")
    print("runtime %.1f s" % (time.time() - T0))
    return nstable


if __name__ == "__main__":
    name = sys.argv[1] if len(sys.argv) > 1 else 'G7'
    kap = NETWORKS[name]['k']
    if kap is None:
        kap = sys.argv[2].split(',')
    main(name, kap)
    sys.exit(1 if FAIL else 0)
