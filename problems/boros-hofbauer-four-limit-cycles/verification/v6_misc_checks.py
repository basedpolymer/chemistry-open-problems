"""
v6_misc_checks.py -- VERIFIER's checks of factual statements in paper.tex (exact, sympy):
  M1 G_k = kappa_k a^(2k) c^(2k) g_k with kappa = (-1/8, 1/1152, -1/442368, 1/530841600), g_k in Z[a,c,s]
     primitive, irreducible over Q, degrees (2k, 2k, 3k-1) in (a, c, s), term counts 9, 48, 136, 296;
  M2 the g_1, g_2 printed in the appendix of paper.tex coincide with the certificate files;
  M3 the mass-action form of the Corollary: with kappa1 = d-1, kappa2 = 1, kappa3 = K^(-a), the system
     x' = kappa2/y - kappa3 x^a y^b, y' = -kappa1 - kappa2/y + d kappa3 x^a y^b has the equilibrium (K, 1), and
     x = K X, y = Y, t = K tau maps it to (BH);
  M4 the explicit parameters of Theorem C: tr J < 0 < det J, tr^2 < 4 det, distance to tr = 0 (exact).
Environment: Python 3.12, sympy 1.14.  Runtime: ~1 min (factorisations).
"""
import os
import re
import sys
import json
import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "..")
results = []
a, c, s = sp.symbols('a c s')


def report(name, ok, info=""):
    results.append(bool(ok))
    print(("PASS" if ok else "FAIL") + "  " + name + (("  [" + info + "]") if info else ""), flush=True)


def loadG(k):
    with open(os.path.join(ROOT, "certificates", "G_%d.json" % k)) as fh:
        data = json.load(fh)
    return sp.Add(*[sp.Rational(t[3]) * a**t[0] * c**t[1] * s**t[2] for t in data["terms"]])


def main():
    kap = [sp.Rational(-1, 8), sp.Rational(1, 1152), sp.Rational(-1, 442368), sp.Rational(1, 530841600)]
    gs = {}
    for k in range(1, 5):
        G = sp.Poly(loadG(k), a, c, s)
        g = sp.Poly(sp.expand(G.as_expr() / (kap[k - 1] * a**(2 * k) * c**(2 * k))), a, c, s)
        integral = all(cf.is_integer for cf in g.coeffs())
        prim = integral and sp.gcd_list([int(x) for x in g.coeffs()]) == 1
        fl = sp.factor_list(g.as_expr())
        irred = len(fl[1]) == 1 and fl[1][0][1] == 1
        degs = (g.degree(a), g.degree(c), g.degree(s))
        ok = prim and irred and degs == (2 * k, 2 * k, 3 * k - 1) and len(g.terms()) == [9, 48, 136, 296][k - 1]
        report("M1 g_%d: integral primitive, irreducible over Q, degrees %s, %d terms" % (k, degs, len(g.terms())), ok,
               "factor_list constant %s" % fl[0])
        gs[k] = g.as_expr()
    tex = open(os.path.join(ROOT, "paper.tex"), encoding="utf-8").read()
    for k in (1, 2):
        m = re.search(r"\\texttt\{g\\_%d\\ =\\ (.*)\}\s*$" % k, tex, re.M)
        txt = m.group(1).replace("\\^{}", "^").replace("\\ ", " ").replace("^", "**")
        gp = sp.sympify(txt, locals={'a': a, 'c': c, 's': s})
        report("M2 appendix g_%d of paper.tex == certificate" % k, sp.expand(gp - gs[k]) == 0)
    # M3
    x, y, X, Y, K, b, d = sp.symbols('x y X Y K b d', positive=True)
    k1, k2, k3 = d - 1, 1, K**(-a)
    fx = k2 / y - k3 * x**a * y**b
    fy = -k1 - k2 / y + d * k3 * x**a * y**b
    eq = [sp.simplify(fx.subs({x: K, y: 1})), sp.simplify(fy.subs({x: K, y: 1}))]
    # dX/dtau = (dx/dt)(dt/dtau)/K = fx ; dY/dtau = K fy
    FX = sp.simplify(sp.powsimp(sp.expand_power_base(fx.subs({x: K * X, y: Y}), force=True), force=True))
    FY = sp.simplify(sp.powsimp(sp.expand_power_base((K * fy).subs({x: K * X, y: Y}), force=True), force=True))
    okX = sp.simplify(FX - (1 / Y - X**a * Y**b)) == 0
    okY = sp.simplify(FY - K * (-(d - 1) - 1 / Y + d * X**a * Y**b)) == 0
    report("M3 Corollary: (K,1) is an equilibrium and x=KX, y=Y, t=K tau gives (BH)", eq == [0, 0] and okX and okY)
    # M4
    with open(os.path.join(ROOT, "certificates", "explicit_params.json")) as fh:
        par = json.load(fh)
    av, bv, dv, Kv = (sp.Rational(par[k]) for k in ("a", "b", "d", "K"))
    tr = -av + Kv * (1 + bv * dv)
    det = Kv * av * (dv - 1)
    report("M4 Theorem C parameters: tr J < 0 < det J, tr^2 < 4 det, in the domain",
           tr < 0 and det > 0 and tr**2 < 4 * det and av > 0 and bv > -1 and dv > 1 and 1 + bv * dv > 0,
           "tr J = %s, K - a/(1+bd) = %s" % (sp.N(tr, 6), sp.N(Kv - av / (1 + bv * dv), 6)))
    print("ALL PASS" if all(results) else "SOME CHECKS FAILED")
    return all(results)


if __name__ == "__main__":
    sys.exit(0 if main() else 1)
