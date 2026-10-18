"""
v5_msolve_uniqueness.py -- VERIFIER's independent check of Theorem B (uniqueness of the triple weak focus)
with a different tool: msolve (Groebner bases + certified real root isolation, exact over Q).

System: g_1 = g_2 = g_3 = 0 in (a, c, s), where G_k = (rational) * a^(2k) c^(2k) * g_k (the monomial factors
are positive on the domain).  {g1=g2=g3=0} has a positive-dimensional complex component (msolve reports
positive dimension even after saturating by a*c*s*rho; it lies in {s = -1}), so we saturate by
a*c*s*rho*(1+s) (rho = 1-(c-1)s; all factors > 0 on the domain) with an extra variable z:
z*a*c*s*rho*(1+s) - 1 = 0.  The resulting system is zero-dimensional.
msolve returns isolating boxes for ALL real solutions; we then keep the ones in the domain
a > 0, c > 0, s > 0, rho > 0 and compare with the Krawczyk box.
Environment: msolve 0.10.1, Python 3.12.  Runtime: seconds to minutes.
"""
import os
import sys
import json
import subprocess
import tempfile
from fractions import Fraction
from math import gcd
from functools import reduce

HERE = os.path.dirname(os.path.abspath(__file__))
CERT = os.path.join(HERE, "..", "certificates")
results = []


def report(name, ok, info=""):
    results.append(bool(ok))
    print(("PASS" if ok else "FAIL") + "  " + name + (("  [" + info + "]") if info else ""), flush=True)


def load(k):
    with open(os.path.join(CERT, "G_%d.json" % k)) as fh:
        return {(t[0], t[1], t[2]): Fraction(t[3]) for t in json.load(fh)["terms"]}


def primitive_g(P):
    ma = min(e[0] for e in P)
    mc = min(e[1] for e in P)
    ms = min(e[2] for e in P)
    Q = {(e[0] - ma, e[1] - mc, e[2] - ms): cf for e, cf in P.items()}
    den = reduce(lambda x, y: x * y // gcd(x, y), (cf.denominator for cf in Q.values()), 1)
    Z = {e: int(cf * den) for e, cf in Q.items()}
    g = reduce(gcd, (abs(v) for v in Z.values()))
    return {e: v // g for e, v in Z.items()}, (ma, mc, ms)


def to_str(P):
    terms = []
    for (i, j, k), cf in sorted(P.items()):
        mon = "*".join(["a^%d" % i] * (i > 0) + ["c^%d" % j] * (j > 0) + ["s^%d" % k] * (k > 0))
        terms.append("%+d%s" % (cf, ("*" + mon) if mon else ""))
    return "".join(terms).lstrip("+")


def parse_msolve(txt):
    """Parse the list of real solutions printed by msolve: [0, [1, [ [[lo,hi],[lo,hi],...], ... ]]]:"""
    txt = txt.strip().rstrip(':')
    import ast
    data = ast.literal_eval(txt.replace('\n', ''))
    return data


def main():
    gs = []
    for k in (1, 2, 3):
        g, mon = primitive_g(load(k))
        gs.append(g)
        print("     g_%d: %d terms, removed monomial a^%d c^%d s^%d" % (k, len(g), *mon))
    polys = [to_str(g) for g in gs]
    # z*a*c*s*rho*(1+s) - 1, expanded (msolve's parser takes no parentheses):
    polys.append("-a*c^2*s^3*z-a*c^2*s^2*z+a*c*s^3*z+2*a*c*s^2*z+a*c*s*z-1")
    import sympy as sp
    z_, a_, c_, s_ = sp.symbols("z a c s")
    assert sp.expand(z_ * a_ * c_ * s_ * (1 - (c_ - 1) * s_) * (1 + s_) - 1 -
                     sp.sympify(polys[-1].replace("^", "**"))) == 0
    with tempfile.TemporaryDirectory() as td:
        fin, fout = os.path.join(td, "in.ms"), os.path.join(td, "out.ms")
        with open(fin, "w") as fh:
            fh.write("z,a,c,s\n0\n" + ",\n".join(polys) + "\n")
        r = subprocess.run(["msolve", "-f", fin, "-o", fout, "-t", "3"], capture_output=True, text=True)
        out = open(fout).read()
    if out.strip().startswith("[1, 4, -1"):
        report("msolve: system is zero-dimensional", False, "positive-dimensional")
        return False
    # msolve real-root output is a nested list of rational endpoints written as numbers / 2^k fractions
    txt = out.strip().rstrip(':')
    import re
    txt = re.sub(r"(\d+)\s*\^\s*(\d+)", lambda mo: str(int(mo.group(1)) ** int(mo.group(2))), txt)
    # evaluate with Fractions safely: tokens are integers, '/', '[', ']', ',', '-'
    tokens = re.findall(r"-?\d+|[\[\],/]", txt)
    pos = 0

    def parse():
        nonlocal pos
        tok = tokens[pos]
        if tok == '[':
            pos += 1
            lst = []
            while tokens[pos] != ']':
                lst.append(parse())
                if tokens[pos] == ',':
                    pos += 1
            pos += 1
            return lst
        pos += 1
        val = Fraction(int(tok))
        if pos < len(tokens) and tokens[pos] == '/':
            pos += 1
            val = val / int(tokens[pos])
            pos += 1
        return val
    data = parse()
    # expected structure: [0, [1, [ sol1, sol2, ... ]]] with sol = [[lo,hi] per variable]
    sols = data[1][1]
    print("     msolve: %d real solutions of the saturated system (variables z, a, c, s)" % len(sols))
    dom = []
    for sol in sols:
        (zl, zh), (al, ah), (cl, ch), (sl, sh) = sol
        # rho = 1-(c-1)s : bound over the box (monotone pieces; just test corners)
        rhos = [1 - (cv - 1) * sv for cv in (cl, ch) for sv in (sl, sh)]
        indom = al > 0 and cl > 0 and sl > 0 and min(rhos) > 0
        maybe = ah > 0 and ch > 0 and sh > 0 and max(rhos) > 0
        print("       a in [%.12f, %.12f], c in [%.12f, %.12f], s in [%.12f, %.12f]  in domain: %s%s"
              % (float(al), float(ah), float(cl), float(ch), float(sl), float(sh), indom,
                 "" if (indom or not maybe) else "  (UNDECIDED)"))
        if maybe and not indom:
            report("no undecided solution box", False)
        if indom:
            dom.append(sol)
    with open(os.path.join(CERT, "box.json")) as fh:
        box = json.load(fh)
    m = [Fraction(v) for v in box["centre"]]
    rad = Fraction(box["radius"])
    ok = len(dom) == 1
    if ok:
        (_, _), (al, ah), (cl, ch), (sl, sh) = dom[0]
        # the box of msolve and the Krawczyk box must intersect (both contain the same isolated zero)
        ok = (al <= m[0] + rad and ah >= m[0] - rad and cl <= m[1] + rad and ch >= m[1] - rad
              and sl <= m[2] + rad and sh >= m[2] - rad)
    report("Theorem B via msolve: exactly one real solution in the domain, and it is the certified zero", ok,
           "%d solution(s) in the domain" % len(dom))
    print("ALL PASS" if all(results) else "SOME CHECKS FAILED")
    return all(results)


if __name__ == "__main__":
    sys.exit(0 if main() else 1)
