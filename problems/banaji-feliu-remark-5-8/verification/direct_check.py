#!/usr/bin/env python3
"""Independent verification (verifier's own code; does not use the Gale-duality reduction of the manuscript).

For each test network, written directly as a list of reactions (source complex -> product complex), with exact
rational rate constants, this script
  1. checks the hypotheses of Banaji-Feliu Theorem 5.5 from the reaction list: n species, n+2 reactions,
     rank Gamma = n, n+2 distinct sources of molecularity <= 2 (quadratic = sources only, BF Sect. 2),
     nonnegative products, rank [A|1] = n+1 (source rank n = r, i.e. not source deficient);
  2. builds the mass-action vector field g(x) = sum_r k_r x^{a_r} (b_r - a_r) and solves g = 0, t*x1*...*xn = 1
     with msolve (exact Groebner basis + certified real root isolation), i.e. all real equilibria with nonzero
     coordinates, counting the solution boxes that lie in the open positive orthant (every box is checked to be
     either strictly positive or to have a coordinate < 0);
  3. Singular: computes vdim of <g, t*x1*...*xn - 1> (number of solutions in (C^*)^n with multiplicity) and checks
     that <g, det Dg, t*x1*...*xn - 1> = <1>, i.e. NO equilibrium in (C^*)^n is degenerate (full rank, so this is
     exactly nondegeneracy in the sense of BF Sect. 2);
  4. PARI/GP: Sturm count of the positive roots of the eliminant of x1 (cross-check of the count in 2);
  5. prints the positive equilibria (Newton-refined, 30 digits) and det Dg there.
A positive equilibrium x gives the positive flux vector k*x^A in ker Gamma, so the network is dynamically nontrivial.

Environment: WSL Fedora, Python 3 + sympy + mpmath, msolve 0.10.1, Singular 4.4.1, PARI/GP 2.17 (single thread).
Runtime: about 1 minute.  Usage: python direct_check.py      Prints PASS/FAIL per claim; exit code 1 on failure.
"""
import os
import re
import subprocess
import sys
import tempfile
from fractions import Fraction as F

import mpmath as mp
import sympy as sp

FAILS = 0


def check(cond, msg):
    global FAILS
    print(("PASS " if cond else "FAIL ") + msg, flush=True)
    FAILS += 0 if cond else 1


# ---------------------------------------------------------------- networks (written from the reaction lists)
def chain_network(n, a, b):
    """0 -> sum a_i X_i ; 2X1 -> 2X1 + sum b_i X_i ; X_i + X_{i+1} -> X_{i+1} ; X_n -> 0 (species 0..n-1)."""
    R = [({}, {i: a[i] for i in range(n)})]
    prod = {i: b[i] for i in range(n)}
    prod[0] += 2
    R.append(({0: 2}, prod))
    for i in range(n - 1):
        R.append(({i: 1, i + 1: 1}, {i + 1: 1}))
    R.append(({n - 1: 1}, {}))
    return R


def vec(d, n):
    return [d.get(i, 0) for i in range(n)]


# ---------------------------------------------------------------- external tools
def run(cmd, inp=None):
    r = subprocess.run(cmd, input=inp, capture_output=True, text=True, timeout=1800)
    if r.returncode != 0:
        raise RuntimeError(f"{cmd} failed: {r.stderr[-2000:]}")
    return r.stdout


def msolve_real(polys_str, names, tmp):
    fin, fout = os.path.join(tmp, "in.ms"), os.path.join(tmp, "out.ms")
    with open(fin, "w") as fh:
        fh.write(",".join(names) + "\n0\n" + ",\n".join(polys_str) + "\n")
    run(["msolve", "-t", "1", "-f", fin, "-o", fout])
    txt = open(fout).read().strip().rstrip(":")
    txt = re.sub(r"(-?\d+)\s*/\s*2\^(\d+)", r"F(\1,2**\2)", txt)
    data = eval(txt, {"F": F, "__builtins__": {}})          # our own msolve output: nested lists of rationals
    if data == [-1]:
        return []
    assert data[0] == 0, "positive-dimensional solution set"
    boxes = data[1][1]
    return [[(F(lo), F(hi)) for lo, hi in box] for box in boxes]


def singular_checks(polys_str, names, tmp):
    n = len(names)
    xs = ",".join(names)
    prodx = "*".join(names)
    lines = [f"ring r = 0, (t,{xs}), dp;", "option(redSB);"]
    for i, p in enumerate(polys_str):
        lines.append(f"poly g{i+1} = {p};")
    gl = ",".join(f"g{i+1}" for i in range(n))
    lines += [f"ideal I = {gl}, t*{prodx} - 1;", "ideal G = std(I);", 'print("VDIM " + string(vdim(G)));',
              f"matrix J[{n}][{n}];"]
    for i in range(n):
        for j in range(n):
            lines.append(f"J[{i+1},{j+1}] = diff(g{i+1}, {names[j]});")
    lines += ["ideal K = I, det(J);", "ideal GK = std(K);",
              'if (reduce(1, GK) == 0) { print("NONDEG yes"); } else { print("NONDEG no"); }',
              f"ideal E = eliminate(I, t*{'*'.join(names[1:]) if n > 1 else '1'});",
              'print("ELIM " + string(size(E)));', 'print("POLY " + string(E[1]));', "quit;"]
    out = run(["Singular", "-q", "--no-rc"], "\n".join(lines) + "\n")
    vd = int(re.search(r"VDIM (-?\d+)", out).group(1))
    nondeg = "NONDEG yes" in out
    nE = int(re.search(r"ELIM (\d+)", out).group(1))
    poly = re.search(r"POLY (.*)", out).group(1).strip()
    return vd, nondeg, nE, poly


def pari_sturm(poly, var):
    prog = (f"p = {poly};\n"
            f"print(polsturm(p, [0, +oo]), \" \", poldegree(p, {var}), \" \", poldegree(gcd(p, deriv(p, {var})), {var}),"
            f" \" \", subst(p, {var}, 0) != 0);\nquit;\n")
    out = run(["gp", "-q", "-f", "--default", "parisizemax=1000000000"], prog).split()
    return int(out[0]), int(out[1]), int(out[2]), out[3] == "1"


# ---------------------------------------------------------------- main certification routine
def certify(label, n, reactions, kappa, expected):
    print(f"--- {label} ---")
    m = len(reactions)
    names = [f"x{i+1}" for i in range(n)]
    xs = sp.symbols(names, positive=True)
    A = sp.Matrix([vec(s, n) for s, _ in reactions])                    # exponent matrix (m x n)
    B = sp.Matrix([vec(p, n) for _, p in reactions])
    Gam = (B - A).T                                                      # stoichiometric matrix (n x m)
    def cplx(d):
        return " + ".join((f"{c}" if c > 1 else "") + f"X{i+1}" for i, c in sorted(d.items()) if c) or "0"
    print("reactions:", "; ".join(cplx(s) + " -> " + cplx(p) for s, p in reactions))
    print("rate constants:", [str(k) for k in kappa])
    check(m == n + 2, f"(n, m) = ({n}, {m}) with m = n+2")
    srcs = [tuple(A.row(j)) for j in range(m)]
    check(len(set(srcs)) == m and max(sum(s) for s in srcs) <= 2, "n+2 distinct sources, all of molecularity <= 2 (quadratic)")
    check(all(c >= 0 for c in B), f"products nonnegative (max product molecularity {max(sum(B.row(j)) for j in range(m))})")
    check(all(tuple(B.row(j)) != srcs[j] for j in range(m)), "each reaction has product != source")
    check(Gam.rank() == n, "rank Gamma = n (full rank, single stoichiometric class R^n_{>0})")
    check(A.row_join(sp.ones(m, 1)).rank() == n + 1, "rank [A|1] = n+1, source rank n = r (not source deficient)")
    # mass-action vector field with exact rational constants
    rate = [sp.Rational(kappa[j]) * sp.prod([xs[i] ** A[j, i] for i in range(n)]) for j in range(m)]
    g = [sp.expand(sum(Gam[i, j] * rate[j] for j in range(m))) for i in range(n)]
    gi = []
    for p in g:
        P = sp.Poly(p, *xs)
        den = sp.ilcm(*[sp.fraction(c)[1] for c in P.coeffs()])
        gi.append(sp.Poly(sp.expand(p * den), *xs))
    pstr = [str(P.as_expr()).replace("**", "^") for P in gi]
    with tempfile.TemporaryDirectory() as tmp:
        boxes = msolve_real(pstr + ["t*" + "*".join(names) + "-1"], names + ["t"], tmp)
        boxes = [b[:n] for b in boxes]                                   # drop the t coordinate
        pos = [b for b in boxes if all(lo > 0 for lo, _ in b)]
        nonpos = [b for b in boxes if any(hi < 0 for _, hi in b)]
        check(len(pos) + len(nonpos) == len(boxes), f"msolve: every real-solution box is strictly positive or has a negative coordinate ({len(boxes)} real solutions in (R*)^n)")
        check(len(pos) == expected, f"msolve: exactly {len(pos)} positive equilibria (expected {expected})")
        vd, nondeg, nE, poly = singular_checks(pstr, names, tmp)
        check(nondeg, f"Singular: <g, det Dg, t*x1*...*xn-1> = <1>: no degenerate equilibrium in (C*)^n ({vd} complex solutions in (C*)^n)")
        nst, deg, gdeg, nz0 = pari_sturm(poly, "x1")
        check(nE == 1 and nz0 and gdeg == 0 and nst == expected,
              f"PARI: eliminant in x1 (degree {deg}) squarefree, nonzero at 0, Sturm count of positive roots = {nst}")
    # numerical display (not part of the certificate)
    mp.mp.dps = 50
    Jsym = sp.Matrix(g).jacobian(xs)
    gf = sp.lambdify(xs, g, "mpmath")
    Jf = sp.lambdify(xs, Jsym, "mpmath")
    for b in pos:
        x0 = [mp.mpf(lo.numerator) / lo.denominator / 2 + mp.mpf(hi.numerator) / hi.denominator / 2 for lo, hi in b]
        if n > 1:
            sol = mp.findroot(lambda *v: gf(*v), x0)
            x = [sol[i] for i in range(n)]
        else:
            x = [mp.findroot(lambda v: gf(v)[0], x0[0])]
        d = mp.det(mp.matrix(Jf(*x)))
        print("   x =", [mp.nstr(t, 12) for t in x], " det Dg =", mp.nstr(d, 6))
    return len(pos)


if __name__ == "__main__":
    # (A) n = 1: 0 -> X (k=1), 2X -> 3X (k=1), X -> 0 (k=3):  x' = 1 + x^2 - 3x
    certify("n = 1: 0 -> X, 2X -> 3X, X -> 0", 1, chain_network(1, [1], [1]), [1, 1, 3], 2)
    # (B) small n = 3 example of the manuscript, kappa_1 = 1/1600 chosen by the verifier, all others 1
    certify("small n = 3 example (1,3,16)/(14,3,1), k(2X1 -> ...) = 1/1600", 3,
            chain_network(3, [1, 3, 16], [14, 3, 1]), [1, F(1, 1600), 1, 1, 1], 4)
    # negative control: same network, k = 1/100 (outside the window): fewer equilibria expected
    print("(negative control, a FAIL on the count is expected and NOT counted:)")
    saved = FAILS
    npos = certify("negative control: same network, k(2X1 -> ...) = 1/100", 3,
                   chain_network(3, [1, 3, 16], [14, 3, 1]), [1, F(1, 100), 1, 1, 1], 4)
    FAILS = saved
    check(npos < 4, f"negative control has {npos} < 4 positive equilibria (the test discriminates)")
    # (C) R = 64 family with EXACTLY the rate constants of the manuscript's proof: kappa_1 = exp(c),
    #     c = M_1 - L/2, M_1 = L(1 + 2 sum_m (-1)^m m) = -nL (n odd), (n+1)L (n even), so kappa_1 = 2^(-6n-3) or 2^(6n+3).
    for n in range(1, 11):
        k1 = F(1, 2 ** (6 * n + 3)) if n % 2 else F(2 ** (6 * n + 3))
        certify(f"R = 64 family, n = {n}, kappa_1 = {k1}", n,
                chain_network(n, [64 ** j for j in range(1, n + 1)], [1] * n), [1, k1] + [1] * n, n + 1)
    # (D) second (3,5,3) example of the manuscript: species X1=X, X2=Y, X3=Z
    #     Z -> 78Y, 2Z -> Z, 2Y -> X+Y, X -> 0, X+Y -> 24X+20Z ; k(2Z -> Z) = 9971888893/778, others 1
    R2 = [({2: 1}, {1: 78}), ({2: 2}, {2: 1}), ({1: 2}, {0: 1, 1: 1}), ({0: 1}, {}), ({0: 1, 1: 1}, {0: 24, 2: 20})]
    certify("second (3,5,3) example, k(2Z -> Z) = 9971888893/778", 3, R2, [1, F(9971888893, 778), 1, 1, 1], 4)
    print("ALL PASS" if FAILS == 0 else f"{FAILS} FAILURES")
    sys.exit(1 if FAILS else 0)
