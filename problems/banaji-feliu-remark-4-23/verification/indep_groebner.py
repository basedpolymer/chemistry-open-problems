#!/usr/bin/env python3
"""Independent verification (verifier's own code) of the family N_n, Banaji-Feliu Remark 4.23.

Written from the reaction table of paper.tex only ("Table 1" = label tab:reactions, numbered Table 3 in the
PDF), without reusing the solver's scripts.
Method differs from the solver's: the equilibria on the stoichiometric class are NOT computed from the
monomial curve x = C o lambda^q.  Instead, for each test case the full polynomial system

        Gamma v(x) = 0  (all n rows),   z . x - T = 0,    z = computed left kernel vector of Gamma,

is saturated by x1*...*xn in Singular (positive equilibria lie in the torus), a lex Groebner basis is
computed, the shape lemma {g(x1), x_i - h_i(x1)} is checked, and the positive real roots of g are counted
exactly with Sturm sequences (sympy).  Nondegeneracy of each equilibrium is checked in two ways:
  (a) e_{n-1}(J) = sum of the principal (n-1)-minors of the Jacobian J != 0  (since J R^n is inside im Gamma,
      which is J-invariant with codimension 1, charpoly(J)(t) = t * charpoly(J|im Gamma)(t), so
      e_{n-1}(J) = det(J|im Gamma));
  (b) the stacked matrix [Jacobian of Gamma v ; z^T] has rank n (simple solution of the system).

Rate constants are built with the recipe of the proof (generalized-Vandermonde kernel), for
lambda_r = r and for other rational lambdas; the count is then checked from scratch.

Environment: WSL Fedora, Python 3 + sympy 1.14, Singular 4.4.1.  Runtime: about 4 s, 75 MB (n <= 6).
Prints PASS/FAIL for each claim and exits 1 on any failure.
Usage: python verification/indep_groebner.py
"""
import subprocess
import sys
import tempfile
from fractions import Fraction as Fr
from itertools import combinations

import sympy as sp

FAILS = 0


def check(cond, msg):
    global FAILS
    print(("PASS " if cond else "FAIL ") + msg, flush=True)
    if not cond:
        FAILS += 1


# ---------------------------------------------------------------- network from Table 1 of paper.tex
def table_network(n):
    """Reactions (source, product) as dicts {species index (1-based): stoichiometry}, from Table 1."""
    eps = lambda i: (-1) ** ((i - 1) // 2)
    s = -eps(n)
    R = []
    # R1: s=+1: X1 -> 2X1 + Xn ; s=-1: X1 -> Xn
    R.append(({1: 1}, {1: 2, n: 1} if s == 1 else {n: 1}))
    for j in range(2, n + 1):
        if j % 2 == 0:
            k = j // 2
            c = (-1) ** (k + 1) * s
            src = {2 * k - 1: 2, 2 * k: 1}
            prod = {2 * k - 1: 1, 2 * k: 2} if c == 1 else {2 * k - 1: 3}
        else:
            k = (j - 1) // 2
            d = (-1) ** k * s
            src = {2 * k: 2, 2 * k + 1: 1}
            prod = {2 * k: 3, 2 * k + 1: 2} if d == 1 else {2 * k: 1}
        R.append((src, prod))
    return R


def matrices(n, R):
    A = sp.zeros(n, n)      # rows = sources
    G = sp.zeros(n, n)      # columns = reaction vectors
    for j, (src, prod) in enumerate(R):
        for i in range(1, n + 1):
            A[j, i - 1] = src.get(i, 0)
            G[i - 1, j] = prod.get(i, 0) - src.get(i, 0)
    return A, G


def complex_str(c):
    t = [(str(v) if v != 1 else "") + f"X{i}" for i, v in sorted(c.items()) if v]
    return " + ".join(t) if t else "0"


# ---------------------------------------------------------------- Singular: saturated lex GB
def singular_gb(n, kappa, A, G, z, T):
    xs = [f"x{i}" for i in range(1, n + 1)]
    v = []
    for j in range(n):
        mon = "*".join(f"x{i+1}^{A[j, i]}" for i in range(n) if A[j, i])
        v.append(f"({kappa[j]})*{mon}")
    eqs = []
    for i in range(n):
        terms = [f"({G[i, j]})*{v[j]}" for j in range(n) if G[i, j] != 0]
        eqs.append(" + ".join(terms) if terms else "0")
    cons = " + ".join(f"({z[i]})*x{i+1}" for i in range(n)) + f" - ({T})"
    order = ",".join(reversed(xs))          # x_n > ... > x_1 in lp: eliminant in x1
    prod = "*".join(xs)
    code = f"""
LIB "elim.lib";
ring r = 0, ({order}), lp;
short = 0;
option(redSB);
ideal I = {", ".join(eqs)}, {cons};
def S = sat(I, {prod});
ideal Js;
if (typeof(S) == "list") {{ Js = S[1]; }} else {{ Js = S; }}
Js = std(Js);
print("VDIM"); print(vdim(Js));
print("GBSTART");
int k;
for (k = 1; k <= size(Js); k++) {{ print(string(Js[k]) + ";"); }}
print("GBEND");
quit;
"""
    with tempfile.NamedTemporaryFile("w", suffix=".sing", delete=False) as f:
        f.write(code)
        fname = f.name
    out = subprocess.run(["Singular", "-q", "--no-rc", fname], capture_output=True, text=True,
                         timeout=1800).stdout
    lines = out.split("\n")
    vdim = int(lines[lines.index("VDIM") + 1].strip())
    i0, i1 = lines.index("GBSTART"), lines.index("GBEND")
    gb = "".join(lines[i0 + 1:i1]).split(";")
    gb = [g.strip() for g in gb if g.strip()]
    return vdim, gb


# ---------------------------------------------------------------- one test case
def run_case(n, lam, label, R=None):
    print(f"--- n = {n}, {label}: lambda = {[str(l) for l in lam]} ---")
    R = table_network(n) if R is None else R
    for j, (a, b) in enumerate(R, start=1):
        print(f"  R{j}: {complex_str(a)} -> {complex_str(b)}")
    A, G = matrices(n, R)
    check(all(v >= 0 for (_, b) in R for v in b.values()), "products nonnegative")
    check(max(sum(a.values()) for a, _ in R) == 3, "source molecularity (max over sources) = 3")
    check(len({tuple(sorted(a.items())) for a, _ in R}) == n
          and all(a != b for a, b in R), "n distinct reactions, source != product")
    check(G.rank() == n - 1, "n species, n reactions, rank Gamma = n-1")
    K = G.nullspace()
    k = K[0] / K[0][0]
    check(len(K) == 1 and all(x > 0 for x in k), f"ker Gamma = span of positive vector {list(k)}")
    L = G.T.nullspace()
    z = L[0] / L[0][0]
    check(len(L) == 1, f"left kernel spanned by z = {list(z)}")
    q = A.LUsolve(sp.ones(n, 1))
    check(A.det() != 0 and all(x.is_integer for x in q), f"q = S^-1 1 = {list(q)}")
    # rate constants by the recipe of the proof
    expo = sorted(list(q) + [0])
    M = sp.Matrix([[sp.Rational(l) ** e for e in expo] for l in lam])
    ker = M.nullspace()
    check(len(ker) == 1, "generalized Vandermonde kernel 1-dim")
    c = ker[0]
    pos = {e: i for i, e in enumerate(expo)}
    if (c[pos[q[0]]] > 0) != (z[0] > 0):
        c = -c
    C = [c[pos[q[i]]] / z[i] for i in range(n)]
    T = -c[pos[0]]
    check(all(x > 0 for x in C) and T > 0, "C_i > 0, T > 0 (sign pattern of z matches)")
    kappa = [sp.prod([C[i] ** (-A[j, i]) for i in range(n)]) for j in range(n)]
    print(f"  kappa = {[str(x) for x in kappa]}" if n <= 4 else "  (kappa computed)")
    # Groebner count, no use of the monomial curve
    vdim, gb = singular_gb(n, kappa, A, G, z, T)
    X = sp.symbols(f"x1:{n+1}")
    loc = {f"x{i+1}": X[i] for i in range(n)}
    polys = [sp.Poly(sp.sympify(g, locals=loc), *X) for g in gb]
    uni = [p for p in polys if p.free_symbols <= {X[0]}]
    lin = {}
    for p in polys:
        for i in range(1, n):
            if p.degree(X[i]) == 1 and p.free_symbols <= {X[0], X[i]} and \
                    sp.Poly(p.as_expr(), X[i]).coeffs()[0].is_number:
                lin[i] = p
    shape = len(uni) == 1 and len(lin) == n - 1 and len(polys) == n
    check(shape, f"saturated ideal has shape-lemma lex basis (|GB| = {len(polys)}, vdim = {vdim})")
    g = sp.Poly(uni[0].as_expr(), X[0])
    check(g.degree() == vdim, f"deg g(x1) = vdim = {vdim} (= max q - min q = {max(q) - min(q)})")
    check(sp.gcd(g, g.diff(X[0])).degree() == 0, "g squarefree: all torus solutions simple")
    npos = g.count_roots(0, None)
    check(g.eval(0) != 0 and npos == n, f"g has exactly {npos} positive real roots (Sturm)")
    # the n claimed equilibria: rational points; express x_i from the shape basis
    pts = []
    for l in lam:
        x1 = C[0] * sp.Rational(l) ** q[0]
        pt = [x1]
        for i in range(1, n):
            sol = sp.solve(lin[i].as_expr().subs(X[0], x1), X[i])
            pt.append(sol[0])
        pts.append(pt)
    ok_roots = all(g.eval(p[0]) == 0 for p in pts) and len({p[0] for p in pts}) == n
    ok_pos = all(all(x > 0 for x in p) for p in pts)
    check(ok_roots and ok_pos, "the n positive roots of g give n points with all coordinates > 0")
    ok_curve = all(all(pts[r][i] == C[i] * sp.Rational(lam[r]) ** q[i] for i in range(n)) for r in range(n))
    check(ok_curve, "these points coincide with x(lambda_r) = C o lambda_r^q of the manuscript")
    # direct substitution and nondegeneracy
    v = [kappa[j] * sp.prod([X[i] ** A[j, i] for i in range(n)]) for j in range(n)]
    f = G * sp.Matrix(v)
    Jf = f.jacobian(X)
    for r, p in enumerate(pts):
        sub = dict(zip(X, p))
        eq = all(sp.simplify(e.subs(sub)) == 0 for e in f)
        cl = sum(z[i] * p[i] for i in range(n)) == T
        Jp = Jf.subs(sub)
        e1 = sum(Jp.extract(list(S), list(S)).det() for S in combinations(range(n), n - 1))
        stacked = Jp.col_join(sp.Matrix([list(z)]))
        check(eq and cl and e1 != 0 and stacked.rank() == n,
              f"  equilibrium {r+1}: Gamma v = 0, z.x = T, det(J|imGamma) = e_(n-1)(J) = "
              f"{sp.N(e1, 6)} != 0, rank[J; z] = n")
        if n == 3:
            tr = Jp.trace()
            print(f"    x = {p}, trace = {tr} ({sp.N(tr, 6)}), det(J|imGamma) = {e1}, "
                  f"eigs = {[sp.N(e, 6) for e in sp.Poly(sp.Symbol('t')**2 - tr*sp.Symbol('t') + e1).nroots()]}")


def main():
    cases = []
    for n in (3, 4, 5, 6):
        cases.append((n, [Fr(r) for r in range(1, n + 1)], "lambda_r = r"))
    cases.append((3, [Fr(1, 2), Fr(7, 5), Fr(4)], "other lambdas"))
    cases.append((4, [Fr(1, 3), Fr(1), Fr(5, 2), Fr(3)], "other lambdas"))
    cases.append((5, [Fr(2, 7), Fr(1, 2), Fr(1), Fr(9, 4), Fr(3)], "other lambdas"))
    cases.append((5, [Fr(9, 10), Fr(1), Fr(11, 10), Fr(6, 5), Fr(13, 10)], "clustered lambdas"))
    for n, lam, label in cases:
        run_case(n, lam, label)
    print("ALL PASS" if FAILS == 0 else f"{FAILS} FAILURES")
    sys.exit(1 if FAILS else 0)


if __name__ == "__main__":
    main()
