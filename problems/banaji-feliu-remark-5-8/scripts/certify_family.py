#!/usr/bin/env python3
"""Exact certificate for the quadratic (n, n+2, n) networks N_n(alpha, beta) with n+1 positive nondegenerate
equilibria (answer to Banaji & Feliu, Remark 5.8).

Network N_n(alpha, beta), species X_1..X_n, reactions R_0..R_{n+1} with sources
    0, 2X_1, X_1+X_2, ..., X_{n-1}+X_n, X_n
and reaction vectors (rows of Gamma: row i = alpha_{i+1} e_0 + beta_{i+1} e_1 - e_{i+1}):
    R_0     : 0        -> sum_i alpha_{i+1} X_i
    R_1     : 2X_1     -> 2X_1 + sum_i beta_{i+1} X_i
    R_{i+1} : X_i+X_{i+1} -> X_{i+1}   (1 <= i <= n-1),      R_{n+1}: X_n -> 0.
For each test case the script checks in exact arithmetic (sympy):
  1. sources quadratic and distinct, products nonnegative, rank Gamma = n, rank [S|1] = n+1;
  2. ker Gamma = span(alpha, beta), alpha = (1,0,alpha_2..), beta = (0,1,beta_2..) > 0 elsewhere;
  3. the affine dependency w of the sources is (w_0, 1, -2, 2, -2, ...);
  4. a rational C such that P_C(u) = prod_{w_j>0} v_j^{w_j} - C prod_{w_j<0} v_j^{-w_j}, v = alpha + u beta,
     has exactly n+1 roots in (0, oo), all simple (Sturm count + gcd(P, P') = 1);
  5. numerically (mpmath, 80 digits): each root gives a positive equilibrium x with Gamma v(x) = 0 and det J(x) != 0.
Rate constants: kappa_j = 1 except kappa_1 = C (w_1 = 1), so that sum_j w_j log kappa_j = log C.
Usage: python certify_family.py            (cases built in; runtime ~1 min)
Prints PASS/FAIL, exits 1 on failure.
"""
import sys
import sympy as sp
import mpmath as mp

u = sp.symbols("u")
mp.mp.dps = 80
FAILS = 0


def check(c, msg):
    global FAILS
    print(("PASS " if c else "FAIL ") + msg)
    FAILS += (not c)


def network(n, al, be):
    """al, be: lists of length n (alpha_2..alpha_{n+1}, beta_2..beta_{n+1}). Returns sources S, Gamma."""
    S = [tuple([0] * n), tuple([2] + [0] * (n - 1))]
    for i in range(1, n):
        s = [0] * n; s[i - 1] = 1; s[i] = 1; S.append(tuple(s))
    s = [0] * n; s[n - 1] = 1; S.append(tuple(s))
    G = sp.zeros(n, n + 2)
    for i in range(n):              # species X_{i+1}
        G[i, 0] = al[i]
        G[i, 1] = be[i]
        G[i, i + 2] = -1
    return S, G


def certify(n, al, be, label):
    print(f"--- {label}: n = {n}, alpha = {al}, beta = {be} ---")
    S, G = network(n, al, be)
    m = n + 2
    A = sp.Matrix(S)
    prods = [[A[j, i] + G[i, j] for i in range(n)] for j in range(m)]
    check(len(set(S)) == m and all(sum(s) <= 2 for s in S), "sources: n+2 distinct complexes of molecularity <= 2")
    check(all(c >= 0 for p in prods for c in p), f"products nonnegative (max molecularity {max(sum(p) for p in prods)})")
    check(G.rank() == n, "rank Gamma = n (full rank)")
    check(A.row_join(sp.ones(m, 1)).rank() == n + 1, "rank [S|1] = n+1 (not source deficient)")
    alpha = sp.Matrix([1, 0] + list(al)); beta = sp.Matrix([0, 1] + list(be))
    check(G * alpha == sp.zeros(n, 1) and G * beta == sp.zeros(n, 1), "Gamma alpha = Gamma beta = 0, so ker Gamma = span(alpha, beta)")
    wv = A.row_join(sp.ones(m, 1)).T.nullspace()[0]
    wv = wv / wv[1]
    w = [int(x) for x in wv]
    check(w[1:] == [1] + [2 * (-1) ** (j - 1) for j in range(2, m)], f"affine dependency w = {w}")
    v = [alpha[j] + u * beta[j] for j in range(m)]
    Pp = sp.prod([v[j] ** w[j] for j in range(m) if w[j] > 0])
    Pm = sp.prod([v[j] ** (-w[j]) for j in range(m) if w[j] < 0])
    # critical values of f = log(Pp/Pm): roots of the numerator of f'
    fp = sp.together(sum(w[j] * beta[j] / v[j] for j in range(m) if beta[j] != 0))
    crit = [r for r in sp.Poly(sp.numer(fp), u).real_roots() if r > 0]
    check(len(crit) == n, f"f' has exactly {len(crit)} = n positive zeros")
    vals = [mp.mpf(sp.N(sp.log(Pp.subs(u, r) / Pm.subs(u, r)), 60)) for r in crit]
    maxima = [vals[i] for i in range(n) if i % 2 == 0]      # f' > 0 near 0, so the first critical point is a maximum
    minima = [vals[i] for i in range(n) if i % 2 == 1]
    lo = max(minima) if minima else min(maxima) - 4
    hi = min(maxima)
    check(lo < hi, "every local minimum of f lies below every local maximum")
    logC = (lo + hi) / 2
    C = sp.nsimplify(sp.Float(mp.nstr(mp.exp(logC), 30), 30), rational=True)
    P = sp.Poly(sp.expand(Pp - C * Pm), u)
    npos = P.count_roots(0, sp.oo) - (1 if P.eval(0) == 0 else 0)
    sqf = sp.gcd(P, P.diff(u)).degree() == 0
    check(npos == n + 1 and sqf, f"P_C has exactly n+1 = {npos} positive roots, all simple (C = rational with {len(str(C.p))}+{len(str(C.q))} digits)")
    # equilibria and nondegeneracy (numerical cross-check)
    roots = [r for r in P.real_roots() if r > 0]
    kappa = [mp.mpf(1)] * m
    kappa[1] = mp.mpf(sp.N(C, 60))
    ok = True
    for r in roots:
        uu = mp.mpf(sp.N(r, 60))
        vv = [mp.mpf(sp.N(alpha[j], 60)) + uu * mp.mpf(sp.N(beta[j], 60)) for j in range(m)]
        # solve S log x - log s = log(v/kappa) in least squares (consistent by construction)
        M = mp.matrix([[S[j][i] for i in range(n)] + [-1] for j in range(m)])
        rhs = mp.matrix([mp.log(vv[j] / kappa[j]) for j in range(m)])
        sol = mp.lu_solve(M.T * M, M.T * rhs)
        res = max(abs(x) for x in (M * sol - rhs))
        x = [mp.e ** sol[i] for i in range(n)]
        rates = [kappa[j] * mp.fprod([x[i] ** S[j][i] for i in range(n)]) for j in range(m)]
        Gm = mp.matrix([[int(G[i, j]) for j in range(m)] for i in range(n)])
        F = Gm * mp.matrix(rates)
        J = mp.matrix(n, n)
        for a in range(n):
            for b in range(n):
                J[a, b] = sum(Gm[a, j] * rates[j] * S[j][b] / x[b] for j in range(m))
        d = mp.det(J)
        scale = max(abs(t) for t in rates)
        ok &= res < mp.mpf(10) ** -60 and max(abs(t) for t in F) < mp.mpf(10) ** -55 * scale and abs(d) > 0
    check(ok, f"the {len(roots)} roots give positive equilibria with |Gamma v(x)| < 1e-55 max_j v_j and det J != 0")
    return C


if __name__ == "__main__":
    cases = [
        (3, [1, 3, 16], [14, 3, 1], "small n = 3 example"),
        (1, [1], [1], "n = 1 (0 -> X, 2X -> 3X, X -> 0)"),
    ]
    for n in range(1, 9):
        cases.append((n, [64 ** (j - 1) for j in range(2, n + 2)], [1] * n, f"R = 64 family, n = {n}"))
    for n, al, be, label in cases:
        certify(n, al, be, label)
    print("ALL PASS" if FAILS == 0 else f"{FAILS} FAILURES")
    sys.exit(1 if FAILS else 0)
