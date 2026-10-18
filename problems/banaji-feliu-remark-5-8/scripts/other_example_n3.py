#!/usr/bin/env python3
# Second, unrelated quadratic (3,5,3) example found by a computer search (sources Z, 2Z, 2Y, X, X+Y):
#   Z -> 78Y, 2Z -> Z, 2Y -> X+Y, X -> 0, X+Y -> 24X+20Z.
# Positive equilibria <-> roots u>0 of prod v_j^{w_j} = C with C = prod kappa_j^{w_j} (Gale duality); e.g. kappa_j = 1
# except kappa_2 = 1/C (w_2 = -1). The script finds the extremal rays of ker Gamma and a rational C with 4 simple roots.
"""Exact verification of a quadratic (3,5,3) network with 4 positive nondegenerate equilibria.

Input: sources S (5 triples), Gamma (3x5 integer, columns = reactions). Steps (all exact, sympy rationals):
 1. checks: distinct quadratic sources, products >= 0, rank Gamma = 3, rank [S|1] = 4;
 2. affine dependency w (w^T S = 0, sum w = 0);
 3. extremal rays alpha, beta of ker Gamma cap R^5_{>=0} (rational);
 4. P_C(u) = prod_{w>0} v_j^{w_j} - C prod_{w<0} v_j^{-w_j}, v = alpha + u beta; picks rational C giving 4 simple
    positive roots (Sturm count on (0, oo)) and checks squarefreeness.
"""
import sympy as sp, itertools, sys

u = sp.symbols("u")


def analyse(S, G, verbose=True):
    S = [tuple(s) for s in S]
    Gm = sp.Matrix(G)
    A = sp.Matrix(S)                              # 5x3, rows = sources
    out = {}
    assert len(set(S)) == 5 and all(sum(s) <= 2 for s in S), "sources must be 5 distinct quadratic complexes"
    prods = [tuple(A[j, i] + Gm[i, j] for i in range(3)) for j in range(5)]
    assert all(c >= 0 for p in prods for c in p), "negative product"
    assert Gm.rank() == 3, "rank Gamma != 3"
    M = A.row_join(sp.ones(5, 1))
    assert M.rank() == 4, "rank [S|1] != 4"
    w = M.T.nullspace()[0]
    w = w * sp.ilcm(*[sp.fraction(x)[1] for x in w])
    g = sp.igcd(*[int(x) for x in w]); w = [int(x) // g for x in w]
    K = Gm.nullspace()
    assert len(K) == 2
    k1, k2 = [sp.Matrix(v) for v in K]
    # extremal rays: vectors of K with one zero coordinate that are >= 0
    rays = []
    for j in range(5):
        # combination s k1 + t k2 with coordinate j = 0
        s, t = k2[j], -k1[j]
        if s == 0 and t == 0:
            continue
        v = s * k1 + t * k2
        for sign in (1, -1):
            vv = sign * v
            if all(x >= 0 for x in vv) and any(x > 0 for x in vv):
                q = sp.ilcm(*[sp.fraction(x)[1] for x in vv]); vv = vv * q
                gg = sp.igcd(*[int(x) for x in vv]); vv = vv / gg
                if not any(list(vv) == list(r) for r in rays):
                    rays.append(vv)
    assert len(rays) == 2, f"expected 2 extremal rays, got {len(rays)}"
    al, be = rays
    assert all(al[j] + be[j] > 0 for j in range(5)), "positive cone not 2-dimensional"
    v = [al[j] + u * be[j] for j in range(5)]
    Pp = sp.prod([v[j] ** w[j] for j in range(5) if w[j] > 0])
    Pm = sp.prod([v[j] ** (-w[j]) for j in range(5) if w[j] < 0])
    f = sum(w[j] * sp.log(v[j]) for j in range(5))
    fp = sp.together(sp.diff(f, u))
    num = sp.Poly(sp.numer(fp), u)
    crit = [r for r in sp.real_roots(num) if r > 0]
    out.update(dict(sources=S, products=prods, w=w, alpha=list(al), beta=list(be), crit=crit))
    if verbose:
        print("w =", w, " alpha =", list(al), " beta =", list(be))
        print("critical points of f on (0,oo):", [sp.N(r, 8) for r in crit])
    if len(crit) < 3:
        return out
    # choose C: f = log(Pp/Pm); values at critical points
    R = [sp.N(Pp.subs(u, r) / Pm.subs(u, r), 30) for r in crit]
    lo, hi = sorted(R)[0], sorted(R)[1]
    # 4 roots for log C strictly between the middle critical value and the nearest outer one: try rational midpoints
    cands = sorted(R)
    best = None
    for a_, b_ in zip(cands, cands[1:]):
        Cq = sp.nsimplify((a_ + b_) / 2, rational=True, tolerance=abs(b_ - a_) / 1000)
        P = sp.Poly(sp.expand(Pp - Cq * Pm), u)
        n = P.count_roots(0, sp.oo) - (1 if P.eval(0) == 0 else 0)
        sqf = sp.degree(sp.gcd(P, P.diff(u)), u) == 0
        if verbose:
            print(f"  C = {Cq}: positive roots = {n}, squarefree = {sqf}")
        if n == 4 and sqf:
            best = (Cq, P)
    out["C"] = best
    return out


if __name__ == "__main__":
    S = [(0, 0, 1), (0, 0, 2), (0, 2, 0), (1, 0, 0), (1, 1, 0)]
    G = [[0, 0, 1, -1, 23], [78, 0, -1, 0, -1], [-1, -1, 0, 0, 20]]
    o = analyse(S, G)
    print("products:", o["products"])
    print("4 equilibria:", o.get("C") is not None)
