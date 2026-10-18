#!/usr/bin/env python3
"""Complementary independent checks (verifier's own code), Banaji-Feliu Remark 4.23.

A. Structural lemmas for N_n, n = 2..80, from Table 1 of paper.tex, exact (python-flint):
   rank Gamma = n-1, Gamma 1 = 0, eps^T Gamma = 0, S lower bidiagonal with unit diagonal,
   q = S^-1 1 = (1 + 2(-2)^(i-1))/3, ordering of the exponents used in the proof.
B. Groebner/Sturm count (method of indep_groebner.py) for
   - N_2 (the family at n = 2, which is the authors' example X + 2Y -> 3Y, Y -> X with X <-> Y),
   - the network displayed as "N_4" in paper.tex line 130 before correction: it is shown to have
     rank 4 (so it is not a (4,4,3) network), and the corrected display is checked against the table.
"Table 1" below means the reaction table of paper.tex (label tab:reactions; numbered Table 3 in the PDF).
C. Construction of the proof with non-integer lambdas for n = 2..12 (exact rationals):
   x(lambda_r) are positive equilibria on eps.x = T, det(J|im Gamma) = e_(n-1)(J) != 0, the coefficient
   vector alternates (Descartes), and, for n <= 7, a Sturm count of the positive roots of
   lambda^(-min q) F(lambda) gives exactly n.
Environment: WSL Fedora, Python 3, python-flint 0.9, sympy 1.14, Singular 4.4.1.  Runtime: about 3 min, 110 MB.
Prints PASS/FAIL per claim, exits 1 on failure.
"""
import sys
from itertools import combinations

import sympy as sp
from flint import fmpq, fmpq_mat

import indep_groebner as ig

FAILS = 0


def check(cond, msg):
    global FAILS
    print(("PASS " if cond else "FAIL ") + msg, flush=True)
    if not cond:
        FAILS += 1


def eps(i):
    return (-1) ** ((i - 1) // 2)


def mats(n, R):
    A = [[R[j][0].get(i + 1, 0) for i in range(n)] for j in range(n)]
    G = [[R[j][1].get(i + 1, 0) - R[j][0].get(i + 1, 0) for j in range(n)] for i in range(n)]
    return A, G


def part_a():
    print("=== A. structural lemmas, n = 2..80 ===")
    ok = True
    for n in range(2, 81):
        R = ig.table_network(n)
        A, G = mats(n, R)
        Gm = fmpq_mat(n, n, [G[i][j] for i in range(n) for j in range(n)])
        Am = fmpq_mat(n, n, [A[j][i] for j in range(n) for i in range(n)])
        e = fmpq_mat(1, n, [eps(i) for i in range(1, n + 1)])
        one = fmpq_mat(n, 1, [1] * n)
        bidiag = all(A[j][i] == (1 if i == j else 2 if i == j - 1 else 0) for j in range(n) for i in range(n))
        q = Am.solve(one)
        qv = [q[i, 0] for i in range(n)]
        qf = [fmpq(1 + 2 * (-2) ** (i - 1), 3) for i in range(1, n + 1)]
        neg = [qv[i - 1] for i in range(2, n + 1, 2)]        # q2, q4, ...
        posq = [qv[i - 1] for i in range(1, n + 1, 2)]       # q1, q3, ...
        order = all(neg[k + 1] < neg[k] for k in range(len(neg) - 1)) and neg[0] < 0 < posq[0] and \
            all(posq[k] < posq[k + 1] for k in range(len(posq) - 1))
        prodmol = max(sum(b.values()) for _, b in R)
        okn = (Gm.rank() == n - 1 and all(x == 0 for x in (Gm * one).entries())
               and all(x == 0 for x in (e * Gm).entries()) and bidiag and qv == qf and order
               and max(sum(a.values()) for a, _ in R) == 3 and prodmol <= 5
               and all(v >= 0 for _, b in R for v in b.values()))
        if not okn:
            print(f"  failure at n = {n}")
        ok = ok and okn
    check(ok, "n = 2..80: rank n-1, Gamma 1 = 0, eps^T Gamma = 0, S bidiagonal unit, q formula and order, "
              "source molecularity 3, products nonnegative of molecularity <= 5")


def part_b():
    print("=== B. Groebner/Sturm: n = 2 and the displayed N_4 ===")
    from fractions import Fraction as Fr
    ig.run_case(2, [Fr(1), Fr(3)], "N_2 from Table 1 (authors' example, X <-> Y)")
    R4 = ig.table_network(4)
    print("  N_4 from the definitions / Table 1:",
          "; ".join(f"{ig.complex_str(a)} -> {ig.complex_str(b)}" for a, b in R4))
    disp = [({1: 1}, {4: 1}), ({1: 2, 2: 1}, {1: 3}), ({2: 2, 3: 1}, {2: 3, 3: 2}), ({3: 2, 4: 1}, {3: 1, 4: 2})]
    print("  network displayed as N_4 at line 130 of the original paper.tex:",
          "; ".join(f"{ig.complex_str(a)} -> {ig.complex_str(b)}" for a, b in disp))
    _, Gd = ig.matrices(4, disp)
    print(f"  its stoichiometric matrix has rank {Gd.rank()} and kernel {Gd.nullspace()}")
    check(R4 != disp and Gd.rank() == 4,
          "the originally displayed 'N_4' is NOT N_4: it has rank 4, i.e. it is a (4,4,4) network "
          "with no positive flux vector (display error, corrected in paper.tex)")
    R4fix = [({1: 1}, {1: 2, 4: 1}), ({1: 2, 2: 1}, {1: 1, 2: 2}), ({2: 2, 3: 1}, {2: 1}), ({3: 2, 4: 1}, {3: 3})]
    check(R4 == R4fix, "corrected display of N_4 (X1 -> 2X1+X4, 2X1+X2 -> X1+2X2, 2X2+X3 -> X2, "
                       "2X3+X4 -> 3X3) agrees with Table 1 and the solver's script")


def rpow(x, k):
    return x ** k if k >= 0 else 1 / x ** (-k)


def part_c():
    print("=== C. proof construction with non-integer lambdas, n = 2..12 ===")
    for n in range(2, 13):
        R = ig.table_network(n)
        A, G = mats(n, R)
        q = [(1 + 2 * (-2) ** (i - 1)) // 3 for i in range(1, n + 1)]
        lam = [fmpq(2 * r + 1, 5) for r in range(n)]           # 1/5, 3/5, 1, 7/5, ...
        expo = sorted(q + [0])
        M = fmpq_mat(n, n + 1, [rpow(l, e) for l in lam for e in expo])
        Rr, rk = M.rref()
        # rank n with pivots in the first n columns -> kernel spanned by (-Rr[:, n], 1)
        piv = all(Rr[i, i] == 1 for i in range(n))
        c = [-Rr[i, n] for i in range(n)] + [fmpq(1)]
        alt = rk == n and piv and all(x != 0 for x in c) and all(c[k] * c[k + 1] < 0 for k in range(n))
        pos = {e: k for k, e in enumerate(expo)}
        if c[pos[q[0]]] < 0:
            c = [-x for x in c]
        C = [c[pos[q[i]]] * eps(i + 1) for i in range(n)]
        T = -c[pos[0]]
        sign_ok = all(x > 0 for x in C) and T > 0
        kappa = []
        for j in range(n):
            v = fmpq(1)
            for i in range(n):
                v /= C[i] ** A[j][i]
            kappa.append(v)
        Gm = fmpq_mat(n, n, [G[i][j] for i in range(n) for j in range(n)])
        Am = fmpq_mat(n, n, [A[j][i] for j in range(n) for i in range(n)])
        allok = True
        dets = []
        for l in lam:
            x = [C[i] * rpow(l, q[i]) for i in range(n)]
            v = []
            for j in range(n):
                t = kappa[j]
                for i in range(n):
                    t *= x[i] ** A[j][i]
                v.append(t)
            eq = all(z == 0 for z in (Gm * fmpq_mat(n, 1, v)).entries())
            cl = sum(eps(i + 1) * x[i] for i in range(n)) == T
            Dv = fmpq_mat(n, n, [v[i] if i == k else 0 for i in range(n) for k in range(n)])
            Dx = fmpq_mat(n, n, [1 / x[i] if i == k else 0 for i in range(n) for k in range(n)])
            J = Gm * Dv * Am * Dx
            e1 = fmpq(0)
            for S in combinations(range(n), n - 1):
                sub = fmpq_mat(n - 1, n - 1, [J[a, b] for a in S for b in S])
                e1 += sub.det()
            dets.append(e1)
            allok = allok and eq and cl and all(xx > 0 for xx in x) and e1 != 0
        sturm = ""
        if n <= 7:
            lmb = sp.Symbol("l")
            P = sp.Poly(sum(sp.Rational(int(c[k].p), int(c[k].q)) * lmb ** (expo[k] - expo[0])
                            for k in range(n + 1)), lmb)
            npos = P.count_roots(0, None)
            sturm = f", Sturm: {npos} positive roots"
            allok = allok and npos == n and P.eval(0) != 0
        signs = "".join("+" if d > 0 else "-" for d in dets)
        check(alt and sign_ok and allok,
              f"n = {n}: alternating c, C>0, T>0, n positive equilibria on eps.x = T, all with "
              f"det(J|imGamma) != 0 (signs {signs}){sturm}")


if __name__ == "__main__":
    part_a()
    part_b()
    part_c()
    tot = FAILS + ig.FAILS
    print("ALL PASS" if tot == 0 else f"{tot} FAILURES")
    sys.exit(1 if tot else 0)
