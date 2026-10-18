#!/usr/bin/env python3
"""Exact certificate for the family N_n answering Banaji & Feliu, Remark 4.23.

For each n in the requested range, this script builds the (n, n, n-1) network N_n of the note
(sources X1 and 2X_{i-1} + X_i, source molecularity 3) and checks in exact rational arithmetic:

  1. N_n has n species, n reactions, rank n-1; sources have molecularity <= 3 (and 3 is attained);
     all products are nonnegative.
  2. ker Gamma = span(1,...,1) (dynamically nontrivial) and ker Gamma^T = span(eps),
     eps = (+1,+1,-1,-1,+1,+1,...).
  3. The source matrix A is invertible and q = A^{-1} 1 has entries s_i = (1 + 2(-2)^{i-1})/3,
     pairwise distinct and nonzero.
  4. With lambda_r = r (r = 1..n), rate constants kappa and a total amount T > 0 constructed as in the
     note, the n points x^(r) = C o r^q are positive equilibria on the class eps.x = T
     (Gamma v(x) = 0 exactly, eps.x = T exactly).
  5. Each x^(r) is nondegenerate: the Jacobian restricted to im Gamma has nonzero determinant.
  6. The conservation polynomial F(lambda) has exactly n sign changes, so by Descartes' rule there is
     no other positive equilibrium on that class.

Usage: python certify_family.py [nmin] [nmax]     (default 3 10; runtime a few seconds)
Environment: Python 3.12, python-flint 0.9 (exact fmpq / fmpq_mat arithmetic).
Prints PASS/FAIL per check and exits with status 1 on any failure.
"""
import sys
from flint import fmpq, fmpq_mat

FAILS = 0


def check(cond, msg):
    global FAILS
    print(("PASS " if cond else "FAIL ") + msg)
    if not cond:
        FAILS += 1


def eps(i):
    """Sign of species X_i (1-based) in the conservation law: + + - - + + - - ..."""
    return 1 if ((i - 1) // 2) % 2 == 0 else -1


def network(n):
    """Return (A, B) with A[j] = source of reaction R_{j+1}, B[j] = product (lists of ints)."""
    s = -eps(n)
    A, B = [], []
    # R1: X1 -> X1 + r, r = s e_1 + e_n
    a = [0] * n; a[0] = 1
    b = a[:]; b[0] += s; b[n - 1] += 1
    A.append(a); B.append(b)
    for j in range(2, n + 1):           # R_j has source 2X_{j-1} + X_j
        a = [0] * n; a[j - 2] = 2; a[j - 1] = 1
        b = a[:]
        if j % 2 == 0:                  # j = 2k: gamma = c_k (e_{2k} - e_{2k-1}), c_k = (-1)^{k+1} s
            k = j // 2
            c = (-1) ** (k + 1) * s
            b[j - 1] += c; b[j - 2] -= c
        else:                           # j = 2k+1: gamma = d_k (e_{2k} + e_{2k+1}), d_k = (-1)^k s
            k = (j - 1) // 2
            d = (-1) ** k * s
            b[j - 2] += d; b[j - 1] += d
        A.append(a); B.append(b)
    return A, B


def fmt(v):
    t = [(f"{c}" if c != 1 else "") + f"X{i+1}" for i, c in enumerate(v) if c]
    return " + ".join(t) if t else "0"


def nullspace(M):
    """Basis of the right kernel of an fmpq_mat, as list of column lists."""
    R, rank = M.rref()
    rows, cols = M.nrows(), M.ncols()
    pivots = []
    r = 0
    for c in range(cols):
        if r < rank and R[r, c] != 0:
            pivots.append(c); r += 1
    free = [c for c in range(cols) if c not in pivots]
    basis = []
    for f in free:
        v = [fmpq(0)] * cols
        v[f] = fmpq(1)
        for i, p in enumerate(pivots):
            v[p] = -R[i, f]
        basis.append(v)
    return basis


def rpow(x, e):
    """x**e for fmpq x > 0 and integer e (possibly negative)."""
    return x ** e if e >= 0 else fmpq(1) / (x ** (-e))


def certify(n):
    print(f"--- n = {n} ---")
    A, B = network(n)
    for j in range(n):
        print(f"  R{j+1}: {fmt(A[j])} -> {fmt(B[j])}")
    check(all(x >= 0 for row in B for x in row), "products are nonnegative")
    check(max(sum(a) for a in A) == 3, "maximum source molecularity is 3")
    check(len({tuple(a) for a in A}) == n, "n distinct sources")
    G = fmpq_mat(n, n, [B[j][i] - A[j][i] for i in range(n) for j in range(n)])  # columns = reactions
    check(G.rank() == n - 1, "rank Gamma = n - 1")
    ones = fmpq_mat(n, 1, [1] * n)
    check(all(x == 0 for x in (G * ones).entries()), "Gamma (1,...,1)^T = 0 (positive flux vector)")
    E = fmpq_mat(1, n, [eps(i) for i in range(1, n + 1)])
    check(all(x == 0 for x in (E * G).entries()), "eps^T Gamma = 0 with eps = (+,+,-,-,...)")
    Am = fmpq_mat(n, n, [A[j][i] for j in range(n) for i in range(n)])  # rows = sources
    check(Am.det() != 0, "source matrix A invertible")
    q = Am.solve(ones)
    qv = [q[i, 0] for i in range(n)]
    expected = [fmpq(1 + 2 * (-2) ** (i - 1), 3) for i in range(1, n + 1)]
    check(qv == expected, "q = A^{-1} 1 has entries (1 + 2(-2)^{i-1})/3")
    check(len(set(qv)) == n and all(x != 0 for x in qv), "exponents distinct and nonzero")
    qi = [int(x.p) for x in qv]
    # terms of F(lambda) = sum_i eps_i C_i lambda^{s_i} - T, ordered by exponent; constant has exponent 0
    expo = sorted(qi + [0])
    lam = [fmpq(r) for r in range(1, n + 1)]
    M = fmpq_mat(n, n + 1, [rpow(l, e) for l in lam for e in expo])
    ker = nullspace(M)
    check(len(ker) == 1, "generalized Vandermonde kernel is one-dimensional")
    c = ker[0]
    check(all(x != 0 for x in c) and all(c[k] * c[k + 1] < 0 for k in range(n)),
          "kernel vector strictly alternates in sign")
    pos = {e: k for k, e in enumerate(expo)}
    sgn = 1 if (c[pos[qi[0]]] > 0) == (eps(1) > 0) else -1    # global sign so that coefficient of X1 has sign eps_1
    c = [sgn * x for x in c]
    check(all((c[pos[qi[i]]] > 0) == (eps(i + 1) > 0) for i in range(n)),
          "coefficient signs match eps_i for every species")
    C = [c[pos[qi[i]]] * eps(i + 1) for i in range(n)]      # C_i > 0
    T = -c[pos[0]]
    check(all(x > 0 for x in C) and T > 0, "C_i > 0 and T > 0")
    # rate constants kappa_j = C^{-a_j}  (so that kappa_j x^{a_j} = lambda on the curve x = C o lambda^q)
    kappa = []
    for j in range(n):
        val = fmpq(1)
        for i in range(n):
            if A[j][i]:
                val /= C[i] ** A[j][i]
        kappa.append(val)
    check(all(k > 0 for k in kappa), "kappa_j > 0")
    # equilibria
    for r, l in enumerate(lam, start=1):
        x = [C[i] * rpow(l, qi[i]) for i in range(n)]
        v = []
        for j in range(n):
            val = kappa[j]
            for i in range(n):
                if A[j][i]:
                    val *= x[i] ** A[j][i]
            v.append(val)
        Gv = G * fmpq_mat(n, 1, v)
        ok_eq = all(z == 0 for z in Gv.entries())
        ok_T = sum(eps(i + 1) * x[i] for i in range(n)) == T
        # Jacobian J = Gamma diag(v) A diag(1/x); restriction to im Gamma via basis of im Gamma = eps^perp
        Dv = fmpq_mat(n, n, [v[i] if i == k else 0 for i in range(n) for k in range(n)])
        Dx = fmpq_mat(n, n, [1 / x[i] if i == k else 0 for i in range(n) for k in range(n)])
        J = G * Dv * Am * Dx
        Bas = fmpq_mat(n, n - 1, [G[i, k] for i in range(n) for k in range(n - 1)])  # first n-1 columns span im Gamma
        check(Bas.rank() == n - 1, f"  r={r}: first n-1 reaction vectors span im Gamma")
        BtB = Bas.transpose() * Bas
        Mres = BtB.solve(Bas.transpose() * J * Bas)
        ok_inv = (Bas * Mres == J * Bas)            # J maps im Gamma into im Gamma; Mres is its matrix
        det = Mres.det()
        check(ok_eq and ok_T and ok_inv and det != 0,
              f"  r={r}: x^(r) is a positive equilibrium on eps.x = T, restricted det = {float(det):.3e} != 0")
    # Descartes: F has exactly n sign changes (n+1 nonzero alternating coefficients) -> at most n positive roots
    changes = sum(1 for k in range(n) if c[k] * c[k + 1] < 0)
    check(changes == n, f"F has {changes} sign changes: no other positive equilibrium on this class")
    print(f"  T = {T}" if len(str(T)) < 80 else f"  T = <rational with {len(str(T))} characters>")


if __name__ == "__main__":
    nmin = int(sys.argv[1]) if len(sys.argv) > 1 else 3
    nmax = int(sys.argv[2]) if len(sys.argv) > 2 else 10
    for n in range(nmin, nmax + 1):
        certify(n)
    print("ALL PASS" if FAILS == 0 else f"{FAILS} FAILURES")
    sys.exit(1 if FAILS else 0)
