#!/usr/bin/env python3
"""Novelty check: exact verification (same pipeline as direct_check.py: msolve + Singular + PARI) of the explicit
networks of the earlier, unrefereed public claim

  Anonymous (maintained by I. Pitchford, Evidence Press), "Sharp positive-equilibrium counts in quadratic
  mass-action networks with two more reactions than species", v1.0.0-candidate, 4 Oct 2026,
  doi:10.5281/zenodo.23146115, github.com/ipitchford/quadratic-mass-action-sharpness (README as of 18 Oct 2026),

and of the Phillipson-Rojas system G_eps (Contemp. Math. 605 (2013), Thm 1.6, eps = 1/4) read as the equilibrium
system of the network N_n(alpha, beta) of the manuscript (ratio R = 16 instead of 64).
Only the reaction lists and rate constants are taken from those sources; no code from them is used.
Runtime: a few seconds.  Usage: python prior_art_check.py   (exit code 1 on failure)
"""
import sys
from fractions import Fraction as F

import direct_check as dc
from direct_check import certify, chain_network

# (P1) compact (3,5,3) example of the earlier claim; species X, Y, Z = X1, X2, X3
#   X -> 0 (1), Y -> 6X + 4Z (1), 2X -> 6Y + 4Z (1), 2Z -> 4Y (2/3137), Y + Z -> 0 (1)
P1 = [({0: 1}, {}), ({1: 1}, {0: 6, 2: 4}), ({0: 2}, {1: 6, 2: 4}), ({2: 2}, {1: 4}), ({1: 1, 2: 1}, {})]
certify("earlier claim, compact (3,5,3) example (product molecularity 10)", 3, P1, [1, 1, 1, F(2, 3137), 1], 4)


# (P2) all-n family of the earlier claim: 0 -> sum a_i X_i with a_i = 24*25^(n-i) (i<n), a_n = 25;
#   X1 -> 0; X_i + X_{i+1} -> 2X_i (i<n); 2X_n -> 3X_n; all rates 1 except the last, 5^((-1)^n (2n+1)).
def pitchford(n):
    a = {i: 24 * 25 ** (n - 1 - i) for i in range(n - 1)}
    a[n - 1] = 25
    R = [({}, a), ({0: 1}, {})]
    R += [({i: 1, i + 1: 1}, {i: 2}) for i in range(n - 1)]
    R.append(({n - 1: 2}, {n - 1: 3}))
    e = (-1) ** n * (2 * n + 1)
    k = [1] * (len(R) - 1) + [F(5) ** e]
    return R, k


for n in range(2, 7):
    R, k = pitchford(n)
    certify(f"earlier claim, all-n family, n = {n}", n, R, k, n + 1)

# (P3) Phillipson-Rojas G_{1/4}: x1 x2 = eps + x1^2, x_i x_{i+1} = 1 + eps^(2i-3) x1^2, x_n = 1 + eps^(2n-3) x1^2.
#   Realised by N_n(alpha, beta) with alpha_{i+1} = 16^(i-1), beta_{i+1} = 1 and rate constants
#   k_0 = 1/4, k_1 = 1, k_2 = 1, k_{i+1} = 4^(2i-3) (i = 2..n); dividing equation i by k_{i+1} gives G_{1/4}.
for n in range(2, 7):
    k = [F(1, 4), F(1), F(1)] + [F(4) ** (2 * i - 3) for i in range(2, n + 1)]
    certify(f"Phillipson-Rojas G_1/4 realised as N_n(alpha,beta), R = 16, n = {n}", n,
            chain_network(n, [16 ** (i - 1) for i in range(1, n + 1)], [1] * n), k, n + 1)

for n in range(2, 7):
    k = [F(1, 512), F(1), F(1)] + [F(8) ** (2 * i - 3) for i in range(2, n + 1)]
    certify(f"Phillipson-Rojas G_1/8 realised as N_n(alpha,beta), R = 64, n = {n}", n,
            chain_network(n, [64 ** i for i in range(1, n + 1)], [1] * n), k, n + 1)

# (P4) the R = 64 networks of the manuscript, with the rate constants of its proof (kappa_1 = 2^(-+(6n+3)), others 1),
#   have the same Gale constant c = sum_j w_j log kappa_j as N_n(R=64) with kappa_0 = 2^-9, kappa_1 = kappa_2 = 1,
#   kappa_{i+1} = 8^(2i-3) (i >= 2), whose equilibrium equations (divided by kappa_{i+1}) are Phillipson-Rojas G_{1/8}.
ok = True
for n in range(2, 101):
    w = [None, 1] + [2 * (-1) ** (j - 1) for j in range(2, n + 2)]
    w[0] = -sum(w[1:])
    log2k = [-9, 0, 0] + [3 * (2 * i - 3) for i in range(2, n + 1)]      # log_2 of the G_{1/8} rate constants
    ok &= sum(a * b for a, b in zip(w, log2k)) == (-(6 * n + 3) if n % 2 else 6 * n + 3)
dc.check(ok, "n = 2..100: Gale constant of the G_{1/8} realisation equals that of the proof's rate constants")

print("ALL PASS" if dc.FAILS == 0 else f"{dc.FAILS} FAILURES")
sys.exit(1 if dc.FAILS else 0)
