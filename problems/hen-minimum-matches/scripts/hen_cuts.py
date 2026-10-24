"""Valid inequalities in the match variables y (proofs in result.md), generated in exact arithmetic.

Notation: H_i(<=t) = heat of hot i in intervals <= t, H_i(>=t) likewise; D_j(<=t), D_j(>=t) for cold j.

(K1) cold prefix knapsack, for cold j and interval T:
        sum_i min(H_i(<=T), D_j(<=T)) y_ij >= D_j(<=T)
     (the demand of j in intervals <= T is met by heat that originates in intervals <= T).
(K2) hot suffix knapsack, for hot i and interval T:
        sum_j min(D_j(>=T), H_i(>=T)) y_ij >= H_i(>=T)
     (heat of i in intervals >= T can only be received in intervals >= T).
(K3) aggregated cold prefix ("all cold streams in J"):  for a set J of cold streams and interval T,
        sum_{i} sum_{j in J} min(H_i(<=T), D_j(<=T)) y_ij >= D_J(<=T)  (implied by K1 summed; not used)
(P)  partner count: from any knapsack  sum_e a_e y_e >= b  with 0 <= a_e, y binary:
        sum_e y_e >= p,   p = least number of coefficients whose sum is >= b (largest first).
     More generally the Chvatal-Gomory rounding  sum_e ceil(min(a_e,b)/lam) y_e >= ceil(b/lam), lam > 0.
Each generated row is returned as (name, {pair: coef}, rhs) meaning sum coef*y >= rhs, with Fraction data.
"""
from fractions import Fraction
import math


def knapsack_rows(inst, kinds=("K1", "K2")):
    n, m, k = inst.n, inst.m, inst.k
    Hc = [[sum(inst.QH[i][: t + 1]) for t in range(k)] for i in range(n)]
    Hs = [[sum(inst.QH[i][t:]) for t in range(k)] for i in range(n)]
    Dc = [[sum(inst.QC[j][: t + 1]) for t in range(k)] for j in range(m)]
    Ds = [[sum(inst.QC[j][t:]) for t in range(k)] for j in range(m)]
    rows = []
    if "K1" in kinds:
        for j in range(m):
            for T in range(k):
                if inst.QC[j][T] == 0:
                    continue  # only prefixes ending where demand occurs (others are dominated)
                b = Dc[j][T]
                coef = {(i, j): min(Hc[i][T], b) for i in range(n) if Hc[i][T] > 0}
                rows.append((f"K1_c{j}_t{T}", coef, b))
    if "K2" in kinds:
        for i in range(n):
            for T in range(k):
                if inst.QH[i][T] == 0:
                    continue
                b = Hs[i][T]
                coef = {(i, j): min(Ds[j][T], b) for j in range(m) if Ds[j][T] > 0}
                rows.append((f"K2_h{i}_t{T}", coef, b))
    return rows


def partner_count(coef, b):
    vals = sorted(coef.values(), reverse=True)
    s = Fraction(0)
    for p, v in enumerate(vals, 1):
        s += v
        if s >= b:
            return p
    return None  # knapsack infeasible (cannot happen for a feasible instance)


def count_rows(rows):
    out = []
    for name, coef, b in rows:
        p = partner_count(coef, b)
        assert p is not None
        if p >= 2:
            out.append(("P_" + name, {e: Fraction(1) for e in coef}, Fraction(p)))
    return out


def cg_rows(rows, lams_per_row=3):
    """CG roundings sum ceil(min(a,b)/lam) y >= ceil(b/lam) for lam = b/q, q = 2..lams_per_row+1 (exact)."""
    out = []
    for name, coef, b in rows:
        for q in range(2, lams_per_row + 2):
            lam = b / q
            c2 = {e: Fraction(math.ceil(min(a, b) / lam)) for e, a in coef.items()}
            out.append((f"CG{q}_" + name, c2, Fraction(math.ceil(b / lam))))
    return out


def check_row_valid(inst, row):
    """Re-derive a K1/K2/P/CG row from the instance and confirm it is identical (used by the checker)."""
    raise NotImplementedError
