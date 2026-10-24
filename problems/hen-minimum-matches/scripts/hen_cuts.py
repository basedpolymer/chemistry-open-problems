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


# ---------------------------------------------------------------------------------------------------------------
# Residual-capped bounds.  In every feasible flow, sum_i r_{i,t} = R_{t+1} (sum of the balance equations of
# intervals 0..t), hence 0 <= r_{i,t} <= min(R_{t+1}, H_i(<=t)) =: rc(i,t)   (rc(i,-1) = 0).
# ---------------------------------------------------------------------------------------------------------------

def rcap(inst, i, t):
    """upper bound on the residual of hot i leaving interval t (t = -1 .. k-1)"""
    if t < 0 or t >= inst.k - 1:
        return Fraction(0)
    return min(inst.R[t + 1], sum(inst.QH[i][: t + 1]))


def cap_window(inst, i, s, t):
    """max heat hot i can deliver to cold nodes in intervals s..t:  r_{i,s-1} + H_i[s..t]"""
    return rcap(inst, i, s - 1) + sum(inst.QH[i][s: t + 1])


def U_dp(inst, i, j, s=0, t=None):
    """Exact min cut (= max flow) of the 2-stream network restricted to window s..t with carry caps rc:
    source->hot_u (QH_iu), hot_u->hot_{u+1} (rc(i,u)), hot_u->cold_u (inf), cold_u->sink (QC_ju), plus an entry
    arc source->hot_s of capacity rc(i,s-1).  Any cut is an upper bound on heat i->j in the window (weak duality).
    DP over states: hot_u in S ('in', pay QC_ju) or not ('out', pay QH_iu); in->out transition pays rc(i,u)."""
    if t is None:
        t = inst.k - 1
    INF = None
    # entry: state before s is 'in' (source side); if hot_s is out we pay the entry arc rc(i,s-1)
    best_in = inst.QC[j][s]
    best_out = inst.QH[i][s] + rcap(inst, i, s - 1)
    for u in range(s + 1, t + 1):
        ni = min(best_in, best_out) + inst.QC[j][u]          # out->in is free (arc goes hot_u -> hot_{u+1})
        no = min(best_in + rcap(inst, i, u - 1), best_out) + inst.QH[i][u]
        best_in, best_out = ni, no
    return min(best_in, best_out)


def u_interval_capped(inst, i, j, t):
    return min(inst.QC[j][t], inst.QH[i][t] + rcap(inst, i, t - 1))


def window_rows(inst, cold=True, hot=True, only_pinch_starts=False):
    """(W1) cold j, window [s,T] (QC_j[T] > 0):  sum_i min(D_j[s..T], cap_i(s,T)) y_ij >= D_j[s..T]
       (W2) hot i, window [s,T] (QH_i[s] > 0):  sum_j min(H_i[s..T], D_j[s..T]) y_ij >= H_i[s..T] - rc(i,T)
    Proofs: (W1) j's demand in [s,T] is supplied by heat delivered in [s,T]; hot i delivers there at most
    r_{i,s-1} + H_i[s..T] and at most D_j[s..T] to j, and nothing if y_ij = 0.  (W2) heat of i originating in
    [s,T] is delivered in [s,T] except the part carried below T, which is <= r_{i,T} <= rc(i,T)."""
    n, m, k = inst.n, inst.m, inst.k
    rows = []
    starts = [s for s in range(k) if (not only_pinch_starts) or s == 0 or inst.R[s] == 0]
    if cold:
        for j in range(m):
            for s in starts:
                Dw = Fraction(0)
                for T in range(s, k):
                    Dw += inst.QC[j][T]
                    if inst.QC[j][T] == 0 or Dw == 0:
                        continue
                    coef = {}
                    for i in range(n):
                        c = min(Dw, cap_window(inst, i, s, T))
                        if c > 0:
                            coef[(i, j)] = c
                    rows.append((f"W1_c{j}_{s}_{T}", coef, Dw))
    if hot:
        for i in range(n):
            for s in starts:
                if inst.QH[i][s] == 0:
                    continue
                Hw = Fraction(0)
                for T in range(s, k):
                    Hw += inst.QH[i][T]
                    b = Hw - rcap(inst, i, T)
                    if b <= 0:
                        continue
                    coef = {}
                    for j in range(m):
                        c = min(Hw, sum(inst.QC[j][s: T + 1]))
                        if c > 0:
                            coef[(i, j)] = c
                    rows.append((f"W2_h{i}_{s}_{T}", coef, b))
    # remove dominated duplicates (same support and coefs scaled) cheaply: keep unique by (support, normalized)
    uniq = {}
    for name, coef, b in rows:
        key = tuple(sorted((e, c / b) for e, c in coef.items() if c < b)) + (tuple(sorted(e for e, c in coef.items() if c >= b)),)
        if key not in uniq:
            uniq[key] = (name, coef, b)
    return list(uniq.values())
