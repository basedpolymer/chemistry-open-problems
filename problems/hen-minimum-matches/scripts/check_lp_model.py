#!/usr/bin/env python3
"""Independent re-derivation of every row and bound of an LP file written by hen_lpfile.write_lp_strong.

usage: python3 check_lp_model.py <instance.dat> <model.lp> <K or ->      (K: the 'maxm' row sum y <= K)

Every LP row must be one of the following families (proved valid in result.md, Sec. 3); this script recomputes
the coefficients from the instance data with code independent of the generator and compares exactly:
  hot_i_t   sum_j q_ijt + r_it - r_i(t-1) = QH_it                              (balance)
  cold_j_t  sum_i q_ijt = QC_jt                                                (balance)
  link_i_j  sum_t q_ijt - U_ij y_ij <= 0, U_ij >= max heat i->j (checked: U_ij = value of an explicit cut)
  dl_i_j_t  q_ijt - u_ijt y_ij <= 0, u_ijt = min(QC_jt, QH_it + rc_i(t-1))
  W1_c{j}_{s}_{T}   sum_i min(D_j[s..T], rc_i(s-1) + H_i[s..T]) y_ij >= D_j[s..T]
  W2_h{i}_{s}_{T}   sum_j min(H_i[s..T], D_j[s..T]) y_ij >= H_i[s..T] - rc_i(T)
  P_<W-row>         sum y over the W-row's support >= p, p <= least number of W-coefficients reaching its rhs
  maxm              sum_ij y_ij <= K
with rc_i(t) = min(R_{t+1}, H_i[0..t]) for 0 <= t <= k-2 and rc_i(-1) = rc_i(k-1) = 0, where R_t is the exact
heat cascade computed from the data.  Bounds: 0 <= q_ijt <= u_ijt, 0 <= r_it <= rc_i(t), y binary.
Variables absent from the file must be forced to 0 in every feasible flow (u_ijt = 0, rc_i(t) = 0, or U_ij = 0):
this is checked too.  Prints PASS/FAIL.  Environment: Python >= 3.10, standard library.  Runtime: < 1 min.
"""
import sys, os, re
from fractions import Fraction
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from hen_instance import parse_dat
from check_vipr_model import parse_lp


def main():
    inst = parse_dat(sys.argv[1])
    rows, bounds, binaries, obj = parse_lp(sys.argv[2])
    K = None if sys.argv[3] == "-" else Fraction(sys.argv[3])
    n, m, k = inst.n, inst.m, inst.k
    head = open(sys.argv[2]).readline()
    p = int(re.search(r"10\^-(\d+)", head).group(1))
    S = Fraction(10) ** p
    QH = [[x / S for x in row] for row in inst.QH]
    QC = [[x / S for x in row] for row in inst.QC]
    # exact cascade from the data (not from the file's R lines)
    R = [Fraction(0)]
    for t in range(k):
        R.append(R[-1] + sum(QH[i][t] for i in range(n)) - sum(QC[j][t] for j in range(m)))
    assert R[-1] == 0 and min(R) >= 0
    def H(i, a, b):
        return sum(QH[i][a:b + 1], Fraction(0))
    def Dm(j, a, b):
        return sum(QC[j][a:b + 1], Fraction(0))
    def rc(i, t):
        if t < 0 or t >= k - 1:
            return Fraction(0)
        return min(R[t + 1], H(i, 0, t))
    def ucap(i, j, t):
        return min(QC[j][t], QH[i][t] + rc(i, t - 1))
    def cutU(i, j):
        # explicit cut of the 2-stream carry-capped network: enumerate by DP and reconstruct the cut set,
        # then evaluate its capacity directly (weak duality: capacity >= any flow i->j)
        INF = None
        cost = {}
        # DP with back-pointers
        best = {("in", 0): (QC[j][0], None), ("out", 0): (QH[i][0], None)}
        for u in range(1, k):
            ci = min((best[("in", u - 1)][0], "in"), (best[("out", u - 1)][0], "out"))
            co = min((best[("in", u - 1)][0] + rc(i, u - 1), "in"), (best[("out", u - 1)][0], "out"))
            best[("in", u)] = (ci[0] + QC[j][u], ci[1])
            best[("out", u)] = (co[0] + QH[i][u], co[1])
        st = min(("in", "out"), key=lambda s_: best[(s_, k - 1)][0])
        val = best[(st, k - 1)][0]
        states = [None] * k
        for u in range(k - 1, -1, -1):
            states[u] = st
            st = best[(st, u)][1]
        # direct capacity of the cut S = {source} + {hot_u, cold_u : state in}
        capv = Fraction(0)
        for u in range(k):
            if states[u] == "in":
                capv += QC[j][u]          # cold_u -> sink
                if u + 1 < k and states[u + 1] == "out":
                    capv += rc(i, u)      # hot_u -> hot_{u+1}
            else:
                capv += QH[i][u]          # source -> hot_u
        assert capv == val
        return val
    errs = []
    yv = {(int(a), int(b)) for a, b in (re.match(r"y_(\d+)_(\d+)$", v).groups() for v in binaries)}
    qv = set(); rv = set()
    for v in bounds:
        mo = re.match(r"q_(\d+)_(\d+)_(\d+)$", v)
        if mo:
            i, j, t = map(int, mo.groups()); qv.add((i, j, t))
            if bounds[v] != (0, ucap(i, j, t)):
                errs.append(f"bound {v}")
        mo = re.match(r"r_(\d+)_(\d+)$", v)
        if mo:
            i, t = map(int, mo.groups()); rv.add((i, t))
            if bounds[v] != (0, rc(i, t)):
                errs.append(f"bound {v}")
    U = {}
    for i in range(n):
        for j in range(m):
            U[i, j] = cutU(i, j)
            if U[i, j] > 0 and (i, j) not in yv:
                errs.append(f"missing y_{i}_{j}")
            if U[i, j] == 0 and (i, j) in yv:
                pass  # harmless
            for t in range(k):
                if ucap(i, j, t) > 0 and U[i, j] > 0 and (i, j, t) not in qv:
                    errs.append(f"missing q_{i}_{j}_{t}")
        for t in range(k - 1):
            if rc(i, t) > 0 and (i, t) not in rv:
                errs.append(f"missing r_{i}_{t}")
    if obj and set(obj) != {f"y_{i}_{j}" for (i, j) in yv} and not (len(obj) == 0):
        if any(c != 0 for c in obj.values()):
            if any(c != 1 for c in obj.values()) or set(obj) != {f"y_{i}_{j}" for (i, j) in yv}:
                errs.append("objective is not sum of all y")
    W = {}
    fam = {}
    for name, coef, sense, rhs in rows:
        f = name.split("_")[0]
        fam[f] = fam.get(f, 0) + 1
        if f == "hot":
            _, i, t = name.split("_"); i, t = int(i), int(t)
            exp = {f"q_{i}_{j}_{t}": 1 for j in range(m) if (i, j, t) in qv}
            if (i, t) in rv: exp[f"r_{i}_{t}"] = 1
            if (i, t - 1) in rv: exp[f"r_{i}_{t-1}"] = -1
            ok = sense == "E" and coef == exp and rhs == QH[i][t]
        elif f == "cold":
            _, j, t = name.split("_"); j, t = int(j), int(t)
            exp = {f"q_{i}_{j}_{t}": 1 for i in range(n) if (i, j, t) in qv}
            ok = sense == "E" and coef == exp and rhs == QC[j][t]
        elif f == "link":
            _, i, j = name.split("_"); i, j = int(i), int(j)
            exp = {f"q_{i}_{j}_{t}": 1 for t in range(k) if (i, j, t) in qv}
            exp[f"y_{i}_{j}"] = -U[i, j]
            ok = sense == "L" and rhs == 0 and coef == exp
        elif f == "dl":
            _, i, j, t = name.split("_"); i, j, t = int(i), int(j), int(t)
            ok = sense == "L" and rhs == 0 and coef == {f"q_{i}_{j}_{t}": 1, f"y_{i}_{j}": -ucap(i, j, t)}
        elif f == "W1":
            mo = re.match(r"W1_c(\d+)_(\d+)_(\d+)$", name); j, s, T = map(int, mo.groups())
            b = Dm(j, s, T)
            exp = {}
            for i in range(n):
                c = min(b, rc(i, s - 1) + H(i, s, T))
                if c > 0 and (i, j) in yv: exp[f"y_{i}_{j}"] = c
            ok = sense == "G" and rhs == b and coef == exp
            W[name] = (coef, rhs)
        elif f == "W2":
            mo = re.match(r"W2_h(\d+)_(\d+)_(\d+)$", name); i, s, T = map(int, mo.groups())
            Hw = H(i, s, T); b = Hw - rc(i, T)
            exp = {}
            for j in range(m):
                c = min(Hw, Dm(j, s, T))
                if c > 0 and (i, j) in yv: exp[f"y_{i}_{j}"] = c
            ok = sense == "G" and rhs == b and b > 0 and coef == exp
            W[name] = (coef, rhs)
        elif f == "P":
            base = name[2:]
            ok = base in W
            if ok:
                wc, wb = W[base]
                vals = sorted(wc.values(), reverse=True)
                acc, need = Fraction(0), None
                for cnt, v in enumerate(vals, 1):
                    acc += v
                    if acc >= wb:
                        need = cnt; break
                ok = sense == "G" and need is not None and rhs <= need and set(coef) == set(wc) and all(c == 1 for c in coef.values())
        elif f == "maxm":
            ok = K is not None and sense == "L" and rhs == K and set(coef) == {f"y_{i}_{j}" for (i, j) in yv} \
                and all(c == 1 for c in coef.values())
        else:
            ok = False
        if not ok:
            errs.append(f"row {name} not re-derived")
    print(f"{inst.name}: {len(rows)} rows by family {fam}; {len(qv)} q, {len(rv)} r, {len(yv)} y variables")
    if errs:
        print("FAIL:", "; ".join(errs[:10]), f"({len(errs)} errors)")
        sys.exit(1)
    print("PASS: every row/bound of the LP file re-derived exactly from the instance data")


if __name__ == "__main__":
    main()
