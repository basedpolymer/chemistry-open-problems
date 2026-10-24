"""MILP formulations of the minimum number of matches problem (independent re-implementation).

Transshipment MILP (Papoulias-Grossmann 1983 type, as in Letsios et al. 2018 Sec. 2):
    min  sum_{ij} y_ij
    hot  i, interval t :  r_{i,t-1} + QH[i,t] = sum_j q_{ijt} + r_{it},   r_{i,-1} = r_{i,k-1} = 0
    cold j, interval t :  sum_i q_{ijt} = QC[j,t]
    link              :  sum_t q_{ijt} <= U_ij y_ij
    optional (valid)  :  q_{ijt} <= u_ijt y_ij,  u_ijt = min(QC[j,t], sum_{s<=t} QH[i,s])          ("disagg")
    optional (valid)  :  sum_{s<=t} q_{ijs} <= min(QC_j(<=t), QH_i(<=t)) y_ij                       ("cumul")
    y binary, q, r >= 0.

Transportation MILP (Cerda-Westerberg 1983 type):
    q_{i s j t} >= 0 for s <= t only ("reduced" transportation model),
    sum_{j,t>=s} q_{isjt} = QH[i,s],  sum_{i,s<=t} q_{isjt} = QC[j,t],  sum_{s,t} q_{isjt} <= U_ij y_ij.

Floating point models are for search only; certificates are checked in exact arithmetic elsewhere.
Environment: Python 3.13, pyscipopt 6.3 (SCIP 10.0), highspy 1.15.
"""
from fractions import Fraction


def model_data(inst, scale=None):
    """float data, optionally scaled (heat / scale)."""
    if scale is None:
        scale = 1.0
    QH = [[float(v) / scale for v in row] for row in inst.QH]
    QC = [[float(v) / scale for v in row] for row in inst.QC]
    return QH, QC


def build_transshipment_scip(inst, disagg=True, cumul=False, scale=None, maxmatches=None, name=None,
                             fix=None, extra_cuts=None, Umode="pair"):
    import pyscipopt as ps
    n, m, k = inst.n, inst.m, inst.k
    if scale is None:
        scale = float(max(max(max(r) for r in inst.QH), max(max(r) for r in inst.QC))) / 1000.0
    QH, QC = model_data(inst, scale)
    M = ps.Model(name or inst.name)
    y, q, r = {}, {}, {}
    pairs = [(i, j) for i in range(n) for j in range(m) if inst.U_pair(i, j) > 0]
    for (i, j) in pairs:
        y[i, j] = M.addVar(vtype="B", name=f"y_{i}_{j}", obj=1.0)
    hc = [[inst.hot_cum(i, t) for t in range(k)] for i in range(n)]
    for (i, j) in pairs:
        for t in range(k):
            if inst.QC[j][t] > 0 and hc[i][t] > 0:
                q[i, j, t] = M.addVar(lb=0.0, ub=float(inst.u_interval(i, j, t)) / scale, name=f"q_{i}_{j}_{t}")
    for i in range(n):
        for t in range(k - 1):
            if hc[i][t] > 0:
                r[i, t] = M.addVar(lb=0.0, ub=float(hc[i][t]) / scale, name=f"r_{i}_{t}")
    for i in range(n):
        for t in range(k):
            terms = [(q[i, j, t], 1.0) for j in range(m) if (i, j, t) in q]
            if (i, t) in r:
                terms.append((r[i, t], 1.0))
            if (i, t - 1) in r:
                terms.append((r[i, t - 1], -1.0))
            rhs = QH[i][t]
            if not terms:
                assert rhs == 0
                continue
            M.addCons(ps.quicksum(c * v for v, c in terms) == rhs, name=f"hot_{i}_{t}")
    for j in range(m):
        for t in range(k):
            if inst.QC[j][t] > 0:
                M.addCons(ps.quicksum(q[i, j, t] for i in range(n) if (i, j, t) in q) == QC[j][t], name=f"cold_{j}_{t}")
    for (i, j) in pairs:
        U = float(inst.U_pair(i, j)) / scale
        M.addCons(ps.quicksum(q[i, j, t] for t in range(k) if (i, j, t) in q) - U * y[i, j] <= 0, name=f"link_{i}_{j}")
        if disagg:
            for t in range(k):
                if (i, j, t) in q:
                    M.addCons(q[i, j, t] - float(inst.u_interval(i, j, t)) / scale * y[i, j] <= 0, name=f"dlink_{i}_{j}_{t}")
        if cumul:
            for t in range(k):
                if (i, j, t) in q:
                    cap = min(inst.cold_cum(j, t), hc[i][t])
                    M.addCons(ps.quicksum(q[i, j, s] for s in range(t + 1) if (i, j, s) in q) - float(cap) / scale * y[i, j] <= 0,
                              name=f"clink_{i}_{j}_{t}")
    if maxmatches is not None:
        M.addCons(ps.quicksum(y.values()) <= maxmatches, name="maxmatches")
    if fix:
        for (i, j), v in fix.items():
            if (i, j) in y:
                M.fixVar(y[i, j], v)
    if extra_cuts:
        for cut in extra_cuts:  # cut = (list of ((i,j),coef), sense, rhs)
            expr = ps.quicksum(c * y[p] for p, c in cut[0] if p in y)
            if cut[1] == ">=":
                M.addCons(expr >= cut[2])
            else:
                M.addCons(expr <= cut[2])
    M.setMinimize()
    return M, y, q, r, scale


def build_transportation_scip(inst, scale=None, disagg=False):
    import pyscipopt as ps
    n, m, k = inst.n, inst.m, inst.k
    if scale is None:
        scale = float(max(max(max(r) for r in inst.QH), max(max(r) for r in inst.QC))) / 1000.0
    QH, QC = model_data(inst, scale)
    M = ps.Model(inst.name + "_transp")
    y, q = {}, {}
    pairs = [(i, j) for i in range(n) for j in range(m) if inst.U_pair(i, j) > 0]
    for (i, j) in pairs:
        y[i, j] = M.addVar(vtype="B", name=f"y_{i}_{j}", obj=1.0)
        for s in range(k):
            if inst.QH[i][s] == 0:
                continue
            for t in range(s, k):
                if inst.QC[j][t] > 0:
                    q[i, s, j, t] = M.addVar(lb=0.0, ub=float(min(inst.QH[i][s], inst.QC[j][t])) / scale,
                                             name=f"q_{i}_{s}_{j}_{t}")
    for i in range(n):
        for s in range(k):
            if inst.QH[i][s] > 0:
                M.addCons(ps.quicksum(q[i, s, j, t] for j in range(m) for t in range(s, k) if (i, s, j, t) in q) == QH[i][s])
    for j in range(m):
        for t in range(k):
            if inst.QC[j][t] > 0:
                M.addCons(ps.quicksum(q[i, s, j, t] for i in range(n) for s in range(t + 1) if (i, s, j, t) in q) == QC[j][t])
    for (i, j) in pairs:
        U = float(inst.U_pair(i, j)) / scale
        M.addCons(ps.quicksum(q[i, s, j, t] for s in range(k) for t in range(s, k) if (i, s, j, t) in q) - U * y[i, j] <= 0)
        if disagg:
            for t in range(k):
                lst = [q[i, s, j, t] for s in range(t + 1) if (i, s, j, t) in q]
                if lst:
                    M.addCons(ps.quicksum(lst) - float(inst.u_interval(i, j, t)) / scale * y[i, j] <= 0)
    M.setMinimize()
    return M, y, q, scale


def build_transshipment_highs(inst, disagg=True, cumul=False, scale=None, maxmatches=None):
    """Same transshipment model in HiGHS (highspy), built row by row."""
    import highspy
    import numpy as np
    n, m, k = inst.n, inst.m, inst.k
    if scale is None:
        scale = float(max(max(max(r) for r in inst.QH), max(max(r) for r in inst.QC))) / 1000.0
    QH, QC = model_data(inst, scale)
    h = highspy.Highs()
    inf = highspy.kHighsInf
    col = {}
    lb, ub, cost, integ = [], [], [], []
    def addcol(key, l, u, c, isint):
        col[key] = len(lb)
        lb.append(l); ub.append(u); cost.append(c); integ.append(isint)
    pairs = [(i, j) for i in range(n) for j in range(m) if inst.U_pair(i, j) > 0]
    hc = [[inst.hot_cum(i, t) for t in range(k)] for i in range(n)]
    for (i, j) in pairs:
        addcol(("y", i, j), 0.0, 1.0, 1.0, True)
    for (i, j) in pairs:
        for t in range(k):
            if inst.QC[j][t] > 0 and hc[i][t] > 0:
                addcol(("q", i, j, t), 0.0, float(inst.u_interval(i, j, t)) / scale, 0.0, False)
    for i in range(n):
        for t in range(k - 1):
            if hc[i][t] > 0:
                addcol(("r", i, t), 0.0, float(hc[i][t]) / scale, 0.0, False)
    ncol = len(lb)
    h.setOptionValue("output_flag", True)
    h.addVars(ncol, np.array(lb), np.array(ub))
    h.changeColsCost(ncol, np.arange(ncol, dtype=np.int32), np.array(cost))
    intidx = np.array([c for c in range(ncol) if integ[c]], dtype=np.int32)
    h.changeColsIntegrality(len(intidx), intidx, np.array([highspy.HighsVarType.kInteger] * len(intidx)))
    rows = []
    def addrow(lo, hi, entries):
        idx = np.array([e[0] for e in entries], dtype=np.int32)
        val = np.array([e[1] for e in entries], dtype=np.float64)
        h.addRow(lo, hi, len(idx), idx, val)
    for i in range(n):
        for t in range(k):
            ent = [(col[("q", i, j, t)], 1.0) for j in range(m) if ("q", i, j, t) in col]
            if ("r", i, t) in col:
                ent.append((col[("r", i, t)], 1.0))
            if ("r", i, t - 1) in col:
                ent.append((col[("r", i, t - 1)], -1.0))
            if not ent:
                assert QH[i][t] == 0
                continue
            addrow(QH[i][t], QH[i][t], ent)
    for j in range(m):
        for t in range(k):
            if inst.QC[j][t] > 0:
                ent = [(col[("q", i, j, t)], 1.0) for i in range(n) if ("q", i, j, t) in col]
                addrow(QC[j][t], QC[j][t], ent)
    for (i, j) in pairs:
        ent = [(col[("q", i, j, t)], 1.0) for t in range(k) if ("q", i, j, t) in col]
        ent.append((col[("y", i, j)], -float(inst.U_pair(i, j)) / scale))
        addrow(-inf, 0.0, ent)
        if disagg:
            for t in range(k):
                if ("q", i, j, t) in col:
                    addrow(-inf, 0.0, [(col[("q", i, j, t)], 1.0), (col[("y", i, j)], -float(inst.u_interval(i, j, t)) / scale)])
        if cumul:
            for t in range(k):
                if ("q", i, j, t) in col:
                    cap = min(inst.cold_cum(j, t), hc[i][t])
                    ent = [(col[("q", i, j, s)], 1.0) for s in range(t + 1) if ("q", i, j, s) in col]
                    ent.append((col[("y", i, j)], -float(cap) / scale))
                    addrow(-inf, 0.0, ent)
    if maxmatches is not None:
        addrow(-inf, float(maxmatches), [(col[("y", i, j)], 1.0) for (i, j) in pairs])
    return h, col, scale
