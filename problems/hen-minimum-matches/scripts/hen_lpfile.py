"""Write the transshipment MILP to a CPLEX-LP file with EXACT decimal coefficients (for SCIP exact mode).

All numbers are Fractions whose denominators divide a power of 10 (instance data are decimals; the bounds U_ij,
u_ijt are sums/minima of data), so they are printed exactly as decimal strings.  Heat is scaled by 10^-p (exact).

Constraints (y binary, q, r >= 0):
  hot_i_t  : sum_j q_ijt + r_it - r_i(t-1) = QH[i,t]
  cold_j_t : sum_i q_ijt = QC[j,t]
  link_i_j : sum_t q_ijt - U_ij y_ij <= 0          U_ij = min_{t0=-1..k-1} ( QH_i(<=t0) + QC_j(>t0) )
  dl_i_j_t : q_ijt - u_ijt y_ij <= 0               u_ijt = min( QC[j,t], QH_i(<=t) )          (if disagg)
  cl_i_j_t : sum_{s<=t} q_ijs - min(QC_j(<=t), QH_i(<=t)) y_ij <= 0                          (if cumul)
  maxm     : sum y_ij <= K                                                                     (if maxmatches=K)
Variables q_ijt exist only if QC[j,t] > 0 and QH_i(<=t) > 0; r_it only for t <= k-2 and QH_i(<=t) > 0
(all omitted variables are zero in every feasible flow).
"""
from fractions import Fraction


def exact_decimal(x):
    x = Fraction(x)
    d = x.denominator
    a = b = 0
    while d % 2 == 0:
        d //= 2; a += 1
    while d % 5 == 0:
        d //= 5; b += 1
    assert d == 1, f"not a finite decimal: {x}"
    e = max(a, b)
    num = x.numerator * (10 ** e) // x.denominator
    assert Fraction(num, 10 ** e) == x
    s = str(abs(num)).rjust(e + 1, "0")
    out = s[:-e] + "." + s[-e:] if e > 0 else s
    if e > 0:
        out = out.rstrip("0").rstrip(".")
    return ("-" if num < 0 else "") + out


def U_mincut(inst, i, j):
    """min over t0 in {-1..k-1} of QH_i(<=t0) + QC_j(>t0): an upper bound on heat i->j (one-line proof:
    heat j receives from i in intervals <= t0 originates in intervals <= t0; in intervals > t0 it is <= demand)."""
    k = inst.k
    best = None
    for t0 in range(-1, k):
        val = sum(inst.QH[i][: t0 + 1]) + sum(inst.QC[j][t0 + 1:])
        best = val if best is None else min(best, val)
    return best


def write_lp(inst, path, disagg=True, cumul=False, maxmatches=None, scale_pow=None, objective=True, extra=None):
    n, m, k = inst.n, inst.m, inst.k
    if scale_pow is None:
        mx = max(max(max(r) for r in inst.QH), max(max(r) for r in inst.QC))
        scale_pow = 0
        while mx / Fraction(10) ** scale_pow > 10000:
            scale_pow += 1
    S = Fraction(10) ** scale_pow
    D = exact_decimal
    hc = [[sum(inst.QH[i][: t + 1]) for t in range(k)] for i in range(n)]
    pairs = [(i, j) for i in range(n) for j in range(m) if U_mincut(inst, i, j) > 0]
    qv = {(i, j, t) for (i, j) in pairs for t in range(k) if inst.QC[j][t] > 0 and hc[i][t] > 0}
    rv = {(i, t) for i in range(n) for t in range(k - 1) if hc[i][t] > 0}
    lines = []
    w = lines.append
    w("\\ transshipment MILP, min number of matches, instance " + inst.name + f", heat scaled by 10^-{scale_pow}")
    w("Minimize" if objective else "Minimize")
    if objective:
        w(" obj: " + " + ".join(f"y_{i}_{j}" for (i, j) in pairs))
    else:
        w(" obj: 0 y_%d_%d" % pairs[0])
    w("Subject To")
    def row(name, terms, sense, rhs):
        # terms: list of (coef Fraction, varname)
        parts = []
        for c, v in terms:
            if c == 0:
                continue
            sgn = "-" if c < 0 else "+"
            parts.append(f"{sgn} {D(abs(c))} {v}")
        s = " ".join(parts)
        if s.startswith("+ "):
            s = s[2:]
        # wrap long rows
        out = f" {name}: " + s + f" {sense} {D(rhs)}"
        w(out)
    for i in range(n):
        for t in range(k):
            terms = [(Fraction(1), f"q_{i}_{j}_{t}") for j in range(m) if (i, j, t) in qv]
            if (i, t) in rv:
                terms.append((Fraction(1), f"r_{i}_{t}"))
            if (i, t - 1) in rv:
                terms.append((Fraction(-1), f"r_{i}_{t-1}"))
            if not terms:
                assert inst.QH[i][t] == 0
                continue
            row(f"hot_{i}_{t}", terms, "=", inst.QH[i][t] / S)
    for j in range(m):
        for t in range(k):
            if inst.QC[j][t] > 0:
                terms = [(Fraction(1), f"q_{i}_{j}_{t}") for i in range(n) if (i, j, t) in qv]
                row(f"cold_{j}_{t}", terms, "=", inst.QC[j][t] / S)
    for (i, j) in pairs:
        U = U_mincut(inst, i, j) / S
        terms = [(Fraction(1), f"q_{i}_{j}_{t}") for t in range(k) if (i, j, t) in qv] + [(-U, f"y_{i}_{j}")]
        row(f"link_{i}_{j}", terms, "<=", 0)
        if disagg:
            for t in range(k):
                if (i, j, t) in qv:
                    u = min(inst.QC[j][t], hc[i][t]) / S
                    row(f"dl_{i}_{j}_{t}", [(Fraction(1), f"q_{i}_{j}_{t}"), (-u, f"y_{i}_{j}")], "<=", 0)
        if cumul:
            for t in range(k):
                if (i, j, t) in qv:
                    cap = min(sum(inst.QC[j][: t + 1]), hc[i][t]) / S
                    terms = [(Fraction(1), f"q_{i}_{j}_{s}") for s in range(t + 1) if (i, j, s) in qv] + [(-cap, f"y_{i}_{j}")]
                    row(f"cl_{i}_{j}_{t}", terms, "<=", 0)
    if maxmatches is not None:
        row("maxm", [(Fraction(1), f"y_{i}_{j}") for (i, j) in pairs], "<=", Fraction(maxmatches))
    if extra:
        for name, terms, sense, rhs in extra:
            row(name, terms, sense, rhs)
    w("Bounds")
    for (i, j, t) in sorted(qv):
        w(f" 0 <= q_{i}_{j}_{t} <= {D(min(inst.QC[j][t], hc[i][t]) / S)}")
    for (i, t) in sorted(rv):
        w(f" 0 <= r_{i}_{t} <= {D(hc[i][t] / S)}")
    w("Binaries")
    for (i, j) in pairs:
        w(f" y_{i}_{j}")
    w("End")
    with open(path, "w") as f:
        f.write("\n".join(lines) + "\n")
    return {"pairs": pairs, "nq": len(qv), "nr": len(rv), "scale_pow": scale_pow}


def write_lp_strong(inst, path, rows=(), maxmatches=None, scale_pow=None, objective=True):
    """Strong transshipment MILP (residual-capped bounds, see hen_cuts.py) with extra y-rows, exact decimals.
      hot_i_t, cold_j_t as in write_lp;  link_i_j: sum_t q_ijt - U_dp(i,j) y_ij <= 0;
      dl_i_j_t: q_ijt - u_cap(i,j,t) y_ij <= 0;  bounds 0 <= q_ijt <= u_cap, 0 <= r_it <= rc(i,t);
      rows: sum coef*y >= rhs (W1/W2/P families);  maxm: sum y <= K."""
    from hen_cuts import U_dp, u_interval_capped, rcap
    n, m, k = inst.n, inst.m, inst.k
    if scale_pow is None:
        mx = max(max(max(r) for r in inst.QH), max(max(r) for r in inst.QC))
        scale_pow = 0
        while mx / Fraction(10) ** scale_pow > 10000:
            scale_pow += 1
    S = Fraction(10) ** scale_pow
    D = exact_decimal
    U = {(i, j): U_dp(inst, i, j) for i in range(n) for j in range(m)}
    pairs = [p for p in sorted(U) if U[p] > 0]
    pset = set(pairs)
    uq = {(i, j, t): u_interval_capped(inst, i, j, t) for (i, j) in pairs for t in range(k)}
    qv = sorted(key for key, u in uq.items() if u > 0)
    qset = set(qv)
    rv = sorted((i, t) for i in range(n) for t in range(k - 1) if rcap(inst, i, t) > 0)
    rset = set(rv)
    L = []
    w = L.append
    w("\\ strong transshipment MILP (residual-capped), instance " + inst.name + f", heat scaled by 10^-{scale_pow}")
    w("Minimize")
    w(" obj: " + (" + ".join("y_%d_%d" % p for p in pairs) if objective else "0 y_%d_%d" % pairs[0]))
    w("Subject To")

    def row(name, terms, sense, rhs):
        parts = []
        for c, v in terms:
            if c == 0:
                continue
            parts.append(("- " if c < 0 else "+ ") + f"{D(abs(c))} {v}")
        s = " ".join(parts)
        if s.startswith("+ "):
            s = s[2:]
        w(f" {name}: {s} {sense} {D(rhs)}")
    for i in range(n):
        for t in range(k):
            terms = [(Fraction(1), f"q_{i}_{j}_{t}") for j in range(m) if (i, j, t) in qset]
            if (i, t) in rset:
                terms.append((Fraction(1), f"r_{i}_{t}"))
            if (i, t - 1) in rset:
                terms.append((Fraction(-1), f"r_{i}_{t-1}"))
            if not terms:
                assert inst.QH[i][t] == 0
                continue
            row(f"hot_{i}_{t}", terms, "=", inst.QH[i][t] / S)
    for j in range(m):
        for t in range(k):
            if inst.QC[j][t] > 0:
                terms = [(Fraction(1), f"q_{i}_{j}_{t}") for i in range(n) if (i, j, t) in qset]
                row(f"cold_{j}_{t}", terms, "=", inst.QC[j][t] / S)
    for (i, j) in pairs:
        terms = [(Fraction(1), f"q_{i}_{j}_{t}") for t in range(k) if (i, j, t) in qset] + [(-U[i, j] / S, f"y_{i}_{j}")]
        row(f"link_{i}_{j}", terms, "<=", 0)
        for t in range(k):
            if (i, j, t) in qset:
                row(f"dl_{i}_{j}_{t}", [(Fraction(1), f"q_{i}_{j}_{t}"), (-uq[i, j, t] / S, f"y_{i}_{j}")], "<=", 0)
    for name, coef, b in rows:
        terms = [(c, "y_%d_%d" % e) for e, c in sorted(coef.items()) if e in pset]
        row(name, terms, ">=", b)
    if maxmatches is not None:
        row("maxm", [(Fraction(1), "y_%d_%d" % p) for p in pairs], "<=", Fraction(maxmatches))
    w("Bounds")
    for (i, j, t) in qv:
        w(f" 0 <= q_{i}_{j}_{t} <= {D(uq[i, j, t] / S)}")
    for (i, t) in rv:
        w(f" 0 <= r_{i}_{t} <= {D(rcap(inst, i, t) / S)}")
    w("Binaries")
    for p in pairs:
        w(" y_%d_%d" % p)
    w("End")
    with open(path, "w") as f:
        f.write("\n".join(L) + "\n")
    return {"pairs": pairs, "nq": len(qv), "nr": len(rv), "scale_pow": scale_pow, "nrows": len(rows)}
