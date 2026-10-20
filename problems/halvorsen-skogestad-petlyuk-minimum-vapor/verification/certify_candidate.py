"""SUPERSEDED ATTEMPT (kept for the record; the certificate is krawczyk_certificate.py).  With all column roots
fixed at rounded rationals, three feed-equation rows of column BCDE become linearly dependent on the others with
right-hand sides inconsistent at the 1e-9 level, so no exact solution of the fixed-root system exists; the
nonlinear system with free roots has full rank and is certified by krawczyk_certificate.py.

Exact rational certificate attempt for a permissive-model (TG literal) point with R < V_FTC.

Idea: once every Underwood root is fixed, all TG constraints are LINEAR in the flows and vapours.  We
 1. take the candidate point found by open_n5_probe.py (basic configuration 161, j = 3, all thermal couplings);
 2. choose rational relative volatilities and RATIONAL feed roots theta_1..theta_4 close to the candidate's,
    fix q = 0, and solve the 5 linear equations  sum_p a_p z_p/(a_p - theta_r) = 1 - q,  sum_p z_p = 1  for an
    exactly rational feed z (checked > 0);
 3. fix the roots of all other columns at rationals close to the candidate's;
 4. solve the resulting LP (min R) in floating point with strict margins on every inequality (flows >= delta_f,
    other inequalities >= delta);
 5. rationalize the LP solution, project it EXACTLY (Fraction Gaussian elimination, least-norm correction) onto
    the linear equality system, and check every equality and inequality of the TG model in exact arithmetic;
 6. compute R and H_j = sum_{p>j} a_p z_p/(theta_j - a_p) exactly and check R < max_j H_j.
Writes logs/fixed_root_attempt.json and prints PASS/FAIL.
Environment: WSL Fedora, Python 3, numpy, scipy (HiGHS).  Runtime: seconds.
Usage: python certify_candidate.py logs/candidate_161_3_t0.json delta_f delta"""
import sys, json, re
from fractions import Fraction as Fr
import numpy as np
from scipy.optimize import linprog
from indep_nlp import configurations, name

D = json.load(open(sys.argv[1]))
delta_f = float(sys.argv[2]) if len(sys.argv) > 2 else 1e-7
delta = float(sys.argv[3]) if len(sys.argv) > 3 else 1e-7
cfg = [c for c in configurations(5) if name(c) == D['name']][0]
N = 5; full = (1, N); j = D['j']
cols = cfg['streams']; split = cfg['split']

# ---------- 1-2: rational data with rational feed roots
alpha = [Fr(a).limit_denominator(10 ** 6) for a in D['alpha']]
alpha[-1] = Fr(1)
q = Fr(D['q']).limit_denominator(1000)
def feed_roots_float(al, z, q):
    out = []
    for r in range(N - 1):
        lo, hi = float(al[r + 1]), float(al[r])
        a, b = lo, hi
        for _ in range(200):
            c = 0.5 * (a + b)
            v = sum(float(x) * float(zz) / (float(x) - c) for x, zz in zip(al, z)) - (1 - float(q))
            if v < 0: a = c
            else: b = c
        out.append(0.5 * (a + b))
    return out
th_f = feed_roots_float(alpha, D['z'], q)
theta = [Fr(t).limit_denominator(10 ** 9) for t in th_f]
# solve for z: rows r: sum_p a_p z_p/(a_p - theta_r) = 1 - q ; last row: sum z = 1
def solve_exact(Mrows, rhs):
    n = len(Mrows[0]); M = [list(r) + [b] for r, b in zip(Mrows, rhs)]
    m = len(M)
    for col in range(n):
        piv = next(i for i in range(col, m) if M[i][col] != 0)
        M[col], M[piv] = M[piv], M[col]
        pv = M[col][col]
        M[col] = [v / pv for v in M[col]]
        for i in range(m):
            if i != col and M[i][col] != 0:
                f = M[i][col]; M[i] = [a - f * b for a, b in zip(M[i], M[col])]
    return [M[i][n] for i in range(n)]
rows = [[alpha[p] / (alpha[p] - theta[r]) for p in range(N)] for r in range(N - 1)] + [[Fr(1)] * N]
z = solve_exact(rows, [1 - q] * (N - 1) + [Fr(1)])
assert all(v > 0 for v in z), z
H = [sum(alpha[p] * z[p] / (theta[r] - alpha[p]) for p in range(r + 1, N)) for r in range(N - 1)]
V = max(H)
print('rational data: alpha =', [str(a) for a in alpha], ' q =', q)
print('  z (exact) ~', [float(v) for v in z], ' (candidate z', D['z'], ')')
print('  H (exact) ~', [float(h) for h in H], ' V_FTC ~', float(V))

# ---------- 3: fix all column roots (feed column: the rational feed roots)
sol = D['solution']
def getv(nm, idx):
    key = f"{nm}[{idx}]"
    return sol[key]
TH = {}
ncopy = 0
for c in sorted(cols, key=lambda c: c[0] - c[1]):          # producers (longer labels) first
    k, l = split[c]
    for r in range(l - 1, k + 1):
        if c == full:
            TH[(c, r)] = theta[r - 1]
            continue
        v = getv('th', (c[0], c[1], r))
        # carried-over roots: a root equal (numerically) to a producer's root at the same r is set EXACTLY equal,
        # otherwise the fixed-root linear system becomes nearly singular
        cands = [TH[(P, r)] for kind, P in cfg['prod'][c] if (P, r) in TH]
        close = [t for t in cands if abs(float(t) - v) < 1e-7 * (1 + abs(v))]
        if close:
            TH[(c, r)] = close[0]; ncopy += 1
        else:
            TH[(c, r)] = Fr(v).limit_denominator(10 ** 9)
        assert alpha[r] < TH[(c, r)] < alpha[r - 1]
print('roots set equal to a producer root (carry-over):', ncopy)

# ---------- 4: LP
keys = []
for c in cols:
    i, m = c; k, l = split[c]
    keys += [('d', c, p) for p in range(i, k + 1)] + [('b', c, p) for p in range(l, m + 1)]
    keys += [('VT', c), ('VB', c), ('Y', c)]
ix = {k: n for n, k in enumerate(keys)}
nv = len(keys)
def dvar(c, p):
    k = split[c][0]; return {ix[('d', c, p)]: Fr(1)} if c[0] <= p <= k else {}
def bvar(c, p):
    l = split[c][1]; return {ix[('b', c, p)]: Fr(1)} if l <= p <= c[1] else {}
def add(*es, coef=None):
    out = {}
    for e, s in zip(es, coef or [1] * len(es)):
        for k, v in e.items(): out[k] = out.get(k, 0) + s * v
    return out
def Fexpr(c, p):     # (linear dict, constant)
    if c == full: return {}, z[p - 1]
    e = {}
    for kind, P in cfg['prod'][c]:
        e = add(e, dvar(P, p) if kind == 'top' else bvar(P, p))
    return e, Fr(0)
def Uexpr(c):
    if c == full: return {}, 1 - q
    e = {}
    for kind, P in cfg['prod'][c]:
        e = add(e, {ix[('VT', P)]: Fr(1)} if kind == 'top' else {ix[('VB', P)]: Fr(-1)})
    return e, Fr(0)
EQ, INEQ, EQL = [], [], []            # (dict, const) meaning dict.x + const = 0  /  >= margin
for c in cols:
    i, m = c; k, l = split[c]
    for p in range(i, m + 1):
        fe, fc = Fexpr(c, p)
        EQ.append((add(dvar(c, p), bvar(c, p), fe, coef=[1, 1, -1]), -fc)); EQL.append(('bal', c, p))
        if p <= k: INEQ.append((dvar(c, p), Fr(0), 'flow'))
        if p >= l: INEQ.append((bvar(c, p), Fr(0), 'flow'))
    ue, uc = Uexpr(c)
    EQ.append((add({ix[('VT', c)]: Fr(1), ix[('VB', c)]: Fr(-1)}, ue, coef=[1, -1]), -uc)); EQL.append(('vap', c))
    INEQ += [({ix[('VT', c)]: Fr(1)}, Fr(0), 'v'), ({ix[('VB', c)]: Fr(1)}, Fr(0), 'v')]
    for r in range(l - 1, k + 1):
        t = TH[(c, r)]
        fe = {}; fc = Fr(0)
        for p in range(i, m + 1):
            e, cst = Fexpr(c, p); w = alpha[p - 1] / (alpha[p - 1] - t)
            fe = add(fe, e, coef=[1, w]); fc += w * cst
        if c == full:
            assert fc == uc, 'feed root not exact'
        else:
            EQ.append((add(fe, ue, coef=[1, -1]), fc - uc)); EQL.append(('feed', c, r))
        gT = {}
        for p in range(i, k + 1):
            gT = add(gT, dvar(c, p), coef=[1, alpha[p - 1] / (alpha[p - 1] - t)])
        if l <= r <= k - 1:
            EQ.append((add(gT, {ix[('Y', c)]: Fr(1)}, coef=[1, -1]), Fr(0))); EQL.append(('int', c, r))
        else:
            INEQ.append((add({ix[('Y', c)]: Fr(1)}, gT, coef=[1, -1]), Fr(0), 'u'))
    INEQ.append((add({ix[('VT', c)]: Fr(1)}, {ix[('Y', c)]: Fr(1)}, coef=[1, -1]), Fr(0), 'u'))
    INEQ.append((add({ix[('Y', c)]: Fr(1)}, ue, coef=[1, -1]), -uc, 'u'))
    INEQ.append(({ix[('Y', c)]: Fr(1)}, Fr(0), 'u'))
for p in range(1, N + 1):
    pr = cfg['prod'][(p, p)]
    if len(pr) == 2:
        top = [P for kd, P in pr if kd == 'top'][0]; bot = [P for kd, P in pr if kd == 'bot'][0]
        EQ.append(({ix[('VT', top)]: Fr(1), ix[('VB', bot)]: Fr(-1)}, Fr(0))); EQL.append(('side', p))
Rexpr = {}
for p in range(1, N + 1):
    pr = cfg['prod'][(p, p)]
    if len(pr) == 1 and pr[0][0] == 'bot': Rexpr = add(Rexpr, {ix[('VB', pr[0][1])]: Fr(1)})
def dense(e):
    v = np.zeros(nv)
    for k, x in e.items(): v[k] = float(x)
    return v
A_eq = np.array([dense(e) for e, _ in EQ]); b_eq = np.array([-float(cst) for _, cst in EQ])
A_ub = np.array([-dense(e) for e, _, _ in INEQ])
b_ub = np.array([float(cst) - (delta_f if kind == 'flow' else delta) for _, cst, kind in INEQ])
cvec = dense(Rexpr)
res = linprog(cvec, A_ub=A_ub, b_ub=b_ub, A_eq=A_eq, b_eq=b_eq, bounds=[(None, None)] * nv, method='highs')
print('LP status', res.status, res.message, ' R_LP =', res.fun, ' ratio', res.fun / float(V) if res.status == 0 else None)
if res.status != 0:
    print('FAIL'); sys.exit(1)

# ---------- 5: exact projection onto the equalities (least-norm correction)
x = [Fr(v).limit_denominator(10 ** 12) for v in res.x]
def ev(e, cst): return sum(v * x[k] for k, v in e.items()) + cst
r = [-ev(e, cst) for e, cst in EQ]          # need A dx = r
Aq = [[e.get(k, Fr(0)) for k in range(nv)] for e, _ in EQ]
# least norm: dx = A^T (A A^T)^+ r ; remove dependent rows first
def rank_rows(M):
    rows, basis = [], []
    for idx, row in enumerate(M):
        v = list(row)
        for (pc, bv) in basis:
            if v[pc] != 0:
                f = v[pc] / bv[pc]; v = [a - f * b for a, b in zip(v, bv)]
        pc = next((t for t, a in enumerate(v) if a != 0), None)
        if pc is not None:
            basis.append((pc, v)); rows.append(idx)
    return rows
ind = rank_rows(Aq)
Ai = [Aq[i] for i in ind]; ri = [r[i] for i in ind]
AAt = [[sum(a * b for a, b in zip(Ai[s], Ai[t])) for t in range(len(Ai))] for s in range(len(Ai))]
y = solve_exact(AAt, ri)
dx = [sum(Ai[s][k] * y[s] for s in range(len(Ai))) for k in range(nv)]
x = [a + b for a, b in zip(x, dx)]
maxeq = max(abs(ev(e, cst)) for e, cst in EQ)
bad = [(EQL[n], float(ev(e, cst)), n in ind) for n, (e, cst) in enumerate(EQ) if ev(e, cst) != 0]
print('dependent rows:', [EQL[n] for n in range(len(EQ)) if n not in ind])
print('nonzero equality residuals (row, value, independent?):', bad[:10], len(bad))
mins = min(ev(e, cst) for e, cst, _ in INEQ)
R = ev(Rexpr, Fr(0))
print('exact point: max |equality| =', maxeq, ' (dependent equality rows:', len(EQ) - len(ind), ')')
print('             min inequality slack =', float(mins))
print('             R =', float(R), ' V_FTC =', float(V), ' R/V_FTC =', float(R / V), ' V_FTC - R =', float(V - R))
flows_pos = all(ev(e, cst) > 0 for e, cst, kind in INEQ if kind == 'flow')
ok = (maxeq == 0) and (mins >= 0) and (R < V)
print('all flows strictly positive (assumption (N)):', flows_pos)
out = {'configuration': D['name'], 'j': j, 'model': 'permissive (TG literal, excess vapor allowed)',
       'alpha': [str(a) for a in alpha], 'q': str(q), 'z': [str(v) for v in z],
       'feed_roots': [str(t) for t in theta], 'H': [str(h) for h in H], 'V_FTC': str(V), 'R': str(R),
       'roots': {f'{c}:{rr}': str(v) for (c, rr), v in TH.items()},
       'variables': {f'{k[0]}{k[1]}' + (f'_{k[2]}' if len(k) > 2 else ''): str(v) for k, v in zip(keys, x)}}
json.dump(out, open('logs/fixed_root_attempt.json', 'w'), indent=1)
print('PASS: exact rational TG-permissive point with R < V_FTC' if ok else 'FAIL')
sys.exit(0 if ok else 1)
