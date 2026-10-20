"""Rigorous (interval) certificate that the PERMISSIVE model (Tumbalam Gooty et al.'s Underwood model read
literally: excess vapor allowed in every column) admits, for N = 5, a feasible point with R < V_FTC = max_j H_j.

Data: rational relative volatilities, q = 0, rational feed roots theta_1..theta_4 and the exactly rational feed z
that has them as roots (so H_j and V_FTC are exact rationals).  Configuration: basic 161 of the Shah-Agrawal
space, all thermal couplings,
  ABCDE->ABCD/BCDE; ABCD->ABC/BCD; BCDE->BCD/CDE; ABC->AB/BC; BCD->BC/D; CDE->C/DE; AB->A/B; BC->B/C; DE->D/E.
Method: the TG constraints form an underdetermined system of polynomial-rational equalities E(x) = 0 plus
inequalities G(x) >= 0.  From a refined numerical point (all inequality slacks >= ~1e-6, all flows >= 1e-6) we
 1. polish E(x) = 0 by least-norm Newton steps in 60-digit arithmetic;
 2. choose m pivot variables by QR with column pivoting of the Jacobian, fix the other variables at their (binary,
    exactly representable) values, and Newton-solve the square system for the pivots;
 3. prove existence of an exact zero of the square system in a box of radius 1e-30 around it with the Krawczyk
    operator, evaluated in interval arithmetic (mpmath.iv) with an interval Jacobian from forward-mode automatic
    differentiation over intervals;
 4. bound every inequality G and the objective R over that box in interval arithmetic and check G > 0, R < V_FTC.
Prints PASS/FAIL and exits non-zero on failure.  Writes verification/permissive_counterexample_N5.json.
Environment: WSL Fedora, Python 3, mpmath, numpy.  Runtime: about 1 minute.
Usage: python krawczyk_certificate.py logs/candidate_161_3_t0_refined.json"""
import sys, json, re
from fractions import Fraction as Fr
import numpy as np
import mpmath as mp
from mpmath import iv
from indep_nlp import configurations, name

mp.mp.dps = 60; iv.dps = 60
D = json.load(open(sys.argv[1]))
cfg = [c for c in configurations(5) if name(c) == D['name']][0]
N = 5; full = (1, N); j = D['j']
cols = cfg['streams']; split = cfg['split']

# ---------------- exact rational data (same construction as refine_candidate.py)
alphaQ = [Fr(a).limit_denominator(10 ** 6) for a in D['alpha']]; alphaQ[-1] = Fr(1)
qQ = Fr(D['q']).limit_denominator(1000)
def roots_f(al, z, q):
    out = []
    for r in range(N - 1):
        a, b = float(al[r + 1]), float(al[r])
        for _ in range(200):
            c = 0.5 * (a + b)
            v = sum(float(x) * float(zz) / (float(x) - c) for x, zz in zip(al, z)) - (1 - float(q))
            if v < 0: a = c
            else: b = c
        out.append(0.5 * (a + b))
    return out
thetaQ = [Fr(t).limit_denominator(10 ** 9) for t in roots_f(alphaQ, D['z'], qQ)]
rows = [[alphaQ[p] / (alphaQ[p] - thetaQ[r]) for p in range(N)] for r in range(N - 1)] + [[Fr(1)] * N]
M_ = [list(r) + [b] for r, b in zip(rows, [1 - qQ] * (N - 1) + [Fr(1)])]
for col in range(N):
    piv = next(i for i in range(col, N) if M_[i][col] != 0); M_[col], M_[piv] = M_[piv], M_[col]
    pv = M_[col][col]; M_[col] = [v / pv for v in M_[col]]
    for i in range(N):
        if i != col and M_[i][col] != 0:
            f = M_[i][col]; M_[i] = [a - f * b for a, b in zip(M_[i], M_[col])]
zQ = [M_[i][N] for i in range(N)]
assert all(v > 0 for v in zQ) and sum(zQ) == 1
for r in range(N - 1):
    assert sum(alphaQ[p] * zQ[p] / (alphaQ[p] - thetaQ[r]) for p in range(N)) == 1 - qQ
    assert alphaQ[r + 1] < thetaQ[r] < alphaQ[r]
HQ = [sum(alphaQ[p] * zQ[p] / (thetaQ[r] - alphaQ[p]) for p in range(r + 1, N)) for r in range(N - 1)]
VQ = max(HQ)

# ---------------- numeric types: mp (point) and interval; constants converted per type
def cst(x, kind):
    if kind == 'iv':
        return iv.mpf(x.numerator) / iv.mpf(x.denominator)
    if kind == 'ivd':                    # interval constant wrapped for automatic differentiation
        return Dual(iv.mpf(x.numerator) / iv.mpf(x.denominator), {})
    return mp.mpf(x.numerator) / mp.mpf(x.denominator)

class Dual:
    """forward-mode AD over intervals: value v (iv.mpf) and sparse gradient g {index: iv.mpf}"""
    __slots__ = ('v', 'g')
    def __init__(self, v, g=None): self.v = v; self.g = g or {}
    @staticmethod
    def lift(o): return o if isinstance(o, Dual) else Dual(o, {})
    def __add__(s, o):
        o = Dual.lift(o); g = dict(s.g)
        for k, x in o.g.items(): g[k] = g[k] + x if k in g else x
        return Dual(s.v + o.v, g)
    __radd__ = __add__
    def __neg__(s): return Dual(-s.v, {k: -x for k, x in s.g.items()})
    def __sub__(s, o): return s + (-Dual.lift(o))
    def __rsub__(s, o): return Dual.lift(o) + (-s)
    def __mul__(s, o):
        o = Dual.lift(o); g = {k: x * o.v for k, x in s.g.items()}
        for k, x in o.g.items(): g[k] = g[k] + x * s.v if k in g else x * s.v
        return Dual(s.v * o.v, g)
    __rmul__ = __mul__
    def __truediv__(s, o):
        o = Dual.lift(o); inv = 1 / o.v
        g = {k: x * inv for k, x in s.g.items()}
        for k, x in o.g.items():
            t = -x * s.v * inv * inv
            g[k] = g[k] + t if k in g else t
        return Dual(s.v * inv, g)
    def __rtruediv__(s, o): return Dual.lift(o) / s

# ---------------- unknowns and model
keys = []
for c in cols:
    i, m = c; k, l = split[c]
    keys += [('d', c, p) for p in range(i, k + 1)] + [('b', c, p) for p in range(l, m + 1)]
    keys += [('VT', c), ('VB', c), ('Y', c)]
    if c != full:
        keys += [('th', c, r) for r in range(l - 1, k + 1)]
ix = {k: n for n, k in enumerate(keys)}

def model(x, kind):
    """x: list of numbers (mp, iv or Dual).  Returns (equalities, inequalities, R) in the same type."""
    A = [cst(a, kind) for a in alphaQ]; Z = [cst(v, kind) for v in zQ]; one_q = cst(1 - qQ, kind)
    TH = {}
    for c in cols:
        k, l = split[c]
        for r in range(l - 1, k + 1):
            TH[(c, r)] = cst(thetaQ[r - 1], kind) if c == full else x[ix[('th', c, r)]]
    def d(c, p):
        k = split[c][0]; return x[ix[('d', c, p)]] if c[0] <= p <= k else 0
    def b(c, p):
        l = split[c][1]; return x[ix[('b', c, p)]] if l <= p <= c[1] else 0
    def F(c, p):
        if c == full: return Z[p - 1]
        s = 0
        for kd, P in cfg['prod'][c]: s = s + (d(P, p) if kd == 'top' else b(P, p))
        return s
    def U(c):
        if c == full: return one_q
        s = 0
        for kd, P in cfg['prod'][c]:
            s = s + (x[ix[('VT', P)]] if kd == 'top' else -x[ix[('VB', P)]])
        return s
    E, G = [], []
    for c in cols:
        i, m = c; k, l = split[c]
        VT, VB, Y = x[ix[('VT', c)]], x[ix[('VB', c)]], x[ix[('Y', c)]]
        for p in range(i, m + 1):
            E.append(d(c, p) + b(c, p) - F(c, p))
            if p <= k: G.append(d(c, p))
            if p >= l: G.append(b(c, p))
        Uc = U(c)
        E.append(VT - VB - Uc)
        G += [VT, VB]
        for r in range(l - 1, k + 1):
            t = TH[(c, r)]
            if c != full:
                G += [t - A[r], A[r - 1] - t]
                fe = 0
                for p in range(i, m + 1): fe = fe + A[p - 1] * F(c, p) / (A[p - 1] - t)
                E.append(fe - Uc)
            gT = 0
            for p in range(i, k + 1): gT = gT + A[p - 1] * d(c, p) / (A[p - 1] - t)
            if l <= r <= k - 1: E.append(gT - Y)
            else: G.append(Y - gT)
        G += [VT - Y, Y - Uc, Y]
    for p in range(1, N + 1):
        pr = cfg['prod'][(p, p)]
        if len(pr) == 2:
            top = [P for kd, P in pr if kd == 'top'][0]; bot = [P for kd, P in pr if kd == 'bot'][0]
            E.append(x[ix[('VT', top)]] - x[ix[('VB', bot)]])
    R = 0
    for p in range(1, N + 1):
        pr = cfg['prod'][(p, p)]
        if len(pr) == 1 and pr[0][0] == 'bot': R = R + x[ix[('VB', pr[0][1])]]
    return E, G, R

# ---------------- starting point from the refined candidate
sol = D['solution']
def sv(k):
    if k[0] == 'th': return sol[f"th[{(k[1][0], k[1][1], k[2])}]"]
    if k[0] in ('d', 'b'): return sol[f"{k[0]}[{(k[1][0], k[1][1], k[2])}]"]
    return sol[f"{k[0]}[{k[1]}]"]
x = [mp.mpf(repr(sv(k))) for k in keys]
n = len(x)

def jac_mp(x):
    xd = [Dual(v, {t: mp.mpf(1)}) for t, v in enumerate(x)]
    E, _, _ = model(xd, 'mp')
    J = mp.matrix(len(E), n)
    for a, e in enumerate(E):
        e = Dual.lift(e)
        for t, g in e.g.items(): J[a, t] = g
    return [Dual.lift(e).v for e in E], J

# 1. least-norm Newton polish
for it in range(25):
    Ev, J = jac_mp(x)
    res = max(abs(v) for v in Ev)
    if res < mp.mpf('1e-55'): break
    y = mp.lu_solve(J * J.T, mp.matrix(Ev))
    dx = J.T * y
    x = [x[t] - dx[t] for t in range(n)]
Ev, J = jac_mp(x)
m = len(Ev)
print(f'unknowns {n}, equalities {m}; polished residual {mp.nstr(max(abs(v) for v in Ev), 3)}')
_, Gp, Rp = model(x, 'mp')
print(f'point: min inequality {mp.nstr(min(Gp), 6)}, R = {mp.nstr(Rp, 15)}, V_FTC = {mp.nstr(cst(VQ, "mp"), 15)}')

# 2. pivots by QR with column pivoting (numpy, double precision is enough to choose them)
Jf = np.array([[float(J[a, t]) for t in range(n)] for a in range(m)])
import scipy.linalg as sla
_, Rq, perm = sla.qr(Jf, pivoting=True)
piv = sorted(perm[:m]); free = [t for t in range(n) if t not in piv]
print('condition number of the pivot block ~', np.linalg.cond(Jf[:, piv]))
xfix = [mp.mpf(float(v)) for v in x]        # exactly representable binary values for the free variables
def full_x(xp):
    out = list(xfix)
    for a, t in enumerate(piv): out[t] = xp[a]
    return out
xp = [x[t] for t in piv]
for it in range(30):
    xx = full_x(xp)
    xd = [Dual(v, {t: mp.mpf(1)}) if t in piv else v for t, v in enumerate(xx)]
    E, _, _ = model(xd, 'mp')
    Ev = [Dual.lift(e).v for e in E]
    if max(abs(v) for v in Ev) < mp.mpf('1e-55'): break
    Jp = mp.matrix(m, m)
    for a, e in enumerate(E):
        for t, g in Dual.lift(e).g.items(): Jp[a, piv.index(t)] = g
    dx = mp.lu_solve(Jp, mp.matrix(Ev))
    xp = [xp[a] - dx[a] for a in range(m)]
xx = full_x(xp)
print('square-system residual', mp.nstr(max(abs(v) for v in model(xx, 'mp')[0]), 3))

# 3. Krawczyk test in interval arithmetic
rho = mp.mpf('1e-30')
Xc = [iv.mpf(v) for v in xp]
Xbox = [iv.mpf([v - rho, v + rho]) for v in xp]
def iv_full(xpi):
    out = [iv.mpf(v) for v in xfix]
    for a, t in enumerate(piv): out[t] = xpi[a]
    return out
Ec, _, _ = model(iv_full(Xc), 'iv')                       # f(x~) enclosure
xd = [Dual(v, {t: iv.mpf(1)}) if t in piv else v for t, v in enumerate(iv_full(Xbox))]
Eb, _, _ = model([v if isinstance(v, Dual) else Dual(v, {}) for v in xd], 'ivd')
JX = [[iv.mpf(0)] * m for _ in range(m)]
for a, e in enumerate(Eb):
    for t, g in Dual.lift(e).g.items(): JX[a][piv.index(t)] = g
# preconditioner Y ~ inverse of midpoint Jacobian (point matrix, any choice is valid)
Jmid = mp.matrix(m, m)
for a in range(m):
    for b_ in range(m): Jmid[a, b_] = JX[a][b_].mid
Yp = mp.inverse(Jmid)
Yi = [[iv.mpf(Yp[a, b_]) for b_ in range(m)] for a in range(m)]
inside = True
for a in range(m):
    s = Xc[a] - sum((Yi[a][b_] * Ec[b_] for b_ in range(m)), iv.mpf(0))
    for b_ in range(m):
        coef = (iv.mpf(1) if a == b_ else iv.mpf(0)) - sum((Yi[a][c_] * JX[c_][b_] for c_ in range(m)), iv.mpf(0))
        s += coef * (Xbox[b_] - Xc[b_])
    if not (s.a > Xbox[a].a and s.b < Xbox[a].b):
        inside = False
print('Krawczyk K(X) inside int(X):', inside)

# 4. inequalities and objective over the box
Eb2, Gb, Rb = model(iv_full(Xbox), 'iv')
gmin = min(g.a for g in Gb)
Vi = cst(VQ, 'iv')
print(f'over the box: min inequality lower bound {mp.nstr(gmin, 6)}; R in [{mp.nstr(Rb.a, 15)}, {mp.nstr(Rb.b, 15)}]; '
      f'V_FTC in [{mp.nstr(Vi.a, 15)}, {mp.nstr(Vi.b, 15)}]')
ok = inside and gmin > 0 and Rb.b < Vi.a
out = {'claim': 'exact feasible point of the permissive (TG literal) model with R < V_FTC, N = 5, j = %d' % j,
       'configuration': D['name'], 'q': str(qQ), 'alpha': [str(a) for a in alphaQ], 'z': [str(v) for v in zQ],
       'feed_roots': [str(t) for t in thetaQ], 'H': [str(h) for h in HQ], 'V_FTC': str(VQ),
       'R_enclosure': [mp.nstr(Rb.a, 30), mp.nstr(Rb.b, 30)], 'box_radius': '1e-30',
       'center': {f'{k[0]}{k[1]}' + (f'_{k[2]}' if len(k) > 2 else ''): mp.nstr(v, 40) for k, v in zip(keys, xx)},
       'fixed_variables': [f'{keys[t][0]}{keys[t][1]}' + (f'_{keys[t][2]}' if len(keys[t]) > 2 else '') for t in free]}
json.dump(out, open('permissive_counterexample_N5.json', 'w'), indent=1)
print('PASS: rigorous permissive-model point with R < V_FTC' if ok else 'FAIL')
sys.exit(0 if ok else 1)
