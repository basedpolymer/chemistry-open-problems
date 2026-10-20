"""Re-examines a candidate point dumped by open_n5_probe.py (permissive model, R < V_FTC) in high precision:
re-evaluates every constraint of the TG model from the stored variables, prints slacks, excesses, the
light-cut decomposition at theta_j, and then polishes the point by a Newton projection in 50-digit arithmetic
onto the equality constraints (inequalities kept with their slack), reporting the final residuals.
Usage: python analyze_candidate.py logs/candidate_161_3_t0.json
Environment: WSL Fedora, Python 3, mpmath."""
import sys, json, re
import mpmath as mp
from indep_nlp import configurations, name

mp.mp.dps = 50
D = json.load(open(sys.argv[1]))
cfg = [c for c in configurations(5) if name(c) == D['name']][0]
N = 5
alpha = [mp.mpf(str(a)) for a in D['alpha']]; z = [mp.mpf(str(x)) for x in D['z']]; q = mp.mpf(str(D['q']))
A = lambda p: alpha[p - 1]
full = (1, N)

# ---- feed roots in high precision
def bis(f, lo, hi):
    a, b = mp.mpf(lo), mp.mpf(hi)
    for _ in range(400):
        c = (a + b) / 2
        if f(c) > 0: a = c
        else: b = c
    return (a + b) / 2
th = [bis(lambda t: -(sum(A(p) * z[p - 1] / (A(p) - t) for p in range(1, N + 1)) - (1 - q)), A(r + 1), A(r))
      for r in range(1, N)]
H = [sum(A(p) * z[p - 1] / (th[j - 1] - A(p)) for p in range(j + 1, N + 1)) for j in range(1, N)]
j = D['j']

# ---- unpack variables: names like d[(1, 5, 2)] in pyomo's flattened form
var = {}
for k, v in D['solution'].items():
    m = re.match(r'(\w+)\[(.*)\]', k)
    nm, idx = m.group(1), m.group(2)
    idx = tuple(int(x) for x in re.findall(r'-?\d+', idx))
    var[(nm,) + idx] = mp.mpf(repr(v)) if v is not None else None
cols = cfg['streams']
split = cfg['split']
# unknown vector: all d, b, VT, VB, Y, theta (non-feed columns)
keys = []
for c in cols:
    i, m = c; k, l = split[c]
    for p in range(i, k + 1): keys.append(('d',) + c + (p,))
    for p in range(l, m + 1): keys.append(('b',) + c + (p,))
    keys += [('VT',) + c, ('VB',) + c, ('Y',) + c]
    if c != full:
        for r in range(l - 1, k + 1): keys.append(('th',) + c + (r,))
x0 = [var[k] for k in keys]

def unpack(x):
    X = dict(zip(keys, x))
    for c in cols:
        if c == full:
            k, l = split[c]
            for r in range(l - 1, k + 1): X[('th',) + c + (r,)] = th[r - 1]
    return X

def model(X):
    """returns (equalities list, inequalities list g >= 0, R)"""
    def d(c, p):
        k = split[c][0]; return X[('d',) + c + (p,)] if c[0] <= p <= k else 0
    def b(c, p):
        l = split[c][1]; return X[('b',) + c + (p,)] if l <= p <= c[1] else 0
    def F(c, p):
        if c == full: return z[p - 1]
        return sum((d(P, p) if kind == 'top' else b(P, p)) for kind, P in cfg['prod'][c])
    def U(c):
        if c == full: return 1 - q
        u = 0
        for kind, P in cfg['prod'][c]:
            u += X[('VT',) + P] if kind == 'top' else -X[('VB',) + P]
        return u
    eq, ineq, labels_eq, labels_in = [], [], [], []
    for c in cols:
        i, m = c; k, l = split[c]
        for p in range(i, m + 1):
            eq.append(d(c, p) + b(c, p) - F(c, p)); labels_eq.append(('bal', c, p))
            if p <= k: ineq.append(d(c, p)); labels_in.append(('d>=0', c, p))
            if p >= l: ineq.append(b(c, p)); labels_in.append(('b>=0', c, p))
        eq.append(X[('VT',) + c] - X[('VB',) + c] - U(c)); labels_eq.append(('vap', c))
        ineq += [X[('VT',) + c], X[('VB',) + c]]; labels_in += [('VT>=0', c), ('VB>=0', c)]
        for r in range(l - 1, k + 1):
            t = X[('th',) + c + (r,)]
            if c != full:
                ineq += [t - A(r + 1), A(r) - t]; labels_in += [('th>lo', c, r), ('th<hi', c, r)]
                eq.append(sum(A(p) * F(c, p) / (A(p) - t) for p in range(i, m + 1)) - U(c)); labels_eq.append(('feed', c, r))
            gT = sum(A(p) * d(c, p) / (A(p) - t) for p in range(i, k + 1))
            if l <= r <= k - 1:
                eq.append(gT - X[('Y',) + c]); labels_eq.append(('int', c, r))
            else:
                ineq.append(X[('Y',) + c] - gT); labels_in.append(('bnd', c, r))
        ineq += [X[('VT',) + c] - X[('Y',) + c], X[('Y',) + c] - U(c), X[('Y',) + c]]
        labels_in += [('Y<=VT', c), ('YB>=0', c), ('Y>=0', c)]
    for p in range(1, N + 1):
        pr = cfg['prod'][(p, p)]
        if len(pr) == 2:
            top = [P for kd, P in pr if kd == 'top'][0]; bot = [P for kd, P in pr if kd == 'bot'][0]
            eq.append(X[('VT',) + top] - X[('VB',) + bot]); labels_eq.append(('sidedraw', p))
    R = 0
    for p in range(1, N + 1):
        pr = cfg['prod'][(p, p)]
        if len(pr) == 1 and pr[0][0] == 'bot': R += X[('VB',) + pr[0][1]]
    return eq, ineq, R, labels_eq, labels_in, d, b, U

X = unpack(x0)
eq, ineq, R, le, li, d, b, U = model(X)
print('config:', D['name'], ' j =', j)
print('alpha', [mp.nstr(a, 8) for a in alpha], 'z', [mp.nstr(v, 6) for v in z], 'q', q)
print('H =', [mp.nstr(h, 10) for h in H], ' V_FTC =', mp.nstr(max(H), 12), ' R =', mp.nstr(R, 12),
      ' R/V =', mp.nstr(R / max(H), 12))
print('max |equality residual| =', mp.nstr(max(abs(e) for e in eq), 5),
      ' min inequality =', mp.nstr(min(ineq), 5))
small = sorted(zip(ineq, li), key=lambda t: t[0])[:12]
print('smallest inequality slacks:', [(mp.nstr(v, 4), lab) for v, lab in small])
for c in cols:
    print('  col', c, 'split', split[c], ' VT', mp.nstr(X[('VT',) + c], 8), ' VB', mp.nstr(X[('VB',) + c], 8),
          ' excess VT-Y', mp.nstr(X[('VT',) + c] - X[('Y',) + c], 6),
          ' d', [mp.nstr(d(c, p), 5) for p in range(c[0], split[c][0] + 1)],
          ' b', [mp.nstr(b(c, p), 5) for p in range(split[c][1], c[1] + 1)])
# light-cut decomposition at theta_j
tj = th[j - 1]
terms = []
for c in cols:
    i, m = c; k, l = split[c]
    if i <= j < m and k <= j:
        lam = X[('VT',) + c] - sum(A(p) * d(c, p) / (A(p) - tj) for p in range(i, k + 1))
        terms.append((c, lam))
print('light exits Lambda(theta_j):', [(c, mp.nstr(v, 8)) for c, v in terms], ' sum', mp.nstr(sum(v for _, v in terms), 8))
print('R - H_j =', mp.nstr(R - H[j - 1], 8), '(all-TC: no condensers on mixed/heavy nodes, no reboilers on light nodes)')

# ---- Newton projection onto the equalities (least-norm steps), 50 digits
x = list(x0)
for it in range(30):
    X = unpack(x); eq = model(X)[0]
    res = max(abs(e) for e in eq)
    if res < mp.mpf('1e-45'): break
    n = len(x); h = mp.mpf('1e-30')
    J = mp.matrix(len(eq), n)
    for jj in range(n):
        xp = list(x); xp[jj] += h
        ep = model(unpack(xp))[0]
        for ii in range(len(eq)): J[ii, jj] = (ep[ii] - eq[ii]) / h
    e = mp.matrix(eq)
    # least-norm step: dx = -J^T (J J^T)^-1 e
    JJt = J * J.T
    y = mp.lu_solve(JJt, e)
    dx = J.T * y
    x = [x[i] - dx[i] for i in range(n)]
X = unpack(x); eq, ineq, R, le, li, d, b, U = model(X)
print('after Newton polish: max |eq| =', mp.nstr(max(abs(e) for e in eq), 5), ' min ineq =', mp.nstr(min(ineq), 6),
      ' R/V_FTC =', mp.nstr(R / max(H), 15), ' R - V_FTC =', mp.nstr(R - max(H), 10))
json.dump({'keys': [list(map(str, k)) for k in keys], 'x': [mp.nstr(v, 45) for v in x]},
          open(sys.argv[1].replace('.json', '_polished.json'), 'w'), indent=0)
