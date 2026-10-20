"""Verifier's independent minimum-vapour solver for regular-column configurations (Shah-Agrawal space).

Written from the Underwood equations and the TG conventions, sharing no code with the solvers' scripts:
own enumeration of configurations, own NLP (pyomo + ipopt, multistart), strict and permissive models,
a-posteriori feasibility check of every returned point, and an independent floating-point evaluation of
the light cut identity  R - H_j = R(light) + Q(mixed) + Q(heavy) + sum_{light exits} Lam(theta_j)  at the
returned points.

Usage: python indep_nlp.py enum            (counts 3/8, 18/152)
       python indep_nlp.py run N case nstart model   (model = strict | permissive)
Environment: WSL Fedora, Python 3, pyomo, ipopt.  Numerical evidence only (local solver).
"""
import sys, itertools, json, random, math
import pyomo.environ as pe

# ------------------------------------------------------------------ enumeration
def configurations(N):
    full = (1, N)
    subs = [(i, m) for i in range(1, N + 1) for m in range(i + 1, N + 1) if (i, m) != full]
    out = []
    for r in range(len(subs) + 1):
        for S in itertools.combinations(subs, r):
            streams = set(S) | {full}
            split = {}
            valid = True
            for (i, m) in streams:
                k = max([kk for kk in range(i, m) if kk == i or (i, kk) in streams])
                l = min([ll for ll in range(i + 1, m + 1) if ll == m or (ll, m) in streams])
                if l > k + 1:
                    valid = False; break
                split[(i, m)] = (k, l)
            if not valid:
                continue
            prod = {}
            for (i, m), (k, l) in split.items():
                prod.setdefault((i, k), []).append(('top', (i, m)))
                prod.setdefault((l, m), []).append(('bot', (i, m)))
            if any(s not in prod for s in S) or any((p, p) not in prod for p in range(1, N + 1)):
                continue
            single = sorted(s for s in S if len(prod[s]) == 1)
            for choice in itertools.product(['TC', 'HX'], repeat=len(single)):
                out.append({'N': N, 'streams': sorted(streams), 'split': split, 'prod': prod,
                            'hx': dict(zip(single, choice)), 'basic': tuple(sorted(S))})
    return out

def name(cfg):
    L = 'ABCDEFGH'
    nm = lambda s: L[s[0] - 1:s[1]]
    parts = []
    for s in sorted(cfg['streams'], key=lambda s: (s[0] - s[1], s)):
        k, l = cfg['split'][s]
        parts.append(f"{nm(s)}->{nm((s[0], k))}/{nm((l, s[1]))}")
    hx = ','.join(f"{nm(s)}:{t}" for s, t in sorted(cfg['hx'].items()))
    return '; '.join(parts) + (' | ' + hx if hx else '')

# ------------------------------------------------------------------ feed data
def feed_roots(alpha, z, q):
    N = len(alpha); out = []
    for r in range(N - 1):
        lo, hi = alpha[r + 1], alpha[r]
        a, b = lo, hi
        for _ in range(200):
            c = 0.5 * (a + b)
            v = sum(al * zz / (al - c) for al, zz in zip(alpha, z)) - (1 - q)
            if v < 0: a = c
            else: b = c
        out.append(0.5 * (a + b))
    return out

def peaks(alpha, z, q):
    th = feed_roots(alpha, z, q); N = len(alpha)
    H = [sum(alpha[p] * z[p] / (th[j] - alpha[p]) for p in range(j + 1, N)) for j in range(N - 1)]
    return th, H

# ------------------------------------------------------------------ NLP
def build(cfg, alpha, z, q, strict, rng, th_feed, Vref, df=0.0, dm=0.0):
    """df, dm: optional margins (flows >= df; inequalities >= dm), used only by refine_candidate.py."""
    N = cfg['N']; A = lambda p: alpha[p - 1]
    M = pe.ConcreteModel()
    cols = cfg['streams']
    full = (1, N)
    eps = 1e-9
    M.d = pe.Var([(c, p) for c in cols for p in range(c[0], cfg['split'][c][0] + 1)], bounds=(df, None))
    M.b = pe.Var([(c, p) for c in cols for p in range(cfg['split'][c][1], c[1] + 1)], bounds=(df, None))
    M.VT = pe.Var(cols, bounds=(0, None)); M.VB = pe.Var(cols, bounds=(0, None)); M.Y = pe.Var(cols, bounds=(0, None))
    rootidx = [(c, r) for c in cols for r in range(cfg['split'][c][1] - 1, cfg['split'][c][0] + 1)]
    def thb(m, i, mm, r):          # pyomo flattens the index (c, r) = (i, m, r)
        lo, hi = A(r + 1), A(r)
        return (lo + eps * (hi - lo), hi - eps * (hi - lo))
    M.th = pe.Var(rootidx, bounds=thb)
    for (c, r) in rootidx:
        if c == full:
            M.th[c, r].fix(th_feed[r - 1])
        else:
            M.th[c, r].value = rng.uniform(*thb(None, c[0], c[1], r))
    def dget(c, p):
        k = cfg['split'][c][0]
        return M.d[c, p] if c[0] <= p <= k else 0
    def bget(c, p):
        l = cfg['split'][c][1]
        return M.b[c, p] if l <= p <= c[1] else 0
    def F(c, p):
        if c == full:
            return z[p - 1]
        e = 0
        for kind, P in cfg['prod'][c]:
            e = e + (dget(P, p) if kind == 'top' else bget(P, p))
        return e
    def U(c):
        if c == full:
            return 1 - q
        e = 0
        for kind, P in cfg['prod'][c]:
            if kind == 'top':
                e = e + (M.VT[P] if cfg['hx'].get(c, 'TC') == 'TC' else sum(dget(P, p) for p in range(P[0], P[1] + 1)))
            else:
                e = e - (M.VB[P] if cfg['hx'].get(c, 'TC') == 'TC' else 0)
        return e
    M.cons = pe.ConstraintList()
    for c in cols:
        i, m = c; k, l = cfg['split'][c]
        for p in range(i, m + 1):
            M.cons.add(dget(c, p) + bget(c, p) == F(c, p))
        M.cons.add(M.VT[c] - M.VB[c] == U(c))
        for r in range(l - 1, k + 1):
            t = M.th[c, r]
            M.cons.add(sum(A(p) * F(c, p) / (A(p) - t) for p in range(i, m + 1)) == U(c))
            gT = sum(A(p) * M.d[c, p] / (A(p) - t) for p in range(i, k + 1))
            if l <= r <= k - 1:
                M.cons.add(gT == M.Y[c])
            else:
                M.cons.add(gT <= M.Y[c] - dm)
        M.cons.add(M.Y[c] <= M.VT[c] - (dm if not (strict and l <= k - 1) else 0))
        M.cons.add(M.Y[c] - U(c) >= dm)
        if strict and l <= k - 1:
            M.cons.add(M.Y[c] == M.VT[c])
    # pure products produced twice: liquid side draw, vapour passes through
    reb = []
    for p in range(1, N + 1):
        pr = cfg['prod'][(p, p)]
        if len(pr) == 2:
            top = [P for kind, P in pr if kind == 'top'][0]; bot = [P for kind, P in pr if kind == 'bot'][0]
            M.cons.add(M.VT[top] == M.VB[bot])
        elif pr[0][0] == 'bot':
            reb.append(M.VB[pr[0][1]])
    for s, t in cfg['hx'].items():
        if t == 'HX' and cfg['prod'][s][0][0] == 'bot':
            reb.append(M.VB[cfg['prod'][s][0][1]])
    M.obj = pe.Objective(expr=sum(reb))
    # random initial point
    for v in M.d.values(): v.value = rng.uniform(0, 0.5)
    for v in M.b.values(): v.value = rng.uniform(0, 0.5)
    for c in cols:
        M.VT[c].value = rng.uniform(0.3, 2.5) * Vref; M.VB[c].value = rng.uniform(0.3, 2.5) * Vref
        M.Y[c].value = 0.8 * M.VT[c].value
    return M, F, U, dget, bget, reb

def feas_check(cfg, M, F, U, dget, bget, alpha, strict):
    N = cfg['N']; A = lambda p: alpha[p - 1]; v = pe.value; worst = 0.0
    for c in cfg['streams']:
        i, m = c; k, l = cfg['split'][c]
        for p in range(i, m + 1):
            worst = max(worst, abs(v(dget(c, p) + bget(c, p) - F(c, p))))
            if p <= k: worst = max(worst, -v(M.d[c, p]))
            if p >= l: worst = max(worst, -v(M.b[c, p]))
        worst = max(worst, abs(v(M.VT[c] - M.VB[c] - U(c))), -v(M.VT[c]), -v(M.VB[c]))
        for r in range(l - 1, k + 1):
            t = v(M.th[c, r])
            if not (A(r + 1) < t < A(r)): worst = max(worst, 1.0)
            fe = sum(A(p) * v(F(c, p)) / (A(p) - t) for p in range(i, m + 1)) - v(U(c))
            worst = max(worst, abs(fe) / (1 + abs(v(U(c)))))
            gT = sum(A(p) * v(M.d[c, p]) / (A(p) - t) for p in range(i, k + 1))
            gap = gT - v(M.Y[c])
            worst = max(worst, abs(gap) if l <= r <= k - 1 else gap)
        worst = max(worst, v(M.Y[c]) - v(M.VT[c]), -(v(M.Y[c]) - v(U(c))))
        if strict and l <= k - 1:
            worst = max(worst, abs(v(M.Y[c]) - v(M.VT[c])))
    return worst

def cut_identity_residual(cfg, M, alpha, z, q, th_feed, F, U):
    """Independent float evaluation, at the returned point, of the light-cut identity for every j."""
    N = cfg['N']; A = lambda p: alpha[p - 1]; v = pe.value
    R = sum(v(x) for x in reb_terms(cfg, M))
    worst = 0.0
    for j in range(1, N):
        th = th_feed[j - 1]
        H = sum(A(p) * z[p - 1] / (th - A(p)) for p in range(j + 1, N + 1))
        rhs = 0.0
        for c in cfg['streams']:
            i, m = c; k, l = cfg['split'][c]
            mixed = i <= j < m
            if mixed and k <= j:      # light exit: Lam of the top section
                rhs += v(M.VT[c]) - sum(A(p) * v(M.d[c, p]) / (A(p) - th) for p in range(i, k + 1))
        # heat duties: condensers on HX distillate links (duty L_T of producer) and product condensers,
        # reboilers on HX residue links and product reboilers; classify node as light/mixed/heavy
        for s, t in cfg['hx'].items():
            kind, P = cfg['prod'][s][0]
            if t != 'HX': continue
            if kind == 'top':     # condenser at node s, duty L_T(P)
                duty = v(M.VT[P]) - sum(v(M.d[P, p]) for p in range(P[0], cfg['split'][P][0] + 1))
                if s[0] <= j < s[1] or s[0] >= j + 1: rhs += duty      # Q(mixed) + Q(heavy)
            else:                 # reboiler at node s, duty V_B(P)
                if s[1] <= j: rhs += v(M.VB[P])                     # R(light)
        for p in range(1, N + 1):
            pr = cfg['prod'][(p, p)]
            if len(pr) == 1:
                kind, P = pr[0]
                if kind == 'top' and p >= j + 1: rhs += v(M.VT[P])   # total condenser on a heavy product
                if kind == 'bot' and p <= j: rhs += v(M.VB[P])       # reboiler on a light product
        worst = max(worst, abs((R - H) - rhs))
    return worst

def reb_terms(cfg, M):
    N = cfg['N']; out = []
    for p in range(1, N + 1):
        pr = cfg['prod'][(p, p)]
        if len(pr) == 1 and pr[0][0] == 'bot':
            out.append(M.VB[pr[0][1]])
    for s, t in cfg['hx'].items():
        if t == 'HX' and cfg['prod'][s][0][0] == 'bot':
            out.append(M.VB[cfg['prod'][s][0][1]])
    return out

CASES = {
    # (alpha, z, q)
    'T61': ([4, 2, 1], [1 / 3, 1 / 3, 1 / 3], 1.0),                       # Halvorsen Table 6.1
    'P2':  ([47.5774, 25.5263, 15.1071, 1], [0.1358, 0.5169, 0.0923, 0.255], 1.0),   # paper's example, H2 top
    'E8':  ([8, 4, 2, 1], [0.25, 0.25, 0.25, 0.25], 1.0),
    'Q04': ([6, 3, 1.6, 1], [0.15, 0.35, 0.35, 0.15], 0.4),
    'Q13': ([5, 2.2, 1.5, 1], [0.2, 0.3, 0.3, 0.2], 1.3),
    'QN':  ([7, 4, 1.4, 1], [0.1, 0.4, 0.4, 0.1], -0.3),
    # middle peak H_2 dominant, q != 1 (the non-extreme case that matters for family K)
    'M05': ([10, 5, 4.2, 1], [0.15, 0.35, 0.35, 0.15], 0.5),
    'N13': ([14.5262, 5.957, 3.327, 1], [0.1368, 0.2208, 0.3374, 0.305], 1.3),   # H2 within 1% of H3
    'N03': ([5.9208, 4.3871, 3.4747, 1], [0.4139, 0.2266, 0.1566, 0.2029], 0.3),  # H2 within 1% of H1
}

def run(N, case, nstart, model, seed=1):
    alpha, z, q = CASES[case]
    assert len(alpha) == N
    th, H = peaks(alpha, z, q)
    V = max(H)
    strict = (model == 'strict')
    cfgs = configurations(N)
    solver = pe.SolverFactory('ipopt')
    solver.options['tol'] = 1e-10; solver.options['max_iter'] = 3000; solver.options['print_level'] = 0
    rng = random.Random(seed)
    rows = []
    for idx, cfg in enumerate(cfgs):
        best = None; best_id = None
        for s in range(nstart):
            M, F, U, dget, bget, reb = build(cfg, alpha, z, q, strict, rng, th, V)
            try:
                solver.solve(M, tee=False, load_solutions=True)
            except Exception:
                continue
            w = feas_check(cfg, M, F, U, dget, bget, alpha, strict)
            if w > 1e-7:
                continue
            R = pe.value(M.obj)
            cid = cut_identity_residual(cfg, M, alpha, z, q, th, F, U)
            if best is None or R < best:
                best, best_id = R, cid
        rows.append({'idx': idx, 'name': name(cfg), 'R': best, 'ratio': (best / V) if best else None,
                     'cut_identity_residual': best_id})
    return {'N': N, 'case': case, 'model': model, 'alpha': alpha, 'z': z, 'q': q, 'H': H, 'V_FTC': V,
            'rows': rows}

if __name__ == '__main__':
    if sys.argv[1] == 'enum':
        for N in (3, 4, 5):
            C = configurations(N)
            print(N, 'basic', len(set(c['basic'] for c in C)), 'total', len(C))
        sys.exit(0)
    N, case, nstart, model = int(sys.argv[2]), sys.argv[3], int(sys.argv[4]), sys.argv[5]
    out = run(N, case, nstart, model)
    rows = [r for r in out['rows'] if r['R'] is not None]
    mn = min(rows, key=lambda r: r['R'])
    ties = [r for r in rows if r['ratio'] < 1 + 1e-6]
    below = [r for r in rows if r['ratio'] < 1 - 1e-7]
    print(f"N={N} case={case} model={model} H={['%.6f' % h for h in out['H']]} V_FTC={out['V_FTC']:.7f}")
    print(f"  solved {len(rows)}/{len(out['rows'])} configurations; min ratio {mn['ratio']:.10f} ({mn['name']})")
    print(f"  ties (ratio < 1+1e-6): {len(ties)}; below 1-1e-7: {len(below)}; "
          f"max cut-identity residual {max(r['cut_identity_residual'] for r in rows):.2e}")
    for r in sorted(ties, key=lambda r: r['name']):
        print('    tie:', f"{r['ratio']:.10f}", r['name'])
    json.dump(out, open(f'logs/nlp_{N}_{case}_{model}.json', 'w'), indent=1)
