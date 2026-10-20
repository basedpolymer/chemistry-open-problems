"""Independent checks written by the verifier (no code shared with the solvers' scripts).

(a) FTC value max_j H_j for the eight five-component cases of Tumbalam Gooty, Agrawal & Tawarmalani
    (arXiv:2010.12113, Table "gap comparison"), compared with their certified optima.
(b) Closed forms for Halvorsen's ternary Table 6.1 (alpha = [4,2,1], z = 1/3, q = 1): direct split,
    indirect split, Petlyuk.
(c) Exact (Fraction) check of the node-balance form of the cut identity on random networks of
    two-section columns (several columns per label, merges, arbitrary heaters/coolers at every node,
    rational evaluation point x, random node subsets U).
(d) Independent construction of the FTC attainment point (preferred split everywhere, excess only in
    binary columns), 60-digit mpmath, N = 3..7, random alpha, z, q in [-1, 2]; checks every strict-model
    constraint, positivity, side-draw vapour equalities and R = max_j H_j.

Environment: WSL Fedora, Python 3 with mpmath.  Runtime: about 20 s.  Prints PASS/FAIL, exit code 1 on failure.
"""
import random, sys
from fractions import Fraction as Fr
import mpmath as mp

ok_all = True
def report(name, ok):
    global ok_all
    ok_all &= bool(ok)
    print(('PASS' if ok else 'FAIL') + ': ' + name)

# ---------------------------------------------------------------- feed roots and peaks
def bisect_decreasing(f, lo, hi):
    """Zero of a function decreasing from +inf (lo^+) to -inf (hi^-), bisection to full working precision."""
    a, b = mp.mpf(lo), mp.mpf(hi)
    for _ in range(int(mp.mp.prec) + 20):
        c = (a + b) / 2
        if c == a or c == b:
            break
        if f(c) > 0: a = c
        else: b = c
    return (a + b) / 2
def feed_roots(alpha, z, q):
    N = len(alpha)
    f = lambda t: sum(a * zz / (a - t) for a, zz in zip(alpha, z)) - (1 - q)
    roots = []
    for r in range(N - 1):
        lo, hi = alpha[r + 1], alpha[r]
        # f -> -inf at hi^-, +inf at lo^+ ; strictly decreasing
        roots.append(bisect_decreasing(lambda t: -f(t), lo, hi))   # f increases from -inf to +inf
    return roots

def peaks_H(alpha, z, q):
    th = feed_roots(alpha, z, q)
    N = len(alpha)
    H = [sum(alpha[p] * z[p] / (th[j] - alpha[p]) for p in range(j + 1, N)) for j in range(N - 1)]
    P = [sum(alpha[p] * z[p] / (alpha[p] - th[j]) for p in range(j + 1)) for j in range(N - 1)]
    return th, H, P

# ---------------------------------------------------------------- (a) TG cases
mp.mp.dps = 40
tg = {
    'A': ([20, 30, 20, 20, 10], [4.1, 3.6, 2.1, 1.42, 1], 402.7),
    'B': ([5, 15, 25, 20, 35], [7.98, 3.99, 3, 1.25, 1], 272.5),
    'C': ([25, 10, 25, 20, 20], [13.72, 3.92, 3.267, 1.21, 1], 260.0),
    'D': ([42.5, 42.5, 5, 5, 5], [3.3275, 3.025, 1.21, 1.1, 1], 896.4),   # TG print 4 flows; 5th = 5 (sum 100)
    'E': ([30, 30, 5, 5, 30], [1.4641, 1.331, 1.21, 1.1, 1], 695.6),
    'F': ([5, 5, 5, 42.5, 42.5], [3.3275, 1.331, 1.21, 1.1, 1], 929.1),
    'G': ([5, 5, 5, 42.5, 42.5], [17.1875, 6.875, 2.75, 1.1, 1], 902.7),
    'H': ([20, 20, 20, 20, 20], [7.5625, 3.025, 1.21, 1.1, 1], 542.0),
}
okA = True
for name, (F, al, opt) in tg.items():
    al = [mp.mpf(str(a)) for a in al]; F = [mp.mpf(str(x)) for x in F]
    th, H, P = peaks_H(al, F, mp.mpf(1))       # q = 1 (saturated liquid), total feed 100
    v = max(H)
    good = abs(v - opt) <= 0.051 + 1e-9 * opt   # TG print one decimal (Case C: '260')
    okA &= good
    print(f'   TG case {name}: max_j H_j = {mp.nstr(v, 7)}  (TG optimum {opt})  peaks {[mp.nstr(h, 6) for h in H]}')
report('(a) FTC value equals the published optimum of the eight TG five-component cases (to print accuracy)', okA)

# ---------------------------------------------------------------- (b) Halvorsen Table 6.1 closed forms
al = [mp.mpf(4), mp.mpf(2), mp.mpf(1)]; z = [mp.mpf(1) / 3] * 3
th, H, P = peaks_H(al, z, mp.mpf(1))
pet = max(H)
def binary_vmin(a1, a2, z1, z2, vapfrac):  # binary feed z1+z2 with vapour fraction vapfrac: returns (V_T, V_B)
    U = vapfrac * (z1 + z2)
    t = bisect_decreasing(lambda t: U - a1 * z1 / (a1 - t) - a2 * z2 / (a2 - t), a2, a1)
    VT = a1 * z1 / (a1 - t)
    return VT, VT - U
# Halvorsen's conventional sequences use TG's conventions: a reboiled residue leaves as saturated liquid,
# a condensed distillate as saturated vapour (partial condenser).
ds = P[0] + binary_vmin(al[1], al[2], z[1], z[2], 0)[1]     # C1 A/BC (V_B = V_T at q=1), C2 B/C, liquid feed
is_ = P[1] + binary_vmin(al[0], al[1], z[0], z[1], 1)[1]    # C1 AB/C, C2 A/B, saturated-vapour feed
is_liq = P[1] + binary_vmin(al[0], al[1], z[0], z[1], 0)[1] # same with a total condenser (liquid transfer)
print(f'   (indirect split with a total condenser on AB instead: {mp.nstr(is_liq, 6)})')
print(f'   Petlyuk {mp.nstr(pet, 6)} (Table 6.1: 1.366), direct split {mp.nstr(ds, 6)} (2.072), '
      f'indirect split {mp.nstr(is_, 6)} (2.032)')
report('(b) Halvorsen Table 6.1: Petlyuk, direct and indirect split values',
       abs(pet - 1.366) < 6e-4 and abs(ds - 2.072) < 6e-4 and abs(is_ - 2.032) < 6e-4)

# ---------------------------------------------------------------- (c) exact cut identity on random networks
def random_network(rng, N):
    """Columns: dict id -> {'lab':(i,m), 'k':k, 'l':l, 'dist':node, 'res':node}; nodes are ('c',id) or ('p',p)."""
    cols = {}
    bylabel = {}
    def new_col(lab):
        cid = len(cols)
        cols[cid] = {'lab': lab}
        bylabel.setdefault(lab, []).append(cid)
        return cid
    root = new_col((1, N))
    queue = [root]
    while queue:
        c = queue.pop(0)
        i, m = cols[c]['lab']
        k = rng.randint(i, m - 1)
        l = rng.randint(i + 1, min(k + 1, m))
        cols[c]['k'], cols[c]['l'] = k, l
        for key, lab in (('dist', (i, k)), ('res', (l, m))):
            if lab[0] == lab[1]:
                cols[c][key] = ('p', lab[0])
            else:
                if bylabel.get(lab) and rng.random() < 0.6:
                    tgt = rng.choice(bylabel[lab])
                else:
                    tgt = new_col(lab); queue.append(tgt)
                cols[c][key] = ('c', tgt)
    return cols

def check_identity_once(rng):
    N = rng.randint(3, 6)
    alpha = sorted({Fr(rng.randint(2, 60), rng.randint(1, 5)) for _ in range(3 * N)}, reverse=True)[:N]
    if len(alpha) < N:
        return True
    z = [Fr(rng.randint(1, 20), 7) for _ in range(N)]
    q = Fr(rng.randint(-10, 20), 10)
    cols = random_network(rng, N)
    nodes = [('c', c) for c in cols] + [('p', p) for p in range(1, N + 1)]
    inflow = {n: {'d': [], 'b': []} for n in nodes}
    for c, C in cols.items():
        inflow[C['dist']]['d'].append(c)
        inflow[C['res']]['b'].append(c)
    order = sorted(cols, key=lambda c: -(cols[c]['lab'][1] - cols[c]['lab'][0]))
    d, b, VT, VB, Qr, Qc = {}, {}, {}, {}, {}, {}
    for c in order:
        C = cols[c]; i, m = C['lab']; k, l = C['k'], C['l']
        if c == 0:
            f = {p: z[p - 1] for p in range(1, N + 1)}
        else:
            f = {p: Fr(0) for p in range(i, m + 1)}
            for P in inflow[('c', c)]['d']:
                for p, v in d[P].items(): f[p] += v
            for P in inflow[('c', c)]['b']:
                for p, v in b[P].items(): f[p] += v
        d[c] = {}; b[c] = {}
        for p in range(i, m + 1):
            if p < l:   dp = f[p]
            elif p > k: dp = Fr(0)
            else:       dp = f[p] * Fr(rng.randint(0, 10), 10)
            if p <= k: d[c][p] = dp
            if p >= l: b[c][p] = f[p] - dp
        Qr[('c', c)] = Fr(rng.randint(0, 5), 3) if rng.random() < 0.5 else Fr(0)
        Qc[('c', c)] = Fr(rng.randint(0, 5), 3) if rng.random() < 0.5 else Fr(0)
        U = (1 - q) if c == 0 else Fr(0)
        U += sum(VT[P] for P in inflow[('c', c)]['d']) - sum(VB[P] for P in inflow[('c', c)]['b'])
        U += Qr[('c', c)] - Qc[('c', c)]
        VT[c] = Fr(rng.randint(0, 40), 4)
        VB[c] = VT[c] - U          # any sign: the identity is algebraic
    for p in range(1, N + 1):
        n = ('p', p)
        bal = sum(VT[P] for P in inflow[n]['d']) - sum(VB[P] for P in inflow[n]['b'])
        r0 = Fr(rng.randint(0, 3), 2)
        Qr[n] = r0 + max(Fr(0), -bal); Qc[n] = r0 + max(Fr(0), bal)
    # sections: (bottom node, top node, V, w)
    secs = []
    for c, C in cols.items():
        secs.append((('c', c), C['dist'], VT[c], d[c]))                              # top section
        secs.append((C['res'], ('c', c), VB[c], {p: -v for p, v in b[c].items()}))     # bottom section
    def Lam(sec, x):
        _, _, V, w = sec
        return V - sum(alpha[p - 1] * wp / (alpha[p - 1] - x) for p, wp in w.items())
    ok = True
    for _ in range(4):
        x = Fr(rng.randint(1, 600), rng.randint(7, 13))
        if x in alpha:
            continue
        Uset = {n for n in nodes if rng.random() < 0.5}
        Uc = [n for n in nodes if n not in Uset]
        lhs = sum(Qr[n] for n in Uc) - sum(Qc[n] for n in Uc)
        lamF = (1 - q) - sum(alpha[p] * z[p] / (alpha[p] - x) for p in range(N))
        rhs = sum(Lam(s, x) for s in secs if s[0] in Uc and s[1] in Uset) \
            - sum(Lam(s, x) for s in secs if s[0] in Uset and s[1] in Uc) \
            - (lamF if ('c', 0) in Uc else 0) \
            - sum(alpha[n[1] - 1] * z[n[1] - 1] / (alpha[n[1] - 1] - x) for n in Uc if n[0] == 'p')
        ok &= (lhs == rhs)
    # component conservation sanity: each product node receives exactly z_p
    for p in range(1, N + 1):
        n = ('p', p)
        got = sum(d[P].get(p, 0) for P in inflow[n]['d']) + sum(b[P].get(p, 0) for P in inflow[n]['b'])
        ok &= (got == z[p - 1])
    return ok

rng = random.Random(20261020)
res = [check_identity_once(rng) for _ in range(400)]
report(f'(c) node-balance cut identity, exact rationals, {len(res)} random networks x 4 random (x, U)', all(res))

# ---------------------------------------------------------------- (d) FTC attainment point
mp.mp.dps = 60
def ftc_point(alpha, z, q):
    N = len(alpha)
    th, H, P = peaks_H(alpha, z, q)
    Vmax = max(H)
    col = {}   # (i,m) 0-based inclusive -> dict
    errs = []
    for length in range(N - 1, 0, -1):
        for i in range(0, N - length):
            m = i + length
            if (i, m) == (0, N - 1):
                f = {p: z[p] for p in range(N)}; U = 1 - q
            else:
                f = {p: mp.mpf(0) for p in range(i, m + 1)}; U = mp.mpf(0)
                if (i, m + 1) in col:
                    Pc = col[(i, m + 1)]
                    for p in range(i, m + 1): f[p] += Pc['d'].get(p, 0)
                    U += Pc['VT']
                if (i - 1, m) in col:
                    Pc = col[(i - 1, m)]
                    for p in range(i, m + 1): f[p] += Pc['b'].get(p, 0)
                    U -= Pc['VB']
            # roots theta_i..theta_{m-1} must solve this column's feed equation (carry-over)
            for r in range(i, m):
                errs.append(abs(sum(alpha[p] * f[p] / (alpha[p] - th[r]) for p in f) - U))
            if length >= 2:
                # preferred split: residues of V_T prod(x-theta_r)/prod(x-alpha_p), p in [i,m-1]
                def res_top(p):
                    num = mp.fprod([alpha[p] - th[r] for r in range(i, m)])
                    den = mp.fprod([alpha[p] - alpha[pp] for pp in range(i, m) if pp != p])
                    return num / den
                VT = f[i] * alpha[i] / res_top(i)
                d = {p: VT * res_top(p) / alpha[p] for p in range(i, m)}
                b = {p: f[p] - d.get(p, 0) for p in range(i + 1, m + 1)}
                errs.append(abs(f[i] - d[i]))
                VB = VT - U
                Y = VT
            else:
                r = i
                e = Vmax - H[r]
                d = {i: f[i]}; b = {m: f[m]}
                Y = alpha[i] * f[i] / (alpha[i] - th[r])
                VT = Y + e
                VB = VT - U
                errs.append(abs((Y - U) - alpha[m] * f[m] / (th[r] - alpha[m])))   # Y_B = g_B(theta_r)
            col[(i, m)] = {'d': d, 'b': b, 'VT': VT, 'VB': VB, 'Y': Y, 'U': U, 'f': f}
    # checks
    good = True
    for (i, m), C in col.items():
        good &= all(v > 0 for v in C['d'].values()) and all(v > 0 for v in C['b'].values())
        good &= C['VT'] > 0 and C['VB'] > 0
        good &= C['VT'] - sum(C['d'].values()) > 0 and C['VB'] + sum(C['b'].values()) > 0   # L_T, L_B > 0
        good &= C['Y'] <= C['VT'] * (1 + mp.mpf('1e-40')) and C['Y'] - C['U'] >= -mp.mpf('1e-40')
        # Underwood at every root of the column's range r in [i, m-1]: g_T(theta_r) <= Y, = Y for interior
        for r in range(i, m):
            gT = sum(alpha[p] * v / (alpha[p] - th[r]) for p, v in C['d'].items())
            if i < r < m - 1:
                errs.append(abs(gT - C['Y']))
            good &= gT <= C['Y'] + mp.mpf('1e-40')
        if m - i >= 2:
            errs.append(abs(C['Y'] - C['VT']))   # strict model: no excess with interior roots
    # side draws: product p (1..N-2, 0-based) is the distillate of (p,p+1) and the residue of (p-1,p)
    for p in range(1, N - 1):
        errs.append(abs(col[(p, p + 1)]['VT'] - col[(p - 1, p)]['VB']))
    R = col[(N - 2, N - 1)]['VB']           # single reboiler (product N)
    errs.append(abs(R - Vmax))
    # products pure: component flows reaching products equal z
    return good, max(errs)

okD = True; worst = mp.mpf(0); count = 0
rng = random.Random(7)
for N in range(3, 8):
    for _ in range(40):
        al = sorted([mp.mpf(1)] + [mp.mpf(1) + mp.mpf(rng.uniform(0.05, 12)) for _ in range(N - 1)], reverse=True)
        if len(set(al)) < N: continue
        z = [mp.mpf(rng.uniform(0.02, 1)) for _ in range(N)]; s = sum(z); z = [x / s for x in z]
        q = mp.mpf(rng.uniform(-1, 2))
        g, e = ftc_point(al, z, q)
        okD &= g and e < mp.mpf('1e-40'); worst = max(worst, e); count += 1
print(f'   FTC attainment: {count} instances, max residual {mp.nstr(worst, 3)}')
report('(d) FTC point (preferred split + binary excess) is strict-feasible, positive, side draws balance, R = max H', okD)

print('ALL PASS' if ok_all else 'SOME FAIL')
sys.exit(0 if ok_all else 1)
