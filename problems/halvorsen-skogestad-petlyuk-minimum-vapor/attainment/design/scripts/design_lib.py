#!/usr/bin/env python3
"""Agent D (design): the operating specification F(delta, eta, pi) of the generalized Petlyuk arrangement and its
infinite-stage limit configuration, in high precision (mpmath), with numerical checks of hypotheses (H1)-(H8) of
attainment/design/hypotheses.md.  Library module (used by check_design.py, export_config.py, crosscheck_N3.py).

Conventions (0-based internally): components 0..N-1 (alpha decreasing), columns (i, m) with 0 <= i < m <= N-1,
column (i, m) has top section ('T', i, m) and bottom section ('B', i, m); nodes ('J', i, m) for i < m (junction /
feed stage of the column) and ('P', p) for products (('P', 0) total condenser, ('P', N-1) partial reboiler, other
side-draw stages).  Top end of ('T', i, m) = ('J', i, m-1) or ('P', i) if m-1 == i; bottom end of ('B', i, m) =
('J', i+1, m) or ('P', m) if i+1 == m.  Net flow w of a section is upward (w = d in T, w = -b in B).

Environment: Python 3, mpmath (WSL Fedora sci environment).  Every function is deterministic.
Numerical evidence only: the proofs are in attainment/design/design.md.
"""
import mpmath as mp


class InvalidSpec(Exception):
    """A non-binary column violates the explicit validity condition (deformed split not positive)."""

# ----------------------------------------------------------------------------------------------------------------
# basic algebra


def bisect_root(f, lo, hi, it=None):
    """Root of a continuous f on the OPEN interval (lo, hi) with a sign change; f is never evaluated at lo, hi."""
    if it is None:
        it = int(mp.mp.prec * 1.2) + 20
    a, b = mp.mpf(lo), mp.mpf(hi)
    fa_sign = None
    # find interior points with opposite signs: f -> +-inf at poles, so test near the ends
    for k in range(1, 400):
        eps = (b - a) * mp.mpf(2) ** (-k)
        fl, fh = f(a + eps), f(b - eps)
        if fl * fh < 0:
            a, b = a + eps, b - eps
            fa_sign = mp.sign(fl)
            break
    if fa_sign is None:
        raise ValueError("no sign change found")
    for _ in range(it):
        c = (a + b) / 2
        fc = f(c)
        if fc == 0:
            return c
        if mp.sign(fc) == fa_sign:
            a = c
        else:
            b = c
    return (a + b) / 2


def feed_roots(al, f, U, idx):
    """Roots of Phi(t) = U - sum_{p in idx} al_p f_p/(al_p - t), one in each (al[idx[k+1]], al[idx[k]])."""
    def Phi(t):
        return U - mp.fsum(al[p] * f[p] / (al[p] - t) for p in idx)
    return [bisect_root(Phi, al[idx[k + 1]], al[idx[k]]) for k in range(len(idx) - 1)]


def Psi(sec, al, t, skip=None):
    """Section flux function Psi_sigma(t) = V - sum al_i w_i/(al_i - t) (optionally without component skip)."""
    return sec['V'] - mp.fsum(al[i] * sec['w'][i] / (al[i] - t) for i in range(len(al))
                              if sec['w'][i] != 0 and i != skip)


def G(x, al, t, skip=None):
    return mp.fsum(al[i] * x[i] / (al[i] - t) for i in range(len(al)) if x[i] != 0 and i != skip)


def section_modes(sec, al):
    """Modes of M = V diag(al) - w al^T divided by V: major roots of Psi (one-signed w assumed and checked) and
    trace modes al_k (w_k == 0).  Returns list of dicts sorted by decreasing value."""
    N = len(al)
    sup = [i for i in range(N) if sec['w'][i] != 0]
    signs = set(mp.sign(sec['w'][i]) for i in sup)
    assert len(signs) <= 1, "section net flow not one-signed"
    poles = sorted([al[i] for i in sup], reverse=True)
    f = lambda t: Psi(sec, al, t)
    roots = []
    for k in range(len(poles) - 1):
        roots.append(bisect_root(f, poles[k + 1], poles[k]))
    if sec['kind'] == 'T':
        roots.append(bisect_root(f, mp.mpf(0), poles[-1]))    # Psi_T(0) = L > 0, -> -inf at min pole
    else:
        hi = poles[0] + 2 * mp.fsum(abs(al[i] * sec['w'][i]) for i in sup) / sec['V'] + 1
        roots.append(bisect_root(f, poles[0], hi))
    modes = [{'type': 'major', 'val': r} for r in roots]
    modes += [{'type': 'trace', 'k': i, 'val': al[i]} for i in range(N) if sec['w'][i] == 0]
    modes.sort(key=lambda m: -m['val'])
    for m in modes:
        phi = m['val']
        if m['type'] == 'major':
            e = [sec['w'][i] / (al[i] - phi) for i in range(N)]
            m['left'] = [al[i] / (al[i] - phi) for i in range(N)]
        else:
            k = m['k']
            e = [(sec['w'][i] / (al[i] - phi) if i != k else Psi(sec, al, phi, skip=k) / al[k]) for i in range(N)]
            m['left'] = [mp.mpf(1) if i == k else mp.mpf(0) for i in range(N)]
        s = mp.fsum(e)                       # = L/phi (proved), positive
        m['right'] = [v / s for v in e]      # normalized pinch composition
        m['norm'] = s
    return modes


def lemma1_vector(al, idx, chi):
    """Astra's Lemma 1: unique normalized X on components idx with G_X(chi_j) = 0, chi interlacing."""
    om = []
    for p in idx:
        num = mp.fprod([al[p] - c for c in chi])
        den = al[p] * mp.fprod([al[p] - al[pp] for pp in idx if pp != p])
        om.append(num / den)
    s = mp.fsum(om)
    return {p: o / s for p, o in zip(idx, om)}


# ----------------------------------------------------------------------------------------------------------------
# the specification F(delta, eta, pi)


def global_data(al, z, q):
    N = len(al)
    th = feed_roots(al, z, 1 - q, list(range(N)))
    H = [mp.fsum(al[p] * z[p] / (th[j] - al[p]) for p in range(j + 1, N)) for j in range(N - 1)]
    return th, H


def preferred_split(al, f, U, i, m):
    """Theorem achieve relative to the actual feed (f, U) of column (i, m), m - i >= 2: returns roots, VT0, d0."""
    idx = list(range(i, m + 1))
    th = feed_roots(al, f, U, idx)

    def coef(p):
        num = mp.fprod([al[p] - t for t in th])
        den = mp.fprod([al[p] - al[pp] for pp in range(i, m) if pp != p])
        return num / den / al[p]
    VT0 = f[i] / coef(i)
    d0 = {p: VT0 * coef(p) for p in range(i, m)}
    return th, VT0, d0


def astra_h(al, th, i, m):
    """Astra's redistribution h (attainment.md (3)-(4)) for column (i, m) with local feed roots th (length m-i)."""
    k = m - i + 1
    h = {p: mp.mpf(0) for p in range(i, m + 1)}
    if k == 3:
        return h, None
    mid = th[1:-1]                                       # theta_2..theta_{k-2} (local)
    inner = list(range(i + 1, m))                        # a_2..a_{k-1}
    beta = mp.fprod([al[p] for p in inner]) / mp.fprod(mid)
    for p in inner:
        num = (al[p] - beta) * mp.fprod([al[p] - t for t in mid])
        den = al[p] * mp.fprod([al[p] - al[pp] for pp in inner if pp != p])
        h[p] = num / den
    return h, beta


def build_spec(al, z, q, delta, eta_nb, eta_bin, pattern, R_override=None):
    """Build F(delta, eta, pi).  eta_nb: dict (i,m)->eta or a scalar; eta_bin: dict p->eta or scalar;
    pattern: string of 'B'/'T' of length N-1.  Returns a dict with columns, sections, products, R, global data."""
    N = len(al)
    al = [mp.mpf(a) for a in al]
    z = [mp.mpf(x) for x in z]
    q = mp.mpf(q)
    th, H = global_data(al, z, q)
    Hs = max(H)
    R = Hs + mp.mpf(delta) if R_override is None else mp.mpf(R_override)
    cols = {}
    order = sorted([(i, m) for i in range(N) for m in range(i + 1, N)], key=lambda c: -(c[1] - c[0]))
    zero = lambda: [mp.mpf(0)] * N
    for (i, m) in order:
        if (i, m) == (0, N - 1):
            f, U = list(z), 1 - q
        else:
            f, U = zero(), mp.mpf(0)
            if m + 1 <= N - 1:
                par = cols[(i, m + 1)]
                f = [a + b for a, b in zip(f, par['d'])]
                U += par['VT']
            if i - 1 >= 0:
                par = cols[(i - 1, m)]
                f = [a + b for a, b in zip(f, par['b'])]
                U -= par['VB']
        col = {'i': i, 'm': m, 'f': f, 'U': U}
        if m - i >= 2:
            eta = mp.mpf(eta_nb[(i, m)] if isinstance(eta_nb, dict) else eta_nb)
            thc, VT0, d0 = preferred_split(al, f, U, i, m)
            h, beta = astra_h(al, thc, i, m)
            d = zero()
            for p in range(i, m):
                d[p] = d0[p] + eta * h[p]
            d[i] = f[i]                                   # exact: h_i = 0 and d0_i = f_i (sharp light end)
            b = [f[p] - d[p] for p in range(N)]
            b[i] = mp.mpf(0)                              # exact zero (trace flow of the limit)
            VT = VT0 + eta
            col.update(kind='nb', eta=eta, roots_feed=thc, VT0=VT0, d0=d0, h=h, beta=beta, d=d, b=b, VT=VT,
                       VB=VT - U)
            if not (all(d[p] > 0 for p in range(i, m)) and all(b[p] > 0 for p in range(i + 1, m + 1))):
                raise InvalidSpec(f"validity condition fails at column {(i, m)} (eta = {mp.nstr(eta, 5)})")
        else:
            col.update(kind='bin')
        cols[(i, m)] = col
    # binary columns: vapors from the reboiler upward, then forcing
    VB = R
    for p in range(N - 2, -1, -1):
        c = cols[(p, p + 1)]
        c['VB'] = VB
        c['VT'] = VB + c['U']
        VB = c['VT']                                    # side draw p: V_B[p-1,p] = V_T[p,p+1]
    for p in range(N - 1):
        c = cols[(p, p + 1)]
        eta = mp.mpf(eta_bin[p] if isinstance(eta_bin, dict) else eta_bin)
        c['eta'], c['pi'] = eta, pattern[p]
        fp, fq = c['f'][p], c['f'][p + 1]
        thb = feed_roots(al, c['f'], c['U'], [p, p + 1])[0]
        c['roots_feed'] = [thb]
        c['kappa'] = c['VT'] - al[p] * fp / (al[p] - thb)          # binary margin at eta_p = 0 (sharp)
        d, b = zero(), zero()
        if pattern[p] == 'B':
            d[p] = fp - eta
            b[p], b[p + 1] = eta, fq
        else:
            d[p], d[p + 1] = fp, eta
            b[p + 1] = fq - eta
        c['d'], c['b'] = d, b
    # phase flows and sections
    secs = {}
    for (i, m), c in cols.items():
        D, B = mp.fsum(c['d']), mp.fsum(c['b'])
        c['D'], c['B'] = D, B
        c['LT'], c['LB'] = c['VT'] - D, c['VB'] + B
        top = ('J', i, m - 1) if m - 1 > i else ('P', i)
        bot = ('J', i + 1, m) if i + 1 < m else ('P', m)
        secs[('T', i, m)] = {'kind': 'T', 'col': (i, m), 'V': c['VT'], 'L': c['LT'], 'w': list(c['d']),
                             'bottom': ('J', i, m), 'top': top}
        secs[('B', i, m)] = {'kind': 'B', 'col': (i, m), 'V': c['VB'], 'L': c['LB'], 'w': [-v for v in c['b']],
                             'bottom': bot, 'top': ('J', i, m)}
    prods = {}
    for p in range(N):
        n = zero()
        if p - 1 >= 0:
            n = [a + b for a, b in zip(n, cols[(p - 1, p)]['b'])]
        if p + 1 <= N - 1:
            n = [a + b for a, b in zip(n, cols[(p, p + 1)]['d'])]
        P = mp.fsum(n)
        prods[p] = {'P': P, 'n': n, 'x': [v / P for v in n]}
    return {'N': N, 'al': al, 'z': z, 'q': q, 'th': th, 'H': H, 'Hstar': Hs, 'delta': mp.mpf(delta), 'R': R,
            'cols': cols, 'secs': secs, 'prods': prods, 'pattern': pattern}


# ----------------------------------------------------------------------------------------------------------------
# the limit configuration


def limit_config(S):
    al, N, cols, secs = S['al'], S['N'], S['cols'], S['secs']
    for key, sec in secs.items():
        sec['modes'] = section_modes(sec, al)
        majors = [mo for mo in sec['modes'] if mo['type'] == 'major']
        sec['pinch'] = majors[-1] if sec['kind'] == 'T' else majors[0]   # T: smallest root; B: largest
    X = {}
    for (i, m), c in cols.items():
        T, B = secs[('T', i, m)], secs[('B', i, m)]
        if c['kind'] == 'nb':
            r = [mo['val'] for mo in T['modes'] if mo['type'] == 'major'][0]
            s = [mo['val'] for mo in B['modes'] if mo['type'] == 'major'][-1]
            chi = [r] + list(c['roots_feed'][1:-1]) + [s]
            v = lemma1_vector(al, list(range(i, m + 1)), chi)
            x = [v.get(p, mp.mpf(0)) for p in range(N)]
            c['chi'] = chi
        else:
            pin = B['pinch'] if c['pi'] == 'B' else T['pinch']
            x = list(pin['right'])
        X[('J', i, m)] = x
    for p, pr in S['prods'].items():
        X[('P', p)] = pr['x']
    S['X'] = X
    return S


def coefficients(sec, x):
    """Coefficients c_r = (l_r . x)/(l_r . e_r) of x in the (normalized) right eigenbasis of the section."""
    out = []
    for mo in sec['modes']:
        num = mp.fsum(l * v for l, v in zip(mo['left'], x))
        den = mp.fsum(l * v for l, v in zip(mo['left'], mo['right']))
        out.append(num / den)
    return out


def matM(sec, al):
    N = len(al)
    return mp.matrix([[sec['V'] * al[i] * (1 if i == j else 0) - sec['w'][i] * al[j] for j in range(N)]
                      for i in range(N)])
