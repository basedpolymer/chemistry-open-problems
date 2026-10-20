#!/usr/bin/env python3
"""Agent D (design): numerical checks of hypotheses (H1)-(H8) (attainment/design/hypotheses.md) on random and
nearly critical data, N = 3..7, all forcing patterns tested at random, in high precision (mpmath).

For each instance: build the specification F(delta, eta, pi) and the limit configuration (design_lib.py), then check
(H1) flows and balances, (H2) net-flow signs/supports and the exactness of the middle roots, (H3) Astra's
conditions for every non-binary column, (H4) binary margins and pinned junctions, (H5) spectra (eigen-residuals of
M = V diag(alpha) - w alpha^T, distinct positive modes), the limit relations (closed modes vanish at the right
end, nonzero pivots), (H6) half-profiles in the simplex, (H7) trace signs and closed ends, (H8) R, purities,
pivot signs; plus kappa_p(eta_nb = 0) = e_p + delta and kappa_p(delta) - kappa_p(0) = delta, the flux
bookkeeping identity, and the stage-level bound R >= sum_P max(0, phi_P(theta_j)) at the limit configuration.
Non-binary slack: eta_nb = c * delta, halved until the explicit validity conditions and kappa_p >= delta/2 hold
(the eta-bar(delta) of Theorem 1, found numerically).  Binary forcing: a fixed fraction of the explicit bound.

Usage: python check_design.py [n_random_per_N=12] [dps=60] [seed=1]
Prints PASS/FAIL per check family and a summary; exit status 1 on any FAIL.
Environment: Python 3 + mpmath (WSL Fedora sci env).  Runtime: about 2-6 minutes for the default.
Numerical evidence only (floating point with 60 digits, no interval certificate).
"""
import sys
import random
import mpmath as mp
from design_lib import (build_spec, limit_config, Psi, G, coefficients, matM, global_data, InvalidSpec)

NR, SEED = 12, 1
mp.mp.dps = 60
TOL = mp.mpf(10) ** (-40)


def set_args(argv):
    """Parse command-line arguments (only when run as a script)."""
    global NR, SEED, TOL
    NR = int(argv[1]) if len(argv) > 1 else 12
    mp.mp.dps = int(argv[2]) if len(argv) > 2 else 60
    SEED = int(argv[3]) if len(argv) > 3 else 1
    TOL = mp.mpf(10) ** (-(mp.mp.dps * 2 // 3))


fails = {}
counts = {}


def chk(name, cond, info=""):
    counts[name] = counts.get(name, 0) + 1
    if not cond:
        fails.setdefault(name, []).append(info)
    return cond


def nodename(n):
    return f"J[{n[1]+1},{n[2]+1}]" if n[0] == 'J' else f"P[{n[1]+1}]"


def admissible_eta_bin(S, p):
    """Explicit binary bounds (design.md Prop. 6): eta_p < f_p (B) or < min(f_{p+1}, L_T^0) (T), and the trace
    ordering bound."""
    al, N, c = S['al'], S['N'], S['cols'][(p, p + 1)]
    fp, fq = c['f'][p], c['f'][p + 1]
    VT, VB = c['VT'], c['VB']
    if S['pattern'][p] == 'B':
        bd = [fp]
        if p >= 1:
            Psi_p = VB - al[p + 1] * fq / (al[p - 1] - al[p + 1])
            bd.append((al[p - 1] - al[p]) * Psi_p / al[p])
    else:
        bd = [fq, VT - fp]
        if p + 2 <= N - 1:
            Psi_q = VT - al[p] * fp / (al[p] - al[p + 2])
            bd.append((al[p + 1] - al[p + 2]) * Psi_q / al[p + 1])
    return min(bd)


def nb_valid(S):
    ok = True
    for (i, m), c in S['cols'].items():
        if c['kind'] != 'nb':
            continue
        ok &= all(c['d'][p] > 0 for p in range(i, m)) and all(c['b'][p] > 0 for p in range(i + 1, m + 1))
    return ok


def make_instance(al, z, q, drel, pattern, fr_bin=mp.mpf('0.25'), c_nb=mp.mpf('0.1'), fixed_eta=None):
    """Find eta_nb (halving) with validity and kappa_p >= delta/2, then eta_bin = fr_bin * explicit bound."""
    th, H = global_data([mp.mpf(a) for a in al], [mp.mpf(x) for x in z], mp.mpf(q))
    Hs = max(H)
    delta = drel * Hs
    eta = c_nb * delta if fixed_eta is None else mp.mpf(fixed_eta)
    for _ in range(60):
        try:
            S = build_spec(al, z, q, delta, eta, 0, pattern)
        except InvalidSpec:
            eta /= 2
            continue
        if nb_valid(S) and all(S['cols'][(p, p + 1)]['kappa'] >= delta / 2 for p in range(len(al) - 1)):
            break
        eta /= 2
    else:
        raise RuntimeError("no admissible eta_nb found")
    eb = {p: fr_bin * admissible_eta_bin(S, p) for p in range(len(al) - 1)}
    S = build_spec(al, z, q, delta, eta, eb, pattern)
    return limit_config(S), eta


def check_instance(S, tag):
    global TOL
    TOL = mp.mpf(10) ** (-(mp.mp.dps * 2 // 3))
    al, N, cols, secs, X = S['al'], S['N'], S['cols'], S['secs'], S['X']
    # ---------------- (H1) flows and balances
    for key, s in secs.items():
        chk('H1 flows', s['V'] > 0 and s['L'] > 0, f"{tag} {key}")
    for (i, m), c in cols.items():
        chk('H1 flows', c['D'] > 0 and c['B'] > 0, f"{tag} col {(i, m)}")
        chk('H1 feeds', all(c['f'][p] > 0 for p in range(i, m + 1)) and
            all(c['f'][p] == 0 for p in range(N) if p < i or p > m), f"{tag} feed {(i, m)}")
        chk('H1 balances', all(abs(c['d'][p] + c['b'][p] - c['f'][p]) < TOL for p in range(N)) and
            abs(c['VT'] - c['VB'] - c['U']) < TOL and
            abs((c['LB'] - c['LT']) - (mp.fsum(c['f']) - c['U'])) < TOL, f"{tag} col {(i, m)}")
        # liquid into the junction from the parents equals F - U
        Lin = mp.mpf(0)
        if (i, m) == (0, N - 1):
            Lin = S['q']
        else:
            if m + 1 <= N - 1:
                Lin -= cols[(i, m + 1)]['LT']
            if i - 1 >= 0:
                Lin += cols[(i - 1, m)]['LB']
        chk('H1 balances', abs(Lin - (mp.fsum(c['f']) - c['U'])) < TOL, f"{tag} liquid {(i, m)}")
    for p, pr in S['prods'].items():
        chk('H1 flows', pr['P'] > 0, f"{tag} P{p}")
        Lin = (cols[(p - 1, p)]['LB'] if p >= 1 else 0) - (cols[(p, p + 1)]['LT'] if p <= N - 2 else 0)
        Vb = (cols[(p, p + 1)]['VT'] if p <= N - 2 else 0) - (cols[(p - 1, p)]['VB'] if p >= 1 else 0)
        if p == 0:
            Lin = cols[(0, 1)]['VT'] - cols[(0, 1)]['LT']       # total condenser
        elif p == N - 1:
            Lin = cols[(p - 1, p)]['LB'] - cols[(p - 1, p)]['VB']  # reboiler
        else:
            chk('H1 balances', abs(Vb) < TOL, f"{tag} side-draw vapor {p}")
        chk('H1 balances', abs(Lin - pr['P']) < TOL, f"{tag} product {p}")
    tot = [mp.fsum(S['prods'][p]['n'][k] for p in range(N)) for k in range(N)]
    chk('H1 balances', all(abs(tot[k] - S['z'][k]) < TOL for k in range(N)), f"{tag} global")
    chk('H8 R', abs(S['R'] - S['Hstar'] - S['delta']) < TOL, tag)
    # ---------------- (H2) signs, supports, middle roots
    for key, s in secs.items():
        sgn = 1 if s['kind'] == 'T' else -1
        chk('H2 signs', all(sgn * v >= 0 for v in s['w']), f"{tag} {key}")
    for (i, m), c in cols.items():
        if c['kind'] == 'nb':
            chk('H2 supports', all(c['d'][p] > 0 for p in range(i, m)) and c['d'][m] == 0 and
                abs(c['d'][i] - c['f'][i]) < TOL and all(c['b'][p] > 0 for p in range(i + 1, m + 1)) and
                c['b'][i] == 0, f"{tag} {(i, m)}")
            for k, t in enumerate(c['roots_feed'][1:-1]):
                chk('H2 middle roots global', abs(t - S['th'][i + 1 + k]) < TOL, f"{tag} {(i, m)} {k}")
    # ---------------- (H3) non-binary columns
    for (i, m), c in cols.items():
        if c['kind'] != 'nb':
            continue
        T, B = secs[('T', i, m)], secs[('B', i, m)]
        mT = [mo['val'] for mo in T['modes'] if mo['type'] == 'major']
        mB = [mo['val'] for mo in B['modes'] if mo['type'] == 'major']
        r, u, v, s_ = mT[0], mT[-1], mB[0], mB[-1]
        th = c['roots_feed']
        chk('H3 interlacing', al[i] > r > th[0] > v > al[i + 1] and al[m - 1] > u > th[-1] > s_ > al[m],
            f"{tag} {(i, m)}")
        chk('H3 common middle roots', all(abs(a - b) < TOL for a, b in zip(mT[1:-1], mB[1:-1])) and
            all(abs(a - b) < TOL for a, b in zip(mT[1:-1], th[1:-1])), f"{tag} {(i, m)}")
        chk('H3 Psi at extreme traces', Psi(T, al, al[m]) > 0 and Psi(B, al, al[i]) > 0, f"{tag} {(i, m)}")
        ec = [mo for mo in T['modes'] if mo['type'] == 'trace' and mo['k'] == m][0]['right']
        eu = T['pinch']['right']
        ev = B['pinch']['right']
        ea = [mo for mo in B['modes'] if mo['type'] == 'trace' and mo['k'] == i][0]['right']
        chk('H3 pinch vectors positive',
            all(ec[p] > 0 for p in range(i, m + 1)) and all(eu[p] > 0 for p in range(i, m)) and
            abs(eu[m]) < TOL and all(ev[p] > 0 for p in range(i + 1, m + 1)) and abs(ev[i]) < TOL and
            all(ea[p] > 0 for p in range(i, m + 1)), f"{tag} {(i, m)}")
        x = X[('J', i, m)]
        A_ = x[m] / ec[m]
        B_ = 1 - A_
        D_ = x[i] / ea[i]
        C_ = 1 - D_
        res1 = max(abs(x[p] - A_ * ec[p] - B_ * eu[p]) for p in range(N))
        res2 = max(abs(x[p] - C_ * ev[p] - D_ * ea[p]) for p in range(N))
        chk('H3 cone weights', A_ > 0 and B_ > 0 and C_ > 0 and D_ > 0 and res1 < TOL and res2 < TOL,
            f"{tag} {(i, m)} A,B,C,D={mp.nstr(A_,5)},{mp.nstr(B_,5)},{mp.nstr(C_,5)},{mp.nstr(D_,5)} res {mp.nstr(res1,3)} {mp.nstr(res2,3)}")
        c['ABCD'] = (A_, B_, C_, D_)
    # ---------------- (H4) binary columns
    for p in range(N - 1):
        c = cols[(p, p + 1)]
        T, B = secs[('T', p, p + 1)], secs[('B', p, p + 1)]
        x = X[('J', p, p + 1)]
        chk('H4 margin kappa>0', c['kappa'] > 0, f"{tag} {p}")
        chk('H4 X positive', x[p] > 0 and x[p + 1] > 0 and all(x[k] == 0 for k in range(N) if k not in (p, p + 1)),
            f"{tag} {p}")
        if c['pi'] == 'B':
            phiT = T['pinch']['val']
            mB = [mo['val'] for mo in B['modes'] if mo['type'] == 'major']
            phi1, phi2 = mB
            marg = Psi(B, al, phiT)
            chk('H4 B-pinned margin', marg > c['kappa'] > 0 and al[p + 1] < phi2 < phiT < al[p] < phi1 and
                (p == 0 or phi1 < al[p - 1]), f"{tag} {p}")
            ec = [mo for mo in T['modes'] if mo['type'] == 'trace' and mo['k'] == p + 1][0]['right']
            A_ = x[p + 1] / ec[p + 1]
            B_ = 1 - A_
            res = max(abs(x[k] - A_ * ec[k] - B_ * (1 if k == p else 0)) for k in range(N))
            chk('H4 decomposition', A_ > 0 and B_ > 0 and res < TOL and all(ec[k] > 0 for k in (p, p + 1)),
                f"{tag} {p}")
            c['AB'] = (A_, B_)
        else:
            phiB = B['pinch']['val']
            mT = [mo['val'] for mo in T['modes'] if mo['type'] == 'major']
            phim, phil = mT
            marg = Psi(T, al, phiB)
            chk('H4 T-pinned margin', marg > c['kappa'] > 0 and al[p + 1] < phiB < phim < al[p] and
                phil < al[p + 1] and phil > 0 and (p + 2 > N - 1 or phil > al[p + 2]), f"{tag} {p}")
            ea = [mo for mo in B['modes'] if mo['type'] == 'trace' and mo['k'] == p][0]['right']
            D_ = x[p] / ea[p]
            C_ = 1 - D_
            res = max(abs(x[k] - C_ * (1 if k == p + 1 else 0) - D_ * ea[k]) for k in range(N))
            chk('H4 decomposition', C_ > 0 and D_ > 0 and res < TOL and all(ea[k] > 0 for k in (p, p + 1)),
                f"{tag} {p}")
            c['CD'] = (C_, D_)
    # ---------------- (H5) spectra and limit relations; (H6) half-profiles; (H7) traces; (H8) pivots
    for key, s in secs.items():
        modes = s['modes']
        vals = [mo['val'] for mo in modes]
        chk('H5 distinct positive modes', all(v > 0 for v in vals) and
            all(vals[k] > vals[k + 1] for k in range(len(vals) - 1)), f"{tag} {key}")
        rho = max(vals[k + 1] / vals[k] for k in range(len(vals) - 1))
        s['rho'] = rho
        chk('H5 rho<1', rho < 1, f"{tag} {key}")
        M = matM(s, al)
        er = mp.mpf(0)
        for mo in modes:
            e = mp.matrix(mo['right'])
            l = mp.matrix(mo['left'])
            er = max(er, mp.norm(M * e - s['V'] * mo['val'] * e, 1), mp.norm(M.T * l - s['V'] * mo['val'] * l, 1))
        chk('H5 eigen residuals', er < TOL * 1e6, f"{tag} {key} {mp.nstr(er, 3)}")
        ip = modes.index(s['pinch'])
        xb, xt = X[s['bottom']], X[s['top']]
        cb, ct = coefficients(s, xb), coefficients(s, xt)
        chk('H5 limit relations', all(abs(cb[k]) < TOL for k in range(ip)) and
            all(abs(ct[k]) < TOL for k in range(ip + 1, len(modes))), f"{tag} {key}")
        chk('H8 pivots positive', cb[ip] > 0 and ct[ip] > 0, f"{tag} {key} {mp.nstr(cb[ip],5)} {mp.nstr(ct[ip],5)}")
        s['pivots'] = (cb[ip], ct[ip])
        # pivot signs as Underwood functions at product ends / child junctions
        phi = s['pinch']['val']
        if s['kind'] == 'T' and (s['top'][0] == 'P' or True):
            chk('H8 G-pivot signs', G(xt, al, phi) > 0, f"{tag} {key} top")
        if s['kind'] == 'B':
            chk('H8 G-pivot signs', G(xb, al, phi) < 0, f"{tag} {key} bottom")
        # (H6) half-profiles of the LIMIT relation, in the exact spectral form (closed modes have zero
        # coefficient; iterating M numerically would amplify round-off along the closed modes):
        #   bottom half x_k ~ sum_{r >= r*} c_r (lambda_r/lambda_*)^k e_r,  top half x_j ~ sum_{r <= r*} c_r (lambda_*/lambda_r)^j e_r
        okpos, okone = True, True
        lam = [mo['val'] for mo in modes]
        for half, cs, rng_ in (('bot', cb, range(ip, len(modes))), ('top', ct, range(0, ip + 1))):
            k = 0
            while True:
                if half == 'bot':
                    fac = {r: cs[r] * (lam[r] / lam[ip]) ** k for r in rng_}
                else:
                    fac = {r: cs[r] * (lam[ip] / lam[r]) ** k for r in rng_}
                v = [mp.fsum(fac[r] * modes[r]['right'][i] for r in rng_) for i in range(N)]
                sv = mp.fsum(v)
                okone &= sv > 0
                okpos &= all(vi / sv >= -TOL for vi in v)
                tail = max([abs(fac[r]) for r in rng_ if r != ip] + [mp.mpf(0)])
                if tail < TOL or k > 10 ** 7:
                    break
                k = k + 1 if k < 100 else 2 * k        # every k up to 100, then k = 200, 400, ...
        chk('H6 half-profiles in simplex', okpos, f"{tag} {key}")
        chk('H6 chart normalization positive', okone, f"{tag} {key}")
        # (H7)
        s['traces'] = []
        for k_, mo in enumerate(modes):
            if mo['type'] != 'trace':
                continue
            k = mo['k']
            end = 'bottom' if k_ < ip else 'top'
            J = s[end]
            xj = X[J]
            Pk = Psi(s, al, al[k], skip=k)
            Gk = G(xj, al, al[k], skip=k)
            majors = [i for i in range(N) if xj[i] != 0]
            light_ok = all(k < i for i in majors) if end == 'bottom' else all(k > i for i in majors)
            chk('H7 trace signs', xj[k] == 0 and Pk > 0 and ((Gk < 0) if end == 'bottom' else (Gk > 0)) and light_ok,
                f"{tag} {key} k={k} {end}")
            s['traces'].append({'k': k, 'end': end, 'node': J, 'Psi_k': Pk, 'G_k': Gk, 'kappa': abs(Gk) / Pk})
        # (v) components with nonzero net flow are present at both ends
        chk('H7 majors at both ends', all((xb[i] > 0 and xt[i] > 0) for i in range(N) if s['w'][i] != 0),
            f"{tag} {key}")
    # (iv) every trace node (junction) is the closed end of a section for that component
    for node, x in X.items():
        if node[0] != 'J':
            continue
        for k in range(N):
            if x[k] != 0:
                continue
            closed = [key for key, s in secs.items() for t in s['traces'] if t['k'] == k and t['node'] == node]
            chk('H7 trace nodes closed', len(closed) > 0, f"{tag} {node} {k}")
    # flux bookkeeping identity Phi_c = Psi_T(parent) - Psi_B(parent)
    t0 = mp.mpf('0.123456789')
    for (i, m), c in cols.items():
        if (i, m) == (0, N - 1):
            continue
        lhs = c['U'] - mp.fsum(al[p] * c['f'][p] / (al[p] - t0) for p in range(N))
        rhs = (Psi(secs[('T', i, m + 1)], al, t0) if m + 1 <= N - 1 else 0) - \
              (Psi(secs[('B', i - 1, m)], al, t0) if i >= 1 else 0)
        chk('flux bookkeeping', abs(lhs - rhs) < TOL, f"{tag} {(i, m)}")
    # stage-level bound at the limit configuration
    for j, th in enumerate(S['th']):
        phis = []
        for p, pr in S['prods'].items():
            phis.append(-mp.fsum(al[k] * pr['n'][k] / (al[k] - th) for k in range(N)))
        lb = mp.fsum(max(0, v) for v in phis)
        chk('stage-level bound at limit', S['R'] - lb > -TOL, f"{tag} j={j} R-lb={mp.nstr(S['R']-lb,5)}")
    # product purities
    for p, pr in S['prods'].items():
        des = mp.mpf(0)
        if p >= 1 and S['pattern'][p - 1] == 'B':
            des += cols[(p - 1, p)]['eta']
        if p <= N - 2 and S['pattern'][p] == 'T':
            des += cols[(p, p + 1)]['eta']
        chk('H8 purities', abs((1 - pr['x'][p]) - des / pr['P']) < TOL, f"{tag} P{p}")


def kappa_identities(al, z, q, drel, pattern, eta):
    """kappa_p(eta_nb=0, delta) = e_p + delta; kappa_p(eta_nb, delta) - kappa_p(eta_nb, 0) = delta."""
    global TOL
    TOL = mp.mpf(10) ** (-(mp.mp.dps * 2 // 3))
    th, H = global_data([mp.mpf(a) for a in al], [mp.mpf(x) for x in z], mp.mpf(q))
    Hs = max(H)
    delta = drel * Hs
    S0 = build_spec(al, z, q, delta, 0, 0, pattern)
    ok1 = all(abs(S0['cols'][(p, p + 1)]['kappa'] - (Hs - H[p] + delta)) < TOL for p in range(len(al) - 1))
    Sa = build_spec(al, z, q, delta, eta, 0, pattern)
    Sb = build_spec(al, z, q, 0, eta, 0, pattern)
    ok2 = all(abs(Sa['cols'][(p, p + 1)]['kappa'] - Sb['cols'][(p, p + 1)]['kappa'] - delta) < TOL
              for p in range(len(al) - 1))
    return ok1, ok2


def rand_data(rng, N, kind):
    rat = [mp.mpf(rng.uniform(1.25, 4.0)) for _ in range(N - 1)]
    if kind == 'close':
        rat[rng.randrange(N - 1)] = mp.mpf(rng.uniform(1.03, 1.1))
    al = [mp.mpf(1)]
    for r in rat:
        al.append(al[-1] * r)
    al = al[::-1]
    w = [rng.gammavariate(1.0, 1.0) for _ in range(N)]
    if kind == 'smallz':
        w[rng.randrange(N)] = 1e-3 * sum(w)
    s = sum(w)
    z = [mp.mpf(v) / s for v in w]
    z[-1] = 1 - mp.fsum(z[:-1])
    q = mp.mpf(rng.uniform(-1.0, 2.0))
    if kind == 'qlow':
        q = mp.mpf(-1)
    if kind == 'qhigh':
        q = mp.mpf(2)
    return al, z, q


def equal_peaks_q(al, z):
    """q in [-1, 2] with the two largest peaks equal (bisection on the difference), or None."""
    def diff(q, a, b):
        th, H = global_data(al, z, q)
        return H[a] - H[b]
    N = len(al)
    for a in range(N - 1):
        for b in range(a + 1, N - 1):
            lo, hi = mp.mpf(-1), mp.mpf(2)
            if diff(lo, a, b) * diff(hi, a, b) < 0:
                for _ in range(120):
                    mid = (lo + hi) / 2
                    if diff(lo, a, b) * diff(mid, a, b) <= 0:
                        hi = mid
                    else:
                        lo = mid
                q = (lo + hi) / 2
                th, H = global_data(al, z, q)
                if max(H) - max(H[a], H[b]) < mp.mpf('1e-30'):
                    return q
    return None


def main():
    rng = random.Random(SEED)
    ninst = 0
    worst = {}
    for N in range(3, 8):
        kinds = ['generic', 'close', 'smallz', 'qlow', 'qhigh', 'equal']
        for t in range(NR):
            kind = kinds[t % len(kinds)] if t >= 2 else 'generic'
            al, z, q = rand_data(rng, N, kind)
            if kind == 'equal':
                qq = equal_peaks_q(al, z)
                if qq is None:
                    kind = 'generic'
                else:
                    q = qq
            pattern = 'B' * (N - 1) if t % 3 == 0 else ''.join(rng.choice('BT') for _ in range(N - 1))
            if t % 3 == 1 and N == 3:
                pattern = 'BT'
            for drel in (mp.mpf('0.05'), mp.mpf('0.001')):
                for fr in (mp.mpf('0.25'), mp.mpf('0.9')):
                    tag = f"N={N} t={t} {kind} pi={pattern} drel={mp.nstr(drel,2)} fr={mp.nstr(fr,2)}"
                    try:
                        S, eta = make_instance(al, z, q, drel, pattern, fr_bin=fr)
                    except Exception as ex:  # noqa
                        chk('construction', False, f"{tag}: {ex}")
                        continue
                    chk('construction', True)
                    check_instance(S, tag)
                    ok1, ok2 = kappa_identities(al, z, q, drel, pattern, eta)
                    chk('kappa(eta_nb=0) = e_p + delta', ok1, tag)
                    chk('kappa(delta) - kappa(0) = delta', ok2, tag)
                    ninst += 1
                    for p in range(N - 1):
                        w = worst.setdefault('min kappa/delta', [mp.inf, ''])
                        v = S['cols'][(p, p + 1)]['kappa'] / S['delta']
                        if v < w[0]:
                            worst['min kappa/delta'] = [v, tag]
                    w = worst.setdefault('max rho', [mp.mpf(0), ''])
                    for key, s in S['secs'].items():
                        if s['rho'] > w[0]:
                            worst['max rho'] = [s['rho'], tag + f" {key}"]
                    w = worst.setdefault('eta_nb/delta', [mp.inf, ''])
                    if eta / S['delta'] < w[0]:
                        worst['eta_nb/delta'] = [eta / S['delta'], tag]
            print(f"N={N} instance {t} ({kind}) done", flush=True)
    print()
    allok = True
    for name in sorted(counts):
        nf = len(fails.get(name, []))
        st = 'PASS' if nf == 0 else 'FAIL'
        allok &= nf == 0
        print(f"{st}  {name}: {counts[name] - nf}/{counts[name]}")
        for info in fails.get(name, [])[:5]:
            print("      ", info)
    print(f"instances: {ninst}")
    for k, v in worst.items():
        print(f"{k}: {mp.nstr(v[0], 6)}  ({v[1]})")
    print("ALL PASS" if allok else "SOME FAIL")
    sys.exit(0 if allok else 1)


if __name__ == '__main__':
    set_args(sys.argv)
    main()
