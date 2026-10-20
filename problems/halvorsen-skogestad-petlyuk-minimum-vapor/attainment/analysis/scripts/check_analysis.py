#!/usr/bin/env python3
"""Agent A (analysis): numerical check of the analytic step of the attainment proof (attainment/analysis/analysis.md).

NUMERICAL EVIDENCE ONLY.  The proofs are in analysis.md; this script checks, on concrete data, that the exact
formulation used there behaves as proved:
  * the network map G_n of analysis.md §3 (unknowns: column splits d^c and junction compositions X^c only;
    section relations in the product form (2.2) with generic numerically computed eigen-data, NOT the Underwood
    formulas) vanishes at the limit configuration for n = infinity (A4), with nonzero pivots (A3);
  * DG_inf(v_inf) is invertible (A6, Theorem L) -- smallest singular value printed;
  * the limit half-profiles are nonnegative (A5);
  * Newton from v_inf converges to a zero v_n of G_n for every n tested, also with unequal section lengths
    (Theorem A.5), and |v_n - v_inf| / rho_inf^(n+1) stays bounded (rate rho_inf of Theorem A.5);
  * the reconstructed stage profiles (Lemma A.6) give S > 0 at every stage and all mole fractions > 0, and EVERY
    stage, junction, side-draw, condenser and reboiler balance holds with the actual streams (Lemma A.1);
  * the composite-profile bound of Lemma A.7: max_k min(|u_k - xi^bot_k|, |u_k - xi^top_{n+1-k}|) is small.
The limit configuration and the specification are built with Agent D's library (attainment/design/scripts/
design_lib.py, imported read-only); everything else (eigen-data, network map, Newton, reconstruction, balance
checks) is independent of D, L and C.

Usage (from the project root, Git Bash):
  outils/wsl "cd chemistry-open-problems/problems/halvorsen-skogestad-petlyuk-minimum-vapor/attainment/analysis/scripts && python check_analysis.py [case ...]"
Cases: n3 (alpha = (4,2,1), z = 1/3, q = 1), n3q (q = -0.3), n4 (alpha = (6,3.5,2,1), q = 0.6), n4T (pattern BTB),
n5 (alpha = (8,5,3,1.8,1), q = 0.8).  Default: n3 n3q n4 n4T.
Environment: Python 3, mpmath (WSL Fedora sci environment); single thread.  Working precision per case:
50 + (n_max + 2) log10(1.5 alpha_1/alpha_N) digits (the forward reconstruction of a section amplifies rounding by
(lambda_max/lambda_min)^(n+1)).  Runtime: about 1 min (N = 3) to 15 min (N = 4) per case; N = 5 longer.  Prints PASS/FAIL per claim; exit code 1 on any FAIL.
"""
import os
import sys
import time
import hashlib

import mpmath as mp

HERE = os.path.dirname(os.path.abspath(__file__))
DLIB = os.path.normpath(os.path.join(HERE, '..', '..', 'design', 'scripts'))
sys.path.insert(0, DLIB)
import design_lib as DL  # noqa: E402

mp.mp.dps = 60
FAILS = []


def check(name, ok, info=''):
    print(('PASS ' if ok else 'FAIL ') + name + ('  ' + info if info else ''))
    if not ok:
        FAILS.append(name)


# ---------------------------------------------------------------------------------------------------------------
# network bookkeeping in the variables of analysis.md §1


class Net:
    def __init__(self, S):
        self.S = S
        self.N = S['N']
        self.al = S['al']
        self.cols = sorted(S['cols'].keys(), key=lambda c: (-(c[1] - c[0]), c[0]))   # top-down
        self.secs = S['secs']
        self.seckeys = sorted(self.secs.keys())
        N = self.N
        # base point v_inf and the parametrization of the affine space V
        self.d0 = {c: [mp.mpf(x) for x in S['cols'][c]['d']] for c in self.cols}
        self.X0 = {c: [mp.mpf(x) for x in S['X'][('J',) + c]] for c in self.cols}
        self.dim = 2 * (N - 1) * len(self.cols)

    def unpack(self, t):
        """t in R^dim (displacement) -> (d^c, X^c) with fixed sums."""
        N = self.N
        d, X = {}, {}
        k = 0
        for c in self.cols:
            dd = [t[k + j] for j in range(N - 1)]
            k += N - 1
            dX = [t[k + j] for j in range(N - 1)]
            k += N - 1
            d[c] = [self.d0[c][j] + (dd[j] if j < N - 1 else -mp.fsum(dd)) for j in range(N)]
            X[c] = [self.X0[c][j] + (dX[j] if j < N - 1 else -mp.fsum(dX)) for j in range(N)]
        return d, X

    def derived(self, t):
        """Net flows of every section, unit compositions (analysis.md §1)."""
        N = self.N
        S = self.S
        d, X = self.unpack(t)
        f, b = {}, {}
        for c in self.cols:                      # top-down
            i, m = c
            if c == (0, N - 1):
                fc = [mp.mpf(x) for x in S['z']]
            else:
                fc = [mp.mpf(0)] * N
                if m + 1 <= N - 1:
                    fc = [a + x for a, x in zip(fc, d[(i, m + 1)])]
                if i - 1 >= 0:
                    fc = [a + x for a, x in zip(fc, b[(i - 1, m)])]
            f[c] = fc
            b[c] = [fc[j] - d[c][j] for j in range(N)]
        comp = {}
        for c in self.cols:
            comp[('J',) + c] = X[c]
        for p in range(N):
            n = [mp.mpf(0)] * N
            if p - 1 >= 0:
                n = [a + x for a, x in zip(n, b[(p - 1, p)])]
            if p + 1 <= N - 1:
                n = [a + x for a, x in zip(n, d[(p, p + 1)])]
            P = S['prods'][p]['P']
            comp[('P', p)] = [x / P for x in n]
        w = {}
        for key in self.seckeys:
            kind, i, m = key
            w[key] = list(d[(i, m)]) if kind == 'T' else [-x for x in b[(i, m)]]
        return w, comp, f, b


# ---------------------------------------------------------------------------------------------------------------
# section eigen-data (generic, numerical) with an analytic normalization


def Mmat(V, w, al):
    N = len(al)
    return mp.matrix([[V * al[i] * (1 if i == j else 0) - w[i] * al[j] for j in range(N)] for i in range(N)])


def eigdata(V, w, al, ref=None):
    """Eigenvalues decreasing, right eigenvectors e_r (columns), left eigenvectors l_r with l_r . e_s = delta_rs.
    Normalization (analytic in w near the reference): e_r scaled so that ref_l[r] . e_r = 1 when a reference is
    given, else 1^T e_r = 1 (or max-abs entry 1 if 1^T e_r ~ 0)."""
    N = len(al)
    M = Mmat(V, w, al)
    E, ER = mp.eig(M)
    lam = []
    for j in range(N):
        ev = E[j]
        if isinstance(ev, mp.mpc):
            assert abs(ev.imag) < mp.mpf(10) ** (-mp.mp.dps // 2), 'complex eigenvalue'
            ev = ev.real
        lam.append(ev)
    order = sorted(range(N), key=lambda j: -lam[j])
    lam = [lam[j] for j in order]
    vecs = []
    for j in order:
        v = [ER[i, j].real if isinstance(ER[i, j], mp.mpc) else ER[i, j] for i in range(N)]
        vecs.append(v)
    out_e = []
    for r, v in enumerate(vecs):
        if ref is not None:
            s = mp.fsum(a * x for a, x in zip(ref['l'][r], v))
        else:
            s = mp.fsum(v)
            if abs(s) < mp.mpf('1e-8') * max(abs(x) for x in v):
                s = max(v, key=abs)
        out_e.append([x / s for x in v])
    Emat = mp.matrix([[out_e[r][i] for r in range(N)] for i in range(N)])
    Linv = mp.inverse(Emat)
    out_l = [[Linv[r, i] for i in range(N)] for r in range(N)]
    return {'lam': lam, 'e': out_e, 'l': out_l}


def dot(a, b):
    return mp.fsum(x * y for x, y in zip(a, b))


# ---------------------------------------------------------------------------------------------------------------
# the network map G_n (analysis.md (2.2)) and G_inf


class System:
    def __init__(self, net):
        self.net = net
        S = net.S
        al = net.al
        # limit eigen-data, pinch indices, references for the normalization
        w0, comp0, _, _ = net.derived([mp.mpf(0)] * net.dim)
        self.ref, self.rstar = {}, {}
        for key in net.seckeys:
            sec = net.secs[key]
            ed = eigdata(sec['V'], w0[key], al)
            # normalize the pinch right eigenvector by 1^T e = 1 (analysis.md (A5)); others by 1^T e or max entry
            self.ref[key] = ed
            pv = sec['pinch']['val']
            rs = min(range(net.N), key=lambda r: abs(ed['lam'][r] / sec['V'] - pv))
            assert abs(ed['lam'][rs] / sec['V'] - pv) < mp.mpf('1e-30'), 'pinch not found'
            self.rstar[key] = rs

    def sec_rows(self, key, w, xb, xt, nsig, ed=None):
        sec = self.net.secs[key]
        if ed is None:
            ed = eigdata(sec['V'], w, self.net.al, ref=self.ref[key])
        rs = self.rstar[key]
        lam, l = ed['lam'], ed['l']
        pb, pt = dot(l[rs], xb), dot(l[rs], xt)
        rows = []
        for r in range(self.net.N):
            if r == rs:
                continue
            if r < rs:      # above the pinch: closed at the bottom
                mu = lam[rs] / lam[r]
                t = mu ** (nsig + 1) if nsig is not None else 0
                rows.append(dot(l[r], xb) * pt - t * dot(l[r], xt) * pb)
            else:           # below the pinch: closed at the top
                mu = lam[r] / lam[rs]
                t = mu ** (nsig + 1) if nsig is not None else 0
                rows.append(dot(l[r], xt) * pb - t * dot(l[r], xb) * pt)
        return rows, pb, pt, ed

    def G(self, t, nmap):
        """nmap: dict section -> n_sigma, or None for G_inf."""
        net = self.net
        w, comp, _, _ = net.derived(t)
        out = []
        for key in net.seckeys:
            sec = net.secs[key]
            xb, xt = comp[sec['bottom']], comp[sec['top']]
            rows, _, _, _ = self.sec_rows(key, w[key], xb, xt, None if nmap is None else nmap[key])
            out += rows
        return out

    def jac(self, t, nmap, h=None):
        if h is None:
            h = mp.mpf(10) ** (-(mp.mp.dps // 2 - 5))
        n = len(t)
        cols = []
        for j in range(n):
            tp = list(t)
            tm = list(t)
            tp[j] += h
            tm[j] -= h
            gp, gm = self.G(tp, nmap), self.G(tm, nmap)
            cols.append([(a - b) / (2 * h) for a, b in zip(gp, gm)])
        return mp.matrix([[cols[j][i] for j in range(n)] for i in range(n)])


def newton(sys_, nmap, t0, tol=None, maxit=60):
    """Damped Newton (backtracking on the max-norm residual; a step that leaves the domain of real simple
    eigen-data is halved)."""
    if tol is None:
        tol = mp.mpf(10) ** (-(mp.mp.dps - 10))
    t = list(t0)
    g = sys_.G(t, nmap)
    res = max(abs(x) for x in g)
    for it in range(maxit):
        if res < tol:
            return t, res, it
        J = sys_.jac(t, nmap)
        dt = mp.lu_solve(J, mp.matrix(g))
        lam = mp.mpf(1)
        while True:
            tn = [a - lam * dt[i] for i, a in enumerate(t)]
            try:
                gn = sys_.G(tn, nmap)
                rn = max(abs(x) for x in gn)
            except (AssertionError, ZeroDivisionError):
                rn = mp.inf
            if rn < res or lam < mp.mpf('1e-6'):
                break
            lam /= 2
        t, g, res = tn, gn, rn
    return t, res, maxit


# ---------------------------------------------------------------------------------------------------------------
# reconstruction and full balance check (Lemma A.1)


def reconstruct_and_check(sys_, t, nmap):
    net = sys_.net
    S = net.S
    al = net.al
    N = net.N
    w, comp, f, b = net.derived(t)
    prof = {}
    minS, minx = mp.inf, mp.inf
    worst_cs = mp.mpf(0)
    for key in net.seckeys:
        sec = net.secs[key]
        V, L = sec['V'], sec['L']
        M = Mmat(V, w[key], al)
        x = [comp[sec['bottom']]]
        n = nmap[key]
        for k in range(n + 1):
            u = M * mp.matrix(x[-1])
            s = mp.fsum(u[i] for i in range(N))
            x.append([u[i] / s for i in range(N)])
        # x[n+1] must equal the top-end composition
        xt = comp[sec['top']]
        worst_cs = max(worst_cs, max(abs(a - c) for a, c in zip(x[n + 1], xt)))
        for k in range(n + 1):
            Sk = dot(al, x[k])
            minS = min(minS, Sk)
        for k in range(1, n + 1):
            minx = min(minx, min(x[k]))
        prof[key] = x
    for u_, xc in comp.items():
        minx = min(minx, min(xc))
        if u_ != ('P', 0):
            minS = min(minS, dot(al, xc))
    # unit balances with actual streams
    y = lambda xx: [a * v / dot(al, xx) for a, v in zip(al, xx)]
    worst_unit = mp.mpf(0)
    worst_stage = mp.mpf(0)
    for key in net.seckeys:          # interior stages
        sec = net.secs[key]
        V, L = sec['V'], sec['L']
        x = prof[key]
        n = nmap[key]
        for k in range(1, n + 1):
            yin = y(x[k - 1]) if k > 1 else y(comp[sec['bottom']]) if sec['bottom'] != ('P', 0) else None
            res = [L * x[k + 1][i] + V * yin[i] - L * x[k][i] - V * y(x[k])[i] for i in range(N)]
            worst_stage = max(worst_stage, max(abs(r) for r in res))
    units = [('J',) + c for c in net.cols] + [('P', p) for p in range(N)]
    for u_ in units:
        bal = [mp.mpf(0)] * N
        xu = comp[u_]
        for key in net.seckeys:
            sec = net.secs[key]
            V, L = sec['V'], sec['L']
            x = prof[key]
            n = nmap[key]
            if sec['bottom'] == u_:      # liquid in from stage 1, vapor out
                bal = [bb + L * x[1][i] - V * y(xu)[i] for i, bb in enumerate(bal)]
            if sec['top'] == u_:         # vapor in from stage n, liquid out
                bal = [bb + V * y(x[n])[i] - L * xu[i] for i, bb in enumerate(bal)]
        if u_ == ('J', 0, N - 1):
            bal = [bb + S['z'][i] for i, bb in enumerate(bal)]
        if u_[0] == 'P':
            P = S['prods'][u_[1]]['P']
            bal = [bb - P * xu[i] for i, bb in enumerate(bal)]
        worst_unit = max(worst_unit, max(abs(r) for r in bal))
    # condenser composition rule: liquid = incoming vapor
    keyT12 = ('T', 0, 1)
    cond = max(abs(a - c) for a, c in zip(y(prof[keyT12][nmap[keyT12]]), comp[('P', 0)]))
    return prof, minS, minx, worst_cs, worst_stage, worst_unit, cond


def spec_balances(S):
    """Total and CMO balances at every unit of the specification (A1)."""
    N = S['N']
    secs = S['secs']
    worst_tot, worst_vap = mp.mpf(0), mp.mpf(0)
    duties = {}
    units = [('J', i, m) for i in range(N) for m in range(i + 1, N)] + [('P', p) for p in range(N)]
    for u in units:
        tin, tout, vin, vout = mp.mpf(0), mp.mpf(0), mp.mpf(0), mp.mpf(0)
        for key, sec in secs.items():
            if sec['bottom'] == u:
                tin += sec['L']
                tout += sec['V']
                vout += sec['V']
            if sec['top'] == u:
                tin += sec['V']
                vin += sec['V']
                tout += sec['L']
        if u == ('J', 0, N - 1):
            tin += 1
            vin += 1 - S['q']
        if u[0] == 'P':
            tout += S['prods'][u[1]]['P']
        worst_tot = max(worst_tot, abs(tin - tout))
        Q = vout - vin
        if u == ('P', 0):
            Q = -vin             # total condenser: no vapor out
        duties[u] = Q
    return worst_tot, duties


# ---------------------------------------------------------------------------------------------------------------


CASES = {
    'n3': dict(al=[4, 2, 1], z=[mp.mpf(1) / 3] * 3, q=1, delta='0.05', eta_nb='0.025', eta_bin='0.003',
               pattern='BB', ns=[8, 16, 32, 64]),
    'n3q': dict(al=[4, 2, 1], z=['0.5', '0.2', '0.3'], q='-0.3', delta='0.05', eta_nb='0.02', eta_bin='0.004',
                pattern='TB', ns=[8, 16, 32, 64]),
    'n4': dict(al=[6, '3.5', 2, 1], z=['0.25', '0.3', '0.2', '0.25'], q='0.6', delta='0.05', eta_nb='0.01',
               eta_bin='0.003', pattern='BBB', ns=[8, 16, 32, 48]),
    'n4T': dict(al=[6, '3.5', 2, 1], z=['0.25', '0.3', '0.2', '0.25'], q='0.6', delta='0.05', eta_nb='0.01',
                eta_bin='0.003', pattern='BTB', ns=[8, 16, 32, 48]),
    'n5': dict(al=[8, 5, 3, '1.8', 1], z=['0.2'] * 5, q='0.8', delta='0.05', eta_nb='0.01', eta_bin='0.003',
               pattern='BBBB', ns=[12, 24, 36]),
}


def run_case(name):
    cs = CASES[name]
    t0 = time.time()
    # working precision: forward reconstruction amplifies rounding by (lambda_max/lambda_min)^(n+1) <~ kappa^(n+1)
    kappa = 1.5 * float(cs['al'][0]) / float(cs['al'][-1])
    nmax = max(max(cs['ns']), 2 * cs['ns'][-2])
    mp.mp.dps = 50 + int((nmax + 2) * mp.log10(kappa)) + 1
    print('working precision: %d digits' % mp.mp.dps)
    print('=' * 100)
    print('case', name, {k: v for k, v in cs.items() if k != 'ns'})
    S = DL.build_spec([mp.mpf(a) for a in cs['al']], [mp.mpf(x) for x in cs['z']], mp.mpf(cs['q']),
                      mp.mpf(cs['delta']), mp.mpf(cs['eta_nb']), mp.mpf(cs['eta_bin']), cs['pattern'])
    S = DL.limit_config(S)
    N = S['N']
    print('H =', [mp.nstr(h, 8) for h in S['H']], ' R = H* + delta =', mp.nstr(S['R'], 10))
    # (A1)
    wt, duties = spec_balances(S)
    pos = all(sec['V'] > 0 and sec['L'] > 0 for sec in S['secs'].values()) and \
        all(pr['P'] > 0 for pr in S['prods'].values())
    zero_duty = all(abs(Q) < mp.mpf('1e-40') for u, Q in duties.items() if u not in (('P', 0), ('P', N - 1)))
    Rr = duties[('P', N - 1)]
    check(name + ' (A1) positive flows, total balances, zero duty except condenser/reboiler, R = H*+delta',
          pos and wt < mp.mpf('1e-40') and zero_duty and abs(Rr - S['R']) < mp.mpf('1e-40'),
          'total-balance residual %s, reboiler duty %s' % (mp.nstr(wt, 3), mp.nstr(Rr, 12)))
    net = Net(S)
    sys_ = System(net)
    zero = [mp.mpf(0)] * net.dim
    # (A2) spectra, compositions in the simplex
    w0, comp0, _, _ = net.derived(zero)
    okspec, rho_inf, rho1 = True, mp.mpf(0), mp.mpf(0)
    for key in net.seckeys:
        lam = sys_.ref[key]['lam']
        okspec &= all(l > 0 for l in lam) and all(lam[r] > lam[r + 1] * (1 + mp.mpf('1e-12')) for r in range(N - 1))
        rs = sys_.rstar[key]
        rho1 = max(rho1, max(lam[r + 1] / lam[r] for r in range(N - 1)))
        adj = []
        if rs + 1 < N:
            adj.append(lam[rs + 1] / lam[rs])
        if rs - 1 >= 0:
            adj.append(lam[rs] / lam[rs - 1])
        rho_inf = max(rho_inf, max(adj))
    insimplex = all(min(x) >= 0 and abs(mp.fsum(x) - 1) < mp.mpf('1e-40') for x in comp0.values())
    check(name + ' (A2) N distinct positive eigenvalues per section; unit compositions in the simplex',
          okspec and insimplex, 'rho_inf = %s, max consecutive ratio = %s' % (mp.nstr(rho_inf, 6), mp.nstr(rho1, 6)))
    # (A3), (A4)
    minpiv, maxres = mp.inf, mp.mpf(0)
    for key in net.seckeys:
        sec = net.secs[key]
        rows, pb, pt, _ = sys_.sec_rows(key, w0[key], comp0[sec['bottom']], comp0[sec['top']], None,
                                        ed=sys_.ref[key])
        minpiv = min(minpiv, abs(pb), abs(pt))
        maxres = max(maxres, max(abs(r) for r in rows))
    check(name + ' (A3) nonzero pivots at both ends of every section', minpiv > mp.mpf('1e-6'),
          'min |pivot| = %s' % mp.nstr(minpiv, 4))
    check(name + ' (A4) G_inf(v_inf) = 0', maxres < mp.mpf('1e-45'), 'max |row| = %s' % mp.nstr(maxres, 3))
    # (A5) half-profiles, normalization 1^T e_{r*} = 1 = l_{r*} . e_{r*}
    worst_neg, minsum = mp.mpf(0), mp.inf
    for key in net.seckeys:
        sec = net.secs[key]
        ed = sys_.ref[key]
        rs = sys_.rstar[key]
        es = ed['e'][rs]
        sc = mp.fsum(es)
        e_star = [x / sc for x in es]
        l_star = [x * sc for x in ed['l'][rs]]
        xb, xt = comp0[sec['bottom']], comp0[sec['top']]
        pb, pt = dot(l_star, xb), dot(l_star, xt)
        for k in range(0, 400):
            xi_b = list(e_star)
            xi_t = list(e_star)
            for r in range(N):
                if r > rs:
                    a = dot(ed['l'][r], xb) / pb
                    xi_b = [u + a * (ed['lam'][r] / ed['lam'][rs]) ** k * v for u, v in zip(xi_b, ed['e'][r])]
                if r < rs:
                    bb = dot(ed['l'][r], xt) / pt
                    xi_t = [u + bb * (ed['lam'][rs] / ed['lam'][r]) ** k * v for u, v in zip(xi_t, ed['e'][r])]
            worst_neg = min(worst_neg, min(xi_b), min(xi_t))
            minsum = min(minsum, mp.fsum(xi_b), mp.fsum(xi_t))
    check(name + ' (A5) limit half-profiles nonnegative (k = 0..399), 1^T xi >= c_0 > 0',
          worst_neg > -mp.mpf('1e-40') and minsum > 0,
          'most negative entry %s, c_0 = %s' % (mp.nstr(worst_neg, 3), mp.nstr(minsum, 5)))
    # (A6)
    J = sys_.jac(zero, None)
    U_, Sv, V_ = mp.svd_r(J)
    smin = min(Sv[i] for i in range(len(Sv)))
    check(name + ' (A6) DG_inf(v_inf) invertible (numerically; Theorem L is the proof)', smin > mp.mpf('1e-8'),
          'dim %d, sigma_min = %s, sigma_max = %s' % (net.dim, mp.nstr(smin, 5),
                                                      mp.nstr(max(Sv[i] for i in range(len(Sv))), 5)))
    # Theorem A.5 / A.8 for several n, equal lengths, then unequal lengths
    # continuation: largest n first (start at v_inf), then decreasing n from the previous zero
    runs = [(n, {key: n for key in net.seckeys}) for n in sorted(cs['ns'], reverse=True)]
    nmix = cs['ns'][-2]
    mix = {key: nmix + (7 * j) % (nmix + 1) for j, key in enumerate(net.seckeys)}   # lengths in [nmix, 2 nmix]
    runs.insert(1, ('mixed(min %d)' % nmix, mix))
    tprev = zero
    for lab, nmap in runs:
        t, res, it = newton(sys_, nmap, tprev)
        dist = max(abs(x) for x in t)
        nmin = min(nmap.values())
        ratio = dist / rho_inf ** (nmin + 1)
        prof, minS, minx, wcs, wst, wun, cond = reconstruct_and_check(sys_, t, nmap)
        imp = []
        w, comp, _, _ = net.derived(t)
        for p in range(N):
            imp.append(1 - comp[('P', p)][p])
        ok = res < mp.mpf('1e-45') and minS > 0 and minx > 0 and wcs < mp.mpf('1e-40') and \
            wst < mp.mpf('1e-40') and wun < mp.mpf('1e-40') and cond < mp.mpf('1e-40')
        check('%s n=%s: zero of G_n found; S > 0 and all x > 0 at every stage; every balance holds' % (name, lab),
              ok, 'Newton it %d res %s | |v_n-v_inf| = %s, /rho_inf^(n+1) = %s | min S %s, min x %s | '
              'balances: cross-sec %s stage %s unit %s cond %s | impurities %s' % (
                  it, mp.nstr(res, 2), mp.nstr(dist, 4), mp.nstr(ratio, 4), mp.nstr(minS, 4), mp.nstr(minx, 3),
                  mp.nstr(wcs, 2), mp.nstr(wst, 2), mp.nstr(wun, 2), mp.nstr(cond, 2),
                  [mp.nstr(x, 3) for x in imp]))
        # Lemma A.7 composite bound
        worst = mp.mpf(0)
        for key in net.seckeys:
            sec = net.secs[key]
            ed = sys_.ref[key]
            rs = sys_.rstar[key]
            sc = mp.fsum(ed['e'][rs])
            e_star = [x / sc for x in ed['e'][rs]]
            l_star = [x * sc for x in ed['l'][rs]]
            xb0, xt0 = comp0[sec['bottom']], comp0[sec['top']]
            pb, pt = dot(l_star, xb0), dot(l_star, xt0)
            n = nmap[key]
            x = prof[key]
            for k in range(n + 2):
                xi_b = list(e_star)
                xi_t = list(e_star)
                for r in range(N):
                    if r > rs:
                        a = dot(ed['l'][r], xb0) / pb
                        xi_b = [u + a * (ed['lam'][r] / ed['lam'][rs]) ** k * v for u, v in zip(xi_b, ed['e'][r])]
                    if r < rs:
                        bb = dot(ed['l'][r], xt0) / pt
                        xi_t = [u + bb * (ed['lam'][rs] / ed['lam'][r]) ** (n + 1 - k) * v
                                for u, v in zip(xi_t, ed['e'][r])]
                xb_n = [u / mp.fsum(xi_b) for u in xi_b]
                xt_n = [u / mp.fsum(xi_t) for u in xi_t]
                dmin = min(max(abs(a - c) for a, c in zip(x[k], xb_n)), max(abs(a - c) for a, c in zip(x[k], xt_n)))
                worst = max(worst, dmin)
        print('     Lemma A.7: max_k min(|x_k - xi^bot_k|, |x_k - xi^top_(n+1-k)|) (normalized) = %s' %
              mp.nstr(worst, 4))
        tprev = t
    print('case %s done in %.1f s' % (name, time.time() - t0))


def main():
    cases = sys.argv[1:] or ['n3', 'n3q', 'n4', 'n4T']
    with open(os.path.join(DLIB, 'design_lib.py'), 'rb') as fh:
        print('design_lib.py sha256 =', hashlib.sha256(fh.read()).hexdigest()[:16], ' mp.dps =', mp.mp.dps)
    for c in cases:
        run_case(c)
    print('=' * 100)
    print('ALL PASS' if not FAILS else 'FAILURES: %s' % FAILS)
    sys.exit(1 if FAILS else 0)


if __name__ == '__main__':
    main()
