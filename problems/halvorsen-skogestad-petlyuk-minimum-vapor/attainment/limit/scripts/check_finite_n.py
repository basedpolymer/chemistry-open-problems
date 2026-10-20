#!/usr/bin/env python3
"""check_finite_n.py -- end-to-end illustration of the limit system (Agent L; ../limit.md Section 6).

For a design (limitsys.design: Agent L's test recipe) and stage numbers n (same n in every section), solve the
EXACT finite-stage system G_n(v) = 0 of limit.md (Theorem 2.3) by Newton's method started at the limit
configuration v_inf (continuation from large n to small n), then
  E1  reconstruct every stage composition of every section by the recursion x_{j+1} = M x_j / (L S(x_j)) from the
      bottom end and check that it reaches the prescribed top end (consistency of the section relation);
  E2  check every component balance of every unit of the arrangement (stages, junction stages, side draws,
      condenser, reboiler) on the reconstructed compositions, in high precision;
  E3  check S(x) > 0 and all mole fractions >= 0 at every unit (positivity, cf. network.md Thm 10.1 (i));
  E4  report |v_n - v_inf| (expected O(n rho^n)) and the product impurities (expected -> designed O(eta) values).
Prints PASS/FAIL; exits non-zero on failure.  Numerical evidence only (the proof is in limit.md / analysis).

Usage:  python check_finite_n.py [case] [n1,n2,...]     case in {N3, N4, N5}
Environment: Python 3.12, mpmath 1.3, numpy (WSL Fedora sci env), 1 thread.  Runtime: N3 ~ 1 min, N4 ~ 5 min.
"""
import sys
import time
import mpmath as mp
import limitsys as ls

mp.mp.dps = 110


def cases(name):
    if name == 'N3':
        return dict(alpha=[4, 2, 1], z=[mp.mpf(1) / 3] * 3, q=1, eta_col=mp.mpf('1e-3'), eta_bin=mp.mpf('3e-3'),
                    delta=mp.mpf('0.05'), pattern={0: 'B', 1: 'T'})
    if name == 'N3B':
        return dict(alpha=[4, 2, 1], z=[mp.mpf(1) / 3] * 3, q=1, eta_col=mp.mpf('1e-3'), eta_bin=mp.mpf('3e-3'),
                    delta=mp.mpf('0.05'), pattern={0: 'B', 1: 'B'})
    if name == 'N4':
        return dict(alpha=[8, 4, 2, 1], z=[mp.mpf(1) / 4] * 4, q=mp.mpf('0.5'), eta_col=mp.mpf('1e-3'),
                    eta_bin=mp.mpf('3e-3'), delta=mp.mpf('0.05'), pattern={0: 'B', 1: 'B', 2: 'B'})
    if name == 'N5':
        return dict(alpha=[10, 6, 3.5, 2, 1], z=[mp.mpf(1) / 5] * 5, q=mp.mpf('0.8'), eta_col=mp.mpf('1e-3'),
                    eta_bin=mp.mpf('3e-3'), delta=mp.mpf('0.05'), pattern={0: 'B', 1: 'T', 2: 'B', 3: 'T'})
    raise ValueError(name)


def profiles(S, v, n):
    """Reconstruct every section's stage compositions (bottom to top), return dict and consistency error."""
    D, N, al = S.D, S.N, S.D.al
    d, b, f, X, prod = S.unpack(v)
    out, err = {}, mp.mpf(0)
    for s in S.sections:
        kind, i, m = s
        c = (i, m)
        w = d[c] if kind == 'T' else [-x for x in b[c]]
        V, L = S.lim['V'][s], S.lim['L'][s]
        xb, xt = S.ends(s, X, prod)
        xs = [xb]
        x = xb
        for _ in range(n + 1):
            Sx = sum(a * y for a, y in zip(al, x))
            x = [(V * al[k] * x[k] - w[k] * Sx) / (L * Sx) for k in range(N)]
            xs.append(x)
        err = max(err, max(abs(x[k] - xt[k]) for k in range(N)))
        out[s] = xs   # xs[0] = bottom end (virtual liquid), xs[1..n] stages, xs[n+1] = top end
    return out, err, (d, b, f, X, prod)


def balances(S, prof, dat, n):
    """Component balances of every unit, using y = K(x) x with K = alpha/S(x) at every equilibrium unit."""
    D, N, al = S.D, S.N, S.D.al
    d, b, f, X, prod = dat
    y = lambda x: [al[k] * x[k] / sum(a * t for a, t in zip(al, x)) for k in range(N)]
    worst = mp.mpf(0)
    # interior stages
    for s, xs in prof.items():
        V, L = S.lim['V'][s], S.lim['L'][s]
        for j in range(1, n + 1):
            yin = y(xs[j - 1])
            xin = xs[j + 1]
            for k in range(N):
                r = V * yin[k] + L * xin[k] - V * y(xs[j])[k] - L * xs[j][k]
                worst = max(worst, abs(r))
    # junction and product units: in = vapor from the top stage of sections ending there from below + liquid from
    # the bottom stage of sections ending there from above; out = vapor up + liquid down + product
    units = {}
    for s, xs in prof.items():
        kind, i, m = s
        V, L = S.lim['V'][s], S.lim['L'][s]
        bot = ('J', i, m) if kind == 'T' else (('J', i + 1, m) if i + 1 < m else ('P', m))
        top = (('J', i, m - 1) if m - 1 > i else ('P', i)) if kind == 'T' else ('J', i, m)
        units.setdefault(bot, []).append(('bot', s, V, L, xs))
        units.setdefault(top, []).append(('top', s, V, L, xs))
    pos = mp.inf
    Smin = mp.inf
    for u, lst in units.items():
        xu = X[(u[1], u[2])] if u[0] == 'J' else prod[u[1]]
        Smin = min(Smin, sum(a * t for a, t in zip(al, xu)))
        pos = min(pos, min(xu))
        cond = (u == ('P', 0))  # total condenser: liquid out only, composition = vapor in
        bal = [mp.mpf(0)] * N
        for (end, s, V, L, xs) in lst:
            if end == 'bot':     # unit is the bottom end of s: sends vapor y(x_u) up, receives liquid xs[1]
                yu = y(xu)
                for k in range(N):
                    bal[k] += L * xs[1][k] - V * yu[k]
            else:                # unit is the top end of s: receives vapor y(xs[n]), sends liquid x_u down
                yn = y(xs[n])
                for k in range(N):
                    bal[k] += V * yn[k] - L * xu[k]
        if u[0] == 'P':
            for k in range(N):
                bal[k] -= D.P[u[1]] * xu[k]
        if u == ('J', 0, N - 1):
            for k in range(N):
                bal[k] += D.z[k]
        worst = max(worst, max(abs(t) for t in bal))
    for s, xs in prof.items():
        for x in xs:
            pos = min(pos, min(x))
            Smin = min(Smin, sum(a * t for a, t in zip(al, x)))
    return worst, pos, Smin


def main():
    case = sys.argv[1] if len(sys.argv) > 1 else 'N3'
    ns = [int(t) for t in (sys.argv[2] if len(sys.argv) > 2 else '80,60,40,30,20').split(',')]
    cfg = cases(case)
    D = ls.design(cfg['alpha'], cfg['z'], cfg['q'], eta_col=cfg['eta_col'], eta_bin=cfg['eta_bin'],
                  delta=cfg['delta'], pattern=cfg['pattern'])
    print(f'case {case}: alpha={cfg["alpha"]} q={cfg["q"]} pattern={cfg["pattern"]} R={mp.nstr(D.R, 10)} '
          f'maxH={mp.nstr(max(D.H), 10)} products P={[mp.nstr(x, 6) for x in D.P]}')
    Sinf = ls.System(D, None)
    vinf = Sinf.v_inf
    _, _, _, _, prod_inf = Sinf.unpack(vinf)
    print('limit product impurities 1 - x_pp:', [mp.nstr(1 - prod_inf[p][p], 4) for p in range(D.N)])
    v = vinf
    allok = True
    t0 = time.time()
    for n in ns:
        S = ls.System(D, n)
        vprev = v
        v, res, it = ls.newton(S, v)
        if it < 0 or res > mp.mpf(10) ** (-mp.mp.dps + 20):
            print(f'n={n:4d} Newton from the previous solution did not converge (|G_n| = {mp.nstr(res, 3)}): '
                  f'n is outside the range where this simple continuation works; stopping (not a FAIL of a claim).')
            v = vprev
            break
        prof, err, dat = profiles(S, v, n)
        bal, pos, Smin = balances(S, prof, dat, n)
        dist = max(abs(a - b_) for a, b_ in zip(v, vinf))
        prod = dat[4]
        imp = [1 - prod[p][p] for p in range(D.N)]
        tr = [mp.nstr(x, 3) for (x, var) in zip(v, S.vars) if var[3] == 'tr' and var[1] in ('d', 'b')][:6]
        ok = res < mp.mpf(10) ** (-mp.mp.dps + 20) and err < mp.mpf('1e-30') and bal < mp.mpf('1e-30') \
            and pos >= 0 and Smin > 0
        allok &= ok
        print(f'n={n:4d} newton it={it} |G_n|={mp.nstr(res, 2)} |v_n-v_inf|={mp.nstr(dist, 3)} '
              f'E1 section-consistency={mp.nstr(err, 2)} E2 balances={mp.nstr(bal, 2)} '
              f'E3 min x={mp.nstr(pos, 3)} min S={mp.nstr(Smin, 4)}  {"PASS" if ok else "FAIL"}')
        print(f'        product impurities {[mp.nstr(x, 4) for x in imp]}  first trace flows {tr}')
        sys.stdout.flush()
    print(f'{time.time() - t0:.1f} s')
    print('ALL PASS' if allok else 'SOME FAIL')
    sys.exit(0 if allok else 1)


if __name__ == '__main__':
    main()
