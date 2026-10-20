#!/usr/bin/env python3
"""
check_astra_prefrac.py (Agent S) -- independent re-implementation (from the formulas only, no code shared) of
GPT-6-Astra's Theorem 1 (isolated prefractionator, total/attainment-astra/attainment.md, Sections 2-4) and of the
feed-matching Lemma 1 / cone Lemma 2, for m = 3, 4, 5 local components and several feed qualities.
For each random instance and several slacks eta it checks, in 60-digit arithmetic:
  (a) preferred split d0 > 0 (i < m), b0 > 0 (i > 1), V_T0, V_B0, L_T0, L_B0 > 0;
  (b) the redistribution h: sum h = 0, h_1 = h_m = 0, K(theta_j) = 0 for the middle roots, K > 0 at the extremes;
  (c) root ordering (8): a_1 > r > v > a_2 and a_{m-1} > u > s > a_m, middle roots common to both sections;
  (d) X from Lemma 1 strictly positive; pinch vectors (11) are eigenvectors (12) with the stated eigenvalues;
      e_c, e_a strictly positive, e_u (resp. e_v) positive on its face;
  (e) cone weights A, B, C, D > 0 with X = A e_c + B e_u = C e_v + D e_a;
  (f) profiles (16) for n = k = 40: every composition positive, every cross-section balance (17), feed balance,
      trace bounds (18), and B, C of order eta (B/eta, C/eta reported).
Prints PASS/FAIL. Environment: Python 3, mpmath. Runtime: under a minute.
"""
import sys
import random
import mpmath as mp

mp.mp.dps = 60
random.seed(7)
ok = True


def roots_in(fun, lo, hi):
    eps = (hi - lo) * mp.mpf('1e-40')
    return mp.findroot(fun, (lo + eps, hi - eps), solver='bisect' if False else 'anderson')


def bisect(fun, lo, hi, it=400):
    flo = fun(lo)
    for _ in range(it):
        mid = (lo + hi) / 2
        fm = fun(mid)
        if (fm > 0) == (flo > 0):
            lo, flo = mid, fm
        else:
            hi = mid
    return (lo + hi) / 2


def check(cond, msg):
    global ok
    if not cond:
        ok = False
        print('   FAIL:', msg)


def instance(m, q):
    a = sorted([mp.mpf(random.uniform(0, 2.5)) for _ in range(m)], reverse=True)
    a = [mp.e ** v for v in a]
    a = [v / a[-1] for v in a]
    f = [mp.mpf(random.uniform(0.2, 1.0)) for _ in range(m)]
    s_ = sum(f); f = [v / s_ for v in f]
    return a, f, mp.mpf(q)


def run(m, q, etas):
    a, f, q = instance(m, q)
    U = 1 - q
    feq = lambda t: sum(a[i] * f[i] / (a[i] - t) for i in range(m)) - U
    th = []
    for j in range(m - 1):
        lo, hi = a[j + 1], a[j]
        eps = (hi - lo) * mp.mpf('1e-45')
        th.append(bisect(feq, lo + eps, hi - eps))
    # preferred split (Theorem achieve): a_p d_p = V_T prod_r (a_p - th_r) / prod_{p' < m, p' != p} (a_p - a_p')
    w = []
    for p in range(m - 1):
        num = mp.fprod([a[p] - t for t in th])
        den = mp.fprod([a[p] - a[k] for k in range(m - 1) if k != p])
        w.append(num / den / a[p])
    VT0 = f[0] / w[0]
    d0 = [VT0 * w[p] for p in range(m - 1)] + [mp.mpf(0)]
    b0 = [f[i] - d0[i] for i in range(m)]
    VB0 = VT0 - U
    LT0 = VT0 - sum(d0); LB0 = VB0 + sum(b0)
    check(all(d0[i] > 0 for i in range(m - 1)) and all(b0[i] > 0 for i in range(1, m)), 'preferred split signs')
    check(abs(b0[0]) < mp.mpf('1e-40'), 'b0_1 = 0')
    check(min(VT0, VB0, LT0, LB0) > 0, 'preferred split flows')
    # redistribution h (eq. 3-5)
    mid_th = th[1:m - 2]
    beta = mp.fprod(a[1:m - 1]) / (mp.fprod(mid_th) if mid_th else 1)
    K = lambda t: (t - beta) * mp.fprod([t - x for x in mid_th]) / mp.fprod([t - a[i] for i in range(1, m - 1)])
    h = [mp.mpf(0)] * m
    for i in range(1, m - 1):
        num = (a[i] - beta) * mp.fprod([a[i] - x for x in mid_th])
        den = a[i] * mp.fprod([a[i] - a[k] for k in range(1, m - 1) if k != i])
        h[i] = num / den
    check(abs(sum(h)) < mp.mpf('1e-40'), 'sum h = 0')
    for x in mid_th:
        check(abs(K(x)) < mp.mpf('1e-40'), 'K at middle roots')
    check(K(th[0]) > 0 and K(th[-1]) > 0, 'K > 0 at extreme roots')
    tt = mp.mpf('0.3') * a[-1]
    check(abs(K(tt) - (1 - sum(a[i] * h[i] / (a[i] - tt) for i in range(m)))) < mp.mpf('1e-40'), 'partial fractions (5)')
    for eta in etas:
        eta = mp.mpf(eta)
        d = [d0[i] + eta * h[i] for i in range(m)]
        b = [b0[i] - eta * h[i] for i in range(m)]
        VT, VB, LT, LB = VT0 + eta, VB0 + eta, LT0 + eta, LB0 + eta
        check(all(d[i] > 0 for i in range(m - 1)) and all(b[i] > 0 for i in range(1, m)), 'perturbed flows positive')
        PsiT = lambda t: VT - sum(a[i] * d[i] / (a[i] - t) for i in range(m - 1))
        PsiB = lambda t: VB + sum(a[i] * b[i] / (a[i] - t) for i in range(1, m))
        eps = mp.mpf('1e-45')
        r = bisect(PsiT, a[1] + eps, a[0] - eps)
        u = bisect(PsiT, a[m - 1] + eps, a[m - 2] - eps)
        v = bisect(PsiB, a[1] + eps, a[0] - eps)
        s = bisect(PsiB, a[m - 1] + eps, a[m - 2] - eps)
        check(a[0] > r > v > a[1] and a[m - 2] > u > s > a[m - 1], f'root ordering (8) eta={eta}')
        for x in mid_th:
            check(abs(PsiT(x)) < mp.mpf('1e-40') and abs(PsiB(x)) < mp.mpf('1e-40'), 'middle roots common')
        # Lemma 1
        chi = [r] + mid_th + [s]
        om = [mp.fprod([a[i] - c for c in chi]) / (a[i] * mp.fprod([a[i] - a[k] for k in range(m) if k != i]))
              for i in range(m)]
        X = [o / sum(om) for o in om]
        G = lambda x, t: sum(a[i] * x[i] / (a[i] - t) for i in range(m))
        check(all(xi > 0 for xi in X), 'X > 0')
        for c in chi:
            check(abs(G(X, c)) < mp.mpf('1e-35'), 'G_X(chi) = 0')
        c_, aa = a[m - 1], a[0]
        ec = [c_ * d[i] / (LT * (a[i] - c_)) for i in range(m - 1)]; ec.append(1 - sum(ec))
        eu = [u * d[i] / (LT * (a[i] - u)) for i in range(m - 1)] + [mp.mpf(0)]
        ev = [mp.mpf(0)] + [v * b[i] / (LB * (v - a[i])) for i in range(1, m)]
        ea = [aa * b[i] / (LB * (aa - a[i])) for i in range(1, m)]; ea = [1 - sum(ea)] + ea
        MT = lambda x: [VT * a[i] * x[i] - d[i] * sum(a[k] * x[k] for k in range(m)) for i in range(m)]
        MB = lambda x: [VB * a[i] * x[i] + b[i] * sum(a[k] * x[k] for k in range(m)) for i in range(m)]
        for vec, M, lam, name in [(ec, MT, VT * c_, 'e_c'), (eu, MT, VT * u, 'e_u'), (ev, MB, VB * v, 'e_v'),
                                  (ea, MB, VB * aa, 'e_a')]:
            Mx = M(vec)
            check(max(abs(Mx[i] - lam * vec[i]) for i in range(m)) < mp.mpf('1e-35'), f'eigen {name}')
            check(abs(sum(vec) - 1) < mp.mpf('1e-35'), f'normalized {name}')
        check(all(x > 0 for x in ec) and all(x > 0 for x in ea), 'e_c, e_a > 0')
        check(all(x > 0 for x in eu[:m - 1]) and all(x > 0 for x in ev[1:]), 'face pinches positive')
        # Lemma 2
        Bw = -G(ec, s) / (G(eu, s) - G(ec, s)); Aw = 1 - Bw
        Cw = G(ea, r) / (G(ea, r) - G(ev, r)); Dw = 1 - Cw
        check(min(Aw, Bw, Cw, Dw) > 0, 'cone weights positive')
        check(max(abs(Aw * ec[i] + Bw * eu[i] - X[i]) for i in range(m)) < mp.mpf('1e-30'), 'X = A e_c + B e_u')
        check(max(abs(Cw * ev[i] + Dw * ea[i] - X[i]) for i in range(m)) < mp.mpf('1e-30'), 'X = C e_v + D e_a')
        # profiles (16) and balances (17)
        n = 40
        S = lambda x: sum(a[k] * x[k] for k in range(m))
        y = lambda x: [a[i] * x[i] / S(x) for i in range(m)]
        top = [[(Aw * (c_ / u) ** j * ec[i] + Bw * eu[i]) / (Aw * (c_ / u) ** j + Bw) for i in range(m)] for j in range(n + 1)]
        bot = [[(Cw * ev[i] + Dw * (v / aa) ** j * ea[i]) / (Cw + Dw * (v / aa) ** j) for i in range(m)] for j in range(n + 1)]
        check(all(xx > 0 for row in top + bot for xx in row), 'profiles positive')
        res = mp.mpf(0)
        for j in range(n):
            yt = y(top[j]); res = max(res, max(abs(VT * yt[i] - LT * top[j + 1][i] - d[i]) for i in range(m)))
            yb = y(bot[j + 1]); res = max(res, max(abs(VB * yb[i] - LB * bot[j][i] + b[i]) for i in range(m)))
        yX = y(X); y1 = y(bot[1])
        feedres = max(abs(LT * top[1][i] + VB * y1[i] + f[i] - VT * yX[i] - LB * X[i]) for i in range(m))
        check(res < mp.mpf('1e-35') and feedres < mp.mpf('1e-35'), 'balances (17) and feed stage')
        tr_top = top[n][m - 1]; tr_bot = bot[n][0]
        check(tr_top <= (Aw / Bw) * (c_ / u) ** n * (1 + mp.mpf('1e-30')) and
              tr_bot <= (Dw / Cw) * (v / aa) ** n * (1 + mp.mpf('1e-30')), 'trace bounds (18)')
        print(f'  m={m} q={float(q):+.2f} eta={float(eta):.0e}: B/eta={float(Bw / eta):.4f} C/eta={float(Cw / eta):.4f} '
              f'traces at n=40: top {float(tr_top):.2e}, bottom {float(tr_bot):.2e}; balance res {float(max(res, feedres)):.1e}')


for m in [3, 4, 5]:
    for q in [-0.5, 0.0, 0.5, 1.0, 1.4]:
        run(m, q, ['1e-2', '1e-4', '1e-6'])
print('PASS' if ok else 'FAIL')
sys.exit(0 if ok else 1)
