#!/usr/bin/env python3
"""Agent D: numerical illustration of Proposition 5.4(d) / Theorem 5.5 (design.md): the binary margins satisfy
|kappa_p(eta_nb; 0) - e_p| <= C_1 max eta_c near eta_nb = 0, so eta-bar(delta) = min(eta_1, delta/(2 C_1)) is linear
in delta.  For several data sets (N = 3..6) and eta_nb = eta * (1, ..., 1) with eta = 10^-1 .. 10^-6, prints
max_p |kappa_p(eta; 0) - e_p| / eta (-> the directional derivative, a lower bound for C_1), the explicit validity
bound eta_0 = min over non-binary columns of eta_0(c) at eta = 0 (design.md §3), and, for delta = 10^-2 .. 10^-5, the
largest eta (bisection on a log grid) with validity and min_p kappa_p >= e_p + delta/2, divided by delta.
PASS if, for each data set, that ratio stays bounded below as delta -> 0 (linear scaling), i.e. the ratio at
delta = 1e-5 is at least half the ratio at delta = 1e-3.

Usage: python eta_scaling.py      Environment: Python 3 + mpmath.  Runtime: about 1-2 minutes.
Numerical evidence only.
"""
import sys
import mpmath as mp
from design_lib import build_spec, global_data, InvalidSpec

mp.mp.dps = 40

DATA = [
    ('N3 (4,2,1) z=1/3 q=1', [4, 2, 1], [mp.mpf(1) / 3] * 3, 1),
    ('N4 (8,4,2,1) z=1/4 q=1', [8, 4, 2, 1], [mp.mpf(1) / 4] * 4, 1),
    ('N4 (6,3.5,1.8,1) q=0', [6, 3.5, 1.8, 1], [0.3, 0.2, 0.15, 0.35], 0),
    ('N5 (16..1) q=0.5', [16, 8, 4, 2, 1], [mp.mpf(1) / 5] * 5, mp.mpf('0.5')),
    ('N5 (5,3.6,2.4,1.6,1) q=0.6', [5, 3.6, 2.4, 1.6, 1], [0.15, 0.25, 0.2, 0.3, 0.1], mp.mpf('0.6')),
    ('N6 (7,5,3.4,2.3,1.5,1) q=0.8', [7, 5, 3.4, 2.3, 1.5, 1], [0.2, 0.1, 0.25, 0.15, 0.1, 0.2], mp.mpf('0.8')),
]


def kap(al, z, q, eta, delta=0):
    S = build_spec(al, z, q, delta, eta, 0, 'B' * (len(al) - 1))
    return [S['cols'][(p, p + 1)]['kappa'] for p in range(len(al) - 1)], S


def eta0(S):
    """explicit validity bound min_c eta_0(c) (design.md §3), computed at the preferred split (eta = 0)."""
    best = mp.inf
    for (i, m), c in S['cols'].items():
        if c['kind'] != 'nb':
            continue
        for p in range(i, m + 1):
            hp = c['h'].get(p, 0)
            if hp < 0 and p < m:
                best = min(best, c['d0'][p] / -hp)
            if hp > 0 and p > i:
                best = min(best, (c['f'][p] - c['d0'].get(p, 0)) / hp)
    return best


def largest_eta(al, z, q, delta, e):
    lo, hi = mp.mpf(0), None
    eta = mp.mpf('0.5')
    for _ in range(200):
        try:
            k, S = kap(al, z, q, eta, delta)
            ok = all(kp >= ep + delta / 2 for kp, ep in zip(k, e))
        except InvalidSpec:
            ok = False
        if ok:
            lo = eta
            if hi is None:
                break
        else:
            hi = eta
        eta = eta / 2 if (hi is not None and lo == 0) else (lo + hi) / 2
        if hi is not None and lo > 0 and (hi - lo) / lo < mp.mpf('1e-3'):
            break
    return lo


allok = True
for name, al, z, q in DATA:
    z = [mp.mpf(x) for x in z]
    z[-1] = 1 - mp.fsum(z[:-1])
    th, H = global_data([mp.mpf(a) for a in al], z, mp.mpf(q))
    e = [max(H) - h for h in H]
    k0, S0 = kap(al, z, q, 0)
    print(f"{name}: H = {[mp.nstr(h, 6) for h in H]}, kappa(0;0) - e = {mp.nstr(max(abs(a - b) for a, b in zip(k0, e)), 3)}, "
          f"explicit validity bound eta_0 = {mp.nstr(eta0(S0), 5)}")
    for ex in range(1, 7):
        eta = mp.mpf(10) ** (-ex)
        try:
            k, _ = kap(al, z, q, eta)
            print(f"   eta = 1e-{ex}: max_p |kappa_p - e_p|/eta = {mp.nstr(max(abs(a - b) for a, b in zip(k, e)) / eta, 6)}")
        except InvalidSpec as ex_:
            print(f"   eta = 1e-{ex}: invalid ({ex_})")
    ratios = {}
    for dex in (2, 3, 4, 5):
        delta = mp.mpf(10) ** (-dex)
        et = largest_eta(al, z, q, delta, e)
        ratios[dex] = et / delta
        print(f"   delta = 1e-{dex}: largest admissible common eta_nb ~ {mp.nstr(et, 5)}  (eta/delta = {mp.nstr(et / delta, 5)})")
    ok = ratios[5] >= ratios[3] / 2
    allok &= ok
    print("   PASS" if ok else "   FAIL")
print("ALL PASS" if allok else "SOME FAIL")
sys.exit(0 if allok else 1)
