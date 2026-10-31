#!/usr/bin/env python3
"""check_construction.py -- runs the proof of Theorem 1 of paper.tex on many values of beta, exact rational arithmetic.

For beta = p/r, block size b (minimal with b(r-p) > p, also b+1 and b+5), the first admissible s among p, 2p, 4p, ...,
and also 7 times that s, the script executes the construction and checks:

 C1  counts: m = m_F + m_g + l, n = N + a + m, ab/n = p/r
 C2  T = tau, i.e. Delta_G(v) = delta - m rho + T = 0 for a block vertex v, with delta and rho recomputed from the general
     closed forms of Lemma 5 (not from the partial fractions)
 C3  the conditions of Lemma 8: 0 <= j_i < Q_i, l <= t <= l(l+1)/2
 C4  Delta_G(h) = delta_h - m rho_h + T > 0 for a hub h (recomputed from Lemma 5)
 C5  Delta_G(z) > 0 for every other non-cut vertex type z of W:
     filler K_{2B}, gadget clique vertices, gadget leaves, tree leaves
 C6  the moduli are pairwise coprime and the denominator of tau divides their product

These are exactly the steps of the proof; the size conditions "for s large" are replaced by their verification on the
given numbers.  Runtime: 30 s, one core.  Environment: Fedora 44 under WSL, Python 3.12.
"""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from fractions import Fraction as Fr
from math import gcd
from hubgraph import *

C = Checker()


def verify(c):
    p, r, b, a, N, n, m, B = c['p'], c['r'], c['b'], c['a'], c['N'], c['n'], c['m'], c['B']
    k, Q, j, l, t, T = c['k'], c['Q'], c['j'], c['l'], c['t'], c['T']
    msgs = []
    # C1
    m_F = B*k*(2*B-1)
    ok1 = m == m_F+c['m_g']+l and n == N+a+m and Fr(a*b, n) == Fr(p, r) and m >= 1
    # C2: recompute delta, rho from Lemma 5
    Y = Fr(a, b)
    d = kf_closed(N, a, a+1, Y)-kf_closed(N-1, a, a+1, Y+Fr(1, B))
    rh = rt_hubless_closed(N-1, a, a+1, Y+Fr(1, B))-rt_hubless_closed(N, a, a+1, Y)
    ok2 = d == c['delta'] and rh == c['rho'] and d-m*rh+T == 0
    # C3
    ok3 = all(0 <= ji < q for ji, q in zip(j, Q)) and l <= t <= l*(l+1)//2 and l >= 0
    # C4
    dh = kf_closed(N, a, a+1, Y)-kf_closed(N, a-1, a+b, Y-Fr(1, b))
    rhh = rt_hubless_closed(N, a-1, a+b, Y-Fr(1, b))-rt_hubless_closed(N, a, a+1, Y)
    ok4 = dh == delta_hub_closed(a, b) and rhh == rho_hub_closed(a, b) and dh > 0 and rhh < 0 and dh-m*rhh+T > 0
    # C5
    RtH = rt_hubless_closed(N, a, a+1, Y)
    ok5 = RtH >= Y
    if k > 0:
        cF = 2*B
        TZ = Fr(2*B-1, B)
        m_R = n-1-(cF-1)
        RtR = RtH+T-TZ
        ok5 = ok5 and 1+m_R*Fr(2, cF*(cF-1))+RtR > 0
    for q, ji in zip(Q, j):
        cc = 2*q
        mZ = 2*q-1+ji
        TZ = 2-Fr(1, q)+ji*(1+Fr(1, q))
        m_R = n-1-mZ
        RtR = RtH+T-TZ
        ok5 = ok5 and RtR > 0
        # clique vertices other than root and y
        ok5 = ok5 and 1+ji+Fr(2*ji, cc*(cc-1))-m_R*Fr(2*(ji-1), cc*(cc-1))+RtR > 0
        # leaves: delta_Z = Rt_Z(leaf) > 0 and rho_Z = -r(root, leaf) < 0, so Delta > RtR > 0: only RtR > 0 is needed
    if l:
        ok5 = ok5 and RtH+T-t > 0     # tree leaves: Delta >= Rt_R(w) > 0
    # C6
    ok6 = True
    D = 1
    for i in range(5):
        D *= Q[i]
        for i2 in range(i+1, 5):
            ok6 = ok6 and gcd(Q[i], Q[i2]) == 1
    ok6 = ok6 and (c['tau']*D).denominator == 1
    return ok1, ok2, ok3, ok4, ok5, ok6


betas = [(1, 100), (1, 7), (1, 3), (1, 2), (3, 5), (2, 3), (7, 10), (3, 4), (4, 5), (5, 6), (9, 10), (19, 20), (99, 100)]
table = []
for (p, r) in betas:
    bmin = 2
    while bmin*(r-p) <= p:
        bmin += 1
    for b in (bmin, bmin+1, bmin+5):
        s = p
        while True:                      # doubling: a small admissible s, not necessarily the smallest
            c = construct(p, r, b, s)
            if c is not None:
                break
            s *= 2
        for ss in (s, 7*s):
            c = construct(p, r, b, ss)
            if c is None:
                C.check(False, 'construction failed at p/r=%d/%d b=%d s=%d' % (p, r, b, ss))
                continue
            res = verify(c)
            C.check(all(res), 'beta=%d/%d b=%d s=%d: n=%d a=%d m=%d, C1-C6 = %s' % (p, r, b, ss, c['n'], c['a'], c['m'], res))
            table.append((Fr(p, r), b, ss, c['x'], c['a'], c['n'], c['k'], c['l'], c['t']))
print('table: beta b s x a n k l t')
for row in table:
    print(*row)
C.finish()
