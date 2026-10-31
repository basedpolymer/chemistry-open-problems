#!/usr/bin/env python3
"""Exploration: smallest graphs with EXACT proportion p/r (root w). W = prime-power gadgets + integer fillers + tree."""
import sys, math
from fractions import Fraction as Fr
from sympy import factorint
sys.path.insert(0, '.')
from x12_rootw import dr, gadget
FILL = [(6, 3, 'K4x2'), (15, 5, 'K6x3'), (28, 7, 'K8x4'), (45, 9, 'K10x5'), (66, 11, 'K12x6'), (3,2,'K3... '),]
FILL = [(6, 3, (4, 2)), (15, 5, (6, 3)), (28, 7, (8, 4)), (45, 9, (10, 5)), (66, 11, (12, 6)), (91, 13, (14, 7)), (120, 15, (16, 8))]
def solve_exact(a, b, m):
    d, r = dr(a, b)
    tau = m*r-d
    if tau <= 0: return None
    D = tau.denominator
    fac = factorint(D)
    qs = [p**e for p, e in fac.items()]
    A = tau.numerator
    T_g = Fr(0); m_g = 0; gad = []
    for q in qs:
        c = A*pow(D//q, -1, q) % q
        if c == 0: continue
        k, j, Tq, mq = gadget(q, c)
        T_g += Tq; m_g += mq; gad.append((k, j))
    if m_g > m: return None
    mr = m-m_g; Tr = tau-T_g
    assert Tr.denominator == 1
    Tr = int(Tr)
    if Tr < 0: return None
    # one filler type f times + tree
    if mr <= Tr <= mr*(mr+1)//2 and (mr > 0 or Tr == 0):
        return gad, (0, 0, None), (mr, Tr)
    for (wv, wT, tag) in FILL:
        # want l = mr - f wv, t = Tr - f wT, l <= t <= l(l+1)/2
        # t - l = Tr - mr + f (wv - wT) >= 0
        need = max(0, -(Tr-mr))
        f0 = (need+(wv-wT)-1)//(wv-wT)
        for f in range(f0, f0+3):
            l = mr-f*wv; t = Tr-f*wT
            if l >= 0 and t >= 0 and l <= t <= l*(l+1)//2 and (l > 0 or t == 0):
                return gad, (f, tag, (wv, wT)), (l, t)
    return None
if __name__ == '__main__':
    p = int(sys.argv[1]); r = int(sys.argv[2]); b = int(sys.argv[3]); amax = int(sys.argv[4])
    cnt = 0
    for a in range(2, amax+1):
        if (a*b*r) % p: continue
        n = a*b*r//p
        m = n-a*b-1-a
        if m < 0: continue
        s = solve_exact(a, b, m)
        if s:
            print('p/r', p, r, 'n', n, 'b', b, 'a', a, 'm', m, s, flush=True)
            cnt += 1
            if cnt >= 4: break
