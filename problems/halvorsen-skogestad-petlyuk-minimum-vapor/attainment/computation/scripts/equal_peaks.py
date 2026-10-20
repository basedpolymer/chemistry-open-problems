#!/usr/bin/env python3
"""equal_peaks.py -- Agent C.  Feed compositions with two (N = 3) or three (N = 4) equal highest peaks
H_1 = H_2 (= H_3), alpha = (4,2,1) resp. (8,4,2,1), q = 1 (N = 3: z_B = z_C; N = 4: z_4 = 1/4 fixed).
Prints z to 40 digits (used by cases.py: n3_equal, n4_equal) and checks the peak equalities.  Numerical only.
Environment: Python 3.12, mpmath.  Runtime: seconds."""
import mpmath as mp
from petnet import Problem
mp.mp.dps = 50


def H(al, z, q):
    return Problem(al, z, q).H


f3 = lambda a: (lambda h: h[0] - h[1])(H([4, 2, 1], [a, (1 - a) / 2, (1 - a) / 2], 1))
lo, hi = mp.mpf('0.05'), mp.mpf('0.95')
for _ in range(170):
    m = (lo + hi) / 2
    if (f3(m) > 0) == (f3(lo) > 0):
        lo = m
    else:
        hi = m
a3 = (lo + hi) / 2
z3 = [a3, (1 - a3) / 2, (1 - a3) / 2]
h3 = H([4, 2, 1], z3, 1)
print("N=3 z =", [mp.nstr(x, 40) for x in z3], "H =", [mp.nstr(x, 20) for x in h3])

z4f = mp.mpf('0.25')


def F4(x):
    z = [x[0], x[1], 1 - z4f - x[0] - x[1], z4f]
    h = H([8, 4, 2, 1], z, 1)
    return mp.matrix([h[0] - h[1], h[1] - h[2]])


best = None
for i in range(1, 30):
    for j in range(1, 30):
        a, b = mp.mpf(i) / 40, mp.mpf(j) / 40
        if a + b >= 0.74:
            continue
        r = F4([a, b])
        m = abs(r[0]) + abs(r[1])
        if best is None or m < best[0]:
            best = (m, a, b)
x = mp.matrix([best[1], best[2]])
for it in range(100):
    Fx = F4(x)
    if mp.norm(Fx) < mp.mpf(10) ** -45:
        break
    J = mp.matrix(2, 2)
    h = mp.mpf(10) ** -20
    for k in range(2):
        xp = x.copy(); xp[k] += h
        Fp = F4(xp)
        for r in range(2):
            J[r, k] = (Fp[r] - Fx[r]) / h
    dx = mp.lu_solve(J, -Fx)
    lam = mp.mpf(1)
    while True:
        xn = x + lam * dx
        if xn[0] > 0 and xn[1] > 0 and 1 - z4f - xn[0] - xn[1] > 0 and mp.norm(F4(xn)) < mp.norm(Fx):
            break
        lam /= 2
        if lam < 1e-12:
            break
    x = xn
z4 = [x[0], x[1], 1 - z4f - x[0] - x[1], z4f]
h4 = H([8, 4, 2, 1], z4, 1)
print("grid best", mp.nstr(best[0], 4))
print("N=4 z =", [mp.nstr(v, 40) for v in z4], "H =", [mp.nstr(v, 20) for v in h4])
ok = abs(h3[0] - h3[1]) < 1e-40 and max(abs(h4[0] - h4[1]), abs(h4[1] - h4[2])) < 1e-40 and min(z3 + z4) > 0
print("PASS" if ok else "FAIL")
