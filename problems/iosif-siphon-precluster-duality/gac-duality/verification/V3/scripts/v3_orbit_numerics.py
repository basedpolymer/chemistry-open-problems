"""Verifier V3 (Claude Opus 5.5, effort xhigh), 27 October 2026: numerical illustration (NOT a proof) of A2's
Theorem O: dual of 3X3 -> X1+X2 -> X2+2X3+3X4 -> 3X3 at rates (17, 400, 70000, 1e6), class T = 2.01, started at
(1, 0.99, 0.02).  Integrated in logarithmic coordinates y = log x (so that the very small values reached near the
vertex (T,0,0) are resolved in relative precision), scipy Radau, rtol 1e-11, atol 1e-11 on y.
Reports: minima of log10 c and log10 b per loop (loops = upward crossings of b = 1), the return times, and the
values of V = c exp(-beta b/alpha) at the crossings.  Usage: python v3_orbit_numerics.py [tmax]"""
import sys
import time
import numpy as np
from scipy.integrate import solve_ivp

al, be, ga, de = 17.0, 400.0, 70000.0, 1e6

def g(t, y):
    a, b, c = np.exp(y)
    r1, r2, r3, r4 = al * b, be * b * c, ga * a**3 * c**2, de * c**3
    return [(r1 + r2 - 3 * r3) / a, (-r1 + 2 * r3 + 3 * r4) / b, (-r2 + r3 - 3 * r4) / c]

def up(t, y):
    return y[1]
up.direction = 1

tmax = float(sys.argv[1]) if len(sys.argv) > 1 else 2e5
t0 = time.time()
sol = solve_ivp(g, (0, tmax), np.log([1.0, 0.99, 0.02]), method="Radau", rtol=1e-11, atol=1e-11, events=up,
                dense_output=True)
cr = sol.t_events[0]
print("success", sol.success, "steps", len(sol.t), "tmax", tmax)
print("upward crossings of b = 1:", len(cr))
prev = 0.0
for k, tc in enumerate(cr[:40]):
    ts = np.linspace(prev, tc, 20001)
    Y = sol.sol(ts)
    a, b, c = np.exp(Y)
    print("loop %2d: t = %.6g, period %.6g, min log10 c = %.2f, min log10 b = %.2f, max a = %.4f, sum-T = %.1e" % (
        k, tc, tc - prev, np.log10(c).min(), np.log10(b).min(), a.max(), abs(a + b + c - 2.01).max()))
    prev = tc
yend = sol.y[:, -1]
print("final point (t=%.3g): a=%.6g b=%.6g c=%.6g" % (sol.t[-1], *np.exp(yend)))
print("# runtime %.1f s (numerical illustration only)" % (time.time() - t0))
