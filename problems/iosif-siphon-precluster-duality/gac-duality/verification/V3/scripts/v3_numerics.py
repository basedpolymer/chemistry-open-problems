"""Verifier V3 (Claude Opus 5.5, effort xhigh), 27 October 2026: numerical illustrations (NOT proofs, NOT
certificates) of hand-proved extinction statements, to catch sign slips.  scipy LSODA, rtol 1e-10, atol 1e-14.

  N1  G's D (unit rates) from (15/8, 1/8): v(t) <= (1/8) exp(-721 t/4096) along the run.
  N2  D2's example U, the three regimes kappa2 > 2kappa1, < 2kappa1, = 2kappa1, small class T: the solution started
      near the vertex (T,0,0) in the regime's region approaches the vertex (k2 + k3 decreasing to ~0).
  N3  Q1's example X dual with kappa2 T^2 < kappa1: k2 -> 0; with kappa2 T^2 > kappa1: k2 -> T - sqrt(kappa1/kappa2).
Usage: python v3_numerics.py
"""
import sys
import time
import numpy as np
from scipy.integrate import solve_ivp

fails = 0


def claim(name, ok, msg=""):
    global fails
    print(("PASS" if ok else "FAIL") + ": " + name + ((" -- " + msg) if msg else ""))
    fails += 0 if ok else 1


def run(f, x0, T):
    return solve_ivp(f, (0, T), x0, method="LSODA", rtol=1e-10, atol=1e-14, dense_output=True)


def main():
    t0 = time.time()
    # N1
    s = run(lambda t, x: [x[1] - x[0]**3 * x[1]**2, x[0]**3 * x[1]**2 - x[1]], [15 / 8, 1 / 8], 200)
    ts = np.linspace(0, 200, 2001)
    v = s.sol(ts)[1]
    claim("N1 D from (15/8,1/8): v below (1/8)exp(-721t/4096), v(200) = %.3e" % v[-1],
          np.all(v <= 0.125 * np.exp(-721 * ts / 4096) + 1e-12))

    # N2: dual of U: 2k3 -> 2k1 (c1), k2+k3 -> k1+k3 (c2), 3k1+2k2 -> 3k2+2k3 (c3)
    def fU(c1, c2, c3):
        def f(t, x):
            a, b, c = x
            r1 = c1 * c**2
            r2 = c2 * b * c
            r3 = c3 * a**3 * b**2
            return [2 * r1 + r2 - 3 * r3, -r2 + r3, -2 * r1 + 2 * r3]
        return f
    for (c1, c2, c3), T, start in (((1, 3, 1), 1.0, (1e-3, 1e-2)), ((1, 1, 1), 0.5, (1e-3, 1e-3)),
                                   ((1, 2, 1), 0.5, (1e-3, 2e-3))):
        a0, b0 = start
        x0 = [T - a0 - b0, a0, b0]
        # logarithmic coordinates (relative precision near the vertex), long horizon, monotone decay of k2 + k3
        fu = fU(c1, c2, c3)
        gl = lambda t, y: [fi / xi for fi, xi in zip(fu(t, np.exp(y)), np.exp(y))]
        s = solve_ivp(gl, (0, 1e9), np.log(x0), method="Radau", rtol=1e-11, atol=1e-11, dense_output=True)
        ts = np.concatenate([np.linspace(0, 10, 1001), np.logspace(1, 9, 3001)])
        X = np.exp(s.sol(ts))
        m = X[1] + X[2]
        mono = np.all(np.diff(m[200:]) <= 1e-12 * m[200:-1])
        claim("N2 U kappa=(%g,%g,%g), T=%g, start k2=%g k3=%g (log coordinates, t <= 1e9): k2+k3 decreases "
              "monotonically after t = 2 to %.2e (vertex attracting)" % (c1, c2, c3, T, a0, b0, m[-1]),
              s.success and mono and m[-1] < 1e-5 * (a0 + b0))

    # N3
    def fX(c1, c2):
        return lambda t, x: [c1 * x[1] - c2 * x[0]**2 * x[1], -c1 * x[1] + c2 * x[0]**2 * x[1]]
    s = run(fX(1, 1), [0.5, 0.4], 200)
    claim("N3 X* kappa=(1,1), T=0.9 (< 1): k2 -> 0 (k2(200) = %.2e)" % s.y[1, -1], s.y[1, -1] < 1e-6)
    s = run(fX(1, 1), [0.5, 1.5], 200)
    claim("N3 X* kappa=(1,1), T=2 (> 1): k2 -> T - 1 = 1 (k2(200) = %.6f)" % s.y[1, -1], abs(s.y[1, -1] - 1) < 1e-6)
    print("# runtime %.1f s (numerical evidence only)" % (time.time() - t0))
    sys.exit(1 if fails else 0)


if __name__ == "__main__":
    main()
