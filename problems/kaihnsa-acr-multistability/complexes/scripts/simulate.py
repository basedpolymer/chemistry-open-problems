"""Illustration only (floating point, not part of the proof; runtime about 1 second; numpy + scipy).

Integrates the mass-action ODEs of N1, N2, N3 from two initial conditions each and shows convergence to the two
exponentially stable steady states. PASS/FAIL (distance < 1e-6 at the final time); exit code 1 on failure.
"""
import sys
import time
import numpy as np
from scipy.integrate import solve_ivp

T0 = time.time()


def fN1(t, u):
    x, y, z = u
    return [8 * y * z - 10 * x * z - 14 * x * x,
            5 * z - 4 * y * z + 7 * x * x,
            -z * z + y - 4 * y * z + 10 * x * z]


def fN3(t, u):
    x, y, z = u
    return [6 * y + 2 * y * z - 7 * x * z - x * x,
            z * z - 6 * y - 5 * y * z + 10 * x * z,
            -z * z + 3 * y * z + x - 3 * x * z]


def fN2(t, u):
    x, y, z = u
    return [z * z + 7 * y * z - 5 * x * z - 6 * x * x + 3 * x,
            -2 * y * z - x * y + 3 * x * x,
            -z * z - 7 * y * z + 5 * x * z + 3 * x * x]


CASES = [("N1", fN1, [0.5, 5.0, 0.2], [1.0, 19 / 4, 0.5]), ("N1", fN1, [2.0, 1.0, 5.0], [1.0, 7 / 4, 3.5]),
         ("N2", fN2, [0.5, 2.0, 0.2], [1.0, 1.5, 0.5]), ("N2", fN2, [2.0, 0.2, 4.0], [1.0, 3 / 7, 3.0]),
         ("N3", fN3, [0.5, 0.5, 0.5], [1.0, 1.0, 1.0]), ("N3", fN3, [2.0, 2.0, 3.0], [1.0, 1.5, 2.0])]
FAIL = 0
for name, f, x0, target in CASES:
    sol = solve_ivp(f, (0, 3000), x0, method="LSODA", rtol=1e-11, atol=1e-13)
    end = sol.y[:, -1]
    d = np.max(np.abs(end - np.array(target)))
    ok = d < 1e-6
    FAIL += not ok
    print(("PASS" if ok else "FAIL") + " %s: from %s the trajectory ends at %s (target %s, distance %.2e)"
          % (name, x0, np.round(end, 8), np.round(target, 8), d))
print("runtime %.2f s" % (time.time() - T0))
sys.exit(1 if FAIL else 0)
