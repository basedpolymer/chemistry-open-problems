#!/usr/bin/env python3
"""
simulate_example.py -- illustration only (floating point, not part of any proof).
Integrates the E1 network G7 with kappa = (1, 1, 5, 4, 1, 6):
    X+Y -> 2Y, Y -> 2Z, 0 -> X+Z, X+Z -> 2Z, 2Z -> X+Z, Y+Z -> Y
from two initial conditions and shows convergence to the two exponentially stable steady states
(1, 2, 1) and (1, 10/9, 5/3); both have x = 1 (ACR in X).
Runtime: about 1 second.  Usage: python simulate_example.py
"""
import numpy as np
from scipy.integrate import solve_ivp

k1, k2, k3, k4, k5, k6 = 1, 1, 5, 4, 1, 6


def f(t, u):
    x, y, z = u
    return [-k1 * x * y + k3 - k4 * x * z + k5 * z * z,
            k1 * x * y - k2 * y,
            2 * k2 * y + k3 + k4 * x * z - k5 * z * z - k6 * y * z]


targets = {"(1, 2, 1)": np.array([1, 2, 1]), "(1, 10/9, 5/3)": np.array([1, 10 / 9, 5 / 3])}
ok = True
for u0 in ([1.0, 2.5, 0.8], [1.0, 1.0, 1.8]):
    sol = solve_ivp(f, (0, 2000), u0, method="LSODA", rtol=1e-10, atol=1e-12)
    end = sol.y[:, -1]
    name = min(targets, key=lambda k: np.linalg.norm(end - targets[k]))
    dist = np.linalg.norm(end - targets[name])
    print("initial condition %s -> (%.8f, %.8f, %.8f), distance to %s: %.2e" % (u0, *end, name, dist))
    ok = ok and dist < 1e-6
print("illustration:", "two trajectories converge to the two stable steady states" if ok else "check failed")
