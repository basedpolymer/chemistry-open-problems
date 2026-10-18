#!/usr/bin/env python3
"""Numerical sanity checks (NOT a proof). Verifier's own code.

Environment: WSL Fedora, Python 3, numpy, scipy. Single thread. Runtime: a few minutes.

Scaling. With x = (k3/k2) X, y = (k3/k2) Y, t = T/k3, the subfamily
    xdot = c1 k1 x^2 - k2 x y + c3 k3 y,  ydot = d1 k1 x^2 - k2 x y
becomes  X' = c1 L X^2 - X Y + c3 Y,  Y' = d1 L X^2 - X Y  with the single parameter L = k1/k2.
(The BBH example is c1 = c3 = 1, d1 = 2.)  Positive equilibrium (d1 > c1):
    X* = c3 d1 / (d1 - c1),  Y* = d1 L X*  (from Y'=0: Y = d1 L X).
Hopf curve (from the symbolic check V8): k2 = k1 (2c1 - d1), i.e. L_H = 1/(2c1 - d1) when c1 < d1 < 2c1;
tr J > 0 iff L > L_H.

Part A: example and other d1 >= 2c1 members, L on a log grid, many initial conditions: every bounded
        trajectory must converge to the equilibrium (no other bounded omega-limit set observed).
Part B: Hopf region c1 < d1 < 2c1: for L slightly above L_H a stable limit cycle must appear
        (supercritical Hopf, BBH Theorem 9), and for L <= L_H trajectories converge to the equilibrium,
        consistently with the Dulac bound (Delta(b*) <= 0 iff L <= L_H).
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
import sys
import numpy as np
from scipy.integrate import solve_ivp

fails = 0


def check(cond, msg):
    global fails
    print(("PASS " if cond else "FAIL ") + msg, flush=True)
    fails += (not cond)


def rhs_factory(c1, c3, d1, L):
    def rhs(t, z):
        X, Y = z
        return [c1 * L * X * X - X * Y + c3 * Y, d1 * L * X * X - X * Y]
    return rhs


def escape(t, z):
    return 1e7 - max(abs(z[0]), abs(z[1]))
escape.terminal = True


def run(c1, c3, d1, L, z0, T):
    sol = solve_ivp(rhs_factory(c1, c3, d1, L), (0, T), z0, method="LSODA", rtol=1e-9, atol=1e-12,
                    events=escape, dense_output=False)
    return sol


def classify(c1, c3, d1, L, z0, T=4000.0):
    Xs = c3 * d1 / (d1 - c1)
    Ys = d1 * L * Xs
    sol = run(c1, c3, d1, L, z0, T)
    if sol.status == 1:
        return "escape", sol
    zf = sol.y[:, -1]
    if np.hypot(zf[0] - Xs, (zf[1] - Ys)) < 1e-4 * (1 + np.hypot(Xs, Ys)):
        return "equilibrium", sol
    # not yet converged: extend once
    sol2 = run(c1, c3, d1, L, zf, 10 * T)
    if sol2.status == 1:
        return "escape", sol2
    zf = sol2.y[:, -1]
    if np.hypot(zf[0] - Xs, (zf[1] - Ys)) < 1e-4 * (1 + np.hypot(Xs, Ys)):
        return "equilibrium", sol2
    return "other", sol2


# ---------------- Part A
rng = np.random.default_rng(12345)
for (c1, c3, d1) in [(1, 1, 2), (1, 1, 3), (1, 2, 2), (2, 1, 4), (2, 3, 5), (1, 3, 6)]:
    counts = {"equilibrium": 0, "escape": 0, "other": 0}
    for L in np.logspace(-2.5, 2.5, 21):
        Xs = c3 * d1 / (d1 - c1); Ys = d1 * L * Xs
        for _ in range(12):
            z0 = [Xs * 10 ** rng.uniform(-2, 2), Ys * 10 ** rng.uniform(-2, 2)]
            kind, _ = classify(c1, c3, d1, L, z0)
            counts[kind] += 1
            if kind == "other":
                print(f"   unclassified: c1={c1} c3={c3} d1={d1} L={L:.4g} z0={z0}")
    check(counts["other"] == 0, f"A (c1,c3,d1)=({c1},{c3},{d1}), d1>=2c1: {counts} (no bounded non-equilibrium limit set)")

# ---------------- Part B: Hopf region
for (c1, c3, d1) in [(2, 1, 3), (3, 1, 4), (3, 2, 5)]:
    LH = 1.0 / (2 * c1 - d1)
    Xs = c3 * d1 / (d1 - c1)
    # below/at the Hopf value: convergence (Dulac excludes cycles there)
    for L in [0.9 * LH, 0.98 * LH]:
        Ys = d1 * L * Xs
        kind, _ = classify(c1, c3, d1, L, [Xs * 1.05, Ys * 0.95])
        check(kind == "equilibrium", f"B (c1,c3,d1)=({c1},{c3},{d1}) L={L/LH:.2f} L_H: trajectory near equilibrium converges ({kind})")
    # above: small stable limit cycle; measure its amplitude from the inside and from the outside
    amps = []
    for L in [1.01 * LH, 1.04 * LH]:
        Ys = d1 * L * Xs
        res = []
        for start in (1.001, 1.3):
            sol = solve_ivp(rhs_factory(c1, c3, d1, L), (0, 6000), [Xs * start, Ys], method="LSODA",
                            rtol=1e-10, atol=1e-12, t_eval=np.linspace(5000, 6000, 20001), events=escape)
            if sol.status == 1 or sol.y.shape[1] == 0:
                res.append(np.nan); continue
            res.append(sol.y[0].max() - sol.y[0].min())
        amps.append(res)
        ok = np.all(np.isfinite(res)) and min(res) > 1e-3 and abs(res[0] - res[1]) < 1e-3 * max(res) + 1e-6
        check(ok, f"B (c1,c3,d1)=({c1},{c3},{d1}) L={L/LH:.2f} L_H: periodic attractor, X-amplitude from inside/outside = {res}")
    # supercritical scaling amplitude ~ sqrt(L - L_H): ratio ~ 2
    if all(np.isfinite(a).all() for a in amps):
        ratio = amps[1][0] / amps[0][0]
        check(1.6 < ratio < 2.4, f"B (c1,c3,d1)=({c1},{c3},{d1}) amplitude ratio (L-L_H x4) = {ratio:.3f} (~2 for supercritical Hopf)")

print("ALL PASS" if not fails else f"{fails} FAILURES")
sys.exit(1 if fails else 0)
