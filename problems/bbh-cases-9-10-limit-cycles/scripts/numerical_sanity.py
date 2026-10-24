#!/usr/bin/env python3
"""Floating-point sanity checks. NOT part of the proof (the proof is result.md + certify_lienard.py).

Environment: Python 3.12/3.13 + numpy + scipy (DOP853, quad, brentq). Single thread. Runtime: about 3-5 minutes.
Deterministic (fixed seeds). Prints PASS/FAIL per claim, exits 1 on failure.

N1  Backward integration, in the ORIGINAL mass-action coordinates (no normalisation), from a point close to the
    positive equilibrium, for random integer networks of the families 9a, 9b, 10a and random rate constants:
    the backward orbit never converges to a periodic orbit (it leaves the closed quadrant, escapes, or tends to
    another limit set). A periodic orbit around the (stable) equilibrium would trap it.
N2  Controls (the detector works): for BBH's Hopf networks just past the supercritical Hopf point, forward
    integration converges to a periodic orbit (successive returns to a section agree to 1e-7).
N3  The Liénard comparison of result.md, computed by quadrature from the normalised system: F(a) - F(b) > 0 at
    G-paired points a < 1 < b for random non-Hopf networks and rho, and < 0 somewhere for Hopf networks.
"""
import sys
import time
import numpy as np
from scipy.integrate import DOP853, quad
from scipy.optimize import brentq

fails = 0
T0 = time.time()


def check(cond, msg):
    global fails
    print(("PASS " if cond else "FAIL ") + msg, flush=True)
    fails += (not cond)


SRC = {9: [(2, 0), (1, 1), (0, 1)], 10: [(2, 0), (1, 1), (0, 0)]}


def field(case, c, d, kap):
    (a1, b1), (a2, b2), (a3, b3) = SRC[case]

    def F(z):
        xx, yy = z
        m = (kap[0] * xx ** a1 * yy ** b1, kap[1] * xx ** a2 * yy ** b2, kap[2] * xx ** a3 * yy ** b3)
        return np.array([c[0] * m[0] + c[1] * m[1] + c[2] * m[2], d[0] * m[0] + d[1] * m[1] + d[2] * m[2]])
    return F


def equilibrium(case, c, d, kap):
    u1 = c[1] * d[2] - c[2] * d[1]; u2 = c[2] * d[0] - c[0] * d[2]; u3 = c[0] * d[1] - c[1] * d[0]
    k1, k2, k3 = kap
    if case == 9:
        xb = k3 * u2 / (k2 * u3); yb = k1 * k3 * u2 ** 2 / (k2 ** 2 * u1 * u3)
    else:
        mu = k3 / u3; xb = np.sqrt(mu * u1 / k1); yb = mu * u2 / (k2 * xb)
    return xb, yb


def sections(F, z0, sign, xb, yb, smax=4e4, nmax=400):
    """integrate sign*F/(1+|z|^2/s^2) (same orbits), record crossings of {y = yb, x > xb}."""
    sc = max(xb, yb)

    def G(t, z):
        return sign * F(z) / (1.0 + (z[0] ** 2 + z[1] ** 2) / sc ** 2)
    sol = DOP853(G, 0.0, np.array(z0, float), smax, rtol=1e-11, atol=1e-13 * sc, max_step=0.05)
    xs, zprev, status = [], sol.y.copy(), "time"
    while sol.status == "running":
        sol.step()
        z = sol.y
        if not np.all(np.isfinite(z)):
            status = "nan"; break
        if abs(z[0]) + abs(z[1]) > 1e7 * sc:
            status = "escape"; break
        if min(z[0], z[1]) < 0:
            status = "left quadrant"; break
        if (zprev[1] - yb) * (z[1] - yb) < 0 and z[0] > xb and zprev[0] > xb:
            dd = sol.dense_output()
            tc = brentq(lambda s: dd(s)[1] - yb, sol.t_old, sol.t, xtol=1e-14)
            xs.append(dd(tc)[0])
            if len(xs) >= nmax:
                status = "nmax"; break
        zprev = z.copy()
    return status, xs


def converged(xs, xb, tol=1e-7):
    if len(xs) < 12:
        return False
    tail = np.array(xs[-6:])
    return np.all(np.abs(np.diff(tail)) < tol * tail[-1]) and tail[-1] > xb * (1 + 1e-3)


def random_network(rng, fam):
    while True:
        c1 = int(rng.integers(1, 7)); c3 = int(rng.integers(1, 7)); d1 = int(rng.integers(1, 20))
        if fam == "9a":
            d3 = int(rng.integers(-1, c3)); c, d = (c1, -1, c3), (d1, -1, d3)
            if d1 > c1 and 2 * c1 * c3 - c1 * d3 - c3 * d1 <= 0: return 9, c, d
        elif fam == "9b":
            c, d = (c1, -1, c3), (d1, 0, -1)
            if c3 * d1 >= c1: return 9, c, d
        else:
            d3 = int(rng.integers(0, c3)); c, d = (c1, -1, c3), (d1, -1, d3)
            if d1 > c1 and 2 * c1 * c3 - c1 * d3 - c3 * d1 <= 0: return 10, c, d


# ---------------------------------------------------------------- N1
rng = np.random.default_rng(20261024)
nets = []
for fam in ("9a", "9b", "10a"):
    for _ in range(4):
        nets.append(random_network(rng, fam))
# boundary networks (k = 1): the BBH example and its relatives
nets += [(9, (1, -1, 1), (2, -1, 0)), (9, (1, -1, 1), (1, 0, -1)), (10, (1, -1, 1), (2, -1, 0)), (9, (2, -1, 3), (5, -1, -1))]
ncycle, nrun, stat = 0, 0, {}
for (case, c, d) in nets:
    for _ in range(8):
        kap = tuple(float(np.exp(rng.uniform(-4, 4))) for _ in range(3))
        F = field(case, c, d, kap)
        xb, yb = equilibrium(case, c, d, kap)
        st, xs = sections(F, [xb * (1 + 1e-3), yb], -1.0, xb, yb)
        nrun += 1
        stat[st] = stat.get(st, 0) + 1
        if converged(xs, xb):
            ncycle += 1
            print("   cycle-like backward limit:", case, c, d, kap, xs[-3:])
print("   N1 outcomes:", stat)
check(ncycle == 0, f"N1 {nrun} backward integrations (16 non-Hopf networks x 8 random kappa): no periodic orbit found")

# ---------------------------------------------------------------- N2 controls
ok2 = True
# BBH (4): 2X -> 3X + Y, X + Y -> Y, Y -> 0 : xdot = k1 x^2 - k2 x y, ydot = k1 x^2 - k3 y, Hopf at k1 = k2
for (case, c, d, kap) in [(9, (1, -1, 0), (1, 0, -1), (1.02, 1.0, 1.0)),
                          (9, (2, -1, 1), (3, -1, 0), None),
                          (10, (2, -1, 1), (3, -1, 0), None)]:
    if kap is None:
        # normalised Hopf point: tr J = 0 at rho_H; choose rho slightly below (unstable focus), map back to kappa
        u1 = c[1] * d[2] - c[2] * d[1]; u2 = c[2] * d[0] - c[0] * d[2]; u3 = c[0] * d[1] - c[1] * d[0]
        # find rho_H numerically from the normalised Jacobian (tr J decreases in rho)
        def trJ(r):
            m3 = (lambda X, Y: Y) if case == 9 else (lambda X, Y: 1.0)
            h = 1e-6
            def Fn(X, Y):
                M = (X * X, X * Y, m3(X, Y))
                return np.array([sum(c[i] * [u1, u2, u3][i] * M[i] for i in range(3)),
                                 r * sum(d[i] * [u1, u2, u3][i] * M[i] for i in range(3))])
            J = np.column_stack([(Fn(1 + h, 1) - Fn(1 - h, 1)) / (2 * h), (Fn(1, 1 + h) - Fn(1, 1 - h)) / (2 * h)])
            return J[0, 0] + J[1, 1]
        rH = brentq(trJ, 1e-6, 1e3)
        r = 0.97 * rH
        # normalised system corresponds to kappa = u_i / M_i(xb, yb) with xb = r, yb = 1 (rho = xb/yb), mu = 1
        xb, yb = r, 1.0
        kap = (u1 / xb ** 2, u2 / (xb * yb), u3 / (yb if case == 9 else 1.0))
    F = field(case, c, d, kap)
    xb, yb = equilibrium(case, c, d, kap)
    st, xs = sections(F, [xb * 1.05, yb], +1.0, xb, yb, nmax=600, smax=2e5)
    conv = converged(xs, xb, tol=1e-7)
    print(f"   control {case} {c} {d} kappa={tuple(round(v, 4) for v in kap)}: {st}, last returns {xs[-2:] if xs else None}")
    ok2 &= bool(conv)
check(ok2, "N2 controls: past the supercritical Hopf point of three BBH Hopf networks, forward orbits converge to a cycle")


# ---------------------------------------------------------------- N3 Liénard comparison by quadrature
def lienard_min(case, c, d, rho, npts=40):
    u1 = c[1] * d[2] - c[2] * d[1]; u2 = c[2] * d[0] - c[0] * d[2]; u3 = c[0] * d[1] - c[1] * d[0]
    al1, al2, al3 = c[0] * u1, -c[1] * u2, c[2] * u3
    be1, be2, be3 = d[0] * u1, -d[1] * u2, d[2] * u3
    if case == 9:
        A = lambda s: al1 * s * s; B = lambda s: al3 - al2 * s
        C = lambda s: be1 * s * s; D = lambda s: be3 - be2 * s
        dA = lambda s: 2 * al1 * s; dB = lambda s: -al2; x0 = al3 / al2
    else:
        A = lambda s: al1 * s * s + al3; B = lambda s: -al2 * s
        C = lambda s: be1 * s * s + be3; D = lambda s: -be2 * s
        dA = lambda s: 2 * al1 * s; dB = lambda s: -al2; x0 = 0.0
    f = lambda s: -(dA(s) * B(s) - A(s) * dB(s)) / B(s) ** 2 - rho * D(s) / B(s)
    g = lambda s: -rho * (B(s) * C(s) - A(s) * D(s)) / B(s) ** 2
    Fq = lambda s: quad(f, 1, s, limit=200)[0]
    Gq = lambda s: quad(g, 1, s, limit=200)[0]
    worst = np.inf
    for aa in x0 + (1 - x0) * np.linspace(0.03, 0.995, npts):
        Ga = Gq(aa); hi = 2.0
        while Gq(hi) < Ga: hi *= 2
        bb = brentq(lambda s: Gq(s) - Ga, 1 + 1e-12, hi, xtol=1e-13)
        worst = min(worst, Fq(aa) - Fq(bb))
    return worst


ok3 = True
for (case, c, d) in nets[:12]:
    rho = float(np.exp(rng.uniform(-5, 3)))
    w = lienard_min(case, c, d, rho)
    ok3 &= bool(w > 0)
    print(f"   non-Hopf {case} {c} {d} rho={rho:.3g}: min F(a)-F(b) = {w:.3e}")
check(ok3, "N3 non-Hopf networks: F(a) > F(b) at all sampled G-paired points")
ok3b = True
for (case, c, d) in [(9, (2, -1, 1), (3, -1, 0)), (9, (2, -1, 1), (1, 0, -1)), (10, (2, -1, 1), (3, -1, 0))]:
    w = lienard_min(case, c, d, 0.01)
    ok3b &= bool(w < 0)
    print(f"   Hopf     {case} {c} {d} rho=0.01: min F(a)-F(b) = {w:.3e}")
check(ok3b, "N3 controls: for Hopf networks (small rho) the comparison fails somewhere, as it must")
print(f"runtime {time.time() - T0:.0f} s")
print("ALL PASS" if not fails else f"{fails} FAILURES")
sys.exit(1 if fails else 0)
