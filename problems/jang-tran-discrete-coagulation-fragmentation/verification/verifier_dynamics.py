#!/usr/bin/env python3
"""Independent verifier checks of the time-dependent argument (floating point: EVIDENCE, not proof).

Written by the independent verifier; does not reuse the solvers' code. Different integrators from the
solvers' scripts: classical fixed-step RK4 for the cluster systems, implicit Radau for characteristics.
Initial datum with INFINITE second moment: rho0_j = c j^{-5/2}, c = m / zeta(3/2), at m = 0.9 and m = 1.

  D1  truncated systems (L = 250, 500, 1000, 2000; tail mass put at size L) stay >= 0 and conserve mass;
      their sizes <= 20 converge, as L grows, to the closed lower-triangular system (paper, eq. (closed)),
      which presupposes mass conservation of the limit (any mass loss would change the coagulation loss
      and the number of clusters); mass carried by sizes > L/2 decreases with L.
  D2  characteristic solution F of the auxiliary HJ equation (Radau, label found by bisection): invariant
      region 0 <= P <= m, U > 0, W >= 0, barrier Gamma = tau F_xx > -1, J > 0; comparison F^L >= F at
      several (x, t), with F^L - F decreasing in L (consistent with F^rho = F, paper's Remark).
  D3  m = 1: F(x, t)/tau(x) increases towards 1 as x -> 0 at t = 1 and t = 3 (slope lemma (a)).

Environment: Python 3, numpy, scipy, mpmath. Runtime: about 7-10 minutes, one thread.
Deterministic. Prints PASS/FAIL per check, exits 1 on failure.
"""
import os
for k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[k] = "1"
import sys
import numpy as np
import mpmath as mp
from scipy.integrate import solve_ivp
from scipy.signal import fftconvolve

FAILS = 0


def report(name, ok, info=""):
    global FAILS
    print(("PASS " if ok else "FAIL ") + name + (("  [" + info + "]") if info else ""), flush=True)
    if not ok:
        FAILS += 1


T_END = 3.0
T_CHECK = [1.0, 3.0]


# ------------------------------------------------------------------ cluster systems (RK4)
def rhs_trunc(rho, L, l):
    a = l * rho
    gain = np.zeros(L)
    gain[1:] = 0.5 * fftconvolve(a, a)[:L - 1]
    cum = np.cumsum(a)                       # cum[k-1] = sum_{j<=k} j rho_j
    S = np.zeros(L)
    S[:L - 1] = cum[L - 2::-1]               # S[n-1] = sum_{j<=L-n} j rho_j
    tail = np.concatenate([np.cumsum(rho[::-1])[::-1][1:], [0.0]])
    return gain - a * S - 0.5 * (l - 1) * rho + tail


def rk4(f, y0, T, dt, save):
    y = y0.copy()
    t = 0.0
    out = {}
    nsteps = int(round(T / dt))
    save_steps = {int(round(s / dt)): s for s in save}
    for k in range(1, nsteps + 1):
        k1 = f(t, y)
        k2 = f(t + dt / 2, y + dt / 2 * k1)
        k3 = f(t + dt / 2, y + dt / 2 * k2)
        k4 = f(t + dt, y + dt * k3)
        y = y + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        t = k * dt
        if k in save_steps:
            out[save_steps[k]] = y.copy()
    return out


def heavy_data(m, L):
    c = m / float(mp.zeta(1.5))
    j = np.arange(1, L + 1, dtype=float)
    r = c * j ** -2.5
    tail_mass = c * float(mp.zeta(1.5, L + 1))          # sum_{j>L} j rho0_j (Hurwitz zeta)
    r[-1] += tail_mass / L
    return r


def closed_solution(m, Lc=40):
    c = m / float(mp.zeta(1.5))
    l = np.arange(1, Lc + 1, dtype=float)
    r0 = c * l ** -2.5
    n00 = c * float(mp.zeta(2.5))
    cn = m * l + (l + 1) / 2

    def f(t, rho):
        a = l * rho
        gain = np.zeros(Lc)
        gain[1:] = 0.5 * np.convolve(a, a)[:Lc - 1]
        n0 = m * (1 - m) + np.exp(-t / 2) * (n00 - m * (1 - m))
        partial = np.concatenate([[0.0], np.cumsum(rho)[:-1]])
        return gain - cn * rho + n0 - partial
    return rk4(f, r0, T_END, 1e-3, T_CHECK)


truncated = {}
for m in [0.9, 1.0]:
    ref = closed_solution(m)
    diffs = []
    highs = []
    ok_pos_mass = True
    for L in [250, 500, 1000, 2000]:
        l = np.arange(1, L + 1, dtype=float)
        r0 = heavy_data(m, L)
        dt = min(1e-3, 0.5 / ((m + 0.5) * L + 2))
        sol = rk4(lambda t, y: rhs_trunc(y, L, l), r0, T_END, dt, T_CHECK)
        truncated[(m, L)] = sol
        d = max(np.abs(sol[t][:20] - ref[t][:20]).max() for t in T_CHECK)
        h = max(np.dot(l[L // 2:], sol[t][L // 2:]) for t in T_CHECK)
        mass_err = max(abs(np.dot(l, sol[t]) - m) for t in T_CHECK)
        ok_pos_mass = ok_pos_mass and min(sol[t].min() for t in T_CHECK) > -1e-14 and mass_err < 1e-10
        diffs.append(d)
        highs.append(h)
    report("D1 heavy tail m=%.1f: truncated systems stay >= 0 and conserve mass (L = 250..2000)" % m, ok_pos_mass)
    report("D1 heavy tail m=%.1f: sizes <= 20 converge to the closed (mass-conserving) system as L grows" % m,
           all(diffs[i + 1] < diffs[i] for i in range(3)) and diffs[-1] < diffs[0] / 4,
           "max diff, L=250,500,1000,2000: " + ", ".join("%.1e" % d for d in diffs))
    report("D1 heavy tail m=%.1f: mass in sizes > L/2 decreases with L (no accumulation at the cutoff)" % m,
           all(highs[i + 1] < highs[i] for i in range(3)),
           ", ".join("%.2e" % h for h in highs))


# ------------------------------------------------------------------ characteristics (Radau)
def tau(x):
    return 2.0 * np.tanh(x / 2.0)


def char_rhs(m):
    def f(t, y):
        X, P, U, J, Pi = y
        T = tau(X)
        k = 1.0 - T * T / 4.0
        kap = 1.0 / T
        kap1 = -k / T ** 2
        kap2 = k / (2 * T) + 2 * k * k / T ** 3
        return [P - m - 0.5,
                k * U / T ** 2 - P / T,
                0.5 * (P * P + m - m * m) - U / T,
                Pi,
                -J * (kap2 * U + 2 * kap1 * P) - kap * Pi]
    return f


def init_heavy(m, y):
    c = m / mp.zeta(1.5)
    z = mp.e ** (-mp.mpf(y))
    F0 = c * (mp.zeta(2.5) - mp.polylog(2.5, z))
    F1 = c * mp.polylog(1.5, z)
    F2 = -c * mp.polylog(0.5, z)
    return [float(y), float(F1), float(F0), 1.0, float(F2)]


def run_char(m, y, t, xstop, dense=False):
    """Integrate up to time t, or until X falls below xstop (X is decreasing, so then X(t) < xstop)."""
    def hit(s, Y):
        return Y[0] - xstop
    hit.terminal = True
    return solve_ivp(char_rhs(m), (0.0, t), init_heavy(m, y), method='Radau', rtol=1e-10, atol=1e-13,
                     events=hit, dense_output=dense)


def char_F(m, x, t):
    """Value F(x,t) of the characteristic solution, and invariants along that characteristic."""
    lo, hi = x + t / 2 - 1e-9, x + (m + 0.5) * t + 1e-9
    for _ in range(50):
        mid = 0.5 * (lo + hi)
        s = run_char(m, mid, t, 0.999 * x)
        if s.status == 1 or s.y[0, -1] < x:
            lo = mid
        else:
            hi = mid
    s = run_char(m, hi, t, 1e-14, dense=True)
    ts = np.linspace(0, s.t[-1], 2000)
    X, P, U, J, Pi = s.sol(ts)
    T = tau(X)
    inv = dict(P_lo=P.min(), P_hi=(P - m).max(), U=U.min(), W=((T * (m + P) - 2 * U) / T).min(),
               Gam=(T * Pi / J).min(), J=J.min())
    return s.y[2, -1], s.y[0, -1], inv


for m in [0.9, 1.0]:
    worst = dict(P_lo=1, P_hi=-1, U=1, W=1, Gam=1, J=1)
    comp_ok = True
    mono_ok = True
    details = []
    for t in T_CHECK:
        for x in [0.02, 0.1, 0.5, 2.0]:
            F, Xend, inv = char_F(m, x, t)
            for kk in worst:
                worst[kk] = max(worst[kk], inv[kk]) if kk == 'P_hi' else min(worst[kk], inv[kk])
            gaps = []
            for L in [250, 500, 1000, 2000]:
                lv = np.arange(1, L + 1, dtype=float)
                FL = np.dot(1 - np.exp(-lv * x), truncated[(m, L)][t])
                gaps.append(FL - F)
            comp_ok = comp_ok and min(gaps) > -1e-8
            mono_ok = mono_ok and all(gaps[i + 1] <= gaps[i] + 1e-10 for i in range(3))
            details.append("t=%g x=%g: %.1e..%.1e" % (t, x, gaps[0], gaps[-1]))
    report("D2 heavy tail m=%.1f: 0<=P<=m, U>0, W>=0, Gamma>-1, J>0 along the characteristics used" % m,
           worst['P_lo'] > -1e-9 and worst['P_hi'] < 1e-9 and worst['U'] > 0 and worst['W'] > -1e-9
           and worst['Gam'] > -1 and worst['J'] > 0,
           "min P=%.1e, max P-m=%.1e, min W/tau=%.1e, min Gamma=%.3f, min J=%.2e"
           % (worst['P_lo'], worst['P_hi'], worst['W'], worst['Gam'], worst['J']))
    report("D2 heavy tail m=%.1f: comparison F^L >= F, gap decreasing in L" % m, comp_ok and mono_ok,
           "; ".join(details[:3]) + "; ...")

# ------------------------------------------------------------------ D3 slope at the boundary, m = 1
m = 1.0
for t in T_CHECK:
    vals = []
    for x in [1e-1, 1e-2, 1e-3, 1e-4]:
        F, _, _ = char_F(m, x, t)
        vals.append(F / tau(x))
    report("D3 m=1, t=%g: F(x,t)/tau(x) increases towards 1 as x -> 0" % t,
           all(vals[i + 1] > vals[i] for i in range(3)) and vals[-1] < 1 + 1e-9 and vals[-1] > 0.9,
           "x=1e-1..1e-4: " + ", ".join("%.5f" % v for v in vals))

print("\nSUMMARY:", "ALL PASS" if FAILS == 0 else "%d FAIL" % FAILS)
sys.exit(1 if FAILS else 0)
