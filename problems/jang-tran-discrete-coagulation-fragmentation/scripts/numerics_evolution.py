#!/usr/bin/env python3
"""Floating-point sanity checks of the part-(ii) argument (numerical evidence, NOT part of the proof).

What is checked, for m in {0.5, 0.75, 0.9, 0.97} and two initial data
  (a) monomers rho^0 = m delta_1,
  (b) heavy tail rho^0_j = c j^{-5/2} (finite mass m, INFINITE second moment):
 N1  along characteristics of the Hamilton-Jacobi equation: 0 <= P <= m, U >= 0, W = tau(m+P) - 2U >= 0,
     Gamma = tau F_xx > -1, J = dX/dx0 > 0, exit through x = 0 in finite time (paper: invariant region,
     curvature barrier);
 N2  eta = m - U/tau obeys the quantitative bound of the slope lemma (b) near the exit;
 N3  truncated system (N sizes): mass conserved, rho >= 0; its transform F^N dominates the
     characteristic solution F (comparison lemma); the gap F^N - F shrinks as N grows;
 N4  the closed lower-triangular system agrees with the truncated system N=1600 for small sizes, and the
     truncated system has negligible mass near the cutoff (no numerical gelation);
 N5  long time: the closed system converges to the recursion (1.9) values; the recursion values are
     >= 0 with sum l rho~_l = m and sum rho~_l = m(1-m) (l <= 3000, m = 0.6, 0.75); truncated = closed
     system up to T = 30 for m = 0.6;
 N6  m = 1: invariant region, barrier and m - U/tau -> 0 at the exit (slope lemma (a)).
Small truncations (N of a few hundred) are NOT accurate for m close to 1 at long times: the
equilibrium then has a heavy n^{-5/2} tail up to very large sizes. This is a truncation error.
Run from the problem folder: python scripts/numerics_evolution.py

Environment: Python 3, numpy, scipy, mpmath. Runtime: about 4-6 minutes, one thread (verifier's run: 4 min 5 s).
Deterministic (fixed grids, no randomness). Prints PASS/FAIL per check; exit status 1 on failure.
"""
import sys
import numpy as np
import mpmath as mp
from scipy.integrate import solve_ivp
from scipy.signal import fftconvolve

FAILS = 0


def report(name, ok, info=""):
    global FAILS
    print(("PASS " if ok else "FAIL ") + name + ("  [" + info + "]" if info else ""))
    if not ok:
        FAILS += 1


def tau(x):
    return 2.0 * np.tanh(x / 2.0)


def taup(x):
    return 1.0 / np.cosh(x / 2.0) ** 2


# ------------------------------------------------------------------ initial data
def data_monomer(m):
    F0 = lambda x: m * (1 - np.exp(-x))
    F1 = lambda x: m * np.exp(-x)
    F2 = lambda x: -m * np.exp(-x)
    return F0, F1, F2


def data_heavy(m):
    z32 = float(mp.zeta(1.5))
    c = m / z32
    z52 = float(mp.zeta(2.5))
    F0 = lambda x: c * (z52 - float(mp.polylog(2.5, mp.e ** (-x))))
    F1 = lambda x: c * float(mp.polylog(1.5, mp.e ** (-x)))
    F2 = lambda x: -c * float(mp.polylog(0.5, mp.e ** (-x)))
    return F0, F1, F2, c


# ------------------------------------------------------------------ characteristics
def char_rhs(m):
    def f(t, y):
        X, P, U, J, Pi = y
        T = tau(X)
        Tp = taup(X)
        kap = 1.0 / T
        kp = -Tp / T ** 2
        # kappa'' = -tau''/tau^2 + 2 tau'^2/tau^3, tau'' = -tau tau'/2
        kpp = (T * Tp / 2) / T ** 2 + 2 * Tp ** 2 / T ** 3
        Xd = P - m - 0.5
        Pd = -kp * U - kap * P
        Ud = 0.5 * (P * P + m - m * m) - kap * U
        Jd = Pi
        Pid = -J * (kpp * U + 2 * kp * P) - kap * Pi
        return [Xd, Pd, Ud, Jd, Pid]
    return f


def exit_event(t, y):
    return y[0] - 1e-7


exit_event.terminal = True
exit_event.direction = -1


def integrate_char(m, x0, F0, F1, F2, tmax, dense=False, rtol=1e-11, atol=1e-13):
    y0 = [x0, F1(x0), F0(x0), 1.0, F2(x0)]
    sol = solve_ivp(char_rhs(m), (0, tmax), y0, method='DOP853', rtol=rtol, atol=atol,
                    events=exit_event, dense_output=dense, max_step=0.05)
    return sol


def check_char_invariants(m, F0, F1, F2, label):
    """N1/N2: invariants along a family of characteristics run until exit."""
    x0s = np.concatenate([np.geomspace(1e-3, 0.5, 25), np.linspace(0.6, 8.0, 25)])
    worst = dict(P_lo=np.inf, P_hi=-np.inf, U=np.inf, W=np.inf, Gam=np.inf, J=np.inf)
    exited = True
    eta_ok = True
    alpha = (1 - m) / (2 * m + 1)
    for x0 in x0s:
        sol = integrate_char(m, x0, F0, F1, F2, tmax=3 * x0 + 1, dense=True)
        if sol.status != 1:
            exited = False
        ts = np.linspace(0, sol.t[-1], 4000)
        Y = sol.sol(ts)
        X, P, U, J, Pi = Y
        T = tau(X)
        W = T * (m + P) - 2 * U
        Gam = T * Pi / J
        worst['P_lo'] = min(worst['P_lo'], P.min())
        worst['P_hi'] = max(worst['P_hi'], P.max())
        worst['U'] = min(worst['U'], U.min())
        worst['W'] = min(worst['W'], (W / np.maximum(T, 1e-300)).min())
        worst['Gam'] = min(worst['Gam'], Gam.min())
        worst['J'] = min(worst['J'], J.min())
        # tau-variable variant of the bound of the slope lemma (b) (paper, Lemma 20(b)), from an earlier
        # draft: once theta = tau(X) <= theta1 := min(tau(x0), 1),
        # eta <= theta^alpha (m theta1^-alpha + C theta1^(2-alpha)), C = b/(c1 (2-alpha)),
        # b = m(m+1/2)/4, c1 = (1 - theta1^2/4)/2.
        th1 = min(tau(x0), 1.0)
        b = m * (m + 0.5) / 4
        c1 = (1 - th1 ** 2 / 4) / 2
        C = b / (c1 * (2 - alpha))
        eta = m - U / T
        mask = T <= th1
        bound = T[mask] ** alpha * (m * th1 ** (-alpha) + C * th1 ** (2 - alpha))
        if mask.any() and np.any(eta[mask] > bound + 1e-8):
            eta_ok = False
    tol = 1e-8
    report("N1 %s m=%.2f: 0<=P<=m along characteristics" % (label, m),
           worst['P_lo'] > -tol and worst['P_hi'] < m + tol,
           "min P=%.3e, max P-m=%.3e" % (worst['P_lo'], worst['P_hi'] - m))
    report("N1 %s m=%.2f: U>=0 and W>=0 along characteristics" % (label, m),
           worst['U'] > -tol and worst['W'] > -tol,
           "min U=%.3e, min W/tau=%.3e" % (worst['U'], worst['W']))
    report("N1 %s m=%.2f: barrier Gamma = tau F_xx > -1 and J > 0" % (label, m),
           worst['Gam'] > -1 and worst['J'] > 0,
           "min Gamma=%.4f, min J=%.3e" % (worst['Gam'], worst['J']))
    report("N1 %s m=%.2f: every characteristic exits through x=0 in finite time" % (label, m), exited)
    report("N2 %s m=%.2f: eta = m - U/tau obeys the bound of the slope lemma (b) (tau-variable form) near the exit" % (label, m), eta_ok)


def char_solution_at(m, F0, F1, F2, x, t, lo, hi):
    """F(x,t) of the characteristic solution: bisection on x0 so that X(t; x0) = x."""
    def Xt(x0):
        sol = integrate_char(m, x0, F0, F1, F2, tmax=t, rtol=1e-12, atol=1e-14)
        if sol.status == 1:      # exited before t
            return -1.0, None
        return sol.y[0, -1], sol.y[2, -1]
    a, b = lo, hi
    for _ in range(60):
        c = 0.5 * (a + b)
        Xc, Uc = Xt(c)
        if Xc < x:
            a = c
        else:
            b = c
    Xc, Uc = Xt(0.5 * (a + b))
    return Uc


# ------------------------------------------------------------------ truncated system
def trunc_rhs(N):
    l = np.arange(1, N + 1, dtype=float)

    def f(t, rho):
        a = l * rho
        conv = fftconvolve(a, a)[:N]
        gain = np.zeros(N)
        gain[1:] = 0.5 * conv[:N - 1]
        A = np.cumsum(a)                       # A[n-1] = sum_{k<=n} k rho_k
        # loss_l = l rho_l * sum_{k <= N-l} k rho_k ; for l = N the sum is empty
        S = np.zeros(N)
        S[:N - 1] = A[N - 2::-1]
        loss = a * S
        tail = np.concatenate([np.cumsum(rho[::-1])[::-1][1:], [0.0]])   # sum_{j>l} rho_j
        frag = -0.5 * (l - 1) * rho + tail
        return gain - loss + frag
    return f


def closed_rhs(m, m00, Lc):
    l = np.arange(1, Lc + 1, dtype=float)
    c = m * l + (l + 1) / 2

    def f(t, rho):
        a = l * rho
        conv = np.convolve(a, a)[:Lc]
        gain = np.zeros(Lc)
        gain[1:] = 0.5 * conv[:Lc - 1]
        m0 = m - m * m + np.exp(-t / 2) * (m00 - m + m * m)
        partial = np.concatenate([[0.0], np.cumsum(rho)[:-1]])          # sum_{i<l} rho_i
        return gain - c * rho + m0 - partial
    return f


def recursion(m, L):
    r = np.zeros(L + 1)
    S = 0.0
    for l in range(1, L + 1):
        i = np.arange(1, l)
        conv = np.sum(i * (l - i) * r[1:l] * r[l - 1:0:-1]) if l > 1 else 0.0
        r[l] = (2 * m * (1 - m) + conv - 2 * S) / ((2 * m + 1) * l + 1)
        S += r[l]
    return r[1:]


def check_truncation(m, label, F0, F1, F2, rho0_full_fn, N, times, xs):
    rho0 = rho0_full_fn(N)
    l = np.arange(1, N + 1)
    mass0 = np.dot(l, rho0)
    sol = solve_ivp(trunc_rhs(N), (0, max(times)), rho0, method='LSODA', rtol=1e-11, atol=1e-15,
                    t_eval=sorted(times))
    ok_mass = np.all(np.abs(sol.y.T @ l - mass0) < 1e-9)
    ok_pos = sol.y.min() > -1e-12
    report("N3 %s m=%.2f N=%d: truncated system keeps rho>=0 and mass" % (label, m, N), ok_mass and ok_pos,
           "mass drift %.1e, min rho %.1e" % (np.abs(sol.y.T @ l - mass0).max(), sol.y.min()))
    worst_neg = 0.0
    worst_abs = 0.0
    for k, t in enumerate(sol.t):
        rho = sol.y[:, k]
        for x in xs:
            FN = np.dot(1 - np.exp(-l * x), rho)
            F = char_solution_at(m, F0, F1, F2, x, t, lo=x + t / 2 - 1e-9, hi=x + (m + 0.5) * t + 1e-9)
            worst_neg = min(worst_neg, FN - F)
            worst_abs = max(worst_abs, abs(FN - F))
    report("N3 %s m=%.2f N=%d: comparison F^N >= F (characteristic solution)" % (label, m, N),
           worst_neg > -1e-8, "min(F^N - F)=%.2e, max|F^N - F|=%.2e" % (worst_neg, worst_abs))
    return sol


def check_m_equals_one():
    """m = 1, monomer and heavy-tail data: invariant region, barrier, and m - U/tau -> 0 at the exit."""
    m = 1.0
    for label, (F0, F1, F2) in [("monomer", data_monomer(m)), ("heavy-tail", data_heavy(m)[:3])]:
        x0s = np.concatenate([np.geomspace(0.05, 0.5, 10), np.linspace(0.75, 6.0, 15)])
        ok_inv, ok_exit = True, True
        worst_eta = 0.0
        for x0 in x0s:
            sol = integrate_char(m, x0, F0, F1, F2, tmax=3 * x0 + 1, dense=True)
            ts = np.linspace(0, sol.t[-1], 4000)
            X, P, U, J, Pi = sol.sol(ts)
            T = tau(X)
            W = T * (m + P) - 2 * U
            Gam = T * Pi / J
            if (P.min() < -1e-8 or P.max() > m + 1e-8 or U.min() < -1e-8 or (W / T).min() < -1e-8
                    or Gam.min() <= -1 or J.min() <= 0):
                ok_inv = False
            Xe, Ue = sol.y[0, -1], sol.y[2, -1]
            eta_exit = m - Ue / tau(Xe)
            worst_eta = max(worst_eta, eta_exit)
            if sol.status != 1:
                ok_exit = False
        report("N6 %s m=1: invariant region, barrier Gamma > -1, J > 0 along characteristics" % label, ok_inv)
        report("N6 %s m=1: m - U/tau is small when characteristics reach x = 1e-7" % label,
               ok_exit and worst_eta < 0.05, "max eta at exit %.3e" % worst_eta)


def rho_mono_fn(m):
    def f(N):
        r = np.zeros(N)
        r[0] = m
        return r
    return f


def main():
    ms = [0.5, 0.75, 0.9, 0.97]
    for m in ms:
        F0, F1, F2 = data_monomer(m)
        check_char_invariants(m, F0, F1, F2, "monomer")
        H0, H1, H2, c = data_heavy(m)
        check_char_invariants(m, H0, H1, H2, "heavy-tail")

    # N3: comparison F^N >= F, monomer data (the gap is the truncation error of F^N)
    for m in [0.75, 0.97]:
        F0, F1, F2 = data_monomer(m)
        check_truncation(m, "monomer", F0, F1, F2, rho_mono_fn(m), N=300, times=[0.5, 2.0, 6.0],
                         xs=[0.02, 0.1, 0.5, 2.0])
    # N3b: the gap F^N - F shrinks as N grows (evidence that F is the transform of the solution)
    m = 0.75
    F0, F1, F2 = data_monomer(m)
    t = 2.0
    xs = [0.05, 0.5]
    Fc = [char_solution_at(m, F0, F1, F2, x, t, lo=x + t / 2 - 1e-9, hi=x + (m + 0.5) * t + 1e-9)
          for x in xs]
    gaps = []
    for N in [100, 200, 400]:
        sol = solve_ivp(trunc_rhs(N), (0, t), rho_mono_fn(m)(N), method='LSODA', rtol=1e-11,
                        atol=1e-15, t_eval=[t])
        lN = np.arange(1, N + 1)
        gaps.append(max(np.dot(1 - np.exp(-lN * x), sol.y[:, -1]) - F for x, F in zip(xs, Fc)))
    report("N3b monomer m=0.75, t=2: gap F^N - F >= 0 and decreasing in N (N=100,200,400)",
           min(gaps) > -1e-9 and gaps[0] > gaps[1] > gaps[2] - 1e-10 and gaps[2] < 1e-7,
           "gaps " + ", ".join("%.2e" % g for g in gaps))

    # N4: closed lower-triangular system (paper, Proposition 7) versus truncated system, sizes <= 10
    for m in [0.6, 0.75, 0.9]:
        Lc = 40
        solc = solve_ivp(closed_rhs(m, m, Lc), (0, 2.0), rho_mono_fn(m)(Lc), method='LSODA',
                         rtol=1e-12, atol=1e-16, t_eval=[2.0])
        N = 1600
        sol = solve_ivp(trunc_rhs(N), (0, 2.0), rho_mono_fn(m)(N), method='LSODA', rtol=1e-11,
                        atol=1e-16, t_eval=[2.0])
        lN = np.arange(1, N + 1)
        diff = np.abs(solc.y[:10, -1] - sol.y[:10, -1]).max()
        high = np.dot(lN[N // 2:], sol.y[N // 2:, -1])
        report("N4 monomer m=%.2f, t=2: closed system = truncated system N=1600 (sizes<=10)" % m,
               diff < 1e-9 and high < 1e-6, "max diff %.1e, mass in sizes>800: %.1e" % (diff, high))

    # N5: long time. Closed system -> recursion (1.9); recursion mass ~ m (evidence for part (i));
    # truncated system N=1600 keeps its mass away from the cutoff up to T=30.
    for m in [0.6, 0.75]:
        Lc = 40
        T = 60.0
        solc = solve_ivp(closed_rhs(m, m, Lc), (0, T), rho_mono_fn(m)(Lc), method='LSODA',
                         rtol=1e-12, atol=1e-16, t_eval=[T])
        rec = recursion(m, 3000)
        d5 = np.abs(solc.y[:, -1] - rec[:Lc]).max()
        report("N5 m=%.2f: closed system at T=60 matches recursion (1.9), l<=40" % m, d5 < 1e-10,
               "max diff %.1e" % d5)
        lr = np.arange(1, 3001)
        mass_rec = np.dot(lr, rec)
        m0_rec = rec.sum()
        report("N5 m=%.2f: recursion (1.9): sum_{l<=3000} l rho~_l = m, sum rho~_l = m(1-m), rho~ >= 0" % m,
               abs(mass_rec - m) < 1e-8 and abs(m0_rec - m * (1 - m)) < 1e-10 and rec.min() > -1e-15,
               "mass %.12f, m0 %.12f, min %.1e" % (mass_rec, m0_rec, rec.min()))
    # truncated system at long times (m = 0.6: the stationary tail is light, so L = 1600 is accurate)
    m = 0.6
    N = 1600
    # explicit DOP853 (verifier's change): LSODA switches to its stiff mode here and spends > 1 h on dense
    # 1600 x 1600 Jacobians; DOP853 needs about 30 s and gives the same check (difference ~5e-13)
    sol = solve_ivp(trunc_rhs(N), (0, 30.0), rho_mono_fn(m)(N), method='DOP853', rtol=1e-11,
                    atol=1e-16, t_eval=[5.0, 15.0, 30.0])
    solc = solve_ivp(closed_rhs(m, m, 40), (0, 30.0), rho_mono_fn(m)(40), method='LSODA', rtol=1e-12,
                     atol=1e-16, t_eval=[5.0, 15.0, 30.0])
    lN = np.arange(1, N + 1)
    high = (lN[N // 2:] @ sol.y[N // 2:, :]).max()
    diff = np.abs(sol.y[:30, :] - solc.y[:30, :]).max()
    report("N5 monomer m=0.60, N=1600: truncated = closed system up to T=30 (sizes<=30), no mass near cutoff",
           diff < 1e-9 and high < 1e-8, "max diff %.1e, max mass in sizes>800 %.1e" % (diff, high))

    # m = 1: invariants, barrier, and lambda = U/tau -> 1 at the exit (paper, Lemma 20(a))
    check_m_equals_one()

    # N3 with a heavy-tail initial datum (infinite second moment)
    m = 0.9
    H0, H1, H2, c = data_heavy(m)

    def rho_heavy(N, c=c, m=m):
        j = np.arange(1, N + 1, dtype=float)
        r = c * j ** -2.5
        tail_mass = m - np.dot(j, r)
        r[-1] += tail_mass / N
        return r
    check_truncation(m, "heavy-tail", H0, H1, H2, rho_heavy, N=800, times=[0.5, 2.0],
                     xs=[0.05, 0.3, 1.5])

    print()
    print("SUMMARY:", "ALL PASS" if FAILS == 0 else "%d FAIL" % FAILS)
    sys.exit(1 if FAILS else 0)


if __name__ == "__main__":
    main()
