#!/usr/bin/env python3
"""v08_simulate.py -- verifier's direct simulation of the full Becker-Doring system near the Hopf point.

End-to-end test of the SIGN CONVENTION of the Lyapunov coefficient (common to the solver's closed form and to the
verifier's dense formula): if l1 < 0 (supercritical), then for K = K_c (1 + delta) a small stable cycle must
exist, attracting perturbations from inside and from outside, with n_1 half-amplitude
      a_pred = 2 |q_1| sqrt(-mu / (omega l1)),     mu = Re lambda(K),  q normalised by |q| = 1 (full coordinates),
and for K = K_c (1 - delta) a perturbation of size ~ a_pred(+delta) must decay to the equilibrium.
Inputs: kappa_c and l1 (unit) from the verifier's own certified file v02_dense_j1.csv; everything else (field,
Jacobian, eigenvectors, integration) is recomputed here with numpy/scipy (LSODA, analytic Jacobian).
Usage: python v08_simulate.py M [M ...]
Environment: WSL Fedora, numpy/scipy. Runtime: a few minutes per M.
"""
import sys, os, csv, math, time
import numpy as np
from scipy.integrate import solve_ivp
HERE = os.path.dirname(os.path.abspath(__file__))

def field(n, K):
    M = len(n); N = M - 1
    J = n[:N]*n[0] - n[1:]
    f = np.empty(M)
    f[0] = -J[0] - J.sum() + M*K*n[M - 1]
    f[1:N] = J[:N - 1] - J[1:]
    f[M - 1] = J[N - 1] - K*n[M - 1]
    return f

def jac(n, K):
    M = len(n); N = M - 1
    # dJ_l/dn: J_l = n_l n_1 - n_{l+1}
    dJ = np.zeros((N, M))
    for l in range(N):
        dJ[l, l] += n[0]; dJ[l, 0] += n[l]; dJ[l, l + 1] -= 1
    Jm = np.zeros((M, M))
    Jm[0] = -dJ[0] - dJ.sum(axis=0); Jm[0, M - 1] += M*K
    Jm[1:N] = dJ[:N - 1] - dJ[1:]
    Jm[M - 1] = dJ[N - 1]; Jm[M - 1, M - 1] -= K
    return Jm

def eig_near(K, M, om):
    A = 1 + K
    w, V = np.linalg.eig(jac(A*np.ones(M), K))
    i = int(np.argmin(np.abs(w - 1j*om)))
    q = V[:, i]/np.linalg.norm(V[:, i])
    return w[i], q

def main():
    rows = {int(r["M"]): r for r in csv.DictReader(open(os.path.join(HERE, "v02_dense_j1.csv")))}
    bad = 0
    for M in [int(a) for a in sys.argv[1:]]:
        r = rows[M]
        Kc = float(r["kappa"])/math.sqrt(M); om = float(r["omega"]); l1 = float(r["l1_unit"])
        lam_c, qc = eig_near(Kc, M, om)
        print("M=%d K_c=%.12g omega=%.6g l1(unit)=%.6e  |Re lam(K_c)|=%.1e" % (M, Kc, om, l1, abs(lam_c.real)))
        res = []
        for delta in (0.04, 0.02, 0.01):
            K = Kc*(1 + delta)
            lam, q = eig_near(K, M, om)
            mu = lam.real
            pred = 2*abs(qc[0])*math.sqrt(-mu/(om*l1))
            nbar = (1 + K)*np.ones(M)
            T = 25.0/mu
            amps = []
            for fac in (0.3, 3.0):              # start inside and outside the predicted cycle
                x0 = nbar + fac*pred/(2*abs(q[0]))*2*np.real(q*np.conj(q[0])/abs(q[0]))
                t0 = time.time()
                sol = solve_ivp(lambda t, n: field(n, K), (0, T), x0, method="LSODA", rtol=1e-10, atol=1e-12,
                                jac=lambda t, n: jac(n, K), dense_output=True)
                tt = np.linspace(T - 6*2*math.pi/lam.imag, T, 6000)
                Y = sol.sol(tt)
                amps.append(0.5*(Y[0].max() - Y[0].min()))
            ratio = [a/pred for a in amps]
            res.append((delta, ratio))
            print("   K=K_c(1+%.2f): mu=%.3e  a_pred=%.4e  a_sim(start 0.3x)=%.4e  a_sim(start 3x)=%.4e  ratios %.4f %.4f (%.0fs)"
                  % (delta, mu, pred, amps[0], amps[1], ratio[0], ratio[1], time.time() - t0), flush=True)
        # below threshold: decay
        delta = 0.02
        K = Kc*(1 - delta)
        lam, q = eig_near(K, M, om)
        mu = lam.real
        a0 = 2*abs(qc[0])*math.sqrt(abs(mu)/(om*abs(l1)))
        x0 = (1 + K)*np.ones(M) + a0/(2*abs(q[0]))*2*np.real(q*np.conj(q[0])/abs(q[0]))
        T = 25.0/abs(mu)
        sol = solve_ivp(lambda t, n: field(n, K), (0, T), x0, method="LSODA", rtol=1e-10, atol=1e-12,
                        jac=lambda t, n: jac(n, K), dense_output=True)
        tt = np.linspace(T - 6*2*math.pi/lam.imag, T, 3000)
        Y = sol.sol(tt)
        afin = 0.5*(Y[0].max() - Y[0].min())
        print("   K=K_c(1-%.2f): mu=%.3e  initial n_1 amplitude %.3e -> final %.3e" % (delta, mu, a0, afin))
        # criteria: both starts converge to the same amplitude, ratio -> 1 as delta -> 0, decay below threshold
        conv = all(abs(r_[0] - r_[1]) < 1e-3 for _, r_ in res)
        # amplitude ratio -> 1: deviation decreasing with delta, or already at the measurement level (< 1e-3)
        dev = [abs(r_[0] - 1) for _, r_ in res]
        trend = dev[-1] < 0.05 and (dev[-1] < dev[0] or max(dev) < 1e-3)
        decay = afin < 1e-3*a0
        ok = l1 < 0 and conv and trend and decay
        bad += not ok
        print("   M=%d: same limit from inside/outside %s, ratio -> 1 %s, decay below K_c %s: %s"
              % (M, conv, trend, decay, "PASS" if ok else "FAIL"), flush=True)
    print("SUMMARY:", "PASS" if bad == 0 else "FAIL")
    sys.exit(1 if bad else 0)

if __name__ == "__main__":
    main()
