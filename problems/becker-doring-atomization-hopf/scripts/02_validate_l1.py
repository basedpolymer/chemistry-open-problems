#!/usr/bin/env python3
"""02_validate_l1.py -- validation of the Lyapunov-coefficient machinery (evidence, not part of the proof).

Environment: WSL Fedora, Python 3.12, numpy/scipy; optional: Julia 1.13 + BifurcationKit 0.8 for V3.
Runtime: about 10 minutes (V2 integrates a stiff ODE over ~40 growth times).
Claims (PASS/FAIL):
  V1  planar test: the generic Kuznetsov routine reproduces the Guckenheimer-Holmes coefficient (l1 = 2a/omega
      for |q| = 1) on a planar system with quadratic and cubic terms;
  V2  direct simulation of the full system (2)-(4), M = 30, K = K_c (1 + delta): a small perturbation converges to a
      small stable cycle whose n_1 half-amplitude equals the normal-form prediction 2|q_1| sqrt(-mu/(omega l1)),
      the ratio tending to 1 linearly in delta (supercritical, and l1 quantitatively right);
  V3  BifurcationKit (independent code, automatic differentiation, reduced coordinates) classifies the Hopf points
      M = 15 as SubCritical and M = 16, 17, 30 as SuperCritical (skipped with a warning if julia is missing).
"""
import sys, math, subprocess, shutil, os
import numpy as np, scipy.linalg as sla
from scipy.integrate import solve_ivp
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import bdlib as bd

fails = 0
def report(name, ok, extra=""):
    global fails
    fails += (not ok)
    print("%-3s %s %s" % (name, "PASS" if ok else "FAIL", extra), flush=True)

# V1 -------------------------------------------------------------------------------------------
om = 1.3
a1, a2, a3, b1, b2, b3, c1, c2 = 0.7, -0.4, 0.9, 0.3, 1.1, -0.6, 0.25, -0.15
Aj = np.array([[0, -om], [om, 0.]])
Bf = lambda x, y: np.array([2*a1*x[0]*y[0] + a2*(x[0]*y[1] + x[1]*y[0]) + 2*a3*x[1]*y[1],
                            2*b1*x[0]*y[0] + b2*(x[0]*y[1] + x[1]*y[0]) + 2*b3*x[1]*y[1]])
Cf = lambda x, y, z: np.array([6*c1*x[0]*y[0]*z[0], 6*c2*x[1]*y[1]*z[1]])
w, VL, VR = sla.eig(Aj, left=True, right=True); i = np.argmax(w.imag)
q = VR[:, i]/np.sqrt(np.vdot(VR[:, i], VR[:, i]).real); p = VL[:, i]/np.conj(np.vdot(VL[:, i], q))
w11 = np.linalg.solve(Aj, Bf(q, np.conj(q))); w20 = np.linalg.solve(2j*om*np.eye(2) - Aj, Bf(q, q))
l1 = (np.vdot(p, Cf(q, q, np.conj(q))) - 2*np.vdot(p, Bf(q, w11)) + np.vdot(p, Bf(np.conj(q), w20))).real/(2*om)
fxx, fxy, fyy, fxxx = 2*a1, a2, 2*a3, 6*c1
gxx, gxy, gyy, gyyy = 2*b1, b2, 2*b3, 6*c2
aGH = (fxxx + gyyy)/16 + (fxy*(fxx + fyy) - gxy*(gxx + gyy) - fxx*gxx + fyy*gyy)/(16*om)
report("V1", abs(l1 - 2*aGH/om) < 1e-12, "l1 = %.12f, 2 a_GH/omega = %.12f" % (l1, 2*aGH/om))

# V2 -------------------------------------------------------------------------------------------
M = 30
cr = [c for c in bd.dense_crossings(M) if c[3] == 1][0]
Kc, lamc, _ = bd.dense_hopf(M, cr[0], cr[1], 0)
_, l1c, qc = bd.lyap_dense(Kc, M, lamc)
omc = lamc.imag
ratios = []
for delta in [0.02, 0.01, 0.005]:
    K = Kc*(1 + delta)
    w, VR = sla.eig(bd.Bmat(K, M)); i = np.argmin(np.abs(w - lamc)); mu = w[i].real
    pred = 2*abs(qc[0])*math.sqrt(-mu/(omc*l1c))
    nbar = (1 + K)*np.ones(M)
    T = 40.0/mu
    sol = solve_ivp(lambda t, n: bd.field(n, K), (0, T), nbar + 1e-3*np.real(VR[:, i]), method="LSODA",
                    rtol=1e-10, atol=1e-13, jac=lambda t, n: bd.jac_full(n, K), dense_output=True)
    Y = sol.sol(np.linspace(T - 5*2*math.pi/w[i].imag, T, 4000))
    amp = 0.5*(Y[0].max() - Y[0].min())
    ratios.append(amp/pred)
    print("     delta=%.3f  n_1 half-amplitude: simulation %.6e, normal form %.6e, ratio %.5f" % (delta, amp, pred, amp/pred))
lin = abs((ratios[0] - 1)/(ratios[1] - 1) - 2) < 0.2 and abs((ratios[1] - 1)/(ratios[2] - 1) - 2) < 0.2
report("V2", l1c < 0 and abs(ratios[-1] - 1) < 0.005 and lin, "ratios %s (-> 1, error ~ delta)" % ["%.5f" % r for r in ratios])

# V3 -------------------------------------------------------------------------------------------
if shutil.which("julia") is None:
    print("V3  SKIPPED (julia not found)")
else:
    args = []
    for Mq in [15, 16, 17, 30]:
        c = [c for c in bd.dense_crossings(Mq) if c[3] == 1][0]
        args += [str(Mq), "%.12g" % bd.dense_hopf(Mq, c[0], c[1], 0)[0]]
    out = subprocess.run(["julia", os.path.join(os.path.dirname(os.path.abspath(__file__)), "bk_check.jl")] + args,
                         capture_output=True, text=True, timeout=3000).stdout
    res = [l for l in out.splitlines() if l.startswith("RESULT")]
    for l in res:
        print("     " + l)
    want = {15: "SubCritical", 16: "SuperCritical", 17: "SuperCritical", 30: "SuperCritical"}
    ok3 = len(res) == 4 and all(("M=%d " % m) in l and ("type=" + t) in l for (m, t), l in zip(want.items(), res))
    report("V3", ok3, "BifurcationKit: M=15 SubCritical; M=16,17,30 SuperCritical")
print("SUMMARY:", "PASS" if fails == 0 else "FAIL")
sys.exit(0 if fails == 0 else 1)
