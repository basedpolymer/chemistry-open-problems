#!/usr/bin/env python3
"""v09_po_shooting.py -- verifier's direct computation of the small periodic orbits near the Hopf point
(no Lyapunov-coefficient formula involved): third, formula-free test of the direction of bifurcation.

For given M and amplitude a, solve by Newton (shooting, variational equations) for (x0, T, K):
    phi_T(x0; K)_l = x0_l   (l = 2..M; l = 1 follows from mass conservation),
    l . x0 = (1+K) M(M+1)/2          (the mass level of the constant equilibrium A = 1+K),
    (x0 - A 1) . s = 0               (phase: Poincare section, s = Im q / |Im q|),
    (x0 - A 1) . r = a               (amplitude along r = Re q / |Re q|),
q the Hopf eigenvector at K_c.  The sign of K(a) - K_c gives the side of the bifurcation:
    K(a) > K_c  for small a  <=> supercritical (orbits exist where the equilibrium is unstable),
    K(a) < K_c  <=> subcritical.  Also printed: (K(a) - K_c)/a^2 (should tend to a constant as a -> 0) and the
nontrivial Floquet multiplier closest to 1 (unstable orbit if > 1 within the mass hyperplane).
Inputs: kappa_c from the verifier's certified file v02_dense_j1.csv.
Usage: python v09_po_shooting.py M a1 a2 ...     Environment: WSL Fedora, numpy/scipy. Runtime: ~1 min per point.
"""
import sys, os, csv, math, time
import numpy as np
from scipy.integrate import solve_ivp
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from v08_simulate import field, jac

def flow_var(x0, T, K, M):
    """phi_T(x0), d phi_T / d x0 (monodromy), d phi_T / d K"""
    def rhs(t, y):
        x = y[:M]; Phi = y[M:M + M*M].reshape(M, M); dK = y[M + M*M:]
        Jx = jac(x, K)
        fK = np.zeros(M); fK[0] = M*x[M - 1]; fK[M - 1] = -x[M - 1]       # d f / d K
        return np.concatenate([field(x, K), (Jx@Phi).ravel(), Jx@dK + fK])
    y0 = np.concatenate([x0, np.eye(M).ravel(), np.zeros(M)])
    sol = solve_ivp(rhs, (0, T), y0, method="Radau", rtol=1e-12, atol=1e-14)
    y = sol.y[:, -1]
    return y[:M], y[M:M + M*M].reshape(M, M), y[M + M*M:]

def solve_po(M, Kc, om, q, a, guess=None):
    A = 1 + Kc
    r = np.real(q); r = r/np.linalg.norm(r)
    s = np.imag(q); s = s - (s@r)*r; s = s/np.linalg.norm(s)
    l = np.arange(1, M + 1, dtype=float)
    if guess is None:
        x0 = A*np.ones(M) + a*r; T = 2*math.pi/om; K = Kc
    else:
        x0, T, K = guess
    for it in range(30):
        xT, Phi, dK = flow_var(x0, T, K, M)
        Ak = 1 + K
        F = np.concatenate([(xT - x0)[1:], [l@x0 - Ak*M*(M + 1)/2, (x0 - Ak)@s, (x0 - Ak)@r - a]])
        J = np.zeros((M + 2, M + 2))
        J[:M - 1, :M] = (Phi - np.eye(M))[1:]
        J[:M - 1, M] = field(xT, K)[1:]
        J[:M - 1, M + 1] = dK[1:]
        J[M - 1, :M] = l; J[M - 1, M + 1] = -M*(M + 1)/2
        J[M, :M] = s; J[M, M + 1] = -s.sum()
        J[M + 1, :M] = r; J[M + 1, M + 1] = -r.sum()
        d = np.linalg.solve(J, -F)
        x0 = x0 + d[:M]; T = T + d[M]; K = K + d[M + 1]
        if np.abs(d).max() < 1e-13*max(1.0, np.abs(x0).max()):
            break
    xT, Phi, _ = flow_var(x0, T, K, M)
    res = np.abs(xT - x0).max()
    mult = np.linalg.eigvals(Phi)
    # multipliers: 1 (time shift), 1 (mass direction), the remaining one closest to 1 is the radial one
    mult = sorted(mult, key=lambda z: abs(z - 1))
    return x0, T, K, res, mult

def main():
    M = int(sys.argv[1]); amps = [float(x) for x in sys.argv[2:]]
    rows = {int(r_["M"]): r_ for r_ in csv.DictReader(open(os.path.join(HERE, "v02_dense_j1.csv")))}
    Kc = float(rows[M]["kappa"])/math.sqrt(M); om = float(rows[M]["omega"])
    Jc = jac((1 + Kc)*np.ones(M), Kc)
    w, V = np.linalg.eig(Jc); i = int(np.argmin(np.abs(w - 1j*om))); q = V[:, i]
    print("M=%d K_c=%.15g omega=%.8g (l1 sign from v02: %s)" % (M, Kc, om, rows[M]["l1_unit"]))
    guess = None
    out = []
    for a in amps:
        t0 = time.time()
        x0, T, K, res, mult = solve_po(M, Kc, om, q, a, guess)
        guess = (x0, T, K)
        nm = 0.5*(x0[M - 1] - (1 + K))
        print("  a=%.4f: K-K_c=%+.6e  (K-K_c)/a^2=%+.6e  T=%.10g  periodicity residual %.1e  "
              "3 multipliers closest to 1: %s  (%.0fs)" % (a, K - Kc, (K - Kc)/a**2, T, res,
              ", ".join("%.12f" % abs(z) for z in mult[:3]), time.time() - t0), flush=True)
        out.append((a, K - Kc))
    side = [np.sign(d) for _, d in out]
    print("  side of bifurcation: %s for all amplitudes -> %s"
          % ("K>K_c" if all(s_ > 0 for s_ in side) else ("K<K_c" if all(s_ < 0 for s_ in side) else "mixed"),
             "supercritical" if all(s_ > 0 for s_ in side) else ("subcritical" if all(s_ < 0 for s_ in side) else "?")))

if __name__ == "__main__":
    main()
