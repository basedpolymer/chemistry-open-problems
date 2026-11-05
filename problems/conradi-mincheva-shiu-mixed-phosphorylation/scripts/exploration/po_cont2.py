# Exploration script (floating point; not a certificate). Environment: WSL Fedora, Python 3.12, numpy, scipy, numba.
# Run from scripts/exploration/ (it imports the sibling modules). Typical runtime: seconds to minutes.
"""Faster periodic-orbit continuation (single shooting, pseudo-arclength) in the reduced coordinates.
x0 = x*(mu) + E z0, unknowns u = (z0 (6), T, mu). Floating point exploration."""
import sys
import numpy as np
from fractions import Fraction
from scipy.integrate import solve_ivp
from fastmodel import g9, J9, var_rhs
from model import steady_state, W

ZIDX = [2, 3, 4, 6, 7, 8]
E = np.zeros((9, 6))
for j, i in enumerate(ZIDX):
    E[i, j] = 1
E[0, :] = -1; E[1, :] = [-1, -1, 0, 0, 0, 0]; E[5, :] = [0, 0, 0, -1, 0, -1]
P = np.zeros((6, 9))
for j, i in enumerate(ZIDX):
    P[j, i] = 1


class PO:
    def __init__(self, k, x1, x6, rtol=1e-12, atol=1e-12):
        self.k = np.array(k, float); self.x1 = x1; self.x6 = x6
        self.xs1 = steady_state(self.k, x1, 1.0, x6) - steady_state(self.k, x1, 0.0, x6)
        self.xs0 = steady_state(self.k, x1, 0.0, x6)
        self.rtol = rtol; self.atol = atol

    def xs(self, mu):
        return self.xs0 + mu * self.xs1

    def flow(self, x0, T):
        y0 = np.concatenate([x0, np.eye(9).ravel()])
        s = solve_ivp(var_rhs, [0, T], y0, method='DOP853', rtol=self.rtol, atol=self.atol, args=(self.k,))
        y = s.y[:, -1]
        return y[:9], y[9:].reshape(9, 9)

    def residual(self, u, zref, nref):
        z0, T, mu = u[:6], u[6], u[7]
        x0 = self.xs(mu) + E @ z0
        xT, Phi = self.flow(x0, T)
        F = np.concatenate([P @ (xT - x0), [nref @ (z0 - zref)]])
        J = np.zeros((7, 8))
        J[:6, :6] = P @ Phi @ E - np.eye(6)
        J[:6, 6] = P @ g9(xT, self.k)
        J[:6, 7] = P @ ((Phi - np.eye(9)) @ self.xs1)
        J[6, :6] = nref
        Mz = P @ Phi @ E      # monodromy in z coordinates
        return F, J, Mz

    def newton(self, u, zref, nref, tang=None, u_pred=None, tol=1e-11, maxit=15, fix_mu=False):
        for it in range(maxit):
            F, J, Mz = self.residual(u, zref, nref)
            if tang is not None:
                F = np.concatenate([F, [tang @ (u - u_pred)]]); Jf = np.vstack([J, tang])
                du = np.linalg.solve(Jf, -F)
            else:
                du = np.concatenate([np.linalg.solve(J[:, :7], -F), [0.0]])
            u = u + du
            if np.linalg.norm(du[:7]) < tol * (1 + np.linalg.norm(u[:7])):
                F, J, Mz = self.residual(u, zref, nref)
                return u, J, Mz, True
        return u, J, Mz, False


def start_from_hopf(po, mu, amp_rel):
    xs = po.xs(mu)
    A = P @ J9(xs, po.k) @ E
    ev, R = np.linalg.eig(A); i = np.argmax(ev.imag); q = R[:, i]
    q = q / np.linalg.norm(q)
    z0 = amp_rel * np.linalg.norm(xs) * q.real / np.linalg.norm(q.real)
    return z0, 2 * np.pi / ev[i].imag, q


if __name__ == '__main__':
    k = [float(Fraction(s)) for s in sys.argv[1].split(',')]
    x1 = float(Fraction(sys.argv[2])); x6 = float(Fraction(sys.argv[3])); mu0 = float(sys.argv[4])
    rel = float(sys.argv[5]); amp = float(sys.argv[6]); nsteps = int(sys.argv[7]); ds = float(sys.argv[8])
    tag = sys.argv[9] if len(sys.argv) > 9 else 'branch'
    po = PO(k, x1, x6)
    mu = mu0 * (1 + rel)
    z0, T0, q = start_from_hopf(po, mu, amp)
    zref = np.zeros(6); nref = q.imag / np.linalg.norm(q.imag)
    zref = z0.copy(); nref = (P @ g9(po.xs(mu) + E @ z0, po.k)); nref /= np.linalg.norm(nref)
    u = np.concatenate([z0, [T0, mu]])
    # first: Newton with mu free and amplitude fixed by an extra condition |z0| ~ amp: use arclength-like constraint
    tang = np.zeros(8); tang[:6] = z0 / np.linalg.norm(z0)
    u, J, Mz, ok = po.newton(u, zref, nref, tang=tang, u_pred=u.copy())
    print('start ok', ok, 'mu=%.12g T=%.6f' % (u[7], u[6]))
    t = np.linalg.svd(J)[2][-1]
    if t[:6] @ u[:6] < 0:
        t = -t   # direction of growing amplitude
    out = []
    for step in range(nsteps):
        mults = np.linalg.eigvals(Mz)
        xs = po.xs(u[7]); x0 = xs + E @ u[:6]
        amp_rel = np.linalg.norm(E @ u[:6]) / np.linalg.norm(xs)
        big = np.sort(np.abs(mults))[::-1][:3]
        print('%3d mu=%.10g (rel %+.3e) T=%.6f amp_rel=%.4e min x0 %.3e |mults| %s' % (
            step, u[7], u[7] / mu0 - 1, u[6], amp_rel, x0.min(), np.array2string(big, precision=6)), flush=True)
        out.append(np.concatenate([u, big]))
        np.save('po_%s.npy' % tag, np.array(out))
        u_pred = u + ds * t
        zref = u[:6].copy(); nref = P @ g9(po.xs(u[7]) + E @ u[:6], po.k); nref /= np.linalg.norm(nref)
        un, J, Mzn, ok = po.newton(u_pred.copy(), zref, nref, tang=t, u_pred=u_pred)
        if not ok:
            ds /= 2; print('  halve ds ->', ds, flush=True)
            if ds < 1e-8:
                break
            continue
        tn = np.linalg.svd(J)[2][-1]
        if tn @ t < 0:
            tn = -tn
        u, t, Mz = un, tn, Mzn
    np.save('po_%s.npy' % tag, np.array(out))
