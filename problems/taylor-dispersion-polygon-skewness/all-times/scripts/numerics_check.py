#!/usr/bin/env python3
"""
numerics_check.py -- NOT a certificate.  Floating-point cross-checks of the formulas used by the certificates:
 (1) square: K(tau) from the exact Fourier lattice formula against a Dubiner-Galerkin computation on the half-sector;
 (2) short-time law K ~ <u^3> tau^3/3 + E tau^4/3 with E = <ut>^2 - <u^2>; triangle: <u^3> = 27/19250, E = 9/175, which is
     eq. (14) of Aminian-Camassa-McLaughlin (mu_3 = 3 Pe^3 K ~ (27/19250) Pe^3 t^3 + (9/175) Pe^3 t^4);
 (3) sign of G(a,b) on a logarithmic grid and of H (observed: G > 0 for n >= 4).
Needs ../explore/galerkin.py.  Runtime: about 2 minutes.
"""
import os, sys, math
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'explore'))
from galerkin import Polygon, K_of_tau, H_of_tau


def square_lattice(P=14, Mser=801):
    """û(p,q) = <ut cos(p pi x) cos(q pi y)> on [-1,1]^2 from the double sine-type series of the torsion function."""
    m = np.arange(1, Mser, 2.0)
    s = (-1.0) ** ((m - 1) / 2)

    def Cmat(p):
        # (1/2) int_{-1}^{1} cos(m pi x/2) cos(p pi x) dx
        a = (m / 2 - p) * np.pi
        b = (m / 2 + p) * np.pi
        return 0.5 * (np.sin(a) / a + np.sin(b) / b)
    Cs = {p: Cmat(p) for p in range(0, 2 * P + 1)}
    mm, kk = np.meshgrid(m, m, indexing='ij')
    A = 128 * np.outer(s, s) / (mm * kk * np.pi ** 4 * (mm ** 2 + kk ** 2))
    uh = np.zeros((2 * P + 1, 2 * P + 1))
    for p in range(2 * P + 1):
        Ap = Cs[p] @ A
        for q in range(p, 2 * P + 1):
            uh[p, q] = uh[q, p] = Ap @ Cs[q]
    return uh


def K_square(taus, P=14):
    uh = square_lattice(P)
    xi = [(p, q) for p in range(-P, P + 1) for q in range(-P, P + 1) if (p, q) != (0, 0)]
    pq = np.array(xi)
    lam = np.pi ** 2 * (pq ** 2).sum(axis=1)
    c = uh[np.abs(pq[:, 0]), np.abs(pq[:, 1])]
    g = c / lam
    # T(xi,eta) = û(xi + eta), since u is even and periodic after reflection
    s1 = np.abs(pq[:, 0][:, None] + pq[:, 0][None, :])
    s2 = np.abs(pq[:, 1][:, None] + pq[:, 1][None, :])
    T = uh[s1, s2]
    T[(s1 == 0) & (s2 == 0)] = 0
    out = []
    la = lam[:, None]; lb = lam[None, :]
    for tau in taus:
        d = la - lb
        with np.errstate(divide='ignore', invalid='ignore'):
            E = (np.exp(-lb * tau) - np.exp(-la * tau)) / d
        i = np.abs(d) < 1e-9
        E = np.where(i, tau * np.exp(-la * tau), E)
        I = tau - (1 - np.exp(-la * tau)) / la - (1 - np.exp(-lb * tau)) / lb + E
        out.append(2 * (g @ (T * I) @ g))
    return np.array(out)


taus = np.array([0.015, 0.03, 0.1, 0.3, 1.0, 3.0])
Kl = K_square(taus)
P4 = Polygon(4, 26)
Kg = K_of_tau(P4.lam, P4.ui, P4.Tij, taus)
print("square: tau, K(lattice), K(Galerkin)")
for a, b, c in zip(taus, Kl, Kg):
    print(f"  {a:6.3f}  {b:.6e}  {c:.6e}  rel diff {abs(b-c)/abs(c):.1e}")

print("short-time law: K / (<u^3> tau^3/3 + E tau^4/3) at tau = 1e-4")
for n in [3, 4, 5, 8]:
    P = Polygon(n, 26)
    E = P.c ** 2 - P.u2
    tau = 1e-4
    K = K_of_tau(P.lam, P.ui, P.Tij, np.array([tau]))[0]
    print(f"  n={n}: <u^3>={P.u3:.6e} E={E:.6f} ratio={K/(P.u3*tau**3/3+E*tau**4/3):.4f}")
print(f"  triangle: 27/19250 = {27/19250:.6e}, 9/175 = {9/175:.6f}")

print("sign of G and H on grids")
grid = np.concatenate([[0.0], np.geomspace(1e-5, 3.0, 61)])
for n in [3, 4, 5, 6, 7, 8, 12, 24]:
    P = Polygon(n, 26)
    l = P.lam[1:]; c = P.ui[1:]; T = P.Tij[1:, 1:]
    E = np.exp(-np.outer(grid, l)) * c
    Gm = E @ T @ E.T
    s = np.geomspace(1e-4, 3, 100)
    Hv = H_of_tau(P.lam, P.ui, P.Tij, s)
    print(f"  n={n}: T11 = <u e1^2> = {T[0,0]:+.5f}  #(G<0) = {(Gm<0).sum()} of {Gm.size}  #(H<0) = {(Hv<0).sum()} of {Hv.size}")
