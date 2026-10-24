#!/usr/bin/env python3
"""
numerics_check.py -- NON-RIGOROUS floating-point sanity checks of the rigorous bounds (not part of any proof).

1. For the boundary-layer function H = ut - Y (sg_uniform.py), compare a numerical value of int_P H^2/(n t)
   (high-order harmonic-polynomial approximation, 600x600 Gauss quadrature on the half-sector) with the rigorous
   upper bound Qhat = (vhat + ehat)^2 of Lemma 10, and <u^3> with M(t) = t^4/360 + 19 t^6/7560.
2. Print the numerical S^G and kappa (circumradius 1) next to the asymptotic law S^G ~ (8 sqrt3/15) tan^4(pi/n).
Runtime ~30 s.  Prints OK/NOT-OK (exit code 0 always: these are not certificates).
"""
import math
import numpy as np


def fit_ut(n, K, M=6000):
    h = math.pi / n
    t = math.tan(h)
    S = 1 / math.cos(h)
    s = (1 - np.cos(np.linspace(0, np.pi, M))) / 2
    tau = t * s
    wb = 1 + 1j * tau
    W = np.sqrt(np.gradient(tau))
    A = np.array([np.real((wb / S) ** (n * k)) for k in range(K + 1)]).T
    rhs = tau ** 2 / 2
    a, *_ = np.linalg.lstsq(A * W[:, None], rhs * W, rcond=None)
    return a, S, t


def quad(n, a, S, t, G=600):
    xg, wg = np.polynomial.legendre.leggauss(G)
    xs = (xg + 1) / 2
    wx = wg / 2
    X, Sg = np.meshgrid(xs, xs, indexing='ij')
    WX, WS = np.meshgrid(wx, wx, indexing='ij')
    Wp = X + 1j * t * X * Sg
    Wq = WX * WS * t * X
    r2 = np.abs(Wp) ** 2
    ut = (1 - r2) / 2 + sum(a[k] * np.real((Wp / S) ** (n * k)) for k in range(len(a)))
    Y = (1 - r2) / 2 + t ** 2 / 6
    return ut, Y, Wq


def Qhat(t, K=6):
    R = math.sqrt(1 + t * t)
    om = 1 - t * t / 3
    v = sum(2 * t * t / (k * k * math.pi ** 2) * math.sqrt(t / (math.pi * k) * math.exp(k * math.pi * t / om)) for k in range(1, K + 1))
    e = math.sqrt(1 + R) * (t * t * math.sqrt(2 / (3 * math.pi ** 4 * K ** 3)) +
                            sum(2 * t * t / (k * k * math.pi ** 2) * (k * math.pi * t * t / (3 - t * t) +
                                math.exp(k * math.pi * t / (2 * om)) - 1) for k in range(1, K + 1)))
    return (v + e) ** 2


print(" n   int H^2/(nt) [num]   Qhat bound   <u^3> [num]   M(t)        S^G [num]    (8sqrt3/15)t^4")
for n in [10, 20, 40, 80, 140]:
    a, S, t = fit_ut(n, 30)
    ut, Y, Wq = quad(n, a, S, t)
    area = 2 * n * np.sum(Wq)
    H = ut - Y
    H2 = 2 * n * np.sum(H ** 2 * Wq) / area
    m = 2 * n * np.sum(ut * Wq) / area
    u = ut - m
    u2 = 2 * n * np.sum(u ** 2 * Wq) / area
    u3 = 2 * n * np.sum(u ** 3 * Wq) / area
    M = t ** 4 / 360 + 19 * t ** 6 / 7560
    flag = "OK" if H2 <= Qhat(t) else "NOT-OK"
    print(f"{n:3d}   {H2:.4e}          {Qhat(t):.4e}   {u3:.5e}   {M:.5e}   {u3 / u2 ** 1.5:.5e}  {8 * math.sqrt(3) / 15 * t ** 4:.5e}  {flag}")
