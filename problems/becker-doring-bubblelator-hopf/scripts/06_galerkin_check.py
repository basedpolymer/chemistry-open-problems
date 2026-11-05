"""Independent floating-point check (NOT part of the proof): Fourier-Galerkin continuation of the
bifurcating branch of  delta v' + (1+v)(g_mu * v) = 0  (2 pi-periodic v, g_mu(y) = y^nu exp(-mu y^(beta+1)))
at small amplitudes eps, fit of delta(eps), mu(eps) to second order, and comparison of
eta_2 = (nu+2) mu_2/mu_0 + (beta+1) delta_2/delta_0 and kappa_2 = (nu+1)/(beta+1) mu_2/mu_0 + delta_2/delta_0
with the certified enclosures (from 04_odd_beta3.py and 03_certify_grid.py).
Prints PASS if the agreement is better than 1e-5 (a floating-point consistency test).

Runtime: about 1 minute.  Environment: Python 3, NumPy, SciPy.
"""
import sys
import numpy as np
from scipy.integrate import quad
from scipy.optimize import fsolve

M = 10      # Fourier modes kept


def run(beta, nu, t0):
    def ghat(mu, k):
        f = lambda y: y ** nu * np.exp(-mu * y ** (beta + 1))
        Y = (60 / mu) ** (1 / (beta + 1))
        if k == 0:
            return quad(f, 0, Y, limit=500)[0]
        re = quad(f, 0, Y, weight='cos', wvar=k, limit=500)[0]
        im = -quad(f, 0, Y, weight='sin', wvar=k, limit=500)[0]
        return re + 1j * im

    def resid(x, eps):
        c = np.zeros(M + 1, complex)
        c[0] = x[0]
        c[1] = eps / 2
        for k in range(2, M + 1):
            c[k] = x[1 + 2 * (k - 2)] + 1j * x[2 + 2 * (k - 2)]
        delta, mu = x[-2], x[-1]
        gh = np.array([ghat(mu, k) for k in range(M + 1)])
        n = 64
        full = np.zeros(n, complex); full[:M + 1] = c; full[n - M:] = np.conj(c[1:][::-1])
        gv = np.zeros(n, complex); gv[:M + 1] = gh * c; gv[n - M:] = np.conj((gh * c)[1:][::-1])
        v = (np.fft.ifft(full) * n).real
        w = (np.fft.ifft(gv) * n).real
        prod = np.fft.fft((1 + v) * w) / n
        R = 1j * np.arange(M + 1) * delta * c + prod[:M + 1]
        out = [R[0].real, R[1].real, R[1].imag]
        for k in range(2, M + 1):
            out += [R[k].real, R[k].imag]
        return np.array(out)

    mu0 = 1 / ((beta + 1) * t0 ** (beta + 1))
    delta0 = -ghat(mu0, 1).imag
    x = np.zeros(2 * M + 1); x[-2] = delta0; x[-1] = mu0
    E, D, U = [], [], []
    for eps in [0.02, 0.04, 0.06, 0.08]:
        x = fsolve(resid, x, args=(eps,), xtol=1e-13)
        E.append(eps); D.append(x[-2]); U.append(x[-1])
    E, D, U = map(np.array, (E, D, U))
    d2 = np.polyfit(E ** 2, D, 2)[1]
    m2 = np.polyfit(E ** 2, U, 2)[1]
    return ((nu + 2) * m2 / mu0 + (beta + 1) * d2 / delta0,
            (nu + 1) / (beta + 1) * m2 / mu0 + d2 / delta0)


CASES = [  # beta, nu, t0, certified eta2, certified kappa2 (midpoints)
    (0.0, 1.0, 1.0, -0.0625, -1 / 12),
    (3.0, 0.0, 2.4419679040, -0.43548182, -0.10842900),
    (0.5, 1.0, 1.1701237487, -0.11992936, -0.09765083),
]

ok = True
for beta, nu, t0, e2c, k2c in CASES:
    e2, k2 = run(beta, nu, t0)
    good = abs(e2 - e2c) < 1e-5 and abs(k2 - k2c) < 1e-5
    ok &= good
    print('PASS' if good else 'FAIL', 'beta=%g nu=%g: Galerkin eta2=%.8f kappa2=%.8f, certified %.8f %.8f'
          % (beta, nu, e2, k2, e2c, k2c))
sys.exit(0 if ok else 1)
