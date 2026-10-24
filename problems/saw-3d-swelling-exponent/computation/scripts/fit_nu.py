#!/usr/bin/env python3
"""fit_nu.py -- NUMERICS ONLY: estimate nu from pivot-algorithm data (data/pivot_runs.txt) and from the exact series.

Environment: Python 3 with numpy, scipy.  Runtime < 5 s.
Model:  E_n|omega_n|^2 = A n^{2 nu} (1 + a n^{-Delta} + b n^{-1}),  Delta = 0.528 fixed (Clisby 2010 / Clisby-Dunweg
2016 value of the leading correction exponent), weighted least squares on the pivot means with their batch-means
errors; n_min varied to expose systematic drift.  Also: effective exponents from pairs (n, 2n), and the same fits on
the exact series n <= 36 (no statistical error; drift with n_min is the only error indicator).
Reference: nu = 0.58759700(40), Clisby and Dunweg, Phys. Rev. E 94 (2016) 052102.
"""
import math
import os
import sys

import numpy as np
from scipy.optimize import curve_fit

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
NU_REF = 0.58759700
DELTA = 0.528


def load_pivot(path):
    rows = []
    for line in open(path):
        t = line.split()
        if t and t[0] == "R":
            rows.append((int(t[1]), float(t[4]), float(t[5]), float(t[3]), float(t[6])))
    return sorted(rows)


def model1(n, logA, nu, a):
    return logA + 2 * nu * np.log(n) + np.log1p(a * n ** (-DELTA))


def model2(n, logA, nu, a, b):
    return logA + 2 * nu * np.log(n) + np.log1p(a * n ** (-DELTA) + b / n)


def model0(n, logA, nu):
    return logA + 2 * nu * np.log(n)


def fits(ns, y, sig, label):
    print("  %s" % label)
    print("    n_min  model                 nu            chi2/dof")
    for nmin in sorted(set(ns)):
        sel = ns >= nmin
        for name, f, p0 in (("pure power", model0, [0, 0.59]), ("+ a n^-0.528", model1, [0, 0.59, 0.0]),
                            ("+ a n^-0.528 + b/n", model2, [0, 0.59, 0.0, 0.0])):
            k = len(p0)
            if sel.sum() <= k:
                continue
            try:
                p, cov = curve_fit(f, ns[sel], y[sel], p0=p0, sigma=None if sig is None else sig[sel],
                                   absolute_sigma=sig is not None, maxfev=20000)
            except Exception as e:  # noqa: BLE001
                print("    %6d  %-20s fit failed (%s)" % (nmin, name, e))
                continue
            res = (y[sel] - f(ns[sel], *p))
            dof = sel.sum() - k
            chi2 = float(np.sum((res / sig[sel]) ** 2)) / dof if sig is not None else float("nan")
            err = math.sqrt(cov[1, 1]) if sig is not None else float("nan")
            print("    %6d  %-20s %.5f +- %.5f  %.2f" % (nmin, name, p[1], err, chi2))


def main():
    path = sys.argv[1] if len(sys.argv) > 1 else os.path.join(ROOT, "data", "pivot_runs.txt")
    print("NUMERICS (floating point; evidence only).  Reference nu = %.8f (Clisby-Dunweg 2016)." % NU_REF)
    if os.path.exists(path):
        rows = load_pivot(path)
        ns = np.array([r[0] for r in rows], dtype=float)
        m = np.array([r[1] for r in rows])
        e = np.array([r[2] for r in rows])
        print("\nPivot data (%s):" % os.path.relpath(path, ROOT))
        print("     n       E_n|w_n|^2        rel.err    acc.frac  tau_int  E_n/n^{2 nu_ref}")
        for r in rows:
            print("  %7d  %14.4f  %9.2e  %8.4f  %7.1f  %.5f" % (r[0], r[1], r[2] / r[1], r[3], r[4],
                                                              r[1] / r[0] ** (2 * NU_REF)))
        print("  effective exponents from consecutive sizes: nu_eff = log(E_n2/E_n1) / (2 log(n2/n1))")
        for i in range(len(rows) - 1):
            n1, n2 = rows[i][0], rows[i + 1][0]
            v = math.log(rows[i + 1][1] / rows[i][1]) / (2 * math.log(n2 / n1))
            dv = math.hypot(rows[i][2] / rows[i][1], rows[i + 1][2] / rows[i + 1][1]) / (2 * math.log(n2 / n1))
            print("    (%6d, %6d): %.5f +- %.5f" % (n1, n2, v, dv))
        fits(ns, np.log(m), e / m, "weighted fits of log E_n (sigma = relative error)")
    # exact series
    ref = {}
    for line in open(os.path.join(ROOT, "data", "reference_series.tsv")):
        if line.startswith("#"):
            continue
        t = line.split()
        ref[int(t[0])] = (int(t[1]), int(t[2]))
    print("\nExact series n <= 36 (SBB 2011): fits on even n and on odd n separately (parity oscillation),")
    print("  uniform weights sigma = 1e-6 (chi2 meaningless; look at the drift of nu with n_min).")
    for par in (0, 1):
        ns = np.array([n for n in range(10, 37) if n % 2 == par], dtype=float)
        y = np.array([math.log(ref[int(n)][1] / ref[int(n)][0]) for n in ns])
        print(" parity %d:" % par)
        for nmin in (10, 14, 18, 22, 26):
            sel = ns >= nmin
            for name, f, p0 in (("+ a n^-0.528", model1, [0, 0.59, 0.0]), ("+ a n^-0.528 + b/n", model2, [0, 0.59, 0, 0])):
                if sel.sum() <= len(p0):
                    continue
                p, cov = curve_fit(f, ns[sel], y[sel], p0=p0, maxfev=20000)
                print("    n_min %2d  %-20s nu = %.5f" % (nmin, name, p[1]))


if __name__ == "__main__":
    main()
