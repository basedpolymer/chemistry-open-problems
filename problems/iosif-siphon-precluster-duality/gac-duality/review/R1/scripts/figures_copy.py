#!/usr/bin/env python3
"""Numerical illustration for Theorem O of the manuscript (result 17): the periodic orbit of the dual of
3X3 -> X1+X2 -> X2+2X3+3X4 -> 3X3 at the rates (alpha, beta, gamma, delta) = (17, 400, 70000, 10^6) in the class
a + b + c = T = 201/100.

THIS IS AN ILLUSTRATION, NOT A PROOF.  Theorem O is proved exactly (scripts/certify_duality.py); nothing in the
proof uses this figure.

Method: the field a' = alpha b + beta b c - 3 gamma a^3 c^2, b' = -alpha b + 2 gamma a^3 c^2 + 3 delta c^3,
c' = -beta b c + gamma a^3 c^2 - 3 delta c^3 is integrated in logarithmic coordinates (log a, log b, log c), with the
stiff implicit Radau solver of scipy (rtol 1e-11, atol 1e-12 in the logarithms).  In linear coordinates the solver
wrongly appears to converge to the vertex (T,0,0), because the orbit comes within about 1e-12 of it (verification
report V3).  The conservation law a + b + c = T is monitored (printed drift).

Environment: Python 3.12, numpy, scipy, matplotlib (WSL Fedora 44).  Deterministic (fixed initial point, fixed
tolerances).  Runtime: about 1-3 minutes.
Usage (from the repository root, in WSL):  python problems/iosif-siphon-precluster-duality/scripts/figures.py
R1 copy: writes gac-duality/review/R1/figures/theoremO-orbit.pdf (output folder = sibling of this scripts/ folder) and prints a short report.
"""
import os
import time
import platform

import numpy as np
import scipy
from scipy.integrate import solve_ivp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

AL, BE, GA, DE = 17.0, 400.0, 70000.0, 1.0e6
T = 2.01
EQ = (1.0, 1.0, 0.01)


def rhs(t, y):
    a, b, c = np.exp(y)
    da = AL * b + BE * b * c - 3 * GA * a ** 3 * c ** 2
    db = -AL * b + 2 * GA * a ** 3 * c ** 2 + 3 * DE * c ** 3
    dc = -BE * b * c + GA * a ** 3 * c ** 2 - 3 * DE * c ** 3
    return [da / a, db / b, dc / c]


def jac(t, y):
    a, b, c = np.exp(y)
    # d(f_i/x_i)/d(log x_j) = x_j * d(f_i/x_i)/dx_j
    fa = AL * b / a + BE * b * c / a - 3 * GA * a ** 2 * c ** 2
    fb = -AL + 2 * GA * a ** 3 * c ** 2 / b + 3 * DE * c ** 3 / b
    fc = -BE * b + GA * a ** 3 * c - 3 * DE * c ** 2
    J = np.zeros((3, 3))
    J[0, 0] = a * (-AL * b / a ** 2 - BE * b * c / a ** 2 - 6 * GA * a * c ** 2)
    J[0, 1] = b * (AL / a + BE * c / a)
    J[0, 2] = c * (BE * b / a - 6 * GA * a ** 2 * c)
    J[1, 0] = a * (6 * GA * a ** 2 * c ** 2 / b)
    J[1, 1] = b * (-2 * GA * a ** 3 * c ** 2 / b ** 2 - 3 * DE * c ** 3 / b ** 2)
    J[1, 2] = c * (4 * GA * a ** 3 * c / b + 9 * DE * c ** 2 / b)
    J[2, 0] = a * (3 * GA * a ** 2 * c)
    J[2, 1] = b * (-BE)
    J[2, 2] = c * (GA * a ** 3 - 6 * DE * c)
    del fa, fb, fc
    return J


def main():
    t0 = time.time()
    print(f"# figures.py  Python {platform.python_version()}, numpy {np.__version__}, scipy {scipy.__version__}, "
          f"matplotlib {matplotlib.__version__}")
    # start close to the unstable equilibrium, in the same class
    x0 = np.array([1.0, 1.0 - 1e-4, 0.01 + 1e-4])
    assert abs(x0.sum() - T) < 1e-15
    tend = 1.2e7
    sol = solve_ivp(rhs, (0.0, tend), np.log(x0), method="Radau", jac=jac, rtol=1e-11, atol=1e-12,
                    dense_output=False, max_step=np.inf)
    a, b, c = np.exp(sol.y)
    drift = np.max(np.abs(a + b + c - T))
    print(f"# status {sol.status}, {sol.t.size} steps, max |a+b+c-T| = {drift:.2e}")
    print(f"# min b = {b.min():.3e}, min c = {c.min():.3e}, max a = {a.max():.15f} (T = {T})")
    # crude period estimate from successive maxima of log c after the transient
    lc = np.log10(c)
    late = sol.t > 1e6
    tl, cl = sol.t[late], lc[late]
    peaks = [tl[i] for i in range(1, len(cl) - 1) if cl[i] > cl[i - 1] and cl[i] >= cl[i + 1] and cl[i] > -3]
    if len(peaks) >= 2:
        print(f"# successive maxima of log10 c at t = {[f'{p:.4e}' for p in peaks[:6]]}")
        print(f"# period estimate {np.diff(peaks).mean():.4e}")
    plt.rcParams.update({"pdf.fonttype": 42, "font.size": 8, "font.family": "DejaVu Sans"})
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(7.3, 2.6))
    # panel (a): log-log plane (b, c); the vertex (T,0,0) is at (-inf, -inf)
    ax1.plot(np.log10(b), np.log10(c), color="0.6", lw=0.4, label="whole trajectory")
    ax1.plot(np.log10(b[late]), np.log10(c[late]), color="black", lw=0.9, label="periodic orbit ($t>10^6$)")
    ax1.plot([0.0], [-2.0], marker="o", ms=4, mfc="white", mec="black", ls="none", label="equilibrium")
    ax1.set_xlabel(r"$\log_{10} b$")
    ax1.set_ylabel(r"$\log_{10} c$")
    ax1.set_title("(a) orbit, logarithmic axes", fontsize=8)
    ax1.legend(fontsize=6.5, loc="upper left", frameon=False)
    # panel (b): the spiral leaving the hyperbolic source, linear coordinates
    sp = sol.t < 0.45
    ax2.plot(b[sp] - 1.0, c[sp] - 0.01, color="black", lw=0.5)
    ax2.plot([0.0], [0.0], marker="o", ms=3, mfc="white", mec="black", ls="none")
    ax2.set_xlabel(r"$b-1$")
    ax2.set_ylabel(r"$c-1/100$")
    ax2.set_title(r"(b) spiral, $0 \leq t \leq 0.45$", fontsize=8)
    ax2.ticklabel_format(style="sci", scilimits=(-2, 2))
    # panel (c): one period of the orbit
    win = (sol.t > 2.6e6) & (sol.t < 5.6e6)
    ax3.plot(sol.t[win] / 1e6, np.log10(c[win]), color="black", lw=0.8, label=r"$\log_{10} c$")
    ax3.plot(sol.t[win] / 1e6, np.log10(b[win]), color="0.5", lw=0.8, ls="--", label=r"$\log_{10} b$")
    ax3.set_xlabel(r"$t$ ($10^6$)")
    ax3.set_title(r"(c) one period (about $2.5\cdot 10^6$)", fontsize=8)
    ax3.legend(fontsize=6.5, loc="center right", frameon=False)
    fig.tight_layout()
    here = os.path.dirname(os.path.abspath(__file__))
    out = os.path.join(os.path.dirname(here), "figures")
    os.makedirs(out, exist_ok=True)
    path = os.path.join(out, "theoremO-orbit.pdf")
    fig.savefig(path, metadata={"Creator": "figures.py", "Producer": "matplotlib"})
    print(f"# wrote figures/theoremO-orbit.pdf; runtime {time.time() - t0:.1f} s")


if __name__ == "__main__":
    main()
