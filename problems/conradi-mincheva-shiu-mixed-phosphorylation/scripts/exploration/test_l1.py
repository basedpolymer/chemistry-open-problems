# Exploration script (floating point; not a certificate). Environment: WSL Fedora, Python 3.12, numpy, scipy, numba.
# Run from scripts/exploration/ (it imports the sibling modules). Typical runtime: seconds to minutes.
"""Check Kuznetsov projection formula implementation against the planar Guckenheimer-Holmes formula."""
import numpy as np

rng = np.random.default_rng(1)
for trial in range(5):
    w = 1.3
    c = rng.normal(size=6)  # f = c0 x^2 + c1 xy + c2 y^2, g = c3 x^2 + c4 xy + c5 y^2
    A = np.array([[0, -w], [w, 0]])

    def B(u, v):
        f = 2 * c[0] * u[0] * v[0] + c[1] * (u[0] * v[1] + u[1] * v[0]) + 2 * c[2] * u[1] * v[1]
        gg = 2 * c[3] * u[0] * v[0] + c[4] * (u[0] * v[1] + u[1] * v[0]) + 2 * c[5] * u[1] * v[1]
        return np.array([f, gg])
    ev, R = np.linalg.eig(A)
    i = np.argmax(ev.imag)
    q = R[:, i] / np.linalg.norm(R[:, i])
    evl, L = np.linalg.eig(A.T)
    j = np.argmin(np.abs(evl - np.conj(ev[i])))
    p = L[:, j]
    p = p / np.conj(np.vdot(p, q))
    s1 = np.linalg.solve(A, B(q, np.conj(q)))
    s2 = np.linalg.solve(2j * w * np.eye(2) - A, B(q, q))
    l1 = np.vdot(p, -2 * B(q, s1) + B(np.conj(q), s2)).real / (2 * w)
    fxx, fxy, fyy = 2 * c[0], c[1], 2 * c[2]
    gxx, gxy, gyy = 2 * c[3], c[4], 2 * c[5]
    a = 1 / (16 * w) * (fxy * (fxx + fyy) - gxy * (gxx + gyy) - fxx * gxx + fyy * gyy)
    print("l1 = %.6f   GH a = %.6f   ratio %.6f" % (l1, a, l1 / a))
