"""
illustrate_cycles.py -- NUMERICAL ILLUSTRATION (not a proof) of the four limit cycles.

Environment: Python 3.12, mpmath 1.3.0, python-flint 0.9.0.  Runtime: ~1-3 min.

1. Choose target radii r_1 < r_2 < r_3 < r_4 (in the time-normalised coordinates X = x - 1 on the ray
   Y = 0, X > 0, i.e. on the segment {x = 1 + rho, y = 1 - (a/(1+b)) rho}).  The truncated displacement map
       delta(r) ~ pi [ (tr/omega) r + V1 r^3 + V2 r^5 + V3 r^7 + V4 r^9 ]
   has these roots iff (tr/omega, V1, V2, V3) = V4 * (e4, -e3, e2, -e1) (elementary symmetric functions
   of the r_i^2).  With V_k = G_k/(a^(2k) c^(4k) s^(k-1/2)) we solve G_k(a,c,s) = target_k by Newton from
   the certified zero, then set K = a/(1+bd) + tr/(1+bd).
2. Integrate the ODE with a high-order Taylor method in 60-digit arithmetic (polynomial lifting with
   E = x^a y^b, I = 1/y, J = 1/x) and compute the displacement delta(rho) = P(rho) - rho at test radii.
   Expected sign pattern (from the inside out): -, +, -, +, -  (stable focus, then Gamma0 unstable,
   Gamma1 stable, Gamma2 unstable, Gamma3 stable).
3. Locate the four zeros of delta by bisection.
"""
import sys
import os
import json
import mpmath as mp

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from certify import poly_from_json, ev_mp  # noqa: E402

mp.mp.dps = 60
CERT = os.path.join(HERE, "..", "certificates")
G = [poly_from_json(json.load(open(os.path.join(CERT, "G_%d.json" % k)))) for k in range(1, 5)]
dG = [[G[k].derivative(i) for i in range(3)] for k in range(3)]


def scal(k, x):
    a, c, s = x
    return a**(2 * k) * c**(4 * k) * s**(k - mp.mpf(1) / 2)


def solve_targets(Vt, x0):
    x = list(x0)
    for _ in range(50):
        F = mp.matrix([ev_mp(G[k], x) - Vt[k] * scal(k + 1, x) for k in range(3)])
        J = mp.matrix([[ev_mp(dG[k][i], x) for i in range(3)] for k in range(3)])
        dx = mp.lu_solve(J, F)   # (scal-dependence neglected in J: quasi-Newton, converges since targets tiny)
        x = [x[i] - dx[i] for i in range(3)]
        if mp.norm(dx) < mp.mpf(10)**-55:
            break
    return x


class TaylorODE:
    def __init__(self, a, b, d, K, order=40, h=mp.mpf('0.08')):
        self.a, self.b, self.d, self.K = a, b, d, K
        self.N, self.h = order, h

    def series(self, st):
        a, b, d, K, N = self.a, self.b, self.d, self.K, self.N
        x, y, E, I, J = ([v] for v in st)
        dx, dy = [], []
        I2, J2, W = [], [], []
        for n in range(N):
            dxn = I[n] - E[n]
            dyn = K * ((-(d - 1) if n == 0 else 0) - I[n] + d * E[n])
            dx.append(dxn)
            dy.append(dyn)
            I2.append(mp.fsum(I[i] * I[n - i] for i in range(n + 1)))
            J2.append(mp.fsum(J[i] * J[n - i] for i in range(n + 1)))
            W.append(a * mp.fsum(J[i] * dx[n - i] for i in range(n + 1)) + b * mp.fsum(I[i] * dy[n - i] for i in range(n + 1)))
            dE = mp.fsum(E[i] * W[n - i] for i in range(n + 1))
            dI = -mp.fsum(I2[i] * dy[n - i] for i in range(n + 1))
            dJ = -mp.fsum(J2[i] * dx[n - i] for i in range(n + 1))
            x.append(dxn / (n + 1))
            y.append(dyn / (n + 1))
            E.append(dE / (n + 1))
            I.append(dI / (n + 1))
            J.append(dJ / (n + 1))
        return x, y, E, I, J

    @staticmethod
    def horner(cf, t):
        r = mp.mpf(0)
        for co in reversed(cf):
            r = r * t + co
        return r

    def return_map(self, rho):
        a, c = self.a, 1 + self.b
        x0, y0 = 1 + rho, 1 - (a / c) * rho
        st = [x0, y0, x0**a * y0**self.b, 1 / y0, 1 / x0]
        t = mp.mpf(0)
        while True:
            S = self.series(st)
            ell = [a * (S[0][n] - (1 if n == 0 else 0)) + c * (S[1][n] - (1 if n == 0 else 0)) for n in range(len(S[0]))]
            l0 = ell[0]
            l1 = self.horner(ell, self.h)
            if t > 2 and l0 < 0 <= l1 and self.horner(S[0], self.h) > 1:
                # Newton for the crossing time inside the step
                dl = [n * ell[n] for n in range(1, len(ell))]
                tau = -l0 / (l1 - l0) * self.h
                for _ in range(100):
                    step = self.horner(ell, tau) / self.horner(dl, tau)
                    tau -= step
                    if abs(step) < mp.mpf(10)**-58:
                        break
                return self.horner(S[0], tau) - 1
            st = [self.horner(comp, self.h) for comp in S]
            t += self.h


def main():
    with open(os.path.join(CERT, "box.json")) as fh:
        box = json.load(fh)
    xstar = [mp.mpf(int(v.split('/')[0])) / int(v.split('/')[1]) for v in box["centre"]]
    V4 = ev_mp(G[3], xstar) / scal(4, xstar)
    radii = [mp.mpf(v) for v in ('0.004', '0.008', '0.016', '0.032')]
    q = [r**2 for r in radii]
    e1 = sum(q)
    e2 = sum(q[i] * q[j] for i in range(4) for j in range(i + 1, 4))
    e3 = sum(q[i] * q[j] * q[k] for i in range(4) for j in range(i + 1, 4) for k in range(j + 1, 4))
    e4 = q[0] * q[1] * q[2] * q[3]
    targets = [-V4 * e3, V4 * e2, -V4 * e1]
    trw = V4 * e4                                  # tr/omega
    x = solve_targets(targets, xstar)
    a, c, s = x
    b = c - 1
    d = (1 + s) / (1 - b * s)
    om = a * mp.sqrt(s)
    tr = trw * om
    K = (a + tr) / (1 + b * d)
    print("parameters (60 digits are needed to reproduce):")
    print("  a =", a)
    print("  b =", b)
    print("  d =", d)
    print("  K =", K, "  (K - a/(1+bd) =", mp.nstr(K - a / (1 + b * d), 5), ")")
    print("  trace =", mp.nstr(tr, 5), " omega =", mp.nstr(om, 10))
    print("  G-targets V1..V3 =", [mp.nstr(t, 5) for t in targets], " V4 =", mp.nstr(V4, 8))
    ode = TaylorODE(a, b, d, K)
    test = [mp.mpf(v) for v in ('0.002', '0.006', '0.012', '0.024', '0.045')]
    signs = []
    for r in test:
        dl = ode.return_map(r) - r
        pred = mp.pi * (trw * r + targets[0] * r**3 + targets[1] * r**5 + targets[2] * r**7 + V4 * r**9)
        signs.append(mp.sign(dl))
        print("  rho = %s : delta = %s   (truncated prediction %s)" % (r, mp.nstr(dl, 6), mp.nstr(pred, 6)), flush=True)
    expected = [-1, 1, -1, 1, -1]
    print("sign pattern", signs, "expected", expected, "->", "OK" if signs == expected else "MISMATCH")
    # locate zeros by bisection
    for lo, hi in zip(test[:-1], test[1:]):
        flo = ode.return_map(lo) - lo
        for _ in range(10):
            mid = (lo + hi) / 2
            fm = ode.return_map(mid) - mid
            if mp.sign(fm) == mp.sign(flo):
                lo, flo = mid, fm
            else:
                hi = mid
        print("  limit cycle crossing the ray at rho ~ %s" % mp.nstr((lo + hi) / 2, 8), flush=True)
    return signs == expected


if __name__ == "__main__":
    sys.exit(0 if main() else 1)
