"""
validated_cycles.py -- Theorem C: four limit cycles at EXPLICIT rational parameters, proved by validated
(interval) Taylor integration of the ODE.  Independent of the focal-value computations.

Environment: Python 3.12, python-flint 0.9.0 (Arb).  Runtime: ~1-3 min.  Deterministic.
Usage: python validated_cycles.py        (reads ../certificates/explicit_params.json)

System (BH):  x' = 1/y - x^a y^b,  y' = K(-(d-1) - 1/y + d x^a y^b), with exact rational a, b, d, K.
Section: Sigma = {(1+rho, 1): rho > 0}; on Sigma, y' = K d ((1+rho)^a - 1) > 0 (transversal, upward).

Validated integrator.  Polynomial lift Z = (x, y, E, I, J) with E = x^a y^b, I = 1/y, J = 1/x:
    x' = I - E,  y' = K(-(d-1) - I + d E),  E' = E (a J x' + b I y'),  I' = -I^2 y',  J' = -J^2 x'.
For consistent initial data the lifted solution projects onto the solution of (BH).  One step of size h from an
enclosure Z_k (Arb balls):
  (i)   a priori enclosure B with Z_k + [0,h] F(B) contained in B (Picard-Lindeloef => the solution exists on
        [t_k, t_k+h] and stays in B, for every initial point in Z_k);
  (ii)  Taylor coefficients z_j(Z_k), j < N, by the standard recurrences, evaluated on the ball Z_k;
  (iii) Lagrange remainder z_N(B) tau^N (each component), z_N evaluated on B;
  so Z(t_k + tau) is in sum_{j<N} z_j(Z_k) tau^j + z_N(B) tau^N for every tau in [0,h].
Crossings of the line y = 1 can only occur in steps whose enclosure B meets y = 1; the script requires that in
each such step either (x < 1 and y' < 0 on B) or (x > 1 and y' > 0 on B), so y is strictly monotone there and
the crossings are counted exactly.  The return time is enclosed by bisection on the Taylor model.

Claims checked (PASS/FAIL):
  V1  tr J(1,1) < 0 < det J(1,1) and tr^2 < 4 det (hyperbolic stable focus)  [exact rationals];
  V2  for each test rho_i the first return to Sigma is enclosed, with exactly one crossing of the left
      half-line {x < 1, y = 1} before it, the polar angle around (1,1) (coordinates u = x-1, v = y-1) is
      strictly increasing along the whole arc (u v' - v u' > 0 on every enclosure box), and the displacement
      delta(rho_i) = P(rho_i) - rho_i has the sign  -, +, -, +, -  for rho = 0.002, 0.004, 0.006, 0.012, 0.024.
Mathematical conclusion (notes.md, Sec. "Theorem C"): at least four limit cycles crossing Sigma in
(0.002,0.004), (0.004,0.006), (0.006,0.012), (0.012,0.024); unstable/stable/unstable/stable; (1,1) stable.
"""
import sys
import os
import json
from flint import arb, fmpq, ctx

HERE = os.path.dirname(os.path.abspath(__file__))
CERT = os.path.join(HERE, "..", "certificates")
results = []


def report(name, ok, info=""):
    results.append(bool(ok))
    print(("PASS" if ok else "FAIL") + "  " + name + (("  [" + info + "]") if info else ""), flush=True)


def q(sv):
    """Exact rational from 'p/q', an integer, or a decimal string like '1.25' (exact decimal fraction)."""
    sv = sv.strip()
    if '/' in sv:
        num, _, den = sv.partition('/')
        return fmpq(int(num), int(den))
    if '.' in sv:
        ip, _, fp = sv.partition('.')
        sign = -1 if ip.startswith('-') else 1
        ip = ip.lstrip('-') or '0'
        return sign * fmpq(int(ip + fp), 10**len(fp))
    return fmpq(int(sv))


class System:
    def __init__(self, a, b, d, K, N=36):
        self.aq, self.bq, self.dq, self.Kq = a, b, d, K
        self.a, self.b, self.d, self.K = arb(a), arb(b), arb(d), arb(K)
        self.N = N

    def F(self, Z):
        x, y, E, I, J = Z
        dx = I - E
        dy = self.K * (-(self.d - 1) - I + self.d * E)
        return [dx, dy, E * (self.a * J * dx + self.b * I * dy), -I * I * dy, -J * J * dx]

    def taylor(self, Z, N):
        """Taylor coefficients 0..N of the solution through Z (componentwise lists)."""
        a, b, d, K = self.a, self.b, self.d, self.K
        x, y, E, I, J = ([v] for v in Z)
        dx, dy, I2, J2, W = [], [], [], [], []
        for n in range(N):
            dxn = I[n] - E[n]
            dyn = K * (((1 - d) if n == 0 else arb(0)) - I[n] + d * E[n])
            dx.append(dxn)
            dy.append(dyn)
            I2.append(sum((I[i] * I[n - i] for i in range(n + 1)), arb(0)))
            J2.append(sum((J[i] * J[n - i] for i in range(n + 1)), arb(0)))
            W.append(a * sum((J[i] * dx[n - i] for i in range(n + 1)), arb(0))
                     + b * sum((I[i] * dy[n - i] for i in range(n + 1)), arb(0)))
            dE = sum((E[i] * W[n - i] for i in range(n + 1)), arb(0))
            dI = -sum((I2[i] * dy[n - i] for i in range(n + 1)), arb(0))
            dJ = -sum((J2[i] * dx[n - i] for i in range(n + 1)), arb(0))
            x.append(dxn / (n + 1))
            y.append(dyn / (n + 1))
            E.append(dE / (n + 1))
            I.append(dI / (n + 1))
            J.append(dJ / (n + 1))
        return [x, y, E, I, J]

    def apriori(self, Z, h):
        H = arb(0).union(h)     # ball containing [0, h]
        F0 = self.F(Z)
        # common radius R for all components: B = Z + [-R, R]^5; double R until Z + [0,h] F(B) is in int(B)
        R = 2 * h * max(abs(f.mid()) for f in F0) + arb(10)**-20
        for _ in range(40):
            B = [Z[i] + arb(0, R) for i in range(5)]
            FB = self.F(B)
            Bn = [Z[i] + H * FB[i] for i in range(5)]
            if all(B[i].contains_interior(Bn[i]) for i in range(5)):
                return B
            R = 2 * R
        raise RuntimeError("a priori enclosure failed")

    def step(self, Z, h):
        """Return (coeffs at Z (j < N), remainder coefficient z_N(B), B)."""
        B = self.apriori(Z, h)
        C = self.taylor(Z, self.N - 1)
        RN = [comp[self.N] for comp in self.taylor(B, self.N)]
        return C, RN, B

    @staticmethod
    def evalc(C, RN, tau, N):
        out = []
        for comp, r in zip(C, RN):
            acc = arb(0)
            for co in reversed(comp):
                acc = acc * tau + co
            out.append(acc + r * tau**N)
        return out

    def return_map(self, rho, h0):
        """Validated first return of (1+rho, 1) to Sigma.  Returns (P enclosure, info) or raises."""
        x0 = arb(1 + rho)
        Z = [x0, arb(1), x0 ** self.a, arb(1), 1 / x0]
        h = arb(h0)
        t = arb(0)
        left_runs = 0
        in_left = False
        nsteps = 0
        first = True
        while True:
            C, RN, B = self.step(Z, h)
            nsteps += 1
            Zend = self.evalc(C, RN, h, self.N)
            # angular velocity around (1,1) in the coordinates (u, v) = (x-1, y-1): u v' - v u' > 0 on B
            FB = self.F(B)
            if not ((B[0] - 1) * FB[1] - (B[1] - 1) * FB[0] > 0):
                raise RuntimeError("angular monotonicity not certified")
            yB = B[1]
            meets = yB.contains(arb(1)) or not ((yB < 1) or (yB > 1))
            if meets:
                dyB = self.F(B)[1]
                if B[0] < 1:
                    if not (dyB < 0):
                        raise RuntimeError("left crossing not monotone")
                    if not in_left:
                        left_runs += 1
                        in_left = True
                        if not (Z[1] > 1):
                            raise RuntimeError("left run does not start above y=1")
                elif B[0] > 1:
                    if not (dyB > 0):
                        raise RuntimeError("right crossing not monotone")
                    if in_left:
                        raise RuntimeError("left run not closed")
                    if left_runs == 0 and (first or Z[1] > 1):
                        # start of the orbit: y(t_k) >= 1 and y' > 0 on B, so y > 1 on (t_k, t_k + h]
                        pass
                    elif left_runs == 0:
                        raise RuntimeError("unexpected right-side contact before the half turn")
                    else:
                        if not (Z[1] < 1):
                            raise RuntimeError("right run does not start below y=1")
                        if Zend[1] > 1:
                            # crossing inside this step: bisection on tau
                            lo, hi = arb(0), h
                            for _ in range(140):
                                mid = (lo + hi) / 2
                                ym = self.evalc(C, RN, arb(mid.mid()), self.N)[1]
                                if ym < 1:
                                    lo = arb(mid.mid())
                                elif ym > 1:
                                    hi = arb(mid.mid())
                                else:
                                    break
                            T = lo.union(hi)
                            XT = self.evalc(C, RN, T, self.N)[0]
                            return XT - 1, {"steps": nsteps, "left_runs": left_runs, "t": t + T,
                                            "tau_width": T.rad()}
                        elif Zend[1] < 1:
                            pass  # crossing not yet reached in this step; continue
                        else:
                            # ambiguous end point: retry with half step
                            h = h / 2
                            nsteps -= 1
                            continue
                else:
                    raise RuntimeError("enclosure meets y=1 with undetermined x side")
            else:
                if in_left:
                    # the run of left contacts is over: y must now be below 1 on the whole step
                    if not (yB < 1):
                        raise RuntimeError("left run does not end below y=1")
                    in_left = False
            Z = Zend
            t = t + h
            first = False
            h = arb(h0)
            if nsteps > 2000:
                raise RuntimeError("too many steps")


def main():
    ctx.prec = 320
    with open(os.path.join(CERT, "explicit_params.json")) as fh:
        par = json.load(fh)
    a, b, d, K = (q(par[k]) for k in ("a", "b", "d", "K"))
    tr = -a + K * (1 + b * d)
    det = K * a * (d - 1)
    report("V1 (1,1) is a hyperbolic stable focus: tr J < 0 < det J, tr^2 < 4 det (exact)",
           tr < 0 and det > 0 and tr * tr < 4 * det, "tr J = %.6e" % (float(tr.p) / float(tr.q)))
    report("V1' parameters in the domain a > 0, b > -1, d > 1, 1+bd > 0, K > 0",
           a > 0 and b > -1 and d > 1 and 1 + b * d > 0 and K > 0)
    sysm = System(a, b, d, K, N=int(par.get("order", 36)))
    rhos = [q(r) for r in par["test_rho"]]
    expected = par["expected_signs"]
    h0 = float(par.get("step", 0.05))
    signs = []
    enclosures = []
    for rho, sg in zip(rhos, expected):
        P, info = sysm.return_map(rho, h0)
        delta = P - arb(rho)
        s = 1 if delta > 0 else (-1 if delta < 0 else 0)
        signs.append(s)
        enclosures.append({"rho": str(rho), "delta": delta.str(12), "steps": info["steps"],
                           "return_time": info["t"].str(15)})
        report("V2 rho = %s: first return enclosed, one left crossing, angle increasing, sign(delta) = %+d (expected %+d)"
               % (rho, s, sg), s == sg and info["left_runs"] == 1,
               "delta in %s, %d steps, return time %s" % (delta.str(8), info["steps"], info["t"].str(10)))
    with open(os.path.join(CERT, "validated_cycles_output.json"), "w") as fh:
        json.dump({"params": par, "enclosures": enclosures, "all_pass": all(results)}, fh, indent=1)
    return all(results)


if __name__ == "__main__":
    ok = main()
    print("ALL PASS" if ok else "SOME CHECKS FAILED")
    sys.exit(0 if ok else 1)
