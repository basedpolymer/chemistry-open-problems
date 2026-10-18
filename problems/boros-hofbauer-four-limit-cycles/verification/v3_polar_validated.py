"""
v3_polar_validated.py -- VERIFIER's independent validated integrator for the first-return map of (BH),
written from scratch with a DIFFERENT formulation from the solver's validated_cycles.py.

Formulation (polar coordinates in near-normal linear coordinates, angle as independent variable):
  u = x-1, v = y-1;  xi = u,  eta = kappa (a u + c v)  with c = 1+b and kappa = a fixed binary rational close to
  1/omega (omega^2 = det J(1,1) = K a (d-1)); the map (u,v) -> (xi,eta) is linear and orientation preserving.
  xi = R cos(phi), eta = R sin(phi).  Along (BH):  xi' = f,  eta' = kappa (a f + c g),
        f = 1/y - x^a y^b,   g = K(-(d-1) - 1/y + d x^a y^b),
        dR/dphi = R (cos xi' + sin eta') / (cos eta' - sin xi').
  The section Sigma = {(1+rho, 1)} = {v = 0, u > 0} is the ray phi = phi0 := atan(kappa a).  If the denominator
  Dn = cos eta' - sin xi' (= R dphi/dt) is certified > 0 on every a-priori enclosure, phi is strictly increasing
  along the orbit, the orbit is a reparametrisation of the solution R(phi) of this scalar ODE, and the FIRST
  return to Sigma is the point phi = phi0 + 2 pi, i.e.  P(rho) = R(phi0 + 2pi) cos(phi0).
  (Earlier attempt with plain polar coordinates around (1,1) in the (u,v) plane: correct but badly conditioned,
  the enclosures blew up near the major axis of the elliptic linear flow; this version removes that.)

  Validated step (phi_k -> phi_k + h), Arb ball arithmetic throughout:
    (i)   a-priori enclosure: ball B with R_k + [0,h] * F(B, Phi_k) in int(B), Phi_k ball containing
          [phi_k, phi_k + h]  => the solution exists on the step and stays in B (Dn > 0, x > 0, y > 0 checked);
    (ii)  Taylor coefficients R_0..R_{N-1} at (R_k, phi_k) by Picard iteration in arb_series arithmetic;
    (iii) Lagrange remainder: N-th coefficient computed with ball initial data (B, Phi_k), times h^N.
  No crossing detection is needed and the ODE is scalar (no wrapping effect).

Tasks:
  T1 Theorem C of the solver: delta(rho) = P(rho) - rho at rho = 0.002, 0.004, 0.006, 0.012, 0.024 for the exact
     parameters of ../certificates/explicit_params.json -> signs -, +, -, +, - and comparison with the solver's
     enclosures.
  T2 (evidence for Theorem A, independent of the focal-value code): at a rational point exactly on tr J = 0 within
     1e-40 of p*, delta(rho)/rho^9 should tend to L4 (1+1/s)^4 < 0 as rho -> 0 (s = omega^2/a^2).
Environment: Python 3.12, python-flint 0.9.0 (single thread).  Runtime: a few minutes.
Usage: python v3_polar_validated.py [T1|T2|all]
"""
import sys
import os
import json
import time
from flint import arb, arb_series, fmpq, ctx

HERE = os.path.dirname(os.path.abspath(__file__))
CERT = os.path.join(HERE, "..", "certificates")
results = []


def report(name, ok, info=""):
    results.append(bool(ok))
    print(("PASS" if ok else "FAIL") + "  " + name + (("  [" + info + "]") if info else ""), flush=True)


def qdec(sv):
    sv = sv.strip()
    if '/' in sv:
        n, _, d = sv.partition('/')
        return fmpq(int(n), int(d))
    if '.' in sv:
        ip, _, fp = sv.partition('.')
        sign = -1 if ip.startswith('-') else 1
        ip = ip.lstrip('-') or '0'
        return sign * fmpq(int(ip + fp), 10 ** len(fp))
    return fmpq(int(sv))


class PolarNormal:
    def __init__(self, a, b, d, K, N=30, M=256):
        self.a, self.b, self.d, self.K = arb(a), arb(b), arb(d), arb(K)
        self.c = 1 + self.b
        om = (self.K * self.a * (self.d - 1)).sqrt()
        self.kappa = arb((1 / om).mid())            # exact binary number close to 1/omega
        self.phi0 = (self.kappa * self.a).atan()
        self.cos0 = self.phi0.cos()
        self.N, self.M = N, M

    def _fg(self, R, cs, sn, series):
        xi = R * cs
        eta = R * sn
        u = xi
        v = (eta / self.kappa - self.a * xi) / self.c
        x = 1 + u
        y = 1 + v
        if series:
            x0, y0 = x.coeffs()[0], y.coeffs()[0]
        else:
            x0, y0 = x, y
        if not (x0 > 0 and y0 > 0):
            raise RuntimeError("x or y not certified positive")
        E = (self.a * x.log() + self.b * y.log()).exp()
        iy = y.inv() if series else 1 / y
        f = iy - E
        g = self.K * (-(self.d - 1) - iy + self.d * E)
        dxi = f
        deta = self.kappa * (self.a * f + self.c * g)
        Dn = cs * deta - sn * dxi
        D0 = Dn.coeffs()[0] if series else Dn
        if not (D0 > 0):
            raise RuntimeError("angular speed not certified positive")
        return R * (cs * dxi + sn * deta) / Dn

    def F_point(self, R, ph):
        return self._fg(R, ph.cos(), ph.sin(), False)

    def taylor(self, R0, ph0, L):
        """Coefficients 0..L-1 of R(ph0 + tau), R(ph0) = R0 (balls allowed), by Picard iteration."""
        sn, cs = arb_series([ph0, arb(1)], prec=L).sin_cos()
        rs = arb_series([R0], prec=L)
        for _ in range(L):
            Fs = self._fg(rs, cs, sn, True)
            co = Fs.integral().coeffs()[:L]
            co = co + [arb(0)] * (L - len(co))
            co[0] = R0
            rs = arb_series(co, prec=L)
        out = rs.coeffs()[:L]
        return out + [arb(0)] * (L - len(out))

    def apriori(self, Rk, phk, h):
        Ph = phk.union(phk + h)
        H = arb(0).union(h)
        F0 = abs(self.F_point(Rk, Ph))       # over the whole angular step (dR/dphi vanishes on Sigma)
        w = 2 * h * F0.upper() + 2 * Rk.rad() + arb(10) ** -60
        for _ in range(80):
            B = Rk + w * arb(0, 1)
            Bn = Rk + H * self.F_point(B, Ph)
            if B.contains_interior(Bn):
                return B, Ph
            w = 2 * w
        raise RuntimeError("a-priori enclosure failed")

    def return_map(self, rho):
        N, M = self.N, self.M
        h = 2 * arb.pi() / M
        R = arb(rho) / self.cos0
        maxrem = 0.0
        for k in range(M):
            phk = self.phi0 + k * h
            B, Ph = self.apriori(R, phk, h)
            co = self.taylor(R, phk, N)
            rem = self.taylor(B, Ph, N + 1)[N]
            acc = arb(0)
            for cj in reversed(co):
                acc = acc * h + cj
            remterm = rem * h ** N
            maxrem = max(maxrem, float(abs(remterm).upper()))
            R = acc + remterm
        return R * self.cos0, maxrem


def T1(N=32, M=1024):
    ctx.prec = 320
    ctx.cap = 64          # arb_series length cap (the default 10 would silently truncate the series)
    with open(os.path.join(CERT, "explicit_params.json")) as fh:
        par = json.load(fh)
    a, b, d, K = (qdec(par[k]) for k in ("a", "b", "d", "K"))
    with open(os.path.join(CERT, "validated_cycles_output.json")) as fh:
        sol = json.load(fh)
    solver = {e["rho"]: e["delta"] for e in sol["enclosures"]}
    P = PolarNormal(a, b, d, K, N=N, M=M)
    expected = [-1, 1, -1, 1, -1]
    for rs, sg in zip(par["test_rho"], expected):
        rho = qdec(rs)
        t0 = time.time()
        Pr, maxrem = P.return_map(rho)
        delta = Pr - arb(rho)
        s = 1 if delta > 0 else (-1 if delta < 0 else 0)
        sv = arb(solver[str(rho)].replace('[', '').replace(']', '').split('+/-')[0].strip())
        rel = abs(delta - sv) / abs(sv)
        report("T1 rho=%s: sign(delta) = %+d (expected %+d)" % (rs, s, sg), s == sg,
               "delta in %s ; solver %s ; |ours-solver|/|solver| = %s ; max remainder/step %.1e ; %.0fs"
               % (delta.str(14), solver[str(rho)], rel.str(3), maxrem, time.time() - t0))


def T2(N=32, M=1024):
    ctx.prec = 400
    ctx.cap = 64
    with open(os.path.join(CERT, "box.json")) as fh:
        box = json.load(fh)
    ma, mc, ms = (qdec(v) for v in box["centre"])
    b = mc - 1
    rho_ = 1 - (mc - 1) * ms
    d = (1 + ms) / rho_
    K = ma * rho_ / mc
    assert K == ma / (1 + b * d)        # trace exactly zero
    P = PolarNormal(ma, b, d, K, N=N, M=M)
    s = arb(ms)
    pred = arb("-2.5801264994944244863565347627e-6") * (1 + 1 / s) ** 4
    print("     prediction eta_9 = L4 (1+1/s)^4 = %s" % pred.str(12))
    rows = []
    for rs in ("0.0005", "0.001", "0.002", "0.004"):
        rho = qdec(rs)
        Pr, maxrem = P.return_map(rho)
        delta = Pr - arb(rho)
        ratio = delta / arb(rho) ** 9
        rows.append((rs, delta, ratio))
        print("     rho = %-6s delta in %s   delta/rho^9 in %s  (max rem/step %.1e)"
              % (rs, delta.str(10), ratio.str(10), maxrem), flush=True)
    r1, r2 = rows[0][2], rows[1][2]
    extrap = 2 * r1 - r2          # Richardson, assuming delta/rho^9 = eta9 + eta10 rho + O(rho^2)
    ok = all(row[1] < 0 for row in rows) and abs(extrap / pred - 1) < arb("1e-3")
    report("T2 near p* (tr J = 0 exactly): delta < 0 at all test radii and delta/rho^9 -> L4 (1+1/s)^4", ok,
           "Richardson(0.0005, 0.001) = %s vs prediction %s" % (extrap.str(10), pred.str(10)))


def main():
    which = sys.argv[1] if len(sys.argv) > 1 else "all"
    if which in ("T1", "all"):
        T1()
    if which in ("T2", "all"):
        T2()
    print("ALL PASS" if all(results) else "SOME CHECKS FAILED")
    return all(results)


if __name__ == "__main__":
    sys.exit(0 if main() else 1)
