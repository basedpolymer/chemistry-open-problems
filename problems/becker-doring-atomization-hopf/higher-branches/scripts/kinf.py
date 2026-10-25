#!/usr/bin/env python3
"""kinf.py -- the limit K -> infinity: unstable eigenvalues of the limit matrix C(M) = lim B/A and the onsets M*_j.

Usage:  python kinf.py MMAX JMAX
Environment: Python 3.13, python-flint 0.9 (Arb), sympy.  One thread.  Runtime: about 1 min for MMAX = 1200.

Background (result.md, Section 2).  B = A C + E with E independent of K, C the matrix with rows
   (Cv)_1 = -(M+2) v_1 - sum_{l=2}^{N} v_l + M v_M,  (Cv)_l = v_{l-1} - v_l (2 <= l <= N),  (Cv)_M = v_1 + v_N - v_M.
nu is an eigenvalue of C iff phi = 1 + nu is a root of
   G(phi) = -phi^(M-2) (phi^3 + M phi^2 - 2M phi + M) + (M+1) phi - M      (G(1) = G'(1) = 0),
and R := G/(phi - 1)^2 is a polynomial of degree M - 1 with integer coefficients.
Claims (PASS/FAIL):
  K1  (sympy, symbolic M via P = phi^M) lim_{A->inf} F(phi)/A^3 = -phi^3 G(phi) once the (A phi)^-M terms are
      dropped, F = PV's characteristic function (71)-(76); and -phi^2 G = F_inf := phi^M (M (phi-1)^2 + phi^3)
      - phi^2 ((M+1) phi - M);  det(nu I - C(M)) = -G(1+nu)/nu checked exactly for M = 3..12.
  K2  for every integer 3 <= M <= MMAX: R(M) has no root on Re phi = 1, and the number of its roots in Re phi > 1
      (argument principle for G on the half-disc {Re phi > 1, |phi - 1| < 4}, indented to the left around the
      double zero phi = 1 by a half-disc of radius 0.02/M^2 on which R != 0 is certified; the half-disc contains every
      root with Re phi >= 1: Lemma 2.2 of result.md) equals 2 u(M); u(M) = number of eigenvalue pairs of C(M) in
      Re nu > 0.
  K3  u(M) is nondecreasing on [3, MMAX], and it jumps exactly at the integers M_j = 15, 44, 87, 142, ... (printed),
      with u(M) = #{j : M >= M_j}.
  K4  for j <= JMAX: Krawczyk enclosure of the unique zero (M*_j, eta_j) of G(1 + i eta; M) (real M, phi^M =
      exp(M log phi)) near the float guess, with M_j - 1 < M*_j < M_j, and d Re(phi)/dM > 0 there (the j-th pair
      crosses Re nu = 0 from left to right as M increases).
"""
import sys, os, time, math
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from flint import arb, acb, ctx, fmpz_poly, acb_poly, arb_mat
import mpmath as mp

fails = 0
def report(name, ok, extra=""):
    global fails
    fails += (not ok)
    print("%-4s %s %s" % (name, "PASS" if ok else "FAIL", extra), flush=True)

# ------------------------------------------------------------------------------------------------ K1
def check_K1():
    import sympy as sp
    A, phi, M, P, nu = sp.symbols("A phi M P nu")
    K = A - 1
    # PV (71)-(75); P = phi^M, the (A phi)^-M terms set to 0
    f1 = (phi - 1) + P * A
    f2 = (1 - A * phi)
    G1 = phi * (1 - A * phi) * (K * phi + (K * M + 1) * (phi - 1))
    G2 = (1 - A * phi) * (A * M * (phi - 1) + K * phi + A * phi * (phi - 1))
    H1 = (phi - 1) * (K + (K * M + 1) * (1 - A * phi))
    F = f1 * H1 - (G1 - P * G2) * f2
    lead = sp.limit(sp.expand(F) / A ** 3, A, sp.oo)
    Gs = -P / phi ** 2 * (phi ** 3 + M * phi ** 2 - 2 * M * phi + M) + (M + 1) * phi - M
    ok1 = sp.simplify(lead + phi ** 3 * Gs) == 0
    Finf = P * (M * (phi - 1) ** 2 + phi ** 3) - phi ** 2 * ((M + 1) * phi - M)
    ok2 = sp.simplify(-phi ** 2 * Gs - Finf) == 0
    ok3 = True
    for m in range(3, 13):
        Cm = sp.zeros(m, m)
        Cm[0, 0] = -(m + 2)
        for l in range(1, m - 1):
            Cm[0, l] = -1
        Cm[0, m - 1] = m
        for l in range(1, m - 1):
            Cm[l, l - 1] = 1; Cm[l, l] = -1
        Cm[m - 1, 0] += 1; Cm[m - 1, m - 2] += 1; Cm[m - 1, m - 1] += -1
        cp = (nu * sp.eye(m) - Cm).det()
        g = (-(1 + nu) ** (m - 2) * ((1 + nu) ** 3 + m * (1 + nu) ** 2 - 2 * m * (1 + nu) + m) + (m + 1) * (1 + nu) - m)
        ok3 &= sp.expand(cp * nu + g) == 0
    report("K1", ok1 and ok2 and ok3, "lim F/A^3 = -phi^3 G: %s; -phi^2 G = F_inf: %s; det(nu - C) = -G(1+nu)/nu (M=3..12): %s"
           % (ok1, ok2, ok3))

# ------------------------------------------------------------------------------------------------ K2
def Rpoly(M):
    c = [0] * (M + 2)          # G as integer polynomial, coefficient of phi^k at index k
    c[M + 1] -= 1; c[M] -= M; c[M - 1] += 2 * M; c[M - 2] -= M; c[1] += M + 1; c[0] -= M
    G = fmpz_poly(c)
    R, rem = divmod(G, fmpz_poly([1, -2, 1]))
    assert rem == 0 and R.degree() == M - 1
    return R

def Gval(M, z):
    """G(z) = -z^(M-2) (z^3 + M z^2 - 2M z + M) + (M+1) z - M and G'(z), closed form (acb balls)"""
    q = z ** 3 + M * z ** 2 - 2 * M * z + M
    p = z ** (M - 3)
    G = -p * z * q + (M + 1) * z - M
    dG = -p * ((M - 2) * q + z * (3 * z ** 2 + 2 * M * z - 2 * M)) + (M + 1)
    return G, dG

def winding(M, pts, maxdepth=40):
    """number of zeros of G inside the closed polygon pts (ccw); each edge is subdivided until the mean-value
    enclosure of G over the edge excludes 0; returns (count, nseg) or (None, info)"""
    total = 0.0; nseg = 0
    for a, b in zip(pts[:-1], pts[1:]):
        stack = [(acb(a), acb(b), 0)]
        while stack:
            z0, z1, d = stack.pop()
            zc = (z0 + z1) / 2
            ball = z0.union(z1)
            enc = Gval(M, zc)[0] + Gval(M, ball)[1] * (ball - zc)
            if enc.contains(0) or not enc.is_finite():
                if d >= maxdepth:
                    return None, ("zero not excluded", complex(z0.mid()))
                zm = acb(zc.mid())
                stack.append((zm, z1, d + 1)); stack.append((z0, zm, d + 1))
                continue
            v0 = Gval(M, z0)[0]; v1 = Gval(M, z1)[0]
            a0 = math.atan2(float(v0.imag.mid()), float(v0.real.mid()))
            a1 = math.atan2(float(v1.imag.mid()), float(v1.real.mid()))
            da = a1 - a0
            while da > math.pi: da -= 2 * math.pi
            while da < -math.pi: da += 2 * math.pi
            for v in (v0, v1):          # endpoint values precise enough for the principal difference
                if not (v.abs_lower() > 0) or float(arb(acb(arb(0, v.real.rad()), arb(0, v.imag.rad())).abs_upper())) > 1e-6 * float(v.abs_lower()):
                    return None, ("imprecise endpoint", complex(z0.mid()))
            total += da; nseg += 1
    w = total / (2 * math.pi); n = round(w)
    if abs(w - n) > 1e-3:
        return None, ("non-integer winding", w)
    return n, nseg

RHO = 1e-3

def contour(H=4, rho=RHO, n_arc=64, n_line=64, n_small=8):
    """ccw boundary of {Re phi > 1, |phi - 1| < H} united with the small left half-disc {|phi - 1| < rho, Re phi <= 1}"""
    pts = []
    for k in range(n_arc + 1):        # arc from 1 - iH to 1 + iH through 1 + H
        th = -math.pi / 2 + math.pi * k / n_arc
        pts.append(complex(1 + H * math.cos(th), H * math.sin(th)))
    pts[0] = complex(1, -H); pts[-1] = complex(1, H)
    for k in range(1, n_line + 1):    # line Re phi = 1 downwards to 1 + i rho
        pts.append(complex(1, H - (H - rho) * k / n_line))
    for k in range(1, n_small):       # small semicircle to the left, from 1 + i rho to 1 - i rho
        th = math.pi / 2 + math.pi * k / n_small
        pts.append(complex(1 + rho * math.cos(th), rho * math.sin(th)))
    pts.append(complex(1, -rho))
    for k in range(1, n_line + 1):    # line downwards to 1 - iH
        pts.append(complex(1, -rho - (H - rho) * k / n_line))
    pts[-1] = complex(1, -H)
    return pts

def small_disc_ok(M, rho=RHO):
    """R = G/(phi-1)^2 has no zero in |phi - 1| <= rho (R(1) + R'(disc) disc excludes 0); R from the exact
    integer polynomial"""
    R = Rpoly(M)
    Rc = acb_poly([acb(int(c)) for c in R.coeffs()])
    disc = acb(arb(1, rho), arb(0, rho))
    enc = Rc(acb(1)) + Rc.derivative()(disc) * (disc - 1)
    return not enc.contains(0)

def u_count(M):
    """2 u(M) = number of zeros of R in Re phi > 1 (all such zeros have |phi - 1| < 4, Lemma 2.2)"""
    rho = 0.02 / M ** 2
    if not small_disc_ok(M, rho):
        return None, "small disc"
    n, info = winding(M, contour(rho=rho))
    if n is None:
        return None, info
    return n - 2, info          # minus the double zero phi = 1

# ------------------------------------------------------------------------------------------------ K4
def Gline(M, eta):
    """G(1 + i eta; M) for real M (balls), with d/dM and d/deta"""
    phi = 1 + acb(0, 1) * eta
    lp = phi.log()
    pM2 = ((M - 2) * lp).exp()
    q = phi ** 3 + M * phi ** 2 - 2 * M * phi + M
    G = -pM2 * q + (M + 1) * phi - M
    dGdphi = -pM2 * ((M - 2) / phi * q + 3 * phi ** 2 + 2 * M * phi - 2 * M) + (M + 1)
    dGdM = -lp * pM2 * q - pM2 * (phi - 1) ** 2 + phi - 1
    return G, dGdM, dGdphi * acb(0, 1), dGdphi

def krawczyk2(Mm, em, rM, re):
    box = [arb(Mm, rM), arb(em, re)]
    def fj(M, e):
        G, gM, ge, _ = Gline(M, e)
        return [G.real, G.imag], [[gM.real, ge.real], [gM.imag, ge.imag]]
    f0, J0 = fj(arb(Mm), arb(em))
    _, JX = fj(box[0], box[1])
    Y = arb_mat([[c.mid() for c in r] for r in J0]).inv()
    Y = arb_mat([[Y[i, k].mid() for k in range(2)] for i in range(2)])
    I2 = arb_mat([[1, 0], [0, 1]])
    dev = arb_mat([[box[0] - Mm], [box[1] - em]])
    Kx = arb_mat([[arb(Mm)], [arb(em)]]) - Y * arb_mat([[f0[0]], [f0[1]]]) + (I2 - Y * arb_mat(JX)) * dev
    if box[0].contains_interior(Kx[0, 0]) and box[1].contains_interior(Kx[1, 0]):
        return [Kx[0, 0], Kx[1, 0]]
    return None

def main():
    Mmax = int(sys.argv[1]) if len(sys.argv) > 1 else 1200
    Jmax = int(sys.argv[2]) if len(sys.argv) > 2 else 13
    ctx.prec = 128
    t0 = time.time()
    check_K1()
    # K2, K3
    us = {}; okK2 = True; maxseg = 0
    for M in range(3, Mmax + 1):
        n, info = u_count(M)
        if n is None or n % 2:
            okK2 = False; print("K2 FAIL at M=%d: %s" % (M, info), flush=True); continue
        us[M] = n // 2; maxseg = max(maxseg, info)
    report("K2", okK2, "u(M) computed for M = 3..%d (max %d segments, %.0fs)" % (Mmax, maxseg, time.time() - t0))
    jumps = [M for M in range(4, Mmax + 1) if us.get(M) != us.get(M - 1)]
    mono = all(us[M] >= us[M - 1] for M in range(4, Mmax + 1) if M in us and M - 1 in us)
    unit = all(us[M] - us[M - 1] in (0, 1) for M in range(4, Mmax + 1) if M in us and M - 1 in us)
    okK3 = okK2 and mono and unit and us[3] == 0
    report("K3", okK3, "u(3) = %d, u nondecreasing with unit jumps: %s; onsets M_j = %s" % (us.get(3), mono and unit, jumps))
    # K4
    ctx.prec = 192
    mp.mp.dps = 40
    def Gm(M, e):
        p = 1 + 1j * e
        return -mp.exp((M - 2) * mp.log(p)) * (p ** 3 + M * p ** 2 - 2 * M * p + M) + (M + 1) * p - M
    okK4 = True
    for j in range(1, min(Jmax, len(jumps)) + 1):
        Mj = jumps[j - 1]
        # float guess: M0 in (Mj - 1, Mj); eta from a scan of |G| / |phi^M| on the line Re phi = 1 at M0, then Newton
        best = None
        mp.mp.dps = 80
        for M0 in (Mj - 0.5, Mj - 0.1, Mj - 0.9, Mj - 0.02, Mj - 0.98):
            M0 = mp.mpf(M0)
            scan = None
            for i in range(4001):
                zz = 2 * mp.pi * j - 3 + 3.5 * i / 4000
                e = zz / M0
                v = abs(Gm(M0, e)) / abs(mp.exp(M0 * mp.log(1 + 1j * e)))
                if scan is None or v < scan[0]:
                    scan = (v, e)
            try:
                r = mp.findroot(lambda M, e: [mp.re(Gm(M, e)), mp.im(Gm(M, e))], (M0, scan[1]), verify=False)
                r = mp.findroot(lambda M, e: [mp.re(Gm(M, e)), mp.im(Gm(M, e))], (r[0], r[1]), verify=False)
                if Mj - 1 < r[0] < Mj and abs(Gm(r[0], r[1])) < mp.mpf(10) ** -50:
                    best = r; break
            except Exception:
                pass
        if best is None:
            okK4 = False; print("K4 j=%d: no float zero found" % j); continue
        Mm, em = arb(mp.nstr(best[0], 50)), arb(mp.nstr(best[1], 50))
        enc = krawczyk2(Mm, em, arb("1e-25"), arb("1e-25"))
        if enc is None:
            okK4 = False; print("K4 j=%d: Krawczyk failed" % j); continue
        G, gM, ge, gphi = Gline(enc[0], enc[1])
        dphidM = -gM / gphi                          # d phi/dM along the root, phi = 1 + i eta + ...
        tr = dphidM.real
        okj = bool(tr > 0) and bool(enc[0] > Mj - 1) and bool(enc[0] < Mj)
        okK4 &= okj
        print("K4   j=%d: M*_j in %s, eta_j in %s, z = M eta = %s, M_j = %d, M_j - M*_j = %.6f, dRe(phi)/dM > %.3g %s"
              % (j, enc[0].str(20), enc[1].str(15), (enc[0] * enc[1]).str(10), Mj, float((Mj - enc[0]).mid()),
                 float(tr.lower()), "PASS" if okj else "FAIL"), flush=True)
    report("K4", okK4, "onsets enclosed for j = 1..%d" % min(Jmax, len(jumps)))
    print("SUMMARY:", "PASS" if fails == 0 else "FAIL", "(%.0fs)" % (time.time() - t0))
    sys.exit(0 if fails == 0 else 1)

if __name__ == "__main__":
    main()
