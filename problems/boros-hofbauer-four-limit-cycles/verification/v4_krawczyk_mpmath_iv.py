"""
v4_krawczyk_mpmath_iv.py -- VERIFIER's re-certification of the common zero of (G_1, G_2, G_3) with a
DIFFERENT interval library (mpmath.iv, outward-rounded interval arithmetic) and a box of the verifier's own
choice (own Newton centre, radius 1e-30), Horner-type nested evaluation, own derivative code.

The polynomials G_k are read from ../certificates/G_k.json; v1_focal_real_exact.py proves that they are
exactly the focal-value polynomials (w * G_k(a,c,w^2) = Vhat_k).
Checks:
  K1  Krawczyk: K(X) in int(X)  => unique zero mu* of (G1,G2,G3) in X  (and all of J(X) nonsingular);
  K2  domain conditions on X (a>0, c>0, s>0, rho=1-(c-1)s>0, d>1, 1+bd>0);
  K3  G_4 < 0 on X (so L_4(mu*) < 0) and the enclosure of L_4 = pi G_4 /(a^8 c^16 s^(7/2));
  K4  det d(G1,G2,G3)/d(a,b,d) excludes 0 on X;
  K5  the solver's box (certificates/box.json) and the verifier's box intersect (same zero).
Environment: Python 3.12, mpmath 1.3.0.  Runtime: ~10 s.
"""
import os
import sys
import json
from fractions import Fraction
import mpmath as mp
from mpmath import iv

HERE = os.path.dirname(os.path.abspath(__file__))
CERT = os.path.join(HERE, "..", "certificates")
results = []


def report(name, ok, info=""):
    results.append(bool(ok))
    print(("PASS" if ok else "FAIL") + "  " + name + (("  [" + info + "]") if info else ""), flush=True)


def load(k):
    with open(os.path.join(CERT, "G_%d.json" % k)) as fh:
        return {(t[0], t[1], t[2]): Fraction(t[3]) for t in json.load(fh)["terms"]}


def deriv(P, i):
    out = {}
    for e, cf in P.items():
        if e[i] > 0:
            e2 = list(e)
            e2[i] -= 1
            out[tuple(e2)] = out.get(tuple(e2), 0) + cf * e[i]
    return out


def horner3(P, x, conv):
    """Nested evaluation: sum_i a^i (sum_j c^j (sum_k cf s^k)).  conv maps Fraction -> number type."""
    by_i = {}
    for (i, j, k), cf in P.items():
        by_i.setdefault(i, {}).setdefault(j, {})[k] = cf
    zero = conv(Fraction(0))

    def h1(dct, var):
        deg = max(dct)
        acc = zero
        for e in range(deg, -1, -1):
            acc = acc * var + dct.get(e, zero)
        return acc
    inner_i = {}
    for i, dj in by_i.items():
        inner_j = {}
        for j, dk in dj.items():
            inner_j[j] = h1({k: conv(cf) for k, cf in dk.items()}, x[2])
        inner_i[i] = h1(inner_j, x[1])
    return h1(inner_i, x[0])


def fmt(z, n=30):
    return "[%s, %s]" % (mp.nstr(mp.mpf(z.a.a), n), mp.nstr(mp.mpf(z.b.b), n))


def iv_of(fr):
    return iv.mpf(fr.numerator) / iv.mpf(fr.denominator)


def main():
    G = [load(k) for k in range(1, 5)]
    dG = [[deriv(G[k], i) for i in range(3)] for k in range(3)]
    # ---- own Newton centre (150 digits) from the paper's rounded point
    mp.mp.dps = 150
    a0, b0, d0 = mp.mpf('1.01282'), mp.mpf('0.65463'), mp.mpf('3.28862')
    x = [a0, 1 + b0, (d0 - 1) / (1 + b0 * d0)]
    convf = lambda fr: mp.mpf(fr.numerator) / fr.denominator
    for _ in range(100):
        F = mp.matrix([horner3(G[k], x, convf) for k in range(3)])
        J = mp.matrix([[horner3(dG[k][i], x, convf) for i in range(3)] for k in range(3)])
        dx = mp.lu_solve(J, F)
        x = [x[i] - dx[i] for i in range(3)]
        if mp.norm(dx) < mp.mpf(10) ** -140:
            break
    J = mp.matrix([[horner3(dG[k][i], x, convf) for i in range(3)] for k in range(3)])
    Cm = J ** -1
    # ---- interval part (mpmath.iv), 200 digits
    iv.dps = 200
    rad = iv.mpf(10) ** -30
    m = [iv.mpf(mp.nstr(xi, 120)) for xi in x]          # thin (or tiny) intervals around decimal centres
    m = [iv.mpf(mi.mid) for mi in m]                      # exact binary centres
    X = [iv.mpf([mi.a - rad.b, mi.b + rad.b]) for mi in m]
    C = [[iv.mpf(Cm[i, j]) for j in range(3)] for i in range(3)]
    Gm = [horner3(G[k], m, iv_of) for k in range(3)]
    JX = [[horner3(dG[k][i], X, iv_of) for i in range(3)] for k in range(3)]
    Kx = []
    for i in range(3):
        acc = m[i] - sum((C[i][j] * Gm[j] for j in range(3)), iv.mpf(0))
        for j in range(3):
            Mij = (iv.mpf(1) if i == j else iv.mpf(0)) - sum((C[i][k] * JX[k][j] for k in range(3)), iv.mpf(0))
            acc = acc + Mij * (X[j] - m[j])
        Kx.append(acc)
    inside = all(Kx[i].a > X[i].a and Kx[i].b < X[i].b for i in range(3))
    width = max(float(max(abs(Kx[i].a - m[i].a), abs(Kx[i].b - m[i].b))) for i in range(3))
    report("K1 Krawczyk (mpmath.iv, radius 1e-30): K(X) in int(X)", inside, "max |K_i - m_i| <= %.3e" % width)
    a, c, s = X
    rho = 1 - (c - 1) * s
    b = c - 1
    d = (1 + s) / rho
    ok = a.a > 0 and c.a > 0 and s.a > 0 and rho.a > 0 and (d - 1).a > 0 and (1 + b * d).a > 0
    report("K2 domain on X", ok, "a in %s, b in %s, d in %s" % (fmt(a), fmt(b), fmt(d)))
    G4 = horner3(G[3], X, iv_of)
    L4 = iv.pi * G4 / (a**8 * c**16 * s**3 * iv.sqrt(s))
    report("K3 G_4 < 0 on X", G4.b < 0, "L_4 in %s" % fmt(L4))
    detJ = (JX[0][0] * (JX[1][1] * JX[2][2] - JX[1][2] * JX[2][1])
            - JX[0][1] * (JX[1][0] * JX[2][2] - JX[1][2] * JX[2][0])
            + JX[0][2] * (JX[1][0] * JX[2][1] - JX[1][1] * JX[2][0]))
    dsdd = (1 + b) / (1 + b * d) ** 2
    det_abd = detJ * dsdd
    report("K4 det d(G1,G2,G3)/d(a,b,d) excludes 0 on X", det_abd.b < 0 or det_abd.a > 0,
           "det in %s" % fmt(det_abd))
    with open(os.path.join(CERT, "box.json")) as fh:
        box = json.load(fh)
    sc = [Fraction(v) for v in box["centre"]]
    sr = Fraction(box["radius"])
    inter = all((iv_of(sc[i] - sr).a <= X[i].b) and (iv_of(sc[i] + sr).b >= X[i].a) for i in range(3))
    # the verifier's zero enclosure K(X) must lie in the solver's box (both contain the unique zero)
    inK = all((iv_of(sc[i] - sr).b < Kx[i].a) and (Kx[i].b < iv_of(sc[i] + sr).a) for i in range(3))
    report("K5 verifier's enclosure K(X) of the zero lies inside the solver's Krawczyk box", inter and inK)
    print("ALL PASS" if all(results) else "SOME CHECKS FAILED")
    return all(results)


if __name__ == "__main__":
    sys.exit(0 if main() else 1)
