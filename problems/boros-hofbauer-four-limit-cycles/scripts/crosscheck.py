"""
crosscheck.py -- consistency checks of the exact focal values (not part of the interval certificate,
but they validate the code that produces the polynomials used by the certificate).

Checks (each prints PASS/FAIL; exit code 1 on any FAIL):
  1. V~_1 equals the closed form of Boros-Hofbauer (JDDE 2024, Sec. 5.1) after the substitution
     b = c-1, d = (1+w^2)/(1-(c-1)w^2):  L_1 = pi a (2a^2 d + a(1+bd) - (1+bd)^2) / (8 a sqrt((d-1)(1+bd)^3)),
     i.e. V~_1 = a^2 c^2 w [2a^2(1+w^2)rho + a c rho - c^2]/8 with rho = 1-(c-1)w^2  (exact identity).
  2. L_k := pi V~_k/(a^(2k) c^(4k) w^(2k)) agrees (60 digits) with the independent floating-point
     implementation focal_mp.py (generic code: Taylor expansion + numerical normalisation, no
     reparametrisation) at several parameter points, k = 1..4.
  3. Paper's checkpoints on the curve {L1 = L2 = 0}: at (a,d) = (1, 165/49) (b = 7/11): L2 = 0 exactly and
     L3 < 0; at (a,d) = ((1+sqrt(3961))/60, 3) (b = a - 1/3): L2 = 0 (to 50 digits) and L3 > 0.
Environment: Python 3.12, python-flint 0.9.0, mpmath 1.3.0. Runtime: a few seconds.
"""
import sys
import os
import mpmath as mp
from flint import fmpq

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from focal_exact import lyapunov_exact, CTX  # noqa: E402
from focal_mp import focal_values, bh_system  # noqa: E402

results = []


def report(name, ok, info=""):
    results.append(ok)
    print(("PASS" if ok else "FAIL") + "  " + name + ("  " + info if info else ""))


def evalpoly(V, a, c, w):
    s = mp.mpf(0)
    for e, cf in V.to_dict().items():
        s += mp.mpf(int(cf.p)) / int(cf.q) * a**int(e[0]) * c**int(e[1]) * w**int(e[2])
    return s


def main():
    mp.mp.dps = 60
    Vt = lyapunov_exact(4, verbose=False)
    a, c, w = CTX.gens()
    rho = 1 + w**2 - c * w**2
    target = a**2 * c**2 * w * (2 * a**2 * (1 + w**2) * rho + a * c * rho - c**2) / 8
    report("V~_1 == closed form of the paper (exact polynomial identity)", Vt[0] == target)
    # structure: V~_k = w * (polynomial in a, c, w^2)
    for k, V in enumerate(Vt, 1):
        odd = all(int(e[2]) % 2 == 1 for e in V.to_dict())
        report("V~_%d is w times a polynomial in (a, c, w^2)" % k, odd)
    # 2. agreement with the generic floating-point implementation
    pts = [('1.2', '0.3', '2.5'), ('0.7', '-0.2', '4.0'), ('1.01282', '0.65463', '3.28862'), ('2.5', '1.7', '1.3')]
    for (sa, sb, sd) in pts:
        A, B, Dd = mp.mpf(sa), mp.mpf(sb), mp.mpf(sd)
        Vs, om, tr, _ = focal_values(*bh_system(A, B, Dd), 4, dps=80)
        C = 1 + B
        W = om / A
        worst = mp.mpf(0)
        for k in range(1, 5):
            Lex = mp.pi * evalpoly(Vt[k - 1], A, C, W) / (A**(2 * k) * C**(4 * k) * W**(2 * k))
            Lnum = mp.pi * Vs[k - 1].real
            worst = max(worst, abs(Lex - Lnum) / max(abs(Lnum), mp.mpf(10)**-40))
        report("exact vs generic float implementation at (a,b,d)=(%s,%s,%s)" % (sa, sb, sd),
               worst < mp.mpf(10)**-40, "max rel. diff %s" % mp.nstr(worst, 3))
    # 3. paper checkpoints. Point 1 is rational in (a, c, w^2): a = 1, b = 7/11, d = 165/49, w^2 = 58/77.
    a1, b1, d1 = fmpq(1), fmpq(7, 11), fmpq(165, 49)
    c1 = 1 + b1
    w2 = (d1 - 1) / (1 + b1 * d1)
    report("checkpoint 1: w^2 = 58/77", w2 == fmpq(58, 77))
    # V~_k / w is a polynomial in a, c, w^2: evaluate exactly in Q

    def eval_over_w(V, av, cv, w2v):
        s = fmpq(0)
        for e, cf in V.to_dict().items():
            s += cf * av**int(e[0]) * cv**int(e[1]) * w2v**((int(e[2]) - 1) // 2)
        return s
    v1 = eval_over_w(Vt[0], a1, c1, w2)
    v2 = eval_over_w(Vt[1], a1, c1, w2)
    v3 = eval_over_w(Vt[2], a1, c1, w2)
    report("checkpoint 1 (a,d)=(1,165/49): L1 = 0 exactly", v1 == 0)
    report("checkpoint 1: L2 = 0 exactly", v2 == 0)
    report("checkpoint 1: L3 < 0 (exact rational sign)", v3 < 0, "V~_3/w = %s ~ %.6e" % (v3, float(v3.p) / float(v3.q)))
    # Point 2: a = (1+sqrt(3961))/60, d = 3, b = a - 1/3 (from the L1 = 0 closed form)
    A = (1 + mp.sqrt(3961)) / 60
    Dd = mp.mpf(3)
    B = (-2 + A * (1 + mp.sqrt(1 + 8 * Dd))) / (2 * Dd)
    C = 1 + B
    W = mp.sqrt((Dd - 1) / (1 + B * Dd))
    l1 = evalpoly(Vt[0], A, C, W)
    l2 = evalpoly(Vt[1], A, C, W)
    l3 = evalpoly(Vt[2], A, C, W)
    report("checkpoint 2 (a,d)=((1+sqrt 3961)/60, 3): |V~_1|,|V~_2| < 1e-50", abs(l1) < 1e-50 and abs(l2) < 1e-50,
           "(%s, %s)" % (mp.nstr(l1, 3), mp.nstr(l2, 3)))
    report("checkpoint 2: L3 > 0", l3 > 0, "V~_3 = %s" % mp.nstr(l3, 6))
    return all(results)


if __name__ == "__main__":
    ok = main()
    print("ALL PASS" if ok else "SOME CHECKS FAILED")
    sys.exit(0 if ok else 1)
