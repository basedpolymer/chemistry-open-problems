#!/usr/bin/env python3
"""v06_limit_algebra.py -- verifier's re-derivation of the Landau coefficient of the formal limit equation (R).

(R):  G/(kappa(kappa+G)) + int_0^1 s G(T-s) ds = 0.   Multiple scales (re-derived by hand by the verifier):
   A' = sigma kt A + ell |A|^2 A,  sigma = 2/(kappa^3 E'(i w)),  ell = b/(kappa^4 E'(i w)),
   b = 2/D(2 i w) + 4/D(0) - 3,  D(z) = 1 + kappa^2 E(z),  E(z) = int_0^1 s e^{-z s} ds.
Checks (PASS/FAIL):
  L1  closed forms of E(it), E(2it), E'(it), D(0) at the Hopf points (e^{it} = -(1+it)/R, kappa^2 = R-1) agree
      with direct numerical quadrature of the definitions at t_1, t_2, t_3 (mpmath, 30 digits);
  L2  Re(b conj E'(it)) = -(7R^2-34R+23)/(R(R+1)(5R-3)) identically (sympy, rational parametrisation);
  L3  L_j = Re(ell)/t_j equals the paper's closed form and its values for j = 1..5;
  L4  sign: 7R^2-34R+23 > 0 for R >= sqrt(1+t_1^2), t_1 > 3.92 proved by interval evaluation of tan t - t on
      [3.92, 3.92] and monotonicity on (pi, 3pi/2).
Environment: WSL Fedora, sympy, mpmath. Runtime: a few seconds.
"""
import sys
import sympy as sp
import mpmath as mp

fails = 0
def rep(name, ok, s=""):
    global fails
    fails += not ok
    print("%-3s %s %s" % (name, "PASS" if ok else "FAIL", s))

mp.mp.dps = 30
def tj(j):
    c = (2*j - mp.mpf(1)/2)*mp.pi
    return mp.findroot(lambda u: mp.tan(u) - u, c - 1/c)

# L1
worst = 0
for j in (1, 2, 3):
    t = tj(j); R = mp.sqrt(1 + t**2); k2 = R - 1
    assert mp.cos(t) < 0
    Eq = lambda z: mp.quad(lambda s: s*mp.exp(-z*s), [0, 1])
    dEq = lambda z: mp.quad(lambda s: -s**2*mp.exp(-z*s), [0, 1])
    vals = [(Eq(1j*t), -1/k2), (Eq(2j*t), (1 - 1j*t)/(2*R**2)), (dEq(1j*t), 1/R - 1j*(R + 1)/(t*R*(R - 1))),
            (1 + k2*Eq(0), (R + 1)/2), (mp.exp(1j*t), -(1 + 1j*t)/R)]
    worst = max(worst, max(abs(a - b) for a, b in vals))
rep("L1", worst < 1e-25, "max deviation %.1e" % worst)

# L2 symbolic
s = sp.symbols("s", positive=True)
t = (s**2 - 1)/(2*s); R = (s**2 + 1)/(2*s)
I = sp.I
Ep = 1/R - I*(R + 1)/(t*R*(R - 1))
E2 = (1 - I*t)/(2*R**2)
k2 = R - 1
D2 = 1 + k2*E2; D0 = (R + 1)/2
b = 2/D2 + 4/D0 - 3
expr = sp.re(sp.expand_complex(b*sp.conjugate(Ep)))
target = -(7*R**2 - 34*R + 23)/(R*(R + 1)*(5*R - 3))
diff = sp.simplify(sp.together(expr - target))
rep("L2", diff == 0, "Re(b conj E') - target simplifies to %s" % diff)

# L3 numerical L_j
ok3 = True
paper = {1: -0.0093114, 2: -0.0079410, 3: -0.0037777, 4: -0.0021661, 5: -0.0013977}
for j in range(1, 6):
    t = tj(j); R = mp.sqrt(1 + t**2); k = mp.sqrt(R - 1)
    Ep = 1/R - 1j*(R + 1)/(t*R*(R - 1)); E2 = (1 - 1j*t)/(2*R**2)
    bb = 2/(1 + k**2*E2) + 4/((R + 1)/2) - 3
    ell = bb/(k**4*Ep)
    Lj = ell.real/t
    closed = -R*(R - 1)*(7*R**2 - 34*R + 23)/(t*(R + 1)*(5*R - 3)*((R - 1)**3 + R + 1))
    ok3 &= abs(Lj - closed) < 1e-20 and abs(Lj - paper[j]) < 1e-7
    print("     j=%d t_j=%.10f kappa_j^0=%.8f  L_j=%.10f (closed form %.10f, paper %.7f)" % (j, t, k, Lj, closed, paper[j]))
rep("L3", ok3)

# L4
iv = mp.iv
iv.dps = 30
x = iv.mpf("3.92")
val = iv.tan(x) - x
R0 = iv.sqrt(1 + x**2)
thr = (17 + 8*iv.sqrt(2))/7
ok4 = bool(val.b < 0) and bool(R0.a > thr.b) and bool(iv.pi.a < x.a) and bool(x.b < 1.5*iv.pi.a)
rep("L4", ok4, "tan(3.92)-3.92 in %s ; sqrt(1+3.92^2) in %s > (17+8 sqrt2)/7 in %s" % (val, R0, thr))
print("SUMMARY:", "PASS" if fails == 0 else "FAIL")
sys.exit(1 if fails else 0)
