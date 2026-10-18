#!/usr/bin/env python3
"""05_limit_coefficient.py -- the M -> infinity limit problem and its Landau (first Lyapunov) coefficients.

Environment: WSL Fedora, Python 3.12, sympy, mpmath, python-flint.  Runtime: about 1 minute.

Limit problem (derived formally in notes.md, Sec. 5), in the Lagrangian clock T:
        G(T) / (kappa (kappa + G(T)))  +  int_0^1 s G(T - s) ds  =  0 .                       (R)
Linearization: Delta(zeta) = 1/kappa^2 + E(zeta),  E(zeta) = int_0^1 s e^{-zeta s} ds, equivalent to the paper's
Q(z;kappa) = e^z (1 + z^2/kappa^2) - (1+z).  Hopf points: zeta = +-i t_j, tan t_j = t_j, cos t_j < 0,
kappa_j^2 = sqrt(1+t_j^2) - 1 (paper (89)-(91)).
Multiple scales give the amplitude equation  A' = sigma kt A + ell |A|^2 A  (kappa = kappa_j + kt) with
        sigma = 2 / (kappa^3 E'(it)),   ell = b / (kappa^4 E'(it)),   b = 2/D(2it) + 4/D(0) - 3,
        D(zeta) = 1 + kappa^2 E(zeta).
Claims (PASS/FAIL):
  C1  with R = sqrt(1+t^2) and e^{it} = -(1+it)/R (consequence of the Hopf conditions), exact identities
        E(it) = -1/kappa^2,  E'(it) = 1/R - i (R+1)/(t R (R-1)),  E(2it) = (1-it)/(2R^2),  D(0) = (R+1)/2;
  C2  Re(b conj(E'(it))) = -(7R^2 - 34R + 23) / (R (R+1) (5R-3))   (exact, rational parametrization);
  C3  the first root t_1 of tan t = t, cos t < 0 satisfies t_1 > 3.92, hence R_j >= R_1 > (17+8 sqrt 2)/7 for all j,
      hence 7R^2-34R+23 > 0, Re ell_j < 0 and Re sigma_j > 0 for EVERY j: all Hopf bifurcations of (R) are
      supercritical;
  C4  numerical values L_j = Re(ell_j)/t_j (prediction for lim l1(q_1 = 1)/M) and a direct simulation of (R):
      cycle amplitude matches the amplitude equation.
"""
import sys
import sympy as sp
import mpmath as mp
import numpy as np
from flint import arb, ctx

fails = 0
def report(name, ok, extra=""):
    global fails
    fails += (not ok)
    print("%-4s %s %s" % (name, "PASS" if ok else "FAIL", extra), flush=True)

# ---------------------------------------------------------------- C1: identities at the Hopf point
t, R = sp.symbols("t R", positive=True)
I = sp.I
kap2 = R - 1                               # kappa^2 = sqrt(1+t^2) - 1
em = -(1 - I*t)/R                          # e^{-it}  (conjugate of e^{it} = -(1+it)/R)
def Ez(z, ez_minus):                       # E(z) = (1 - e^{-z}(1+z))/z^2
    return (1 - ez_minus*(1 + z))/z**2
red = lambda e: sp.simplify(sp.expand(e).subs(t**2, R**2 - 1))
E_it = Ez(I*t, em)
c1a = sp.simplify(sp.together(E_it + 1/kap2).subs(t, sp.sqrt(R**2 - 1))) == 0
# E'(z) = e^{-z}/z - 2 E(z)/z
Ep_it = em/(I*t) - 2*E_it/(I*t)
Ep_claim = 1/R - I*(R + 1)/(t*R*(R - 1))
c1b = sp.simplify((Ep_it - Ep_claim).subs(t, sp.sqrt(R**2 - 1))) == 0
E_2it = Ez(2*I*t, em**2)
c1c = sp.simplify((E_2it - (1 - I*t)/(2*R**2)).subs(t, sp.sqrt(R**2 - 1))) == 0
report("C1", c1a and c1b and c1c, "(E(it)=-1/kappa^2: %s, E'(it): %s, E(2it): %s)" % (c1a, c1b, c1c))
# sanity: E'(z) formula vs direct integral
mp.mp.dps = 40
zz = mp.mpc(0.3, 4.1)
d1 = -mp.quad(lambda s: s**2*mp.e**(-zz*s), [0, 1]); d2 = mp.e**(-zz)/zz - 2*((1 - mp.e**(-zz)*(1 + zz))/zz**2)/zz
report("C1'", abs(d1 - d2) < 1e-25, "(E'(z) = e^-z/z - 2E(z)/z checked at z=0.3+4.1i)")

# ---------------------------------------------------------------- C2: exact real part
s = sp.symbols("s", positive=True)           # t = (s^2-1)/(2s), R = (s^2+1)/(2s), s > 1
sub = {t: (s**2 - 1)/(2*s), R: (s**2 + 1)/(2*s)}
D2 = 1 + kap2*(1 - I*t)/(2*R**2)
D0 = (R + 1)/2
b = 2/D2 + 4/D0 - 3
X = sp.together((b*sp.conjugate(Ep_claim)).subs(sub))
num, den = sp.fraction(X)
num = sp.expand(num*sp.conjugate(den)); den = sp.expand(den*sp.conjugate(den))
ReX = sp.re(num)/den
claim = (-(7*R**2 - 34*R + 23)/(R*(R + 1)*(5*R - 3))).subs(sub)
c2 = sp.simplify(sp.together(ReX - claim)) == 0
report("C2", c2, "Re(b conj E') = -(7R^2-34R+23)/(R(R+1)(5R-3))")

# ---------------------------------------------------------------- C3: sign for every j
ctx.prec = 200
x = arb("3.92")
tanx = (x - arb.pi()).tan()                  # tan is pi-periodic
c3a = bool(tanx - x < 0)                     # tan t - t < 0 at 3.92, increasing on (pi, 3pi/2) -> t_1 > 3.92
Rstar = (17 + 8*arb(2).sqrt())/7             # larger root of 7R^2 - 34R + 23
c3b = bool((1 + x*x).sqrt() > Rstar)
report("C3", c3a and c3b, "t_1 > 3.92 and sqrt(1+3.92^2) = %s > R* = %s; so 7R_j^2-34R_j+23 > 0 for all j"
       % ((1 + x*x).sqrt().str(8), Rstar.str(8)))

# ---------------------------------------------------------------- C4: values and simulation of (R)
mp.mp.dps = 30
def tj(j):
    c = (2*j - mp.mpf(1)/2)*mp.pi
    return mp.findroot(lambda u: mp.tan(u) - u, c - 1/c)
def coeffs(j):
    tt = tj(j); RR = mp.sqrt(1 + tt**2); kap = mp.sqrt(RR - 1)
    Ep = 1/RR - 1j*(RR + 1)/(tt*RR*(RR - 1))
    D2v = 1 + (RR - 1)*(1 - 1j*tt)/(2*RR**2)
    bb = 2/D2v + 8/(RR + 1) - 3
    ell = bb/(kap**4*Ep); sig = 2/(kap**3*Ep)
    L = -RR*(RR - 1)*(7*RR**2 - 34*RR + 23)/(tt*(RR + 1)*(5*RR - 3)*((RR - 1)**3 + RR + 1))
    return tt, RR, kap, ell, sig, L
print("   j      t_j            kappa_j^0        Re sigma      Re ell        L_j = Re ell / t_j   (closed form)")
okL = True
for j in range(1, 11):
    tt, RR, kap, ell, sig, L = coeffs(j)
    okL &= abs(ell.real/tt - L) < 1e-20 and mp.cos(tt) < 0 and abs(mp.tan(tt) - tt) < 1e-20
    print("  %2d  %.10f  %.10f  %+.6e  %+.6e  %+.10f  %+.10f" % (j, tt, kap, sig.real, ell.real, ell.real/tt, L))
report("C4a", okL, "closed form of L_j agrees with Re ell/t_j; cos t_j < 0")

# direct simulation of (R): G(T) = Nk(W(T)), W(T) = -kappa^2 int_0^1 s G(T-s) ds, Nk(W) = W/(1 - W/kappa) is the
# inverse of G -> kappa^2 G/(kappa(kappa+G)).  The integral is evaluated by the trapezoidal rule on a grid of step
# h = 1/nd; the weight of s = 0 vanishes, so the scheme is explicit (second order, no auxiliary ODE variables).
try:
    from numba import njit
except ImportError:                      # pure python fallback (slow)
    njit = lambda f: f

@njit
def _march(G, kap, h, nd, n):
    k2 = kap*kap
    for i in range(nd, n + nd):
        acc = 0.0
        for k in range(1, nd + 1):
            w = 0.5 if k == nd else 1.0
            acc += w*(k*h)*G[i + 1 - k]
        W = -k2*h*acc
        G[i + 1] = W/(1 - W/kap)
    return 0

def simulate(kap, h=2e-3, Tend=8000.0, amp0=1e-3, freq=4.5):
    nd = int(round(1/h)); n = int(Tend/h)
    G = np.zeros(n + nd + 1)
    u = np.arange(-nd, 1)*h
    G[:nd + 1] = amp0*np.cos(freq*u)
    _march(G, kap, h, nd, n)
    return G[nd:], h

tt, RR, kap0, ell, sig, L = coeffs(1)
kap0 = float(kap0)
res = []
for kt in [0.004, 0.002]:
    Gs, h = simulate(kap0 + kt, h=2e-3, Tend=8000.0, amp0=1e-3, freq=float(tt))
    tail = Gs[-int(20/h):]
    amp_sim = 0.5*(tail.max() - tail.min()); assert abs(tail.mean()) < 0.5   # stays on the physical branch
    A2 = -float(sig.real)*kt/float(ell.real)
    amp_pred = 2*np.sqrt(A2)
    res.append(amp_sim/amp_pred)
    print("   simulation of (R), kappa = kappa_1^0 + %.3f : cycle half-amplitude %.5f, predicted 2|A| = %.5f, ratio %.4f"
          % (kt, amp_sim, amp_pred, amp_sim/amp_pred))
okS = all(abs(r - 1) < 0.1 for r in res) and abs(res[1] - 1) < abs(res[0] - 1)
report("C4b", okS, "(stable cycle; amplitude ratio -> 1 as kappa -> kappa_1^0)")
print("SUMMARY:", "PASS" if fails == 0 else "FAIL")
sys.exit(0 if fails == 0 else 1)
