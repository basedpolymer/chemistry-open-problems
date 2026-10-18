#!/usr/bin/env python3
"""Independent verifier check V4: the paper's own validation case (Appendix B, Table B.6, M1 = 4, M2 = 1),
re-implemented from the verifier's reading of Table B.6 (journal p. 48).

Exact univariate computation in u = w1/w0 (sympy rational functions, each exchanger solved with sympy.solve).
Claims (PASS/FAIL):
  B1  corrected CV (economic, prices of Table B.6) == -(2/w0) dJ/du identically in u.
  B2  dJ/du has exactly one root in (0,1); it is in [0.71411, 0.71413] (paper 'Optimized' u = 0.7141).
  B3  the paper's CV (B.1, published recursion, printing typos corrected) has a root u_P in (0.716, 0.718), and
      dJ/du != 0 on the isolating interval of u_P  (an exact M = 4 counterexample to Theorem 3 on the paper's data).
  B4  same with all prices = 1 (Theorem 2): the paper's CV root differs from the root of dQtot/du.
Also prints (no PASS/FAIL, exploratory): J, T_end at both points, and roots of variants of the printed (B.1)
(literal misprints) and of an LMTD model, in an attempt to reproduce the paper's reported 0.7152.
Environment: WSL Fedora, Python 3, sympy 1.14, mpmath. Runtime: < 1 minute.
"""
import sys
import sympy as sp
import mpmath as mp

FAIL = 0
def report(name, ok):
    global FAIL
    print(("PASS " if ok else "FAIL ") + name, flush=True)
    if not ok:
        FAIL += 1

u = sp.Symbol('u')
R = sp.Rational
w0, T0 = 100, 130
UA1 = [5, 7, 10, 12]; wh1 = [50, 30, 15, 25]; Th1 = [190, 203, 220, 235]
UA2, wh2, Th2 = 11, 70, 225
p1 = [R(1), R(6, 5), R(13, 10), R(3, 2)]; p2 = R(3, 2)

def branch(w, UA, wh, Th):
    th, Q = [sp.Integer(0)], []
    for i in range(len(UA)):
        t, to, q = sp.symbols('t to q')
        h = Th[i] - T0
        sol = sp.solve([sp.Eq(q, w*(t - th[-1])), sp.Eq(q, wh[i]*(h - to)),
                        sp.Eq(q, R(UA[i], 2)*((h - t) + (to - th[-1])))], [t, to, q], dict=True)[0]
        th.append(sp.cancel(sol[t])); Q.append(sp.cancel(sol[q]))
    return th, Q

th1, Q1 = branch(w0*u, UA1, wh1, Th1)
th2, Q2 = branch(w0*(1 - u), [UA2], [wh2], [Th2])
hh1 = [T - T0 for T in Th1]; hh2 = Th2 - T0

def recs(th, hh):
    c = [(th[i+1] - th[i])/(hh[i] - th[i]) for i in range(len(hh))]
    ac, ap = [], []
    for i in range(len(hh)):
        ac.append(c[i]*(th[i+1] + th[i] - sum(ac)))
        ap.append(c[i]*(th[i+1] + th[i] - (ap[-1] if ap else 0)))
    return ac, ap

ac1, ap1 = recs(th1, hh1)
a12 = th2[1]**2/hh2
J = -(sum(p*q for p, q in zip(p1, Q1)) + p2*Q2[0])
dJ = sp.cancel(sp.diff(J, u))
cv_corr = sp.cancel(sum(p*a for p, a in zip(p1, ac1)) - p2*a12)
cv_pap = sp.cancel(sum(p*a for p, a in zip(p1, ap1)) - p2*a12)
report("B1 corrected CV == -(2/w0) dJ/du identically", sp.cancel(cv_corr + R(2, w0)*dJ) == 0)

def roots01(expr):
    num, den = sp.fraction(sp.cancel(expr))
    P = sp.Poly(num, u)
    out = []
    for (lo, hi), m in P.intervals(eps=R(1, 10**11)):
        if hi > 0 and lo < 1:
            assert sp.Poly(den, u).eval(lo) != 0 or sp.Poly(den, u).eval(hi) != 0
            out.append((lo, hi))
    return out, P

rJ, PJ = roots01(dJ)
report(f"B2 dJ/du has exactly one root in (0,1), in [0.71411,0.71413]: {[(float(a), float(b)) for a, b in rJ]}",
       len(rJ) == 1 and R(71411, 10**5) < rJ[0][0] and rJ[0][1] < R(71413, 10**5))
rP, PP = roots01(cv_pap)
print("   paper-CV roots in (0,1):", [(float(a), float(b)) for a, b in rP])
inside = [(a, b) for a, b in rP if R(716, 1000) < a and b < R(718, 1000)]
ok = len(inside) == 1
if ok:
    a, b = inside[0]
    # dJ/du has no root in [a,b]: its numerator has no root there (sign check + root count)
    nroots = PJ.count_roots(a, b)
    ok = nroots == 0 and a > rJ[0][1]
report("B3 paper CV vanishes at u_P in (0.716, 0.718) where dJ/du != 0 (exact root counting)", ok)

def Tend(uv):
    return T0 + float(((w0*uv)*th1[-1] + w0*(1 - uv)*th2[-1]).subs(u, uv))/w0
uo = (rJ[0][0] + rJ[0][1])/2
for name, uv in [("optimum", uo)] + [("paper-CV root", (a + b)/2) for a, b in inside]:
    print(f"   {name}: u = {float(uv):.10f}  J = {float(J.subs(u, uv)):.6f}  Tend = {Tend(uv):.6f}")

# B4 prices all 1 (Theorem 2)
Qt = sum(Q1) + Q2[0]
dQ = sp.cancel(sp.diff(Qt, u))
cvP1 = sp.cancel(sum(ap1) - a12)
rq, _ = roots01(dQ); rp1, _ = roots01(cvP1)
print("   prices 1: dQ/du roots", [float(a) for a, b in rq], " paper-CV roots", [float(a) for a, b in rp1])
report("B4 prices 1: no root of the paper CV in (0,1) coincides with the root of dQtot/du (isolating intervals disjoint)",
       all(b < a2 or b2 < a for a, b in rp1 for a2, b2 in rq) and len(rq) >= 1)

# ---------- exploratory: attempts to reproduce the paper's u = 0.7152 ----------
print("\nExploratory (not a certificate):")
t = th1
c1, c2, c3, c4 = [(t[i+1] - t[i])/(hh1[i] - t[i]) for i in range(4)]
# literal printed (B.1): p11 th_h1^2/th_h1; p31-term inner denominator (th_h3 - th_2); p41-term inner th_h1^2/th_h1
a1L = sp.Rational(hh1[0]**2, hh1[0])
a2L = (t[2] - t[1])*(t[2] + t[1] - t[1]**2/hh1[0])/(hh1[1] - t[1])   # printed correctly in (B.1)
a2x = (t[2] - t[1])*(t[2] + t[1] - t[1]**2/hh1[0])
a3L = (t[3] - t[2])*(t[3] + t[2] - a2x/(hh1[2] - t[2]))/(hh1[2] - t[2])
a2y = (t[2] - t[1])*(t[2] + t[1] - hh1[0])/(hh1[1] - t[1])
a3y = (t[3] - t[2])*(t[3] + t[2] - a2y)/(hh1[2] - t[2])
a4L = (t[4] - t[3])*(t[4] + t[3] - a3y)/(hh1[3] - t[3])
cvL = p1[0]*a1L + p1[1]*a2L + p1[2]*a3L + p1[3]*a4L - p2*a12
rL, _ = roots01(cvL)
print("   literal printed (B.1), all roots in (0,1) (exact isolation):", [round(float(a), 6) for a, b in rL])
print("   at the paper's reported u = 0.7152: Tend =", round(Tend(R(7152, 10000)), 4),
      " paper CV =", round(float(cv_pap.subs(u, R(7152, 10000))), 4),
      " corrected CV =", round(float(cv_corr.subs(u, R(7152, 10000))), 4))

# LMTD (counter-current) model, numerical
mp.mp.dps = 30
def lmtd_branch(w, UA, wh, hh):
    th, Q = mp.mpf(0), []
    for i in range(len(UA)):
        N = mp.mpf(UA[i])/w; Rr = w/mp.mpf(wh[i])
        E = mp.e**(-N*(1 - Rr))
        eps = (1 - E)/(1 - Rr*E)
        q = w*eps*(hh[i] - th)
        th = th + eps*(hh[i] - th); Q.append(q)
    return Q
def J_lmtd(uv):
    Qa = lmtd_branch(w0*uv, UA1, wh1, hh1); Qb = lmtd_branch(w0*(1 - uv), [UA2], [wh2], [hh2])
    return -(sum(p*q for p, q in zip([1, 1.2, 1.3, 1.5], Qa)) + 1.5*Qb[0])
uL = mp.findroot(lambda x: mp.diff(J_lmtd, x), 0.7)
print(f"   LMTD model optimum u = {float(uL):.6f}, J = {float(J_lmtd(uL)):.4f}")
print("ALL PASS" if FAIL == 0 else f"{FAIL} FAIL")
sys.exit(1 if FAIL else 0)
