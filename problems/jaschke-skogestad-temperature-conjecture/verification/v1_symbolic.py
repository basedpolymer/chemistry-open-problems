#!/usr/bin/env python3
"""Independent verifier check V1 (written from scratch, no code shared with the solver's scripts).

Model of ONE branch of the Jaschke-Skogestad network (Comput. Chem. Eng. 70 (2014) 35-49, eqs (9), (10), (12)):
for exchanger i (shifted temperatures theta = T - T0, theta_0 = 0)
    Q_i = w (th_i - th_{i-1})                                  (9)
    Q_i = wh_i (thh_i - tho_i)                                 (10)  [dropped when wh_i = oo, tho_i = thh_i]
    Q_i = UA_i/2 * ((thh_i - th_i) + (tho_i - th_{i-1}))       (12)  counter-current AMTD
The three equations are solved with sympy.solve (no hand-derived closed form), then dQ_i/dw is obtained by sympy.diff.

Claims checked (PASS/FAIL):
  C1(M, hot)  2 dQ_i/dw == a_i^corr  for every i <= M, as an identity of rational functions in all parameters, where
              a_i^corr = c_i (th_i + th_{i-1} - sum_{l<i} a_l^corr),  c_i = (th_i - th_{i-1})/(thh_i - th_{i-1}).
  C2(M, hot)  the published recursion (18) a_i^P = c_i (th_i + th_{i-1} - a_{i-1}^P) satisfies the same identity for i<=2
              and FAILS for i = 3 (if M >= 3), and the branch total sum a^P != 2 d(sum Q)/dw for M >= 3.
  C3          single exchanger: dJ/dw = -p th_1^2/(2 thh_1)  (factor -1/2 relative to the paper's (16)/(A.7)).
Fully symbolic for: finite hot capacities M = 1..3, constant hot temperatures M = 1..4.
Environment: WSL Fedora, Python 3, sympy 1.14. Runtime: about 1-3 minutes.
"""
import sys, time
import sympy as sp

FAIL = 0
def report(name, ok):
    global FAIL
    print(("PASS " if ok else "FAIL ") + name, flush=True)
    if not ok:
        FAIL += 1

def branch(M, finite_hot):
    w = sp.Symbol('w', positive=True)
    UA = sp.symbols(f'UA1:{M+1}', positive=True)
    thh = sp.symbols(f'thh1:{M+1}')
    wh = sp.symbols(f'wh1:{M+1}', positive=True) if finite_hot else [None]*M
    th = [sp.Integer(0)]
    Q = []
    for i in range(M):
        t, to, q = sp.symbols('t to q')
        eqs = [sp.Eq(q, w*(t - th[-1])),
               sp.Eq(q, UA[i]/2*((thh[i] - t) + ((to if finite_hot else thh[i]) - th[-1])))]
        unk = [t, q]
        if finite_hot:
            eqs.append(sp.Eq(q, wh[i]*(thh[i] - to)))
            unk.append(to)
        sol = sp.solve(eqs, unk, dict=True)
        assert len(sol) == 1
        th.append(sp.factor(sol[0][t]))
        Q.append(sp.factor(sol[0][q]))
    return w, th, thh, Q

def is_zero(expr):
    return sp.cancel(sp.together(expr)) == 0

def run(M, finite_hot):
    t0 = time.time()
    w, th, thh, Q = branch(M, finite_hot)
    dQ = [sp.diff(q, w) for q in Q]
    c = [(th[i+1] - th[i])/(thh[i] - th[i]) for i in range(M)]
    a_corr, a_pap = [], []
    for i in range(M):
        a_corr.append(c[i]*(th[i+1] + th[i] - sum(a_corr)))
        a_pap.append(c[i]*(th[i+1] + th[i] - (a_pap[-1] if a_pap else 0)))
    tag = f"M={M}, {'finite wh' if finite_hot else 'wh=oo'}"
    ok = all(is_zero(2*dQ[i] - a_corr[i]) for i in range(M))
    report(f"C1 corrected recursion 2 dQ_i/dw == a_i for all i ({tag})", ok)
    ok12 = all(is_zero(2*dQ[i] - a_pap[i]) for i in range(min(M, 2)))
    report(f"C2a published recursion correct for i<=2 ({tag})", ok12)
    if M >= 3:
        bad3 = not is_zero(2*dQ[2] - a_pap[2])
        report(f"C2b published recursion wrong at i=3 ({tag})", bad3)
        badtot = not is_zero(2*sum(dQ) - sum(a_pap))
        report(f"C2c published branch total sum a^P != 2 dQ_tot/dw ({tag})", badtot)
        d3 = sp.cancel(sp.together(a_pap[2] - a_corr[2] - c[2]*a_corr[0]))
        report(f"C2d a3^P - a3 == c3*a1 ({tag})", d3 == 0)
    print(f"   [{tag}: {time.time()-t0:.1f} s]", flush=True)

# C3: single exchanger, finite wh, the factor -1/2
w, th, thh, Q = branch(1, True)
p = sp.Symbol('p')
J = -p*Q[0]
report("C3 single exchanger: dJ/dw == -p th1^2/(2 thh1) (paper (16)/(A.7) says +p th1^2/thh1)",
       is_zero(sp.diff(J, w) + p*th[1]**2/(2*thh[0])))
# the paper's (A.5) theta = thh/(1 - w(1/wh + 2/UA)) versus the true solution
UA1, wh1 = sp.Symbol('UA1', positive=True), sp.Symbol('wh1', positive=True)
true_th = 2*thh[0]/(1 + w*(1/wh1 + 2/UA1))
report("C3b model solution th1 == 2 thh1/(1 + w(1/wh1 + 2/UA1))  [paper's (A.5) differs]", is_zero(th[1] - true_th))

for M in (1, 2, 3):
    run(M, True)
for M in (1, 2, 3, 4):
    run(M, False)

print("ALL PASS" if FAIL == 0 else f"{FAIL} FAIL")
sys.exit(1 if FAIL else 0)
