#!/usr/bin/env python3
"""
reproduce_appendixB.py -- the paper's own validation case (Appendix B, Table B.6: M_1 = 4, M_2 = 1, prices).

Environment : WSL Fedora, Python 3 + sympy (1.12+).  Runtime: ~10-60 s.
Usage       : python -I reproduce_appendixB.py

The cost is J = -sum p_ij Q_ij (eq. 14) with the data of Table B.6, AMTD model, split u = w_1/w_0.
Everything is an exact rational function of u.  Root locations are certified by exact real-root isolation
(sympy Poly.intervals, rational endpoints) of the numerators, after checking that the denominators have no
root in (0,1).

Claims (PASS/FAIL, exit 1 on failure):
  B1  identity: (T*_Je,1 - T*_Je,2) * w0 + 2 dJ/du == 0 exactly (corrected economic Jaschke temperature
      = -2 x marginal cost; so its zero is exactly the stationary point);
  B2  dJ/du has exactly one root u_opt in (0,1); it lies in an isolating interval contained in
      [0.71411, 0.71413] (paper, Table B.7, "Optimized": 0.7141);
  B3  the paper's controlled variable (B.1) (recursion with a_{i-1}, typos of the printed (B.1) corrected)
      has a root u_P certified in [0.71670, 0.71671], disjoint from u_opt; its only other root in (0,1)
      (u ~ 0.0028) lies where the AMTD efficiency of exchanger 11 exceeds 1 (unphysical region);
      hence controlling the paper's T_J differences to zero does NOT give the optimum of the AMTD model.
Note: the paper reports u = 0.7152 for "Equal Jaschke Temp."; we could not reproduce this number from (B.1)
(neither with nor without the printing typos); the discrepancy with the AMTD optimum, which the paper
attributes to "flatness of the optimum and numerical round-off errors", is in fact (at least in part) caused
by the incorrect recursion.
"""
import sys
import sympy as sp

FAILED = []


def report(name, ok):
    print(f"{'PASS' if ok else 'FAIL'}  {name}")
    if not ok:
        FAILED.append(name)


R = sp.Rational
u = sp.Symbol("u")
w0, T0 = 100, 130
L1 = dict(wh=[50, 30, 15, 25], Th=[190, 203, 220, 235], UA=[5, 7, 10, 12], p=[1, R(12, 10), R(13, 10), R(15, 10)])
L2 = dict(wh=[70], Th=[225], UA=[11], p=[R(15, 10)])


def branch(w, L):
    th, Q = [sp.Integer(0)], []
    for wh, Th, UA in zip(L["wh"], L["Th"], L["UA"]):
        thh, tp = Th - T0, th[-1]
        tt, tho, q = sp.symbols("tt tho q")
        sol = sp.solve([q - w * (tt - tp), q - wh * (thh - tho), q - sp.Integer(UA) / 2 * ((thh - tt) + (tho - tp))],
                       [tt, tho, q], dict=True)[0]
        th.append(sp.cancel(sol[tt]))
        Q.append(sp.cancel(sol[q]))
    return th, Q


def TJe(th, L, mode):
    a = []
    for i in range(1, len(th)):
        c = (th[i] - th[i - 1]) / (L["Th"][i - 1] - T0 - th[i - 1])
        sub = sum(a) if mode == "corrected" else (a[-1] if a else 0)
        a.append(c * (th[i] + th[i - 1] - sub))
    return sum(p * x for p, x in zip(L["p"], a))


def unique_root_in_01(expr, label):
    num, den = sp.fraction(sp.cancel(sp.together(expr)))
    P, Dn = sp.Poly(num, u), sp.Poly(den, u)
    den_roots = [iv for iv in Dn.intervals() if iv[0][1] > 0 and iv[0][0] < 1]
    ivs = [iv for iv in P.intervals(eps=R(1, 10**9)) if iv[0][1] > 0 and iv[0][0] < 1]
    print(f"      {label}: numerator degree {P.degree()}, roots in (0,1): "
          + ", ".join(f"[{float(a):.10f}, {float(b):.10f}]" for (a, b), m in ivs)
          + f"; denominator roots near (0,1): {len(den_roots)}")
    return ivs, den_roots


def main():
    th1, Q1 = branch(w0 * u, L1)
    th2, Q2 = branch(w0 * (1 - u), L2)
    J = -(sum(p * q for p, q in zip(L1["p"], Q1)) + sum(p * q for p, q in zip(L2["p"], Q2)))
    dJ = sp.diff(J, u)
    c_cor = TJe(th1, L1, "corrected") - TJe(th2, L2, "corrected")
    c_pap = TJe(th1, L1, "paper") - TJe(th2, L2, "paper")
    report("B1 corrected CV * w0 + 2 dJ/du == 0 identically", sp.cancel(sp.together(c_cor * w0 + 2 * dJ)) == 0)

    ivo, dro = unique_root_in_01(dJ, "dJ/du")
    okB2 = len(ivo) == 1 and not dro and ivo[0][0][0] >= R(71411, 100000) and ivo[0][0][1] <= R(71413, 100000)
    report("B2 unique optimum u_opt in (0,1), certified in [0.71411, 0.71413]", okB2)
    ivp, drp = unique_root_in_01(c_pap, "paper CV (B.1)")
    main_root = [iv for iv in ivp if iv[0][0] >= R(71670, 100000) and iv[0][1] <= R(71671, 100000)]
    others = [iv for iv in ivp if iv not in main_root]
    # any other root must lie where the AMTD model is unphysical: efficiency c_11 = 2/(1 + w1 k_11) > 1,
    # i.e. w1 k_11 < 1 with k_11 = 2/UA_11 + 1/wh_11 (c_11 is decreasing in u, so check the right end).
    k11 = R(2, L1["UA"][0]) + R(1, L1["wh"][0])
    others_unphysical = all(w0 * iv[0][1] * k11 < 1 for iv in others)
    okB3 = (len(main_root) == 1 and not drp and okB2 and main_root[0][0][0] > ivo[0][0][1] and others_unphysical)
    report("B3 paper CV root u_P certified in [0.71670, 0.71671], disjoint from u_opt; its only other root in (0,1) "
           "lies where c_11 > 1 (unphysical AMTD region)", okB3)
    ivp = main_root

    if okB2 and okB3:
        uo = (ivo[0][0][0] + ivo[0][0][1]) / 2
        up = (ivp[0][0][0] + ivp[0][0][1]) / 2
        Tend = T0 + (sum(Q1) + sum(Q2)) / w0
        print(f"      u_opt ~ {float(uo):.8f}: J = {float(J.subs(u, uo)):.6f} $/s, Tend = {float(Tend.subs(u, uo)):.6f} C")
        print(f"      u_P   ~ {float(up):.8f}: J = {float(J.subs(u, up)):.6f} $/s, Tend = {float(Tend.subs(u, up)):.6f} C")
        print("      paper Table B.7: optimized u = 0.7141, J = -3.9388 (x1000), Tend = 158.7738; "
              "equal Jaschke temp. u = 0.7152, Tend = 158.7741")
    if FAILED:
        print(f"FAIL: {len(FAILED)} claim(s) failed")
        sys.exit(1)
    print("ALL PASS")


if __name__ == "__main__":
    main()
