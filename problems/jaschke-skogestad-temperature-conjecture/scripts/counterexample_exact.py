#!/usr/bin/env python3
"""
counterexample_exact.py -- exact rational counterexample to the Jaschke-Skogestad conjecture (M_1 = 3, M_2 = 1).

Environment : any Python 3.8+ (standard library only: fractions); mpmath optional (only for the
              illustrative floating-point part at the end, clearly marked NUMERICAL).
Runtime     : < 1 s.
Usage       : python -I counterexample_exact.py [path/to/certificates/counterexample.json]

Everything in the claims E1-E5 is exact rational arithmetic.  Derivatives are computed by forward-mode
automatic differentiation with exact Fraction "dual numbers" applied to the raw model equations
(each exchanger's 3x3 linear system is solved by Gaussian elimination), i.e. WITHOUT using any of the
closed-form marginal-cost formulas that are being tested.

Network (paper's Fig. A.13): cold feed w0 split into branch 1 (3 exchangers in series) and branch 2
(1 exchanger); AMTD heat transfer; objective: maximise total heat recovery Q = sum Q_ij (all prices 1,
Theorem 2), i.e. minimise J = -Q.  Temperatures are shifted: theta = T - T0 (T0 = 20 C below).

Claims (PASS/FAIL, exit code 1 on failure):
  E1  the operating point is physically meaningful (every exchanger heats the cold stream, cold outlet
      below hot inlet, hot outlet above cold inlet, i.e. positive driving force at both terminals,
      efficiencies in (0,1), UA and heat capacities positive);
  E2  the paper's controlled variable c = T_J,1 - T_J,2 (eq. 17-18, recursion with a_{i-1}) is EXACTLY 0;
  E3  the exact gradient dQ/dw1 (w2 = w0 - w1) is NOT zero; it equals (T*_J,1 - T*_J,2)/2 with the corrected
      cumulative recursion, and equals -c_3 a_1 / 2;
  E4  a nearby rational split w1' gives STRICTLY larger total heat recovery (exact), so the point where the
      paper's controlled variable is zero is not optimal;
  E5  single-branch form of Theorem 2: for branch 1 alone, dJ_1/dw1 = -(1/2) T*_J,1 exactly, whereas
      T_J,1 (paper) is not proportional to it (dJ_1/dw1 + T_J,1/2 != 0).
"""
import json
import sys
from fractions import Fraction as F

FAILED = []


def report(name, ok):
    print(f"{'PASS' if ok else 'FAIL'}  {name}")
    if not ok:
        FAILED.append(name)


class D:
    """Dual number v + d*eps over the rationals (exact forward-mode derivative)."""
    __slots__ = ("v", "d")

    def __init__(self, v, d=0):
        self.v, self.d = F(v), F(d)

    @staticmethod
    def c(x):
        return x if isinstance(x, D) else D(x, 0)

    def __add__(s, o):
        o = D.c(o)
        return D(s.v + o.v, s.d + o.d)

    __radd__ = __add__

    def __sub__(s, o):
        o = D.c(o)
        return D(s.v - o.v, s.d - o.d)

    def __rsub__(s, o):
        return D.c(o) - s

    def __neg__(s):
        return D(-s.v, -s.d)

    def __mul__(s, o):
        o = D.c(o)
        return D(s.v * o.v, s.v * o.d + s.d * o.v)

    __rmul__ = __mul__

    def __truediv__(s, o):
        o = D.c(o)
        return D(s.v / o.v, (s.d * o.v - s.v * o.d) / (o.v * o.v))

    def __rtruediv__(s, o):
        return D.c(o) / s


def gauss3(A, b):
    """Solve a 3x3 system with entries in D (exact), partial pivoting on nonzero values."""
    A = [row[:] + [bb] for row, bb in zip(A, b)]
    n = 3
    for k in range(n):
        p = next(r for r in range(k, n) if A[r][k].v != 0)
        A[k], A[p] = A[p], A[k]
        for r in range(n):
            if r != k:
                f = A[r][k] / A[k][k]
                A[r] = [A[r][j] - f * A[k][j] for j in range(n + 1)]
    return [A[i][n] / A[i][i] for i in range(n)]


def branch(w, UA, wh, thh):
    """Raw AMTD model of one branch. w is a D (value, dw/dw1). Returns lists theta, tho, Q (as D)."""
    theta, tho_l, Q = [D(0)], [], []
    for U, W, h in zip(UA, wh, thh):
        tp = theta[-1]
        # unknowns (t, tho, q):  q = w (t - tp);  q = W (h - tho);  q = U/2 [(h - t) + (tho - tp)]
        A = [[-w, D(0), D(1)], [D(0), D(W), D(1)], [D(F(U, 2)), D(-F(U, 2)), D(1)]]
        b = [-w * tp, D(W * h), D(F(U, 2) * h) - F(U, 2) * tp]
        t, tho, q = gauss3(A, b)
        theta.append(t)
        tho_l.append(tho)
        Q.append(q)
    return theta, tho_l, Q


def jaschke(theta, thh, mode):
    """Jaschke temperature from temperatures only (Fractions). mode 'paper' (a_{i-1}) or 'corrected' (sum)."""
    a = []
    for i in range(1, len(theta)):
        c = (theta[i] - theta[i - 1]) / (thh[i - 1] - theta[i - 1])
        sub = sum(a, F(0)) if mode == "corrected" else (a[-1] if a else F(0))
        a.append(c * (theta[i] + theta[i - 1] - sub))
    return sum(a, F(0)), a


# ---------------- data (all exact rationals; units kW/K and K) ----------------
T0 = F(20)                       # cold feed temperature [C] (only for display)
w0 = F(100)                      # cold feed heat capacity rate [kW/K]
w1 = F(50)                       # branch-1 heat capacity rate at the operating point
UA1 = [F(56), F(34), F(28)]      # branch 1, exchangers 1..3
wh1 = [F(88), F(51), F(78)]
thh1 = [F(48), F(52), F(60)]     # hot inlet temperatures, shifted (T_h = 68, 72, 80 C)
c12 = F(4, 5)                    # target efficiency of the single exchanger on branch 2
wh12 = F(50)


def build():
    w2 = w0 - w1
    th1v, _, _ = branch(D(w1, 1), UA1, wh1, thh1)
    TJ1_paper, a_pap = jaschke([x.v for x in th1v], thh1, "paper")
    # Branch 2 is tuned so that the paper's controlled variable vanishes exactly:
    #   T_J,2 = theta12^2/thh12 = c12^2 * thh12 := TJ1_paper  ->  thh12 = TJ1_paper / c12^2,
    #   and UA12 from c12 = 2/(1 + w2 (2/UA12 + 1/wh12)).
    thh12 = TJ1_paper / c12 ** 2
    inv = ((2 / c12 - 1) / w2 - 1 / wh12) / 2           # = 1/UA12
    UA12 = 1 / inv
    return w2, thh12, UA12


def evaluate(w1val, w2, thh12, UA12, dw1=1):
    th1, tho1, Q1 = branch(D(w1val, dw1), UA1, wh1, thh1)
    th2, tho2, Q2 = branch(D(w0 - w1val, -dw1), [UA12], [wh12], [thh12])
    return th1, tho1, Q1, th2, tho2, Q2


def main():
    w2, thh12, UA12 = build()
    th1, tho1, Q1, th2, tho2, Q2 = evaluate(w1, w2, thh12, UA12)

    # E1: physical meaning
    ok = UA12 > 0 and w2 > 0 and thh12 > 0
    ratios = []
    for thv, thov, hh in ((th1, tho1, thh1), (th2, tho2, [thh12])):
        for i in range(len(hh)):
            tin, tout, ho = thv[i].v, thv[i + 1].v, thov[i].v
            c = (tout - tin) / (hh[i] - tin)
            ok &= (hh[i] > tout > tin) and (ho > tin) and (hh[i] > ho) and (0 < c < 1)
            r = (hh[i] - tout) / (ho - tin)
            ratios.append(max(r, 1 / r))
    report("E1 operating point physically meaningful (positive terminal driving forces, 0<c<1)", ok)
    print(f"      max ratio of terminal temperature differences = {float(max(ratios)):.4f} (AMTD = LMTD when 1)")

    TJ1p, a1p = jaschke([x.v for x in th1], thh1, "paper")
    TJ1c, a1c = jaschke([x.v for x in th1], thh1, "corrected")
    TJ2p, _ = jaschke([x.v for x in th2], [thh12], "paper")
    TJ2c, _ = jaschke([x.v for x in th2], [thh12], "corrected")
    report("E2 paper's controlled variable T_J,1 - T_J,2 == 0 exactly", TJ1p - TJ2p == 0)

    Qtot = sum((q for q in Q1 + Q2), D(0))
    grad = Qtot.d                               # exact dQ/dw1 with w2 = w0 - w1
    c3 = (th1[3].v - th1[2].v) / (thh1[2] - th1[2].v)
    report("E3 exact gradient dQ/dw1 != 0", grad != 0)
    report("E3 dQ/dw1 == (T*_J,1 - T*_J,2)/2 (corrected recursion)", grad == (TJ1c - TJ2c) / 2)
    report("E3 dQ/dw1 == -c_3 a_1 / 2", grad == -c3 * a1c[0] / 2)

    # E4: exhibit a better rational split
    delta = F(1, 2) if grad < 0 else F(-1, 2)
    w1b = w1 - abs(delta) if grad < 0 else w1 + abs(delta)
    _, _, Q1b, _, _, Q2b = evaluate(w1b, w2, thh12, UA12)
    Qb = sum((q.v for q in Q1b + Q2b), F(0))
    report(f"E4 split w1' = {w1b} gives strictly more heat than w1 = {w1}", Qb > Qtot.v)

    # E5: single-branch statement of Theorem 2 (J_1 = -sum Q_i1)
    dJ1 = -sum((q.d for q in Q1), F(0))          # with dw1 = 1
    report("E5 dJ_1/dw1 == -(1/2) T*_J,1 (corrected) exactly", dJ1 == -TJ1c / 2)
    report("E5 dJ_1/dw1 + (1/2) T_J,1 (paper) != 0", dJ1 + TJ1p / 2 != 0)

    cert = {
        "description": "Exact counterexample, M1=3, M2=1, AMTD, maximise total heat (prices 1). Shifted temperatures theta=T-T0.",
        "T0_C": str(T0), "w0": str(w0), "w1": str(w1), "w2": str(w2),
        "branch1": {"UA": [str(x) for x in UA1], "wh": [str(x) for x in wh1], "theta_h": [str(x) for x in thh1]},
        "branch2": {"UA": str(UA12), "wh": str(wh12), "theta_h": str(thh12)},
        "theta_branch1": [str(x.v) for x in th1], "theta_hot_out_branch1": [str(x.v) for x in tho1],
        "theta_branch2": [str(x.v) for x in th2], "theta_hot_out_branch2": [str(x.v) for x in tho2],
        "TJ1_paper": str(TJ1p), "TJ2": str(TJ2p), "TJ1_corrected": str(TJ1c),
        "a_paper_branch1": [str(x) for x in a1p], "a_corrected_branch1": [str(x) for x in a1c],
        "dQtot_dw1_exact": str(grad), "Qtot_at_w1": str(Qtot.v), "w1_better": str(w1b), "Qtot_at_w1_better": str(Qb),
    }
    print("\n--- exact data ---")
    for k in ("w1", "w2"):
        print(f"{k} = {cert[k]}")
    print("branch 2: UA12 =", UA12, " wh12 =", wh12, " theta_h12 =", thh12, f"(~{float(thh12):.6f})")
    print("theta (branch 1) ~", [f"{float(x.v):.6f}" for x in th1])
    print("T_J,1 paper     =", f"{float(TJ1p):.10f}", "  T_J,2 =", f"{float(TJ2p):.10f}")
    print("T*_J,1 corrected=", f"{float(TJ1c):.10f}")
    print("dQ/dw1 exact    =", grad, f"(~{float(grad):.6e} K)")
    print("Q(w1) =", f"{float(Qtot.v):.10f}", " Q(w1') =", f"{float(Qb):.10f}", " gain =", f"{float(Qb - Qtot.v):.3e} kW")

    # NUMERICAL illustration only (not part of the certificate): locate the true optimum.
    try:
        import mpmath as mp
        mp.mp.dps = 40

        def gradf(x):
            x = F(mp.nstr(x, 45))
            _, _, a, _, _, b = evaluate(x, w2, thh12, UA12)
            return mp.mpf(sum((q.d for q in a + b), F(0)).numerator) / sum((q.d for q in a + b), F(0)).denominator

        wopt = mp.findroot(gradf, (mp.mpf(30), mp.mpf(50)), solver="anderson")
        _, _, a, _, _, b = evaluate(F(mp.nstr(wopt, 45)), w2, thh12, UA12)
        Qopt = sum((q.v for q in a + b), F(0))
        print(f"NUMERICAL: optimum w1* ~ {mp.nstr(wopt, 12)} (u* ~ {mp.nstr(wopt/100, 8)}), "
              f"Q* ~ {float(Qopt):.8f}; loss at paper point ~ {float(Qopt - Qtot.v):.4e} kW "
              f"({float((Qopt - Qtot.v)/Qopt)*100:.4f} %)")
        cert["NUMERICAL_w1_opt"] = mp.nstr(wopt, 20)
    except ImportError:
        pass

    if len(sys.argv) > 1:
        with open(sys.argv[1], "w") as fh:
            json.dump(cert, fh, indent=1)
        print("certificate written to", sys.argv[1])
    if FAILED:
        print(f"FAIL: {len(FAILED)} claim(s) failed")
        sys.exit(1)
    print("ALL PASS")


if __name__ == "__main__":
    main()
