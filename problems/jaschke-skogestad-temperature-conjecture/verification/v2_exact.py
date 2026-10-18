#!/usr/bin/env python3
"""Independent verifier check V2 (written from scratch; shares no code with the solver's scripts).

Exact forward-mode differentiation with dual numbers over Fractions. Each exchanger's three model equations
(9), (10), (12) of Jaschke & Skogestad (2014) are solved as a 3x3 linear system in (th_i, tho_i, Q_i) by Gaussian
elimination on dual numbers, i.e. the affine 'Lemma 1' of the solver is NOT used.

Claims (PASS/FAIL):
  R1  random exact tests, finite hot capacities, M = 1..8 (40 draws each): 2 dQ_i/dw == a_i^corr exactly for all i.
  R2  the published recursion (18) is exact for i <= 2 and its branch total differs for every draw with M >= 3.
  X1  solver's counterexample (certificates/counterexample.json), re-derived from the raw parameters only:
      physical operating point, paper CV == 0 exactly, dQtot/dw1 == -4972391424/1518467647, Qtot(99/2) > Qtot(50).
  X2  constant-hot-temperature (wh = oo) counterexample built here from scratch for M1 = 3 (the case the paper
      reports as proved): paper CV == 0 and gradient != 0.
  X3  M1 = 4 and M1 = 5 counterexamples built here from scratch (all parameters positive, physical point).
Environment: WSL Fedora, Python 3 standard library only. Runtime: a few seconds.
"""
import sys, json, random
from fractions import Fraction as F

FAIL = 0
def report(name, ok):
    global FAIL
    print(("PASS " if ok else "FAIL ") + name, flush=True)
    if not ok:
        FAIL += 1

class D:
    """dual number v + d*eps over Fractions"""
    __slots__ = ("v", "d")
    def __init__(s, v, d=0):
        s.v, s.d = F(v), F(d)
    @staticmethod
    def c(x):
        return x if isinstance(x, D) else D(x)
    def __add__(s, o): o = D.c(o); return D(s.v + o.v, s.d + o.d)
    __radd__ = __add__
    def __sub__(s, o): o = D.c(o); return D(s.v - o.v, s.d - o.d)
    def __rsub__(s, o): return D.c(o) - s
    def __mul__(s, o): o = D.c(o); return D(s.v*o.v, s.v*o.d + s.d*o.v)
    __rmul__ = __mul__
    def __truediv__(s, o):
        o = D.c(o)
        return D(s.v/o.v, (s.d*o.v - s.v*o.d)/(o.v*o.v))
    def __rtruediv__(s, o): return D.c(o)/s
    def __neg__(s): return D(-s.v, -s.d)
    def iszero(s): return s.v == 0

def gauss(A, b):
    n = len(b)
    A = [row[:] + [b[k]] for k, row in enumerate(A)]
    for col in range(n):
        piv = next(r for r in range(col, n) if not A[r][col].iszero())
        A[col], A[piv] = A[piv], A[col]
        for r in range(n):
            if r != col and not A[r][col].iszero():
                f = A[r][col]/A[col][col]
                A[r] = [A[r][k] - f*A[col][k] for k in range(n+1)]
    return [A[k][n]/A[k][k] for k in range(n)]

def solve_branch(w, UA, wh, thh):
    """w: D (derivative seed 1 for the branch flow). wh[i] None means constant hot temperature.
    returns lists th (len M+1), tho (len M), Q (len M) of dual numbers"""
    th, tho, Q = [D(0)], [], []
    for i in range(len(UA)):
        U, h = F(UA[i]), F(thh[i])
        if wh[i] is None:
            # unknowns (t, q):  q - w t = -w th_prev ;  q + U/2 t = U/2 (2h - th_prev)
            A = [[-w, D(1)], [D(U/2), D(1)]]
            b = [-w*th[-1], (U/2)*(2*h - th[-1])]
            t, q = gauss(A, b)
            to = D(h)
        else:
            W = F(wh[i])
            # unknowns (t, to, q): q - w t = -w th_prev ; q + W to = W h ; q + U/2 t - U/2 to = U/2 (h - th_prev)
            A = [[-w, D(0), D(1)], [D(0), D(W), D(1)], [D(U/2), D(-U/2), D(1)]]
            b = [-w*th[-1], D(W*h), (U/2)*(h - th[-1])]
            t, to, q = gauss(A, b)
        th.append(t); tho.append(to); Q.append(q)
    return th, tho, Q

def recursions(th, thh):
    M = len(thh)
    c = [(th[i+1] - th[i])/(F(thh[i]) - th[i]) for i in range(M)]
    ac, ap = [], []
    for i in range(M):
        ac.append(c[i]*(th[i+1] + th[i] - sum(ac)))
        ap.append(c[i]*(th[i+1] + th[i] - (ap[-1] if ap else 0)))
    return c, ac, ap

def rnd(lo, hi, den=7):
    return F(random.randint(lo*den, hi*den), den)

# ---------------- R1, R2 ----------------
random.seed(20261018)
okR1, okR2a, okR2b = True, True, True
for M in range(1, 9):
    for _ in range(40):
        UA = [rnd(5, 80) for _ in range(M)]
        wh = [rnd(5, 120) if random.random() < 0.8 else None for _ in range(M)]
        thh = [rnd(-20, 150) for _ in range(M)]
        w = D(rnd(5, 120), 1)
        th, tho, Q = solve_branch(w, UA, wh, thh)
        thv = [x.v for x in th]
        if any(F(thh[i]) == thv[i] for i in range(M)):
            continue
        c, ac, ap = recursions(thv, thh)
        okR1 &= all(2*Q[i].d == ac[i] for i in range(M))
        okR2a &= all(2*Q[i].d == ap[i] for i in range(min(M, 2)))
        if M >= 3:
            okR2b &= (sum(ap) != 2*sum(q.d for q in Q))
report("R1 2 dQ_i/dw == corrected a_i exactly, M=1..8, 40 random draws each (finite and infinite wh mixed)", okR1)
report("R2a published recursion exact for i<=2 in every draw", okR2a)
report("R2b published branch total != 2 dQ_tot/dw in every draw with M>=3", okR2b)

# ---------------- X1: solver's counterexample ----------------
def network(w1, data):
    """two branches; returns (Qtot, dQtot/dw1 with w2 = w0 - w1), paper CV, corrected CV, branch data"""
    w0 = F(data["w0"])
    b1, b2 = data["branch1"], data["branch2"]
    th1, tho1, Q1 = solve_branch(D(w1, 1), b1["UA"], b1["wh"], b1["theta_h"])
    th2, tho2, Q2 = solve_branch(D(w0 - F(w1), -1), [b2["UA"]], [b2["wh"]], [b2["theta_h"]])
    Qtot = sum(Q1, D(0)) + sum(Q2, D(0))
    _, ac1, ap1 = recursions([x.v for x in th1], b1["theta_h"])
    _, ac2, ap2 = recursions([x.v for x in th2], [b2["theta_h"]])
    return Qtot, sum(ap1) - sum(ap2), sum(ac1) - sum(ac2), (th1, tho1, Q1, th2, tho2, Q2)

def physical(thh, th, tho):
    ok = True
    for i in range(len(thh)):
        h = F(thh[i])
        ok &= th[i+1].v > th[i].v            # cold stream heated
        ok &= th[i+1].v < h                  # cold outlet below hot inlet
        ok &= tho[i].v > th[i].v             # hot outlet above cold inlet (counter-current)
        ok &= tho[i].v <= h                  # hot stream cooled (or constant)
        c = (th[i+1].v - th[i].v)/(h - th[i].v)
        ok &= 0 < c < 1
    return ok

path = sys.argv[1] if len(sys.argv) > 1 else "../certificates/counterexample.json"
data = json.load(open(path))
data_raw = {"w0": data["w0"],
            "branch1": {k: data["branch1"][k] for k in ("UA", "wh", "theta_h")},
            "branch2": {k: data["branch2"][k] for k in ("UA", "wh", "theta_h")}}
Qt, cvP, cvC, (th1, tho1, Q1, th2, tho2, Q2) = network(F(data["w1"]), data_raw)
report("X1a counterexample operating point is physical (both branches)",
       physical(data_raw["branch1"]["theta_h"], th1, tho1) and physical([data_raw["branch2"]["theta_h"]], th2, tho2)
       and F(data["w1"]) > 0 and F(data["w0"]) - F(data["w1"]) > 0)
report("X1b paper's controlled variable T_J1^P - T_J2^P == 0 exactly", cvP == 0)
report("X1c dQtot/dw1 == -4972391424/1518467647 (recomputed by dual numbers)", Qt.d == F(-4972391424, 1518467647))
report("X1d dQtot/dw1 == -(1/2)*... consistency: dQtot/dw1 == (T_J1^corr - T_J2^corr)/2", Qt.d == cvC/2)
Qb, _, _, _ = network(F(99, 2), data_raw)
report("X1e Qtot(w1=99/2) > Qtot(w1=50) exactly", Qb.v > Qt.v)
print(f"   gradient = {float(Qt.d):.6f}, Qtot gain at 99/2 = {float(Qb.v - Qt.v):.6f} kW")
report("X1f stored temperatures in the json agree with the recomputation",
       [str(x.v) for x in th1] == data["theta_branch1"] and str(Qt.v) == data["Qtot_at_w1"])

# ---------------- X2, X3: counterexamples constructed here ----------------
def construct(UA, wh, thh, w1, w0, UA2, wh2):
    """branch 2 = one exchanger; choose its hot inlet so that the paper's CV is exactly zero."""
    th1, tho1, Q1 = solve_branch(D(w1, 1), UA, wh, thh)
    _, ac1, ap1 = recursions([x.v for x in th1], thh)
    TP1 = sum(ap1)
    # single exchanger: c = theta/thh is independent of thh; T_J2 = c^2 thh
    th2, _, _ = solve_branch(D(w0 - w1), [UA2], [wh2], [1])
    c2 = th2[1].v
    thh2 = TP1 / c2**2
    data = {"w0": w0, "branch1": {"UA": UA, "wh": wh, "theta_h": thh},
            "branch2": {"UA": UA2, "wh": wh2, "theta_h": thh2}}
    return data

for (label, UA, wh, thh) in [
        ("X2 M1=3, constant hot temperatures (wh=oo)", [30, 25, 20], [None]*3, [40, 55, 70]),
        ("X3a M1=4, finite wh", [40, 30, 25, 20], [60, 45, 80, 50], [45, 50, 62, 75]),
        ("X3b M1=5, finite wh", [40, 30, 25, 20, 30], [60, 45, 80, 50, 70], [45, 50, 62, 75, 90])]:
    w0, w1 = F(100), F(45)
    d = construct([F(x) for x in UA], [None if x is None else F(x) for x in wh], [F(x) for x in thh], w1, w0, F(120), F(60))
    Qt, cvP, cvC, (th1, tho1, Q1, th2, tho2, Q2) = network(w1, d)
    phys = physical(d["branch1"]["theta_h"], th1, tho1) and physical([d["branch2"]["theta_h"]], th2, tho2)
    report(f"{label}: physical, paper CV == 0, exact gradient != 0", phys and cvP == 0 and Qt.d != 0)
    print(f"   theta_h(branch 2) = {float(d['branch2']['theta_h']):.4f}, dQtot/dw1 = {Qt.d} ~ {float(Qt.d):.5f}")

print("ALL PASS" if FAIL == 0 else f"{FAIL} FAIL")
sys.exit(1 if FAIL else 0)
