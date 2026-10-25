#!/usr/bin/env python3
"""Independent classification of the case-9 and case-10 networks of Banaji-Boros-Hofbauer (BBH 2024).

Written by the verifier from the source paper alone (arXiv:2304.02303v1, Section 4.1 and Theorem 9 = Theorem 3
of the published version), before reading the solver's scripts.

Setting (BBH (13)-(16)): reaction i is a_i X + b_i Y -> (a_i + c_i) X + (b_i + d_i) Y, rate kappa_i x^a_i y^b_i.
  case 9 : sources (2,0), (1,1), (0,1)   (2X, X+Y, Y; positively oriented)
  case 10: sources (2,0), (1,1), (0,0)   (2X, X+Y, 0; positively oriented)
Targets are arbitrary: c_i >= -a_i, d_i >= -b_i, (c_i, d_i) != (0, 0).

For every network in the box c_i, d_i <= N the script computes, in exact integer arithmetic:
  u = c x d; the class (trivial reaction / rank < 2 / dynamically trivial / saddle / u > 0);
  for u > 0 the trace numerators P = sum a_i c_i u_i and Q = sum b_i d_i u_i, so that at the unique positive
  equilibrium tr J = mu (P / x + Q / y), mu > 0, and (x, y) ranges over the whole open quadrant as kappa varies.
Checks (PASS/FAIL):
  C1  the trace formula, symbolically (sympy), against the Jacobian of the mass-action vector field;
  C2  "tr J = 0 for some kappa with det J > 0" (P Q < 0) agrees with BBH's literal Theorem 9 conditions;
  C3  P = Q = 0 never occurs for u > 0 (no network with tr J identically zero);
  C4  every non-Hopf network with u > 0 lies in exactly one of the verifier's families (derived by hand):
        case 9 : F1 generalised LVA (c1>=1, c2=-1, c3=0, d1=0, d2>=1, d3=-1), tr J > 0 for all kappa;
                 F2a c1>=1, c2=-1, c3>=1, d1>=1, d2=0, d3=-1, c1 <= c3 d1;
                 F2b c1>=1, c2=-1, c3>=1, d1>c1, d2=-1, -1<=d3<c3, c3 d1 + c1 d3 >= 2 c1 c3;
                 F3  c1 in {-2,-1,0}, c2=-1 (u > 0);
                 F4  c1 in {-2,-1}, c2>=0, d2=-1 (u > 0);
        case 10: d2 = -1 always, and G1 c1>=1, c2=-1 (d1/c1 + d3/c3 >= 2); G2 c1<=0, c2=-1; G3 c1<0, c2>=0;
      and conversely each family member is non-Hopf with u > 0;
  C5  outside F1, every non-Hopf network with u > 0 has P <= 0 and Q <= 0 (tr J < 0 for all kappa).
Usage: python3 classify_independent.py 8     (box c_i, d_i <= 8; log classify_independent.log)
       python3 classify_independent.py T20   (all targets with coefficients 0..20, the solver's box, 1.7e8
                                              networks; log classify_independent_T20.log)
Runtime (12 processes): about 10 s for N = 8, about 4 minutes for T20. Environment: WSL Fedora, Python 3, sympy.
"""
import itertools
import sys
from fractions import Fraction as Fr

import sympy as sp

ARG = sys.argv[1] if len(sys.argv) > 1 else '8'
TARGET_BOX = ARG.startswith('T')          # 'T20': all targets with coefficients 0..20 (the solver's box)
N = int(ARG[1:]) if TARGET_BOX else int(ARG)
T = N
SOURCES = {9: [(2, 0), (1, 1), (0, 1)], 10: [(2, 0), (1, 1), (0, 0)]}
fails = 0


def report(ok, msg):
    global fails
    print(("PASS " if ok else "FAIL ") + msg)
    if not ok:
        fails += 1


def cross(c, d):
    return (c[1] * d[2] - c[2] * d[1], c[2] * d[0] - c[0] * d[2], c[0] * d[1] - c[1] * d[0])


# ---------- C1: trace formula ----------
x, y, mu = sp.symbols('x y mu', positive=True)
cs = sp.symbols('c1:4')
ds = sp.symbols('d1:4')
for case, src in SOURCES.items():
    u = cross(cs, ds)
    # choose kappa so that (x, y) is an equilibrium: kappa_i x^a y^b = mu u_i
    kap = [mu * u[i] / (x ** src[i][0] * y ** src[i][1]) for i in range(3)]
    X, Y = sp.symbols('X Y', positive=True)
    f = sum(cs[i] * kap[i] * X ** src[i][0] * Y ** src[i][1] for i in range(3))
    g = sum(ds[i] * kap[i] * X ** src[i][0] * Y ** src[i][1] for i in range(3))
    eq_ok = sp.simplify(f.subs({X: x, Y: y})) == 0 and sp.simplify(g.subs({X: x, Y: y})) == 0
    tr = (sp.diff(f, X) + sp.diff(g, Y)).subs({X: x, Y: y})
    P = sum(src[i][0] * cs[i] * u[i] for i in range(3))
    Q = sum(src[i][1] * ds[i] * u[i] for i in range(3))
    report(eq_ok and sp.simplify(tr - mu * (P / x + Q / y)) == 0,
           f"C1 case {case}: (x,y) is an equilibrium and tr J = mu (P/x + Q/y)")
    det = sp.Matrix([[sp.diff(f, X), sp.diff(f, Y)], [sp.diff(g, X), sp.diff(g, Y)]]).det().subs({X: x, Y: y})
    orient = (src[1][0] - src[0][0]) * (src[2][1] - src[0][1]) - (src[1][1] - src[0][1]) * (src[2][0] - src[0][0])
    report(sp.simplify(det - mu ** 2 * u[0] * u[1] * u[2] * orient / (x * y)) == 0 and orient > 0,
           f"C1 case {case}: det J = mu^2 u1 u2 u3 * {orient} / (x y), orientation positive")


# ---------- BBH Theorem 9 conditions, transcribed literally ----------
def bbh_hopf(case, c, d):
    c1, c2, c3 = c
    d1, d2, d3 = d
    if case == 9:
        if c1 > 0 and c2 == -1 and c3 > 0 and d1 > 0 and d2 >= -1 and d3 >= -1:
            return Fr(1, 2) * (Fr(d3, c3) + Fr(d1, c1)) < Fr(d2, c2) < Fr(d1, c1)
        if c1 > 0 and c2 == -1 and c3 == 0 and d1 > 0 and d2 >= -1 and d3 == -1:
            return Fr(d2, c2) < Fr(d1, c1)
        return False
    if c1 > 0 and c2 == -1 and c3 > 0 and d1 > 0 and d2 == -1 and d3 >= 0:
        return Fr(1, 2) * (Fr(d3, c3) + Fr(d1, c1)) < Fr(d2, c2) < Fr(d1, c1)
    return False


def family(case, c, d):
    c1, c2, c3 = c
    d1, d2, d3 = d
    fams = []
    if case == 9:
        if c1 >= 1 and c2 == -1 and c3 == 0 and d1 == 0 and d2 >= 1 and d3 == -1:
            fams.append('F1')
        if c1 >= 1 and c2 == -1 and c3 >= 1 and d1 >= 1 and d2 == 0 and d3 == -1 and c1 <= c3 * d1:
            fams.append('F2a')
        if (c1 >= 1 and c2 == -1 and c3 >= 1 and d1 > c1 and d2 == -1 and -1 <= d3 < c3
                and c3 * d1 + c1 * d3 >= 2 * c1 * c3):
            fams.append('F2b')
        u = cross(c, d)
        upos = all(t > 0 for t in u)
        if c1 in (-2, -1, 0) and c2 == -1 and upos:
            fams.append('F3')
        if c1 in (-2, -1) and c2 >= 0 and d2 == -1 and upos:
            fams.append('F4')
    else:
        u = cross(c, d)
        upos = all(t > 0 for t in u)
        if c1 >= 1 and c2 == -1 and d2 == -1 and upos and Fr(d1, c1) + Fr(d3, c3) >= 2:
            fams.append('G1')
        if c1 <= 0 and c2 == -1 and d2 == -1 and upos:
            fams.append('G2')
        if c1 < 0 and c2 >= 0 and d2 == -1 and upos:
            fams.append('G3')
    return fams


def work(args):
    """classify all networks with a fixed c1 (one process); returns counters"""
    case, c1, ranges = args
    src = SOURCES[case]
    counts, famcount, examples = {}, {}, {}
    mism = c3fail = c4fail = c5fail = 0
    for c2, c3, d1, d2, d3 in itertools.product(*ranges[1:]):
        c, d = (c1, c2, c3), (d1, d2, d3)
        if any(c[i] == 0 and d[i] == 0 for i in range(3)):
            if TARGET_BOX:
                continue            # target = source is not a reaction (as in the solver's box)
            cls = 'trivial reaction'
        else:
            u = cross(c, d)
            if u == (0, 0, 0):
                cls = 'rank<2'
            elif all(t > 0 for t in u):
                cls = 'u>0'
            elif all(t < 0 for t in u):
                cls = 'saddle (u<0)'
            else:
                cls = 'dynamically trivial'
        hopf = False
        fams = family(case, c, d)
        if cls == 'u>0':
            u = cross(c, d)
            P = sum(src[i][0] * c[i] * u[i] for i in range(3))
            Q = sum(src[i][1] * d[i] * u[i] for i in range(3))
            if P == 0 and Q == 0:
                c3fail += 1
            hopf = P * Q < 0
            if not hopf:
                if len(fams) != 1:
                    c4fail += 1
                else:
                    famcount[fams[0]] = famcount.get(fams[0], 0) + 1
                    examples.setdefault(fams[0], (c, d))
                    if fams[0] != 'F1' and not (P <= 0 and Q <= 0):
                        c5fail += 1
                    if fams[0] == 'F1' and not (P >= 0 and Q >= 0 and P + Q > 0):
                        c5fail += 1
                cls = 'u>0, no Hopf'
            else:
                cls = 'u>0, Hopf'
                if fams:
                    c4fail += 1
        elif fams:
            c4fail += 1
        if hopf != bbh_hopf(case, c, d):
            mism += 1
        counts[cls] = counts.get(cls, 0) + 1
    return counts, famcount, examples, mism, c3fail, c4fail, c5fail


if __name__ == '__main__':
    import multiprocessing as mp
    import time
    t0 = time.time()
    pool = mp.Pool(min(12, mp.cpu_count()))
    for case, src in SOURCES.items():
        if TARGET_BOX:     # targets with both coefficients in 0..T (the solver's box)
            ranges = [range(-a, T - a + 1) for (a, b) in src] + [range(-b, T - b + 1) for (a, b) in src]
            label = f"targets <= {T}"
        else:
            ranges = [range(-a, N + 1) for (a, b) in src] + [range(-b, N + 1) for (a, b) in src]
            label = f"c_i, d_i <= {N}"
        res = pool.map(work, [(case, c1, ranges) for c1 in ranges[0]])
        counts, famcount, examples = {}, {}, {}
        mism = c3fail = c4fail = c5fail = 0
        for (cn, fc, ex, a1, a2, a3, a4) in res:
            for kk, v in cn.items():
                counts[kk] = counts.get(kk, 0) + v
            for kk, v in fc.items():
                famcount[kk] = famcount.get(kk, 0) + v
            for kk, v in ex.items():
                examples.setdefault(kk, v)
            mism += a1; c3fail += a2; c4fail += a3; c5fail += a4
        total = sum(counts.values())
        print(f"case {case}, box {label}: {total} networks; " + ", ".join(f"{k}: {v}" for k, v in sorted(counts.items())))
        print(f"  non-Hopf u>0 by family: {dict(sorted(famcount.items()))}")
        print(f"  first member of each family (c, d): {dict(sorted(examples.items()))}")
        report(mism == 0, f"C2 case {case}: P*Q<0 criterion == BBH Theorem 9 literal conditions on the whole box")
        report(c3fail == 0, f"C3 case {case}: no u>0 network with P = Q = 0")
        report(c4fail == 0, f"C4 case {case}: non-Hopf u>0 networks = disjoint union of the verifier's families")
        report(c5fail == 0, f"C5 case {case}: trace sign constant (F1: > 0; all others: < 0) for every kappa")
    print(f"runtime {time.time() - t0:.0f} s")
    print("ALL PASS" if fails == 0 else f"{fails} FAIL")
    sys.exit(1 if fails else 0)
