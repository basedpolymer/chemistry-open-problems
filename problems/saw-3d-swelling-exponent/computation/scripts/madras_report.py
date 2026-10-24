#!/usr/bin/env python3
"""madras_report.py -- summarise the lower-bound solver's requests from the outputs of
src/cond_prefix.cpp (CL3, per-prefix (I2), exact), src/madras_exact.cpp (CL1 exact sanity check) and
src/madras_mc.cpp (CL1, CL2, CL5: HEURISTIC Monte Carlo).

Environment: Python 3.8+ standard library.  Runtime < 1 s.
Usage: python3 madras_report.py [--prefix FILE] [--exact FILE] [--mc FILE ...]
"""
import math
import sys
from fractions import Fraction as F

DIRS = ["+x", "-x", "+y", "-y", "+z", "-z"]


def prefix_section(path):
    print("## CL3 per-prefix (I2): E_n[x1(w_n) | first k steps = u] >= 1 for every k-step prefix u (first step +e1)")
    print("   exact (src/cond_prefix.cpp, no symmetry reduction).  Columns: n, #prefixes, #failing, fraction of walks in")
    print("   failing prefixes, min conditional mean (exact) and its prefix")
    rows = {}
    for line in open(path):
        t = line.split()
        if not t or t[0] != "C":
            continue
        k, n = int(t[1]), int(t[2])
        rows[(k, n)] = (int(t[3]), int(t[4]), int(t[5]), int(t[6]), F(int(t[7]), int(t[8])), t[9] if len(t) > 9 else "")
    K = max(k for k, _ in rows)
    N = max(n for _, n in rows)
    for k in range(1, K + 1):
        fails = [n for n in range(k, N + 1) if rows[(k, n)][1] > 0]
        if not fails:
            print("[CL3-prefix k=%d] HOLDS for all n in [%d, %d]" % (k, k, N))
        else:
            n0 = fails[0]
            r = rows[(k, n0)]
            print("[CL3-prefix k=%d] FAILS at n = %d (exact): prefix %s has conditional mean %s = %.4f; failing n: %s" % (
                k, n0, " ".join(DIRS[int(c)] for c in r[5]), r[4], float(r[4]), fails))
        for n in range(max(k, 2), N + 1):
            r = rows[(k, n)]
            print("     k=%d n=%2d prefixes=%d failing=%d walk-fraction=%.4f min=%.4f at %s" % (
                k, n, r[0], r[1], r[2] / r[3] if r[3] else 0, float(r[4]), " ".join(DIRS[int(c)] for c in r[5])))


def exact_section(path):
    print("\n## CL1 exact sanity check (src/madras_exact.cpp, all SAWs, near set ||w(n)||_inf <= m)")
    print("   n  m   near          kappa=E[K|near]   P(K=2 | near, span>=2)   self-test")
    st = "?"
    for line in open(path):
        if "self-test" in line:
            st = "PASS" if "PASS" in line else "FAIL"
        t = line.split()
        if not t or t[0] != "X":
            continue
        n, m, near, sK = int(t[1]), int(t[2]), int(t[3]), int(t[4])
        hist2 = list(map(int, t[14:22]))
        near2 = int(t[13])
        p2 = hist2[1] / near2 if near2 else float("nan")
        print("   %2d %d %12d   %.5f            %.5f                  %s" % (n, m, near, sK / near if near else 0, p2, st))


def mc_section(paths):
    print("\n## CL1 / CL2 / CL5 Monte Carlo (HEURISTIC; src/madras_mc.cpp; errors = jackknife over 64 blocks)")
    for p in paths:
        lines = open(p).read().splitlines()
        head = [l for l in lines if l.startswith("# madras_mc")]
        print("  " + (head[0] if head else p))
        ms, ks = [], []
        print("    m   P(near)   n_K   kappa +- err     cand    ren +- err      saw +- err     frac(span<=m)")
        for l in lines:
            t = l.split()
            if t and t[0] == "M" and len(t) > 3:
                m, nK, pn, kap, ek = int(t[1]), float(t[2]), float(t[3]), float(t[4]), float(t[5])
                print("    %3d  %.5f  %6.0f  %.3f +- %.3f  %7.2f  %.3f +- %.3f  %.3f +- %.3f  %.4f" % (
                    m, pn, nK, kap, ek, float(t[6]), float(t[8]), float(t[9]), float(t[10]), float(t[11]), float(t[12])))
                if nK >= 50 and float(t[12]) < 0.01:
                    ms.append(m)
                    ks.append((kap, ek))
        if len(ms) >= 3:
            xs = [math.log(m) for m in ms]
            ys = [math.log(k) for k, _ in ks]
            mx, my = sum(xs) / len(xs), sum(ys) / len(ys)
            sl = sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / sum((x - mx) ** 2 for x in xs)
            print("    NUMERICS log kappa vs log m (m in %s): slope %.3f, i.e. kappa ~ m^{%.3f} (alpha ~ %.3f)" %
                  (ms, sl, sl, 1 - sl))
        js, rs = [], []
        for l in lines:
            t = l.split()
            if t and t[0] == "JE":
                j, r, er, rh, erh = int(t[1]), float(t[2]), float(t[3]), float(t[4]), float(t[5])
                print("    CL2 j=%3d  r_N(j) = %.5f +- %.5f   r'_N(j) = %.5f +- %.5f" % (j, r, er, rh, erh))
                if r > 0:
                    js.append(j)
                    rs.append(r + rh)
        if len(js) >= 3:
            sel = [(j, r) for j, r in zip(js, rs) if j >= 4]
            if len(sel) >= 3:
                xs = [math.log(j) for j, _ in sel]
                ys = [math.log(r) for _, r in sel]
                mx, my = sum(xs) / len(xs), sum(ys) / len(ys)
                sl = sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / sum((x - mx) ** 2 for x in xs)
                print("    NUMERICS r_N(j) + r'_N(j) ~ j^{-beta}, beta ~ %.3f (fit over j >= 4)" % (-sl))
        for l in lines:
            if l.startswith("Q ") or l.startswith("P eta"):
                print("    CL5 " + l)


def main():
    a = sys.argv[1:]
    i = 0
    mc = []
    while i < len(a):
        if a[i] == "--prefix":
            prefix_section(a[i + 1])
            i += 2
        elif a[i] == "--exact":
            exact_section(a[i + 1])
            i += 2
        elif a[i] == "--mc":
            i += 1
            while i < len(a) and not a[i].startswith("--"):
                mc.append(a[i])
                i += 1
        else:
            i += 1
    if mc:
        mc_section(mc)


if __name__ == "__main__":
    main()
