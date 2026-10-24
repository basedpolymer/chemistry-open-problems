#!/usr/bin/env python3
"""check_aux.py -- PASS/FAIL validation of the auxiliary programs against brute force / exact values (small n).

Environment: Linux, g++ >= 9 (C++17), Python 3.8+ standard library.  Runtime about 2-4 min (at most 3 threads).
Usage: python3 check_aux.py

  [a] src/end_steps.cpp N=12: sum s_1.s_n equals the Gram sums of src/saw_gram.cpp, and c_n = OEIS A001412.
  [b] src/halfspace_enum.cpp N=9 (bridges by height, renewals, arches, half-space totals) = Python brute force.
  [c] src/cond_prefix.cpp N=9, k<=4 (per-prefix conditional means of x1(w_n)) = Python brute force.
  [d] src/madras_exact.cpp N=8: multiplicity routine (src/madras.hpp) = direct count of preimages (built-in self-test).
  [e] src/pivot.cpp: E_n|w_n|^2 for n = 10, 20 within 4 sigma of the exact values (2e6 attempts each).
  [f] src/madras_mc.cpp N=8: kappa_8(m) for m = 2, 4 within 4 sigma of the exhaustive values of [d].
"""
import os
import subprocess
import sys
import tempfile
from collections import defaultdict
from fractions import Fraction as F

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import brute_force  # noqa: E402

ROOT = os.path.dirname(HERE)
SRC = os.path.join(ROOT, "src")
fails = 0


def report(ok, msg):
    global fails
    print(("PASS " if ok else "FAIL ") + msg)
    fails += 0 if ok else 1


def build(tmp, name, opt="-O3"):
    exe = os.path.join(tmp, name)
    subprocess.run(["g++", opt, "-std=c++17", "-pthread", "-I" + SRC, os.path.join(SRC, name + ".cpp"), "-o", exe],
                   check=True)
    return exe


def run(exe, *args):
    return subprocess.run([exe] + [str(a) for a in args], capture_output=True, text=True).stdout


def main():
    ref = {}
    for line in open(os.path.join(ROOT, "data", "reference_series.tsv")):
        if not line.startswith("#"):
            t = line.split()
            ref[int(t[0])] = (int(t[1]), int(t[2]))
    with tempfile.TemporaryDirectory() as tmp:
        # [a]
        es, sg = build(tmp, "end_steps"), build(tmp, "saw_gram")
        G = {}
        for line in run(sg, 12, 1).splitlines():
            t = line.split()
            if t and t[0] == "G":
                G[(int(t[1]), int(t[2]), int(t[3]))] = int(t[4])
        ok = True
        for line in run(es, 12, 1).splitlines():
            t = line.split()
            if t and t[0] == "E":
                n, v = int(t[1]), list(map(int, t[2:]))
                gram = G[(n, 1, n)] - (G[(n, 1, n - 1)] if n > 1 else 0)
                ok = ok and 6 * (v[0] - v[1]) == gram and 6 * sum(v) == ref[n][0] and v[2] == v[3] == v[4] == v[5]
        report(ok, "[a] end_steps 12: sum s_1.s_n = Gram sums (saw_gram 12), c_n = A001412, transverse symmetry")
        # [b]
        hs = build(tmp, "halfspace_enum")
        Bb, Rb, Ab, Hb = brute_force.halfspace_all(9)
        Bc, Rc, Ac, Hc = {}, {}, {}, {}
        for line in run(hs, 9, 1).splitlines():
            t = line.split()
            if not t or t[0] == "#":
                continue
            if t[0] == "H":
                Hc[int(t[1])] = int(t[2])
            else:
                {"B": Bc, "R": Rc, "A": Ac}[t[0]][(int(t[1]), int(t[2]))] = int(t[3])
        Ab = {k: v for k, v in Ab.items() if k[0] <= 9}
        report(Bb == Bc and Rb == Rc and Ab == Ac and Hb == Hc,
               "[b] halfspace_enum 9 = brute force (bridges by height, renewals at half height, arches, h_n)")
        # [c]
        cp = build(tmp, "cond_prefix")
        N, K = 9, 4
        acc = defaultdict(lambda: [0, 0])
        walk, seen, dirs = [(0, 0, 0), (1, 0, 0)], {(0, 0, 0), (1, 0, 0)}, [0]

        def visit():
            m = len(walk) - 1
            for k in range(1, min(m, K) + 1):
                a = acc[(k, m, tuple(dirs[:k]))]
                a[0] += 1
                a[1] += walk[m][0]
            if m == N:
                return
            e = walk[-1]
            for d, s in enumerate(brute_force.STEPS):
                q = (e[0] + s[0], e[1] + s[1], e[2] + s[2])
                if q in seen:
                    continue
                seen.add(q)
                walk.append(q)
                dirs.append(d)
                visit()
                dirs.pop()
                walk.pop()
                seen.discard(q)
        visit()
        ok = True
        for line in run(cp, N, K).splitlines():
            t = line.split()
            if not t or t[0] != "C":
                continue
            k, n = int(t[1]), int(t[2])
            items = [(F(v[1], v[0]), key) for key, v in acc.items() if key[0] == k and key[1] == n and v[0] > 0]
            ok = ok and F(int(t[7]), int(t[8])) == min(items)[0] and int(t[3]) == len(items)
            ok = ok and int(t[4]) == sum(1 for r, _ in items if r < 1)
        report(ok, "[c] cond_prefix 9 4 = brute force (number of prefixes, failures, minimal conditional mean)")
        # [d]
        me = build(tmp, "madras_exact")
        out = run(me, 8, 8)
        report("self-test (direct image counts, n <= 8): PASS" in out,
               "[d] madras_exact 8: multiplicity routine = direct preimage counts for every near walk, m <= 4")
        exact = {}
        for line in out.splitlines():
            t = line.split()
            if t and t[0] == "X" and t[1] == "8":
                exact[int(t[2])] = F(int(t[4]), int(t[3]))
        # [e]
        pv = build(tmp, "pivot")
        ok = True
        msg = []
        for n in (10, 20):
            t = [l for l in run(pv, n, 2000000, 10000, 3).splitlines() if l.startswith("R ")][0].split()
            mean, err = float(t[4]), float(t[5])
            ex = ref[n][1] / ref[n][0]
            ok = ok and abs(mean - ex) < 4 * err
            msg.append("n=%d: %.4f +- %.4f vs exact %.4f" % (n, mean, err, ex))
        report(ok, "[e] pivot: " + "; ".join(msg))
        # [f]
        mc = build(tmp, "madras_mc")
        ok = True
        msg = []
        for line in run(mc, 8, 3000000, 5, 1, 1000000).splitlines():
            t = line.split()
            if t and t[0] == "M" and int(t[1]) in (2, 4):
                m, kap, err = int(t[1]), float(t[4]), float(t[5])
                ok = ok and abs(kap - float(exact[m])) < 4 * err + 1e-9
                msg.append("m=%d: %.4f +- %.4f vs exact %.5f" % (m, kap, err, float(exact[m])))
        report(ok and len(msg) == 2, "[f] madras_mc N=8: " + "; ".join(msg))
    print("%s" % ("ALL PASS" if fails == 0 else "%d FAIL" % fails))
    sys.exit(1 if fails else 0)


if __name__ == "__main__":
    main()
