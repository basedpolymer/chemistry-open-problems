#!/usr/bin/env python3
"""Checks the interval arithmetic of src/interval.hpp (both rounding variants) against exact rational
arithmetic: for random operands (and adversarial ones: tiny/huge magnitudes, exact ties, cancellation),
the computed interval must contain the exact range of the operation.  Prints PASS/FAIL per variant.
Environment: g++ (C++20), Python 3.13 (fractions).  Runtime: ~20 s.
Usage: python3 test_interval.py [NTESTS]
"""
import sys, os, random, subprocess, tempfile
from fractions import Fraction as F

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, "..", "src")

def build(variant, out):
    flags = ["-O2", "-std=c++20", "-ffp-contract=off", "-fno-fast-math"]
    if variant == "upward":
        flags += ["-frounding-math", "-DHEX_UPWARD"]
    subprocess.run(["g++"] + flags + ["-o", out, os.path.join(SRC, "test_interval.cpp")], check=True)

def rnd(r):
    k = r.random()
    if k < 0.1: return r.choice([0.0, 1.0, -1.0, 0.5, 3.0, 1e-300, -1e-300, 1e300])
    if k < 0.3: return r.choice([-1, 1]) * r.random() * 10 ** r.randint(-30, 30)
    if k < 0.5: return float(r.randint(-1000, 1000)) / 2 ** r.randint(0, 60)
    return r.uniform(-4, 4)

def exact(op, a, b, c, d):
    A = [F(a), F(b)]; B = [F(c), F(d)]
    if op == "add": return A[0] + B[0], A[1] + B[1]
    if op == "sub": return A[0] - B[1], A[1] - B[0]
    if op == "mul":
        p = [x * y for x in A for y in B]; return min(p), max(p)
    if op == "div":
        p = [x / y for x in A for y in B]; return min(p), max(p)

def main():
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 200000
    r = random.Random(20261024)
    tests = []
    for _ in range(n):
        op = r.choice(["add", "sub", "mul", "div", "sqrt"])
        a, b = sorted([rnd(r), rnd(r)]); c, d = sorted([rnd(r), rnd(r)])
        if op == "div" and c <= 0 <= d: c, d = abs(c) + 0.25, abs(d) + 0.5
        if op == "sqrt": a, b = abs(a), abs(a) + abs(b)
        tests.append((op, a, b, c, d))
    inp = "".join("%s %s %s %s %s\n" % (t[0], t[1].hex(), t[2].hex(), t[3].hex(), t[4].hex()) for t in tests)
    allok = True
    with tempfile.TemporaryDirectory() as td:
        for variant in ("nextafter", "upward"):
            exe = os.path.join(td, "ti_" + variant)
            build(variant, exe)
            out = subprocess.run([exe], input=inp, capture_output=True, text=True, check=True).stdout.split("\n")
            bad = 0
            for t, line in zip(tests, out):
                lo, hi = [float.fromhex(x) for x in line.split()]
                op, a, b, c, d = t
                if op == "sqrt":
                    ok = (lo <= 0 or F(lo) ** 2 <= F(a)) and F(hi) ** 2 >= F(b) and lo <= hi
                else:
                    elo, ehi = exact(op, a, b, c, d)
                    ok = (lo == float("-inf") or F(lo) <= elo) and (hi == float("inf") or F(hi) >= ehi)
                if not ok:
                    bad += 1
                    if bad < 5: print("  violation", variant, t, lo, hi)
            print(("PASS" if bad == 0 else "FAIL"), "interval arithmetic variant %s: %d random operations, %d violations" % (variant, len(tests), bad))
            allok = allok and bad == 0
    sys.exit(0 if allok else 1)

if __name__ == "__main__":
    main()
