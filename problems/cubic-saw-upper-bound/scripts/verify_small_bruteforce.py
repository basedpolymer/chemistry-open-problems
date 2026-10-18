#!/usr/bin/env python3
"""Brute-force, independent re-check of a small SAWCERT1 certificate (pure Python, exact integers).

Environment: Python >= 3.8, no third-party package.  Runtime: ~1 min for d=3, K=8 (206 states);
~1 h for K=10 (2723 states).  Usage:  python3 verify_small_bruteforce.py file.cert

For every listed state s and every unit vector e not in F(s) this script scans ALL listed states t and ALL
2^d d! signed permutations g, keeps those with g(F(t)) subset of (F(s) u {0}) - e, and takes the minimal
weight among them.  It then checks Q * sum_e minweight <= P * V(s) exactly.  No path search, no canonical
form, no hashing: this is the hypothesis of Lemma 1 / Lemma 2 of notes.md checked in the most literal way.
Prints PASS/FAIL lines; exit code 0 iff all pass.
"""
import itertools
import struct
import sys
from fractions import Fraction


def main(path):
    data = open(path, "rb").read()
    if data[:8] != b"SAWCERT1":
        print("FAIL bad magic")
        return 1
    d, K = struct.unpack_from("<II", data, 8)
    N, P, Q = struct.unpack_from("<QQQ", data, 16)
    off = 40
    codes = struct.unpack_from("<%dQ" % N, data, off)
    off += 8 * N
    W = struct.unpack_from("<%dQ" % N, data, off)
    bits = 3 if 2 * d <= 8 else 4
    print("certificate d=%d K=%d N=%d claimed P/Q=%d/%d=%.9f" % (d, K, N, P, Q, P / Q))

    def forbidden(code):
        l = (code >> 57) & 63
        p = [0] * d
        F = []
        for i in range(l):
            t = (code >> (bits * i)) & ((1 << bits) - 1)
            ax, neg = t // 2, t % 2
            if ax >= d:
                raise ValueError("bad step")
            p = list(p)
            p[ax] -= (-1 if neg else 1)
            F.append(tuple(p))
        return F

    Fs = [forbidden(c) for c in codes]
    ok = True
    if any(w < 1 for w in W):
        print("FAIL weights must be >= 1")
        ok = False
    else:
        print("PASS all weights >= 1")
    start = [i for i, F in enumerate(Fs) if len(F) == 0]
    print(("PASS" if start else "FAIL") + " a state with empty forbidden set exists")
    ok = ok and bool(start)
    group = []
    for perm in itertools.permutations(range(d)):
        for signs in itertools.product((1, -1), repeat=d):
            group.append((perm, signs))
    # images of every F(t) under every g, as frozensets
    images = []
    for F in Fs:
        imgs = set()
        for perm, signs in group:
            img = []
            for x in F:
                y = [0] * d
                for a in range(d):
                    y[perm[a]] = signs[a] * x[a]
                img.append(tuple(y))
            imgs.add(frozenset(img))
        images.append(list(imgs))
    order = sorted(range(N), key=lambda i: W[i])  # scan children by increasing weight
    units = []
    for a in range(d):
        for s in (1, -1):
            u = [0] * d
            u[a] = s
            units.append(tuple(u))
    worst = Fraction(0)
    viol = 0
    for s in range(N):
        Fset = set(Fs[s])
        total = 0
        for e in units:
            if e in Fset:
                continue
            Wset = {tuple(x[a] - e[a] for a in range(d)) for x in list(Fset) + [tuple([0] * d)]}
            best = None
            for t in order:
                if len(Fs[t]) > len(Wset):
                    continue
                if any(img <= Wset for img in images[t]):
                    best = W[t]
                    break
            if best is None:
                print("FAIL no valid child for state", s, "direction", e)
                return 1
            total += best
        r = Fraction(total, W[s])
        if r > worst:
            worst = r
        if total * Q > P * W[s]:
            viol += 1
    print(("PASS" if viol == 0 else "FAIL") + " Collatz-Wielandt inequality (violations: %d); max ratio = %.12f" % (viol, float(worst)))
    ok = ok and viol == 0
    print(("PASS  => mu(Z^%d) <= %d/%d" % (d, P, Q)) if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1]))
