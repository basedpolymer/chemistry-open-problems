#!/usr/bin/env python3
"""Tamper test for the checker: every corrupted copy of a valid LATCERT1 certificate must be rejected.

Environment: Python 3 (standard library).  Runtime: a few seconds for certificates of <= 10^5 states.
Usage: python3 tamper_test.py <latsaw_check binary> <valid .lcert> <scratch dir>
Corruptions: claimed bound P lowered by one unit (after setting P to the exact ceiling of the max ratio), one
weight lowered by 10 %, a weight set to 0, a child pointer redirected to another state, a child pointer removed,
the empty walk (start state) removed, two codes swapped (order), a step label out of range.
Prints PASS/FAIL per corruption; exit code 0 iff the original passes and every corrupted copy fails.
"""
import os
import struct
import subprocess
import sys
from fractions import Fraction


def run(binary, path):
    r = subprocess.run([binary, path, '2'], capture_output=True, text=True)
    return r.returncode == 0


def main():
    binary, cert, scratch = sys.argv[1:4]
    os.makedirs(scratch, exist_ok=True)
    data = bytearray(open(cert, 'rb').read())
    lat, K = struct.unpack_from('<II', data, 8)
    N, P, Q = struct.unpack_from('<QQQ', data, 16)
    nd = 12 if lat == 1 else 4
    bits = 4 if lat == 1 else 2
    oc, ov, och = 40, 40 + 8 * N, 40 + 16 * N
    codes = list(struct.unpack_from(f'<{N}Q', data, oc))
    V = list(struct.unpack_from(f'<{N}Q', data, ov))
    ch = list(struct.unpack_from(f'<{N * nd}I', data, och))
    ok = True
    # tight claim: P = ceil(Q * max ratio) recomputed from the explicit automaton (ratios as exact fractions)
    best = Fraction(0)
    worst = 0
    for s in range(N):
        tot = sum(V[c] for c in ch[s * nd:(s + 1) * nd] if c != 0xffffffff)
        r = Fraction(tot, V[s])
        if r > best:
            best, worst = r, s
    Pt = -((-best.numerator * Q) // best.denominator)
    base = bytearray(data)
    struct.pack_into('<Q', base, 24, Pt)
    p0 = os.path.join(scratch, 'tamper_base.lcert')
    open(p0, 'wb').write(base)
    good = run(binary, p0)
    print(('PASS' if good else 'FAIL') + f' original with tight claim P = {Pt} accepted')
    ok &= good

    def expect_fail(name, mutate):
        nonlocal ok
        d = bytearray(base)
        mutate(d)
        p = os.path.join(scratch, 'tamper_' + name + '.lcert')
        open(p, 'wb').write(d)
        rejected = not run(binary, p)
        print(('PASS' if rejected else 'FAIL') + f' corrupted copy rejected: {name}')
        ok &= rejected
        os.remove(p)

    expect_fail('P_minus_1', lambda d: struct.pack_into('<Q', d, 24, Pt - 1))
    expect_fail('weight_minus_10pct', lambda d: struct.pack_into('<Q', d, ov + 8 * worst, V[worst] * 9 // 10 if V[worst] > 10 else V[worst]))
    expect_fail('weight_zero', lambda d: struct.pack_into('<Q', d, ov + 8 * (N // 2), 0))
    # redirect a child pointer of the worst state to a state with a strictly larger forbidden set requirement:
    # point it at the longest listed walk (almost never contained) -- containment or inequality must fail
    longest = max(range(N), key=lambda i: (codes[i] >> 57) & 63)
    js = [j for j in range(nd) if ch[worst * nd + j] != 0xffffffff]
    target = None
    for s in range(N):
        for j in range(nd):
            c = ch[s * nd + j]
            if c != 0xffffffff and c != longest and ((codes[c] >> 57) & 63) < ((codes[longest] >> 57) & 63):
                target = (s, j)
                break
        if target:
            break
    if target:
        expect_fail('child_redirected', lambda d: struct.pack_into('<I', d, och + 4 * (target[0] * nd + target[1]), longest))
    expect_fail('child_removed', lambda d: struct.pack_into('<I', d, och + 4 * (worst * nd + js[0]), 0xffffffff))
    start = codes.index(1 << 63)

    def drop_start(d):
        struct.pack_into('<Q', d, oc + 8 * start, (1 << 63) | (1 << 57) | 0)   # becomes a duplicate of the 1-step walk
    expect_fail('start_state_removed', drop_start)
    if N > 3:
        def swap(d):
            struct.pack_into('<Q', d, oc + 8 * 1, codes[2])
            struct.pack_into('<Q', d, oc + 8 * 2, codes[1])
        expect_fail('codes_swapped', swap)
    lmax = max((c >> 57) & 63 for c in codes)
    i = next(i for i in range(N) if (codes[i] >> 57) & 63 == lmax)
    expect_fail('label_out_of_range' if lat == 0 else 'label_out_of_range',
                lambda d: struct.pack_into('<Q', d, oc + 8 * i, codes[i] | (((1 << bits) - 1) << 0)) if lat == 1 else
                struct.pack_into('<Q', d, oc + 8 * i, codes[i] | (1 << (bits * lmax))))
    os.remove(p0)
    print('ALL PASS' if ok else 'SOME FAIL')
    return 0 if ok else 1


if __name__ == '__main__':
    sys.exit(main())
