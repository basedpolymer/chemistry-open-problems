#!/usr/bin/env python3
"""Exact combinatorial cross-checks. Python >=3.10 standard library only.
Runtime: approximately 5 seconds, one thread, <100 MB RAM.
No network. Prints PASS/FAIL and exits nonzero on any failed assertion.
"""
from collections import Counter
from fractions import Fraction
from itertools import product
from pathlib import Path
import json

ROOT = Path(__file__).resolve().parents[1]


def check(ok, label):
    print(('PASS ' if ok else 'FAIL ') + label, flush=True)
    if not ok:
        raise AssertionError(label)


def enumerate_polynomial(charges):
    """Independent exhaustive path enumeration: exponent = sum of site charges squared."""
    out = Counter()
    n = len(charges)
    if not n:
        return Counter({0: 1})
    for mask in range(1 << (n-1)):
        energy = 0
        charge = charges[0]
        for i in range(n-1):
            if mask & (1 << i):
                energy += charge*charge
                charge = 0
            charge += charges[i+1]
        out[energy + charge*charge] += 1
    return out


def mul(a, b, shift=0):
    c = Counter()
    for i, x in a.items():
        for j, y in b.items():
            c[i+j+shift] += x*y
    return c


def recurrence(charges, q):
    z = [Fraction(1)]
    for j in range(1, len(charges)+1):
        charge = 0
        total = Fraction(0)
        for i in range(j-1, -1, -1):
            charge += charges[i]
            total += z[i] * q**(charge*charge)
        z.append(total)
    return z


def main():
    polys = []
    for n in range(1, 9):
        a = [enumerate_polynomial([1]*k) for k in range(n+1)]
        rhs = mul(a[n], a[n])  # a cut exactly at the sign change
        for i in range(n):
            for j in range(n):
                rhs.update(mul(a[i], a[j], (i-j)**2))
        direct = enumerate_polynomial([1]*n + [-1]*n)
        check(rhs == direct, f'diblock polynomial identity for N={2*n}')
        check(sum(rhs.values()) == 2**(2*n-1), f'path count for N={2*n}')
        polys.append({'half_length': n, 'coefficients': dict(sorted(rhs.items()))})
    q = Fraction(4, 5)
    for n in range(1, 9):
        homogeneous = recurrence([1]*n, q)
        for charges in product([-1, 1], repeat=n):
            z = recurrence(charges, q)
            check_value = sum(c*q**k for k,c in enumerate_polynomial(charges).items())
            assert z[n] == check_value
            assert z[n] >= homogeneous[n]
            for cut in range(1, n):
                right = recurrence(charges[cut:], q)[-1]
                assert z[n] >= z[cut]*right
                # Exact no-cut decomposition, independent of the recurrence for z[n].
                nocut = Fraction(0)
                for left in range(cut):
                    for end in range(cut+1, n+1):
                        charge = sum(charges[left:end])
                        nocut += z[left]*q**(charge*charge)*recurrence(charges[end:],q)[-1]
                assert z[n] == z[cut]*right + nocut
        check(True, f'all 2^{n} charge sequences: enumeration, domination, gluing and cut identity')
    path = ROOT/'certificates'/'exact_polynomials.json'
    if path.exists():
        recorded = json.loads(path.read_text())
        check(recorded == json.loads(json.dumps(polys)), 'stored exact integer polynomial certificate')
    else:
        path.write_text(json.dumps(polys, indent=2)+'\n')
        check(True, 'created exact integer polynomial certificate')


if __name__ == '__main__':
    try:
        main()
    except Exception as exc:
        print(f'FAIL {type(exc).__name__}: {exc}', flush=True)
        raise
