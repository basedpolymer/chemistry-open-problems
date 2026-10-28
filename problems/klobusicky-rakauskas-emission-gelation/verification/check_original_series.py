#!/usr/bin/env python3
"""Verifier: original equation (37), exact Taylor arithmetic.

Author: GPT-6-Astra, verification agent.
Environment: Python 3.12+ standard library. Typical runtime: under one minute.
Run from any directory: python verification/check_original_series.py

Compute Taylor coefficients in the original coupled system, including its
actual mass denominator; do not insert the proposed linear mass law while
solving. Only after constructing the coefficients, compare against the closed
formula and the mass/count laws. Rational arithmetic only, no floating point.

At Taylor order r, only clusters assembled from at most r+1 initial k-mers can
occur: this follows inductively from the gain term. Thus truncation above
ell+(k-ell)*(ORDER+1) is exact through order ORDER, including the mass series.
Finite tests supplement the analytic verification; they are not its proof.
"""
from fractions import Fraction as Q
import json
from pathlib import Path

ORDER = 7
PAIRS = [(1, 2), (1, 5), (2, 3), (2, 7), (3, 4), (3, 8),
         (7, 8), (11, 12), (19, 20), (100, 101)]
checks = []


def check(name, value):
    ok = bool(value)
    checks.append({"claim": name, "pass": ok})
    print(("PASS " if ok else "FAIL ") + name, flush=True)
    if not ok:
        raise AssertionError(name)


def conv(a, b, r):
    return sum((a[h]*b[r-h] for h in range(r+1)), Q(0))


def choose(a, n):
    value = Q(1)
    for i in range(n):
        value *= (a-i)/Q(i+1)
    return value


def original_series(ell, k):
    top = ell+(k-ell)*(ORDER+1)
    u = {n: [Q(n == k)]+[Q(0)]*ORDER for n in range(ell+1, top+1)}
    mass = [Q(0)]*(ORDER+1)
    invmass = [Q(0)]*(ORDER+1)
    for r in range(ORDER):
        mass[r] = sum(n*u[n][r] for n in u)
        invmass[r] = (1/mass[0] if r == 0 else
                      -sum(mass[h]*invmass[r-h] for h in range(1, r+1))/mass[0])
        invsquare = [conv(invmass, invmass, h) for h in range(r+1)]
        for n in u:
            gain = [Q(0)]*(r+1)
            for i in range(ell+1, n):
                b = n+ell-i
                for h in range(r+1):
                    gain[h] += i*b*conv(u[i], u[b], h)
            rhs = conv(gain, invsquare, r)-2*n*conv(u[n], invmass, r)
            u[n][r+1] = rhs/Q(r+1)
    mass[ORDER] = sum(n*u[n][ORDER] for n in u)
    return u, mass


for ell, k in PAIRS:
    u, mass = original_series(ell, k)
    d = k-ell
    p = Q(2*k-ell, ell)
    expected = {n: [Q(0)]*(ORDER+1) for n in u}
    for j in range(1, ORDER+2):
        n = ell+d*j
        lead = Q(ell**(j-1)*k**j, j*n*k**(2*(j-1)))*choose(p*j, j-1)
        for r in range(j-1, ORDER+1):
            h = r-j+1
            expected[n][r] = lead*choose(Q(2*n, ell), h)*Q(-ell, k)**h
    check(f"original coupled ODE Taylor coefficients through {ORDER}: ell={ell}, k={k}",
          u == expected)
    check(f"actual mass Taylor series: ell={ell}, k={k}",
          mass == [Q(k), Q(-ell)]+[Q(0)]*(ORDER-1))
    count = [sum(u[n][r] for n in u) for r in range(ORDER+1)]
    check(f"actual count Taylor series: ell={ell}, k={k}",
          count == [Q(1), Q(-1)]+[Q(0)]*(ORDER-1))
    check(f"off-support coefficients vanish: ell={ell}, k={k}",
          all(not any(row) for n, row in u.items() if (n-ell) % d))

report = {"author": "GPT-6-Astra, verifier",
          "arithmetic": "exact fractions; Python standard library",
          "taylor_order": ORDER, "checks": checks,
          "all_pass": all(row["pass"] for row in checks)}
path = Path(__file__).with_name("original-series-results.json")
path.write_text(json.dumps(report, indent=2)+"\n", encoding="utf-8")
print(f"PASS ALL ({len(checks)} checks); report: {path}")
