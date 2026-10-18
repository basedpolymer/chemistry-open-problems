#!/usr/bin/env python3
"""Source-convention and rational-root checks, added by the solver's internal-review sub-agent
(GPT-6-Astra); despite its name it is not the PROTOCOL.md section 4 independent verification.
Environment: Python >=3.10 standard library only; one thread, <100 MB RAM.
Runtime: approximately 5 seconds. No network or floating-point arithmetic.
PASS/FAIL is printed for each claim; any failed claim exits nonzero.
This finite check is supplementary to the independent analytic verification.
"""
from fractions import Fraction as F
from itertools import product
from pathlib import Path
import json
import time

ROOT = Path(__file__).resolve().parents[1]


def check(ok, label):
    print(('PASS ' if ok else 'FAIL ') + label, flush=True)
    if not ok:
        raise AssertionError(label)


def enumerate_original(charges, q):
    """Sum the original i<j pair Hamiltonian, without a block recurrence."""
    n = len(charges)
    if n == 0:
        return F(1), []
    total = F(0)
    cut_weights = [F(0) for _ in range(n-1)]
    for steps in product((0, 1), repeat=n-1):
        positions = [1]
        for step in steps:
            positions.append(positions[-1]+step)
        h = sum(charges[i]*charges[j]
                for i in range(n) for j in range(i+1, n)
                if positions[i] == positions[j])
        site_charge = {}
        for charge, position in zip(charges, positions):
            site_charge[position] = site_charge.get(position, 0)+charge
        squared_charge = sum(charge**2 for charge in site_charge.values())
        assert squared_charge == n+2*h
        reduced_weight = q**n * q**(2*h)
        assert reduced_weight == q**squared_charge
        total += reduced_weight
        for i, step in enumerate(steps):
            if step:
                cut_weights[i] += reduced_weight
    return total, cut_weights


def exponential_bounds(x, degree=256):
    """Rational Taylor bounds: odd truncation <= exp(-x) <= even truncation.

    Taylor's integral remainder has the stated sign for every x>=0;
    decreasing terms throughout the series are not required.
    """
    assert x >= 0 and degree % 2 == 0
    term = total = F(1)
    for k in range(1, degree+1):
        term *= -x/k
        total += term
    upper = total
    lower = total + term*(-x)/(degree+1)
    assert 0 < lower <= upper
    return lower, upper


def rational_series_bounds(beta, cutoff=12):
    lower = upper = F(0)
    for k in range(1, cutoff+1):
        lo, hi = exponential_bounds(beta*k*k/2)
        lower += lo
        upper += hi
    first_upper = exponential_bounds(beta*(cutoff+1)**2/2)[1]
    ratio_upper = exponential_bounds(beta*(2*cutoff+3)/2)[1]
    assert ratio_upper < 1
    return lower, upper+first_upper/(1-ratio_upper)


def main():
    start = time.perf_counter()
    q = F(3, 4)
    homogeneous = [enumerate_original((1,)*n, q)[0] for n in range(8)]
    for n in range(1, 8):
        for charges in product((-1, 1), repeat=n):
            total, cuts = enumerate_original(charges, q)
            assert total >= homogeneous[n]
            for i, cut_weight in enumerate(cuts, 1):
                left = enumerate_original(charges[:i], q)[0]
                right = enumerate_original(charges[i:], q)[0]
                assert cut_weight == left*right
                bound = sum((1/homogeneous[a] for a in range(1, i+1)), F(0))
                bound *= sum((1/homogeneous[b] for b in range(1, n-i+1)), F(0))
                assert total-cut_weight <= cut_weight*bound
        check(True, f'independent original pair energy, normalization, cut odds for all words N={n}')

    # A deliberately coarser bracket than the 50-decimal Arb certificate:
    # acceptance here uses only integer/fraction operations and Taylor signs.
    lo = F('0.69813170079919')
    hi = F('0.69813170079921')
    check(rational_series_bounds(lo)[0] > 1,
          'rational Taylor lower sum at 0.69813170079919 exceeds one')
    check(rational_series_bounds(hi)[1] < 1,
          'rational Taylor upper sum plus tail at 0.69813170079921 is below one')
    parameters = json.loads((ROOT/'certificates'/'parameters.json').read_text())
    check(lo < F(parameters['beta_lower']) < F(parameters['beta_upper']) < hi,
          'independent rational root bracket contains supplied Arb bracket')
    print(f'PASS all independent finite checks; elapsed={time.perf_counter()-start:.2f}s', flush=True)


if __name__ == '__main__':
    try:
        main()
    except Exception as exc:
        print(f'FAIL {type(exc).__name__}: {exc}', flush=True)
        raise
