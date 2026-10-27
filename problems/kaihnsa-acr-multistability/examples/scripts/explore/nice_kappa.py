#!/usr/bin/env python3
"""
nice_kappa.py -- numerical exploration: for a given design network (X+Y->2Y first, Y->W second), scan rate
constants over a grid of simple rationals and list those with two exponentially stable positive steady states,
ranked by simplicity.  Usage: python nice_kappa.py <network-id>
"""
import numpy as np, itertools, sys
from fractions import Fraction

NETS = {
    # coordinator H2-type: X+Y->2Y, Y->2Z, 0->X+Z, X+Z->2Z, 2Z->X+Z, Y+Z->Y   (7 complexes, 6 reactant complexes)
    'H2': [((1, 1, 0), (0, 2, 0)), ((0, 1, 0), (0, 0, 2)), ((0, 0, 0), (1, 0, 1)), ((1, 0, 1), (0, 0, 2)),
           ((0, 0, 2), (1, 0, 1)), ((0, 1, 1), (0, 1, 0))],
    # coordinator brief network: X+Y->2Y, Y->Z, X+Z->2Z, 0->X, 2Z->2X, Y+Z->Y (10 complexes)
    'B': [((1, 1, 0), (0, 2, 0)), ((0, 1, 0), (0, 0, 1)), ((1, 0, 1), (0, 0, 2)), ((0, 0, 0), (1, 0, 0)),
          ((0, 0, 2), (2, 0, 0)), ((0, 1, 1), (0, 1, 0))],
    # H1: X+Y->2Y, Y->2Z, 0->X+Z, X+Z->0, Y+Z->Y, 2Z->X+Z
    'H1': [((1, 1, 0), (0, 2, 0)), ((0, 1, 0), (0, 0, 2)), ((0, 0, 0), (1, 0, 1)), ((1, 0, 1), (0, 0, 0)),
           ((0, 1, 1), (0, 1, 0)), ((0, 0, 2), (1, 0, 1))],
    # H3: X+Y->2Y, Y->2Z, Y+Z->Y, 2Z->X+Z, 0->2Z, X+Z->0
    'H3': [((1, 1, 0), (0, 2, 0)), ((0, 1, 0), (0, 0, 2)), ((0, 1, 1), (0, 1, 0)), ((0, 0, 2), (1, 0, 1)),
           ((0, 0, 0), (0, 0, 2)), ((1, 0, 1), (0, 0, 0))],
}
rx = NETS[sys.argv[1]]
vals = [Fraction(1, 50), Fraction(1, 20), Fraction(1, 10), Fraction(1, 5), Fraction(1, 4), Fraction(1, 2), Fraction(1),
        Fraction(2), Fraction(3), Fraction(4), Fraction(5), Fraction(10), Fraction(20)]


def simplicity(ks):
    return sum(abs(k.numerator) + k.denominator for k in ks)


def analyse(K):
    c = K[1] / K[0]
    PX0 = np.zeros(3); PX1 = np.zeros(2); PZ0 = np.zeros(3); PZ1 = np.zeros(2)
    for j, (a, b) in enumerate(rx):
        v = np.array(b) - np.array(a)
        coef = K[j] * c ** a[0]
        dz = a[2]
        if a[1] == 0:
            PX0[dz] += coef * v[0]; PZ0[dz] += coef * v[2]
        else:
            PX1[dz] += coef * v[0]; PZ1[dz] += coef * v[2]
    R = np.zeros(4)
    for i in range(3):
        for j in range(2):
            R[i + j] += PX0[i] * PZ1[j] - PZ0[i] * PX1[j]
    z = np.roots(R[::-1])
    z = z[(np.abs(z.imag) < 1e-9) & (z.real > 1e-12)].real
    out = []
    for zz in z:
        den = PX1[0] + PX1[1] * zz
        y = -(PX0[0] + PX0[1] * zz + PX0[2] * zz * zz) / den
        if y <= 0:
            continue
        x = np.array([c, y, zz])
        J = np.zeros((3, 3))
        for j, (a, b) in enumerate(rx):
            a_ = np.array(a); v = np.array(b) - a_
            for i in range(3):
                if a_[i] > 0:
                    J[:, i] += v * K[j] * a_[i] * np.prod(x ** (a_ - np.eye(3, dtype=int)[i]))
        ev = np.linalg.eigvals(J)
        out.append((x, ev, np.max(ev.real) < 0))
    return out


hits = []
for ks in itertools.product(vals, repeat=len(rx)):
    K = np.array([float(k) for k in ks])
    out = analyse(K)
    if sum(st for _, _, st in out) >= 2:
        # robustness margin: min |Re ev| / max |ev| over stable ones, and spread of steady states
        marg = min(-np.max(ev.real) / np.max(np.abs(ev)) for _, ev, st in out if st)
        hits.append((simplicity(ks), -marg, ks, out))
hits.sort(key=lambda h: (h[0], h[1]))
print("hits:", len(hits))
for h in hits[:15]:
    print([str(k) for k in h[2]], "margin %.3g" % -h[1],
          [(np.round(x, 4).tolist(), np.round(ev, 3).tolist()) for x, ev, st in h[3]])
