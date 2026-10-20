#!/usr/bin/env python3
"""Exact check of the stage-wise flux identity and of the one-sign-change lemma (Agent N, network.md, Lemma 2/3).

For a section of equilibrium stages (constant relative volatility, constant molar overflow) with vapor V, liquid L
and net upward component flows w (ANY signs), and any theta > 0:

    Psi(theta) := V - sum_p alpha_p w_p / (alpha_p - theta)
               =  L * u_{k+1}(theta) - (V * theta / S_k) * u_k(theta)          for every stage pair (k, k+1),

with u_k(theta) = sum_p alpha_p x_{k,p} / (alpha_p - theta) (liquid of stage k) and S_k = sum_p alpha_p x_{k,p}.
Consequences: along the section (bottom to top) the sign of u_k(theta) changes at most once, from - to + if Psi > 0,
from + to - if Psi < 0, never if Psi = 0.

Profiles are generated stage to stage in exact rational arithmetic (fractions.Fraction): given the liquid x_{k+1}
leaving stage k+1 and w, the vapor leaving stage k is y_k = (L x_{k+1} + w)/V and x_k is its bubble-point liquid
x_{k,p} = (y_{k,p}/alpha_p) / sum_q (y_{k,q}/alpha_q) (rational).  Profiles with a negative mole fraction are
discarded.  Net flows w of mixed signs (reverse flows) are included.

Environment: Python 3 standard library only.  Runtime: about 10 s.  Prints PASS/FAIL per claim, exits 1 on failure.
"""
import random
import sys
from fractions import Fraction as Fr


def u_val(alpha, x, th):
    return sum(a * xi / (a - th) for a, xi in zip(alpha, x))


def run(seed=1, trials=400):
    rng = random.Random(seed)
    n_prof = n_pairs = 0
    ok_id = ok_sign = True
    n_changes_seen = {0: 0, 1: 0}
    neg_w_profiles = 0
    while n_prof < trials:
        N = rng.randint(2, 5)
        alpha = sorted({Fr(rng.randint(11, 120), 10) for _ in range(N)}, reverse=True)
        N = len(alpha)
        if N < 2:
            continue
        # top liquid composition (random rational, positive)
        raw = [Fr(rng.randint(1, 50)) for _ in range(N)]
        s = sum(raw)
        xtop = [r / s for r in raw]
        V = Fr(rng.randint(5, 40), 10)
        L = Fr(rng.randint(1, 40), 10)
        # net upward flows: V - L = sum w ; random signs allowed
        w = [Fr(rng.randint(-20, 20), 40) for _ in range(N - 1)]
        w.append((V - L) - sum(w))
        # stage-to-stage downward from the top stage
        xs = [xtop]
        valid = True
        nst = rng.randint(2, 9)
        for _ in range(nst):
            xk1 = xs[-1]
            y = [(L * xi + wi) / V for xi, wi in zip(xk1, w)]
            if any(yi < 0 for yi in y) or sum(y) == 0:
                valid = False
                break
            t = [yi / a for yi, a in zip(y, alpha)]
            st = sum(t)
            xk = [ti / st for ti in t]
            xs.append(xk)
        if not valid:
            continue
        xs.reverse()  # now xs[0] = bottom stage, xs[-1] = top stage
        n_prof += 1
        if any(wi < 0 for wi in w):
            neg_w_profiles += 1
        # test several theta values, including interval points and feed-root-like points
        thetas = []
        for r in range(N - 1):
            lo, hi = alpha[r + 1], alpha[r]
            for frac in (Fr(1, 7), Fr(1, 2), Fr(5, 6)):
                thetas.append(lo + (hi - lo) * frac)
        thetas.append(alpha[-1] * Fr(1, 3))
        thetas.append(alpha[0] * Fr(3, 2))
        for th in thetas:
            Psi = V - sum(a * wi / (a - th) for a, wi in zip(alpha, w))
            us = [u_val(alpha, x, th) for x in xs]
            for k in range(len(xs) - 1):
                Sk = sum(a * xi for a, xi in zip(alpha, xs[k]))
                rhs = L * us[k + 1] - V * th / Sk * us[k]
                n_pairs += 1
                if rhs != Psi:
                    ok_id = False
            sg = [1 if uu >= 0 else -1 for uu in us]
            changes = [(sg[k], sg[k + 1]) for k in range(len(sg) - 1) if sg[k] != sg[k + 1]]
            if len(changes) > 1:
                ok_sign = False
            if changes:
                c = changes[0]
                if Psi > 0 and c != (-1, 1):
                    ok_sign = False
                if Psi < 0 and c != (1, -1):
                    ok_sign = False
                if Psi == 0:
                    ok_sign = False
            n_changes_seen[min(len(changes), 1)] += 1
    return ok_id, ok_sign, n_prof, n_pairs, n_changes_seen, neg_w_profiles


def main():
    ok_id, ok_sign, n_prof, n_pairs, seen, negw = run()
    print(f"profiles: {n_prof} (with some negative net flow: {negw}); stage pairs x theta tested: {n_pairs}")
    print(f"theta tests without / with a sign change of u: {seen[0]} / {seen[1]}")
    print(("PASS" if ok_id else "FAIL") + "  identity Psi(theta) = L u_{k+1} - (V theta/S_k) u_k (exact)")
    print(("PASS" if ok_sign else "FAIL") + "  at most one sign change of u along a section, direction given by sign(Psi)")
    sys.exit(0 if (ok_id and ok_sign) else 1)


if __name__ == "__main__":
    main()
