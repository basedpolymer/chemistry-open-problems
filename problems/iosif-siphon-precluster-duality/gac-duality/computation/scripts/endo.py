"""Exact tests of endotacticity, strong endotacticity and of a 'drain species' (a rate-independent certificate of
non-persistence), used by tests.py.

Environment: Python 3.12, standard library only; exact rational arithmetic (Fourier-Motzkin elimination over
fractions.Fraction; nullspaces from crn.py).

Definitions (Craciun-Nazarov-Pantea 2013; Gopalkrishnan-Miller-Shiu 2014).  For w in R^n, a reaction y -> y' is
w-essential if w.(y'-y) != 0.
  endotactic: for every w, every w-essential reaction whose source maximises w.y among the sources of w-essential
              reactions has w.(y'-y) < 0.
  strongly endotactic: for every w with w.(y'-y) != 0 for some reaction, every reaction whose source maximises w.y
              among ALL sources has w.(y'-y) <= 0, and at least one of them has w.(y'-y) < 0.
Both are decided exactly: the negations are finite unions of polyhedral cones (enumerated below); strict homogeneous
inequalities are replaced by '>= 1' (valid by scaling), equalities are removed by parametrising their solution space.
  drain species i (of a conservative network): some reaction k0 has source d0*e_i (only species i) and consumes X_i,
  and every reaction producing X_i has a source with X_i-coefficient > d0.  Then, for every rate constants,
  x_i' <= x_i^d0 (C x_i - c0) on each (bounded) compatibility class, so x_i -> 0 from every positive initial point
  with x_i < c0/C: the network is not persistent, whatever the rate constants.
"""
from fractions import Fraction
from itertools import combinations

import ecore


def _fm_feasible(cons, d):
    """cons: list of (a, b) meaning a.t >= b, a in Q^d.  Exact Fourier-Motzkin feasibility."""
    cons = [([Fraction(x) for x in a], Fraction(b)) for a, b in cons]
    for v in range(d):
        pos = [c for c in cons if c[0][v] > 0]
        neg = [c for c in cons if c[0][v] < 0]
        new = [c for c in cons if c[0][v] == 0]
        for ap, bp in pos:
            for an, bn in neg:
                fp, fn = ap[v], -an[v]
                a = [x / fp + y / fn for x, y in zip(ap, an)]
                a[v] = Fraction(0)
                new.append((a, bp / fp + bn / fn))
        # drop trivial duplicates
        seen = set()
        cons = []
        for a, b in new:
            key = (tuple(a), b)
            if key not in seen:
                seen.add(key)
                cons.append((a, b))
        if len(cons) > 4000:
            raise RuntimeError("Fourier-Motzkin blow-up")
    return all(Fraction(0) >= b for a, b in cons)


def _feasible(eqs, ineqs, n):
    """Exists w in R^n with e.w = 0 (e in eqs) and a.w >= b ((a, b) in ineqs)?"""
    if eqs:
        B = ecore.crn.nullspace([list(map(Fraction, e)) for e in eqs], n)  # w = sum t_k B_k
    else:
        B = [[Fraction(int(i == j)) for i in range(n)] for j in range(n)]
    d = len(B)
    if d == 0:
        return all(Fraction(0) >= b for a, b in ineqs)
    cons = []
    for a, b in ineqs:
        cons.append(([sum(Fraction(a[i]) * B[k][i] for i in range(n)) for k in range(d)], b))
    return _fm_feasible(cons, d)


def _vec(rx):
    srcs = [a for a, b in rx]
    deltas = [tuple(bb - aa for aa, bb in zip(a, b)) for a, b in rx]
    return srcs, deltas


def endotactic(n, rx):
    srcs, deltas = _vec(rx)
    ess = [k for k in range(len(rx)) if any(deltas[k])]
    for j in ess:
        others = [k for k in ess if k != j]
        for z in range(len(others) + 1):
            for Z in combinations(others, z):
                eqs = [deltas[k] for k in Z]
                ineqs = [(deltas[j], 1)]
                ineqs += [(tuple(a - b for a, b in zip(srcs[j], srcs[k])), 0) for k in others if k not in Z]
                if _feasible(eqs, ineqs, n):
                    return False
    return True


def strongly_endotactic(n, rx):
    srcs, deltas = _vec(rx)
    r = len(rx)
    if not any(any(d) for d in deltas):
        return True
    # (A) a reaction with globally maximal source points outward
    for j in range(r):
        if not any(deltas[j]):
            continue
        ineqs = [(deltas[j], 1)] + [(tuple(a - b for a, b in zip(srcs[j], srcs[k])), 0) for k in range(r)]
        if _feasible([], ineqs, n):
            return False
    # (B) all reactions with globally maximal source are w-orthogonal, some reaction is not
    S = sorted(set(srcs))
    for t in range(1, len(S) + 1):
        for T in combinations(S, t):
            c0 = T[0]
            eqs = [tuple(a - b for a, b in zip(c, c0)) for c in T[1:]]
            eqs += [deltas[k] for k in range(r) if srcs[k] in T]
            base = [(tuple(a - b for a, b in zip(c0, c)), 1) for c in S if c not in T]
            for l in range(r):
                if srcs[l] in T or not any(deltas[l]):
                    continue
                for sg in (1, -1):
                    if _feasible(eqs, base + [(tuple(sg * x for x in deltas[l]), 1)], n):
                        return False
    return True


def drain_species(n, rx):
    """List of drain species (see module docstring)."""
    out = []
    for i in range(n):
        d0s = [a[i] for a, b in rx if b[i] < a[i] and all(a[t] == 0 for t in range(n) if t != i)]
        if not d0s:
            continue
        d0 = min(d0s)
        if all(a[i] > d0 for a, b in rx if b[i] > a[i]):
            out.append(i)
    return out
