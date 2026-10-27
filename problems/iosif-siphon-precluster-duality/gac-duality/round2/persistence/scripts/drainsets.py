"""Drain sets: a set-valued, rate-independent certificate of non-persistence (agent Q1, round 2).

Environment: Python 3.12+, standard library only; exact rational Fourier-Motzkin elimination (endo._fm_feasible of
agent E, imported, not edited).

Definition (Theorem D of ../notes.md).  Let M be a mass-action network with species 1..p.  A *drain set* is a
nonempty set Z of species together with weights w in R^Z_{>0} and integers d in Z^Z_{>=1} such that
  (D1) for every i in Z some reaction has source exactly d_i e_i (no other species) and strictly decreases
       V = <w, x_Z>, i.e. <w, Delta|_Z> < 0;
  (D2) every reaction that strictly increases V has a source y with sum_{i in Z} y(i)/d_i > 1.
Then, for every choice of rate constants and every compatibility class P, there is delta > 0 such that every
solution in P with V(x(0)) < delta has x_Z(t) -> 0; such solutions exist in every class containing points with x_Z
small (e.g. x_Z = eta, x_{Z^c} = 1).  A drain species of agent E (Prop. E.4) is a drain set with |Z| = 1.

Decision procedure: for each Z (all of whose species have a pure-source reaction), each choice of d (from the
pure-source degrees), and each choice of one candidate draining reaction per species, an exact LP feasibility test
(w_i >= 1; <w, Delta_j|_Z> <= 0 for every reaction with sum y(i)/d_i <= 1; <w, Delta_{j_i}|_Z> <= -1 for the chosen
draining reactions).  Homogeneity makes '>= 1' and '<= -1' equivalent to the strict inequalities.
"""
import os
import sys
from fractions import Fraction
from itertools import combinations, product

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.normpath(os.path.join(HERE, "..", "..", "..", "computation", "scripts")))
from endo import _fm_feasible  # noqa: E402


def pure_degrees(p, rx):
    out = {i: {} for i in range(p)}
    for j, (src, _) in enumerate(rx):
        supp = [i for i in range(p) if src[i]]
        if len(supp) == 1:
            i = supp[0]
            out[i].setdefault(src[i], []).append(j)
    return out


def drain_sets(p, rx, maxsize=None, first_only=False):
    """List of (Z, d, w-feasible draining choice) for all drain sets; with first_only, stop at the first."""
    deltas = [tuple(b - a for a, b in zip(s, t)) for s, t in rx]
    pure = pure_degrees(p, rx)
    cand = [i for i in range(p) if pure[i]]
    found = []
    top = len(cand) if maxsize is None else min(maxsize, len(cand))
    for k in range(1, top + 1):
        for Z in combinations(cand, k):
            dz = [deltas[j] for j in range(len(rx))]
            for d in product(*[sorted(pure[i]) for i in Z]):
                cons_rows = []
                for j, (src, _) in enumerate(rx):
                    dj = [dz[j][i] for i in Z]
                    if not any(dj):
                        continue
                    if sum(Fraction(src[i], di) for i, di in zip(Z, d)) <= 1:
                        cons_rows.append(([-x for x in dj], 0))
                drains_per = []
                ok = True
                for i, di in zip(Z, d):
                    cs = [j for j in pure[i][di] if any(dz[j][t] for t in Z)]
                    if not cs:
                        ok = False
                        break
                    drains_per.append(cs)
                if not ok:
                    continue
                for choice in product(*drains_per):
                    cons = list(cons_rows)
                    for t in range(k):
                        e = [0] * k
                        e[t] = 1
                        cons.append((e, 1))
                    for j in set(choice):
                        cons.append(([-dz[j][i] for i in Z], 1))
                    try:
                        feas = _fm_feasible(cons, k)
                    except RuntimeError:
                        feas = False
                    if feas:
                        found.append((Z, d, choice))
                        if first_only:
                            return found
                        break
                else:
                    continue
                break
    return found


def drain_species_E(p, rx):
    """Agent E's drain species (Prop. E.4), re-implemented for the comparison."""
    out = []
    for i in range(p):
        d0s = [a[i] for a, b in rx if b[i] < a[i] and all(a[t] == 0 for t in range(p) if t != i)]
        if not d0s:
            continue
        d0 = min(d0s)
        if all(a[i] > d0 for a, b in rx if b[i] > a[i]):
            out.append(i)
    return out
