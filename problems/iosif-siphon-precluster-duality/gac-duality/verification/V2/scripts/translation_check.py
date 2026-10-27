"""Verifier V2: do WF-EC networks fit Johnston's translated networks (arXiv:1305.5845, Bull. Math. Biol. 76 (2014),
Def. 4.1, Thm 4.1)?

Fact used (proved in verification.md): a translation keeps the reaction vectors (bijection h1, Def. 4.1(1)), so
Gamma = Ytilde * I (I = incidence matrix of the translated graph).  If the translation has structural deficiency
delta~ = |C~| - l~ - s = 0 then ker Gamma = ker I (dimensions r - s = r - (|C~| - l~)), the cycle space of a directed
graph.  Every minimal-support vector (circuit) of a graph's cycle space is, up to scaling, a signed cycle with all
nonzero entries of absolute value 1.  So: if some circuit of ker Gamma has two nonzero entries of different absolute
value, no translation with delta~ = 0 exists, and Thm 4.1 of Johnston (which needs delta~ = delta~_K = 0) does not
apply.

Checks: (1) the WF-EC network W = {A -> C, 2C -> 2A} (ker Gamma = R(2,1)); (2) the fraction of random WF-EC networks
(own generator) with such a circuit.  Exact.  Environment: WSL Fedora 44, Python 3.12.
"""
import sys, os, time, random
from itertools import combinations
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from v2lib import *  # noqa

T0 = time.time()
FAILS = []


def check(name, cond, info=""):
    print(("PASS " if cond else "FAIL ") + name + ((" -- " + info) if info else ""))
    if not cond:
        FAILS.append(name)


def circuits(M):
    """All minimal-support vectors of ker M (up to scaling), by support enumeration."""
    ncol = len(M[0])
    out, sups = [], []
    for size in range(1, ncol + 1):
        for T in combinations(range(ncol), size):
            if any(set(S) <= set(T) for S in sups):
                continue
            sub = [[row[j] for j in T] for row in M]
            ns = nullspace(sub)
            if len(ns) == 1 and all(x != 0 for x in ns[0]):
                out.append((T, ns[0]))
                sups.append(T)
    return out


def non_graphic(N):
    for T, v in circuits(Gamma(N)):
        if len({abs(x) for x in v}) > 1:
            return (T, v)
    return None


W = parse("A->C, 2C->2A", ["A", "C"])
check("W = {A -> C, 2C -> 2A} is WF-EC (s = 1, m = 2) and satisfies (H3)", is_WFEC(W) and rows_distinct(W))
ns = nullspace(Gamma(W))
check("ker Gamma = R(2,1): not a cycle-space vector, so W has no translation with delta~ = 0",
      len(ns) == 1 and ns[0][0] == 2 * ns[0][1] and non_graphic(W) is not None, str([str(t) for t in ns[0]]))
u = u_of_k(W, [F(3), F(5)])
check("Theorem T on W, k = (3,5): ker A_k = R u, u > 0 (steady states: 3 x_A = 2 * 5 x_C^2 up to the common factor)",
      u is not None and u[0] / u[1] == F(2 * 5, 3), str([str(t) for t in u]) if u else "")
rng = random.Random(77)
tot = bad = 0
example = None
while tot < 300:
    n = rng.randint(3, 5)
    m = rng.randint(2, min(4, n))
    N = gen_wfec(rng, n, m)
    if N is None:
        continue
    tot += 1
    ng = non_graphic(N)
    if ng:
        bad += 1
        if example is None or len(N) < len(example[0]):
            example = (N, ng)
print(f"random WF-EC networks: {tot}; with a non-graphic circuit of ker Gamma (no delta~ = 0 translation): {bad}")
if example:
    print("   smallest such:", example[0], "circuit", example[1][0], [str(t) for t in example[1][1]])
check("some random WF-EC networks admit no delta~ = 0 translation", bad > 0)
print(f"runtime {time.time() - T0:.1f} s; {len(FAILS)} FAIL")
sys.exit(1 if FAILS else 0)
