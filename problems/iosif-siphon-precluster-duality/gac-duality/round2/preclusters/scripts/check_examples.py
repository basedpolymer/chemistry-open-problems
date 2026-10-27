"""Examples for the corrected preclusters (agent Q2, Claude Opus 5.5, xhigh).

Environment: Python 3.12 (WSL Fedora 44; sympy only for the one-block certificates). Exact arithmetic.
Prints PASS/FAIL per claim; exit status 1 on any FAIL.  Runtime: a few seconds.

Claims (ids as in ../notes.md and ../candidates.md):
  X1  chain P <= Pcs <= pi <= pi_plus <= CF_upper on every example, and pi = Pcs v M* (Props 2.1, 2.2).
  X2  B3: P has 2 blocks, Pcs = pi = CF = 1 block; the binary split of B3 has P = 1 block (P not split-invariant).
  X3  N4.12: P = Pcs = singletons, pi = CF = M* = 1 block; with R1 duplicated P is 1 block (P not dup-invariant).
  X4  B2 with Z->X duplicated: pi has 2 blocks, CF has 1 block (pi, Pcs not duplication-invariant; CF != pi).
  X5  K3 (complete digraph on 3 species, WF-FC): pi = 3 blocks (= source classes), CF = 1 block (certificate).
  X6  XY4 = {X+Y->2X, X+Y->2Y, X+2Y->2X+Y, X+2Y->3Y}: kappa_V generic = 1 block, CF = pi = 2 blocks (split witness);
      Z4 = {2Y->2X, 2X->2Y, X+Y->2X, X+Y->2Y}: kappa_gen = CF = 3 blocks, kappa at V = (1,1,2,2) coarser (2 blocks).
  X9  B2dup: pi = Pcs = 2 blocks; after splitting the choice {R1,R2} behind W, pi = Pcs = 1 block (not split-invariant).
  X7  duplication / choice-splitting invariance of the CF bounds on all examples (pulled back), when exact.
  X8  reaction side: M* <= CF and [CF = M*  iff  every doubling pair is fully coupled] on all examples;
      species side (N conservative): M(N) <= CF(N*), equality iff every twin pair is conservation-coupled.
"""
import sys
import time
import random
import platform

import pclust as pc
from pclust import Network, refines

T0 = time.time()
FAIL = 0


def check(cid, cond, msg=""):
    global FAIL
    print(f"{'PASS' if cond else 'FAIL'} {cid} {msg}")
    if not cond:
        FAIL += 1


def v(n, d):
    x = [0] * n
    for i, c in d.items():
        x[i] = c
    return tuple(x)


def net(n, rx, names=None):
    return Network(n, [(v(n, a), v(n, b)) for a, b in rx], names=names)


# networks ---------------------------------------------------------------------------------------
def Bk(k):
    # species: X=0, Y_a = 1..k, Z = k+1
    n = k + 2
    rx = [({0: 1}, {a: 1}) for a in range(1, k + 1)] + [({a: 1}, {k + 1: 1}) for a in range(1, k + 1)]
    rx.append(({k + 1: 1}, {0: 1}))
    return net(n, rx, ["X"] + [f"Y{a}" for a in range(1, k + 1)] + ["Z"])


EX = {}
EX["R"] = net(2, [({0: 1}, {1: 1}), ({1: 1}, {0: 1})])
EX["T"] = net(3, [({0: 1}, {1: 1}), ({1: 1}, {2: 1}), ({2: 1}, {0: 1})])
EX["W"] = net(3, [({0: 1}, {1: 1}), ({0: 1}, {2: 1}), ({1: 1}, {0: 1}), ({2: 1}, {0: 1})])
EX["U"] = net(4, [({3: 1}, {1: 1, 2: 1}), ({2: 1, 3: 1}, {0: 1}), ({0: 1, 1: 1}, {3: 2})])
EX["D"] = net(3, [({0: 1}, {1: 1}), ({0: 1}, {2: 1}), ({1: 1, 2: 1}, {0: 2})])
EX["F"] = net(6, [({0: 1, 2: 1}, {4: 1}), ({4: 1}, {0: 1, 2: 1}), ({4: 1}, {1: 1, 2: 1}),
                  ({1: 1, 3: 1}, {5: 1}), ({5: 1}, {1: 1, 3: 1}), ({5: 1}, {0: 1, 3: 1})],
              ["S0", "S1", "E", "F", "ES0", "FS1"])
EX["N3"] = net(2, [({0: 2}, {0: 1, 1: 1}), ({0: 1, 1: 1}, {0: 2}), ({0: 1, 1: 1}, {1: 2})])
EX["N4.12"] = net(2, [({0: 2}, {0: 1, 1: 1}), ({0: 1, 1: 1}, {0: 2})])
EX["B2"] = Bk(2)
EX["B3"] = Bk(3)
EX["K3"] = net(3, [({0: 1}, {1: 1}), ({0: 1}, {2: 1}), ({1: 1}, {0: 1}), ({1: 1}, {2: 1}),
                   ({2: 1}, {0: 1}), ({2: 1}, {1: 1})], ["X", "Y", "Z"])
EX["XY4"] = net(2, [({0: 1, 1: 1}, {0: 2}), ({0: 1, 1: 1}, {1: 2}), ({0: 1, 1: 2}, {0: 2, 1: 1}),
                    ({0: 1, 1: 2}, {1: 3})], ["X", "Y"])
EX["Y5"] = net(5, [({0: 1}, {1: 1}), ({1: 1}, {2: 1}), ({1: 1}, {3: 1}), ({1: 1}, {4: 1}),
                   ({2: 1}, {0: 1}), ({3: 1}, {0: 1}), ({4: 1}, {0: 1})], ["A", "B", "C", "D", "E"])
EX["Z4"] = net(2, [({1: 2}, {0: 2}), ({0: 2}, {1: 2}), ({0: 1, 1: 1}, {0: 2}), ({0: 1, 1: 1}, {1: 2})])
EX["B3split"] = pc.split_choice(EX["B3"], {1, 2})
EX["B2dup"] = pc.duplicate(EX["B2"], 4)
EX["N4.12dup"] = pc.duplicate(EX["N4.12"], 0)

RNG = random.Random(20261027)


def all_objects(N):
    P = pc.iosif_P(N)
    Pcs = pc.P_cs(N)
    pi = pc.flat_closure(N)
    M = pc.full_coupling(N)
    lo = pc.pi_plus(N)
    up, wit = pc.cf_upper(N, lo, RNG)
    return dict(P=P, Pcs=Pcs, pi=pi, M=M, lo=lo, up=up, wit=wit)


def main():
    print(f"# check_examples.py; Python {platform.python_version()} on {platform.system()}")
    O = {}
    for name, N in EX.items():
        assert N.is_consistent(), name
        o = all_objects(N)
        O[name] = o
        nm = N.names
        r = N.r
        names_r = [f"R{t+1}" for t in range(r)]
        print(f"## {name}: n={N.n} r={r} s={N.rank()}  P={pc.part_str(o['P'], names_r)}  "
              f"Pcs={pc.part_str(o['Pcs'], names_r)}  pi={pc.part_str(o['pi'], names_r)}  "
              f"pi+={pc.part_str(o['lo'], names_r)}  CFup={pc.part_str(o['up'], names_r)}  "
              f"M*={pc.part_str(o['M'], names_r)}")
        ch = (refines(o["P"], o["Pcs"]) and refines(o["Pcs"], o["pi"]) and refines(o["pi"], o["lo"])
              and refines(o["lo"], o["up"]))
        check("X1-chain", ch, name)
        check("X1-join", o["pi"] == pc.join(o["Pcs"], o["M"], r), name)
    # X2
    check("X2", len(O["B3"]["P"]) == 2 and len(O["B3"]["Pcs"]) == 1 and len(O["B3"]["pi"]) == 1
          and len(O["B3"]["up"]) == 1, "B3: P 2 blocks, Pcs = pi = CF = 1 block")
    check("X2", len(O["B3split"]["P"]) == 1, "B3 split binary: P = 1 block")
    # X3
    o = O["N4.12"]
    check("X3", len(o["P"]) == 2 and len(o["Pcs"]) == 2 and len(o["pi"]) == 1 and len(o["up"]) == 1,
          "N4.12: P = Pcs = singletons; pi = CF = 1 block")
    check("X3", len(O["N4.12dup"]["P"]) == 1, "N4.12 with R1 duplicated: P = 1 block")
    # X4
    o = O["B2dup"]
    check("X4", len(o["pi"]) == 2 and len(o["lo"]) == 1 and len(o["up"]) == 1,
          "B2dup: pi 2 blocks, CF = 1 block (lower bound pi+ = upper bound)")
    check("X4", len(O["B2"]["pi"]) == 1, "B2: pi = 1 block")
    # X5
    o = O["K3"]
    check("X5", len(o["pi"]) == 3 and len(o["up"]) == 1, "K3: pi = 3 blocks, CF upper = 1 block")
    try:
        ok, rows = pc.one_block_certificate(EX["K3"], list(o["pi"]))
        check("X5", ok, f"K3: one-block certificate (minor rows {rows})")
        ok, rows = pc.one_block_certificate(EX["B2dup"], list(O["B2dup"]["pi"]))
        check("X5", ok, f"B2dup: one-block certificate (minor rows {rows})")
    except ImportError:
        print("SKIP X5 certificate (no sympy)")
    # X6
    N = EX["XY4"]
    o = O["XY4"]
    V = pc.random_flux(N, RNG)
    kg = pc.kappa(N, V)
    check("X6", len(kg) == 1, f"XY4: kappa at generic V {V} is one block")
    check("X6", len(o["pi"]) == 2 and o["lo"] == o["up"] and len(o["up"]) == 2,
          "XY4: pi = pi+ = CF = 2 blocks (split witness V_a = V_b, V_d = V_e)")
    Vs = [1, 1, 2, 2]
    check("X6", len(pc.kappa(N, Vs)) == 2, "XY4: kappa at V = (1,1,2,2) has 2 blocks")
    # X9: pi and Pcs are not split-invariant (B2dup, split the choice {R1,R2} of X behind W)
    N2 = pc.split_choice(EX["B2dup"], {0, 1})
    check("X9", len(pc.flat_closure(EX["B2dup"])) == 2 and len(pc.flat_closure(N2)) == 1
          and len(pc.P_cs(EX["B2dup"])) == 2 and len(pc.P_cs(N2)) == 1,
          "B2dup: pi = Pcs = 2 blocks; after splitting {R1,R2}: pi = Pcs = 1 block")
    # Z4: kappa at a special V can be COARSER than CF = kappa_gen (thesis clusters jump at special k)
    N = EX["Z4"]
    o = O["Z4"]
    V = pc.random_flux(N, RNG)
    check("X6", len(pc.kappa(N, V)) == 3 and o["lo"] == o["up"] and len(o["up"]) == 3
          and len(pc.kappa(N, [1, 1, 2, 2])) == 2,
          "Z4: kappa_gen = CF = {R1},{R2},{R3,R4}; kappa at V = (1,1,2,2) = {R1,R2},{R3,R4} (coarser)")
    # X7 invariance
    for name, N in EX.items():
        o = O[name]
        if o["lo"] != o["up"]:
            print(f"SKIP X7 {name}: CF not determined by the bounds")
            continue
        for j in range(N.r):
            N2 = pc.duplicate(N, j)
            lo2 = pc.pi_plus(N2)
            up2, _ = pc.cf_upper(N2, lo2, RNG)
            mp = list(range(N.r)) + [j]
            # pull back: a partition of N2's reactions maps to N's reactions via mp (new j' -> j)
            ok = lo2 == up2 and pc.normalize_partition([{mp[x] for x in B} for B in lo2]) == o["lo"]
            check("X7-dup", ok, f"{name} duplicate R{j+1}")
        for c in pc.source_classes(N):
            c = sorted(c)
            for k in range(1, len(c) + 1):
                for A in __import__("itertools").combinations(c, k):
                    if k == len(c) and k == 1 and False:
                        continue
                    N2 = pc.split_choice(N, A)
                    lo2 = pc.pi_plus(N2)
                    up2, _ = pc.cf_upper(N2, lo2, RNG)
                    newt = N.r
                    mp = list(range(N.r)) + [A[0]]
                    img = pc.normalize_partition([{mp[x] for x in B} for B in lo2])
                    ok = lo2 == up2 and img == o["lo"] and any({newt, A[0]} <= B for B in lo2)
                    check("X7-split", ok, f"{name} split {[f'R{a+1}' for a in A]}")
    # X8 duality criteria
    for name, N in EX.items():
        o = O[name]
        allc = all(any({i, j} <= B for B in o["M"]) for i, j in N.doubling_pairs())
        check("X8-reaction", refines(o["M"], o["lo"]) and ((o["up"] == o["M"]) == allc and (o["lo"] == o["M"]) == allc),
              f"{name}: doubling pairs all fully coupled = {allc}")
        if N.is_conservative():
            D = N.dual()
            loD = pc.pi_plus(D)
            upD, _ = pc.cf_upper(D, loD, RNG)
            MN = N.mips_ray()
            twins = D.doubling_pairs()
            allt = all(any({i, j} <= B for B in MN) for i, j in twins)
            check("X8-species", refines(MN, loD) and ((upD == MN) == allt) and ((loD == MN) == allt),
                  f"{name}: twin pairs all conservation-coupled = {allt}; "
                  f"CF(N*) {pc.part_str(upD, N.names)} vs M(N) {pc.part_str(MN, N.names)}")
    print(f"# runtime {time.time() - T0:.1f} s; {FAIL} FAIL")
    sys.exit(1 if FAIL else 0)


if __name__ == "__main__":
    main()
