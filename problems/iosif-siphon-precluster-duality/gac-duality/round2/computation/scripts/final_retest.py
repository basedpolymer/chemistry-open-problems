"""Final re-test of the latest candidates (Q3-Cgen and K5gen, round2/ec/candidates.md 19:35) on every record file
of both rounds.

Environment: Python 3.12 (WSL Fedora 44), standard library only, exact (records computed with crn.py / ecore.py).
Usage (from round2/computation/): python scripts/final_retest.py
Claims (PASS/FAIL, proved by Q3):
  K5gen  every connected network with (H1), (H2) has s >= #clusters - 1 (Petri clusters = components of the
         species-reaction graph with an edge when the species is a reactant);
  Cgen   connected, (H1), (H2), s = #clusters - 1  =>  neither N nor N* has a critical siphon.
Sides: N, and N* when N* is connected with (H1), (H2) (i.e. N has (H2), (H1)); the dual's clusters are its own.
"""
import os
import sys
import time
import platform

import e2lib


def main():
    t0 = time.time()
    print(f"# final_retest.py  Python {platform.python_version()} on {platform.system()}; exact records")
    files = e2lib.E_FILES + [os.path.join(e2lib.DATA_E2, f) for f in ("sample_5_3_5_8_s3.pkl.gz",
                                                                        "sample_6_3_5_8_s4.pkl.gz")]
    recs = e2lib.load(files)
    print(f"# {len(recs)} distinct networks (first-round files and the round-2 molecularity-3 samples)")
    k5 = [0, []]
    cg = [0, []]
    noec = 0
    for rc in recs:
        n, rx = rc["n"], rc["rx"]
        for side in ("N", "D"):
            X = rc[side]
            if not (X["cons"] and X["consist"]):
                continue
            m, xrx = (n, rx) if side == "N" else (len(rx), e2lib.dual_rx(n, rx))
            if not e2lib.petri_connected(m, xrx):
                continue
            ncl = len(X["CLr"])
            k5[0] += 1
            if rc["s"] < ncl - 1:
                k5[1].append((side, n, rx))
            if rc["s"] == ncl - 1:
                cg[0] += 1
                noec += not e2lib.is_EC(xrx)
                other = rc["D"] if side == "N" else rc["N"]
                if X["CS"] or other["CS"]:
                    cg[1].append((side, n, rx))
    ok = not k5[1] and not cg[1]
    print(f"{'PASS' if not k5[1] else 'FAIL'} K5gen: s >= #clusters - 1 on {k5[0]} connected sides with (H1), (H2) "
          f"({len(k5[1])} failures)")
    for f in k5[1][:3]:
        print(f"     failure side {f[0]}: {e2lib.describe(f[1], f[2])}")
    print(f"{'PASS' if not cg[1] else 'FAIL'} Cgen: no critical siphon in the side nor in its dual on {cg[0]} sides "
          f"with s = #clusters - 1 ({noec} of them not EC) ({len(cg[1])} failures)")
    for f in cg[1][:3]:
        print(f"     failure side {f[0]}: {e2lib.describe(f[1], f[2])}")
    print(f"{'PASS' if ok else 'FAIL'}: all claims")
    print(f"# runtime {time.time() - t0:.1f} s")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
