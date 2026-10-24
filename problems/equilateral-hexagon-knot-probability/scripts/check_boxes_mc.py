#!/usr/bin/env python3
"""NUMERICAL consistency check (not part of the proof) of the certifier hexcert on sample boxes.

For each test box D x (0,pi)^3 it runs hexcert and draws M Monte Carlo points (d in D with the exact
density d1^2 s of the (d1,s,t) chart, theta uniform), evaluates the pierce patterns R and L in floating
point (independent Python/numba implementation of the geometry), and checks
      inR/n^3 - 4 se <= freq_R <= outR/n^3 + 4 se   (same for L),
and that certR = inR, certL = inL (all certified components verified knotted).
Prints PASS/FAIL per box; exit code 1 on any failure.
Environment: Python 3.13, numpy, numba.  Runtime: ~1 min.
Usage: python3 check_boxes_mc.py HEXCERT TRIG_TABLE [M]
"""
import sys, subprocess
import numpy as np, numba as nb
sys.path.insert(0, __file__.rsplit("/", 1)[0])
from mc_validate import hexagon, pierce

@nb.njit(cache=True)
def sample_box(d1lo, d1hi, slo, shi, tlo, thi, M, seed):
    np.random.seed(seed)
    V = np.zeros((6, 3)); cR = 0; cL = 0
    for it in range(M):
        u = np.random.random(); d1 = (d1lo**3 + u*(d1hi**3 - d1lo**3))**(1/3)
        u = np.random.random(); s = np.sqrt(slo**2 + u*(shi**2 - slo**2))
        t = tlo + (thi - tlo)*np.random.random()
        d2 = d1*(1 + s)/2; d3 = d1*((1 - s)/2 + t*s)
        t1 = np.pi*np.random.random(); t2 = np.pi*np.random.random(); t3 = np.pi*np.random.random()
        hexagon(d1, d2, d3, t1, t2, t3, V)
        A12 = pierce(V, 3, 4, 0, 1, 2) != 0; B12 = pierce(V, 0, 1, 2, 3, 4) != 0
        A23 = pierce(V, 5, 0, 2, 3, 4) != 0; B23 = pierce(V, 2, 3, 4, 5, 0) != 0
        A13 = pierce(V, 1, 2, 4, 5, 0) != 0; B13 = pierce(V, 4, 5, 0, 1, 2) != 0
        if A12 and A23 and A13: cR += 1
        if B12 and B23 and B13: cL += 1
    return cR, cL

BOXES = [  # L1 a1 L2 a2 L3 a3 Lth
    (9, 224, 8, 208, 8, 152, 8), (9, 224, 8, 208, 8, 152, 9), (8, 112, 7, 104, 7, 76, 8),
    (7, 56, 6, 52, 6, 38, 7), (9, 180, 8, 240, 8, 100, 8), (9, 150, 8, 250, 8, 200, 8),
    (8, 100, 7, 120, 7, 30, 8), (10, 400, 9, 500, 9, 300, 9), (6, 25, 5, 30, 5, 20, 7),
    (5, 0, 4, 0, 4, 0, 6),
    # boxes touching degenerate faces of the chart (t = 0: flat central triangle; d1 = 2; corner d3 = 0)
    (8, 100, 7, 100, 7, 0, 8), (8, 80, 7, 120, 7, 0, 8), (9, 511, 8, 200, 8, 100, 8), (8, 120, 7, 127, 7, 0, 8),
    (8, 90, 7, 0, 7, 64, 8), (9, 300, 8, 255, 8, 255, 9),
]

def main():
    binp, table = sys.argv[1], sys.argv[2]
    M = int(sys.argv[3]) if len(sys.argv) > 3 else 400000
    inp = "".join("%d %d %d %d %d %d %d %d\n" % ((i,) + b) for i, b in enumerate(BOXES))
    out = subprocess.run([binp, table, "1"], input=inp, capture_output=True, text=True, check=True).stdout.split("\n")
    ok = True
    for i, b in enumerate(BOXES):
        f = list(map(int, out[i].split()[:14]))
        Lth, outR, inR, certR, outL, inL, certL = f[1:8]
        n3 = float(8 ** Lth)
        L1, a1, L2, a2, L3, a3, _ = b
        d1lo, d1hi = 2 * a1 / 2 ** L1, 2 * (a1 + 1) / 2 ** L1
        cR, cL = sample_box(d1lo, d1hi, a2 / 2 ** L2, (a2 + 1) / 2 ** L2, a3 / 2 ** L3, (a3 + 1) / 2 ** L3, M, 1000 + i)
        good = True
        for name, c, lo, hi, cert, inn in (("R", cR, inR, outR, certR, inR), ("L", cL, inL, outL, certL, inL)):
            p = c / M; se = np.sqrt(max(p * (1 - p), 1.0 / M) / M)
            g = (lo / n3 - 4 * se <= p <= hi / n3 + 4 * se) and cert == inn
            good = good and g
            print("  box %d %s: certified [%.3e, %.3e]  MC %.3e +- %.1e  cert==in %s" % (i, name, lo / n3, hi / n3, p, se, cert == inn))
        print(("PASS" if good else "FAIL"), "box", i, b)
        ok = ok and good
    sys.exit(0 if ok else 1)

if __name__ == "__main__":
    main()
