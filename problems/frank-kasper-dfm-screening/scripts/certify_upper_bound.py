#!/usr/bin/env python3
"""Exact-arithmetic certificate of an upper bound on the DFM infimum of one
periodic cell topology in a given box.

Input: a relaxed mesh (state.pkl written by run_dfm.py).  Steps:
 1. Round every vertex position and the box matrix to the dyadic grid 2^-K
    (K = 44): this defines an explicit piecewise-flat partition Pi of the torus
    R^3 / L Z^3 with the same cell adjacency (topology) as the foam.
 2. Embeddedness (exact integers): for each cell, choose an apex c_i (its
    rounded floating-point centroid) and check that every tetrahedron
    (c_i, a, b, c) over the oriented boundary triangles of the cell has
    strictly positive volume.  Together with sum_i V_i = det L (exact), the
    cones then tile the torus with degree one, so the cells have disjoint
    interiors and cover the torus: Pi is a genuine partition.
 3. Exact rationals: V_i, first moments M_i, second moments S_i (about the
    apex) of the polyhedral cells; I_i = S_i - |M_i|^2 / V_i.
 4. Area: each triangle area sqrt(|n|^2)/2 is bounded above with an integer
    ceiling square root.
 5. F^3 = 20 A_tot^2 (sum_i I_i) / (27 V_box^3)  (A_tot = total interface area,
    R0 and pi cancel), bounded above exactly; F_upper = a decimal number whose
    cube is checked (exactly) to exceed the bound.
Output: PASS if the exact checks hold, the certified value F_upper, and the
comparison with an optional threshold (e.g. a published value):
   F_upper < threshold  certifies  inf F (this topology, this box) < threshold.

usage: certify_upper_bound.py STATE.pkl [--threshold X] [--digits 9]
Environment: Python 3 integers/fractions only (plus numpy to load the mesh).
Runtime: ~1 min per 10^5 triangles.
"""
import os
import sys
import math
import time
import pickle
import argparse
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import dfm  # noqa: E402
import numpy as np  # noqa: E402
from dfm import energy  # noqa: E402

K = 44
SCALE = 1 << K


def to_int(x):
    return int(round(float(x) * SCALE))


def mesh_to_cert(M):
    """Integer certificate data (positions on the 2^-K grid) from a mesh."""
    Atot_f, V_f, Mo_f, S_f = energy.geometry(M)
    cen = M.ref + Mo_f / V_f[:, None]
    r = lambda a: np.round(np.asarray(a, dtype=float) * SCALE).astype(np.int64)
    return {"K": np.int64(K), "L": r(M.L), "X": r(M.X), "C": r(cen), "tri": M.tri.astype(np.int32),
            "tsh": M.tsh.astype(np.int8), "tface": M.tface.astype(np.int32), "fi": M.fi.astype(np.int32),
            "fj": M.fj.astype(np.int32), "fsj": M.fsj.astype(np.int8), "ncell": np.int64(M.ncell)}


class _CertMesh:
    pass


def certify(M, threshold=None, digits=9, verbose=True):
    """M: a Mesh (rounded here) or a dict of integer certificate arrays."""
    t0 = time.time()
    if not isinstance(M, dict):
        M = mesh_to_cert(M)
    C = M
    assert int(C["K"]) == K
    M = _CertMesh()
    M.tri, M.tsh, M.tface, M.fi, M.fj, M.fsj = C["tri"], C["tsh"], C["tface"], C["fi"], C["fj"], C["fsj"]
    n = int(C["ncell"])
    Li = [[int(C["L"][r, c]) for c in range(3)] for r in range(3)]
    Xi = [[int(v) for v in row] for row in C["X"]]
    Ci = [[int(v) for v in row] for row in C["C"]]

    def shiftvec(s):
        return [Li[0][0] * s[0] + Li[0][1] * s[1] + Li[0][2] * s[2],
                Li[1][0] * s[0] + Li[1][1] * s[1] + Li[1][2] * s[2],
                Li[2][0] * s[0] + Li[2][1] * s[1] + Li[2][2] * s[2]]

    fi, fj, fsj, tface = M.fi, M.fj, M.fsj, M.tface
    face_off = [shiftvec([int(x) for x in fsj[f]]) for f in range(len(fi))]
    V6 = [0] * n            # 6 V_i * SCALE^3
    M24 = [[0, 0, 0] for _ in range(n)]   # 24 M_i * SCALE^4
    S60 = [0] * n           # 60 S_i * SCALE^5
    area4_up = 0            # sum of ceil(sqrt(|n|^2)) = upper bound of 2*area*SCALE^2 ... per triangle
    min_cone = None
    tri, tsh = M.tri, M.tsh
    for t in range(len(tri)):
        P = []
        for k in range(3):
            v = int(tri[t, k])
            sv = shiftvec([int(x) for x in tsh[t, k]])
            P.append((Xi[v][0] + sv[0], Xi[v][1] + sv[1], Xi[v][2] + sv[2]))
        e1 = (P[1][0] - P[0][0], P[1][1] - P[0][1], P[1][2] - P[0][2])
        e2 = (P[2][0] - P[0][0], P[2][1] - P[0][1], P[2][2] - P[0][2])
        nx = e1[1] * e2[2] - e1[2] * e2[1]
        ny = e1[2] * e2[0] - e1[0] * e2[2]
        nz = e1[0] * e2[1] - e1[1] * e2[0]
        nn = nx * nx + ny * ny + nz * nz
        r = math.isqrt(nn)
        if r * r < nn:
            r += 1
        area4_up += r          # |n| <= r ; area = |n|/2 / SCALE^2
        f = int(tface[t])
        for side in (0, 1):
            if side == 0:
                ci = int(fi[f]); o = Ci[ci]; sg = 1
            else:
                ci = int(fj[f]); off = face_off[f]
                o = (Ci[ci][0] + off[0], Ci[ci][1] + off[1], Ci[ci][2] + off[2]); sg = -1
            a = (P[0][0] - o[0], P[0][1] - o[1], P[0][2] - o[2])
            b = (P[1][0] - o[0], P[1][1] - o[1], P[1][2] - o[2])
            c = (P[2][0] - o[0], P[2][1] - o[1], P[2][2] - o[2])
            det = (a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0])
                   + a[2] * (b[0] * c[1] - b[1] * c[0])) * sg
            if min_cone is None or det < min_cone:
                min_cone = det
            V6[ci] += det
            sx = a[0] + b[0] + c[0]; sy = a[1] + b[1] + c[1]; sz = a[2] + b[2] + c[2]
            M24[ci][0] += det * sx; M24[ci][1] += det * sy; M24[ci][2] += det * sz
            q = (a[0] * a[0] + a[1] * a[1] + a[2] * a[2] + b[0] * b[0] + b[1] * b[1] + b[2] * b[2]
                 + c[0] * c[0] + c[1] * c[1] + c[2] * c[2] + a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
                 + b[0] * c[0] + b[1] * c[1] + b[2] * c[2] + c[0] * a[0] + c[1] * a[1] + c[2] * a[2])
            S60[ci] += det * q
    # combinatorial check: each cell boundary is a closed oriented surface
    # (every directed edge of its outward-oriented triangles is matched by the reverse edge)
    from collections import Counter
    dedges = [Counter() for _ in range(n)]
    for t in range(len(tri)):
        f = int(tface[t])
        for side in (0, 1):
            ci = int(fi[f]) if side == 0 else int(fj[f])
            fr = (0, 0, 0) if side == 0 else tuple(int(x) for x in fsj[f])
            vs = [(int(tri[t, k]), tuple(int(tsh[t, k, c]) - fr[c] for c in range(3))) for k in range(3)]
            if side == 1:
                vs = vs[::-1]
            for k in range(3):
                dedges[ci][(vs[k], vs[(k + 1) % 3])] += 1
    closed_ok = True
    for ci in range(n):
        for (a, b), m in dedges[ci].items():
            if m != 1 or dedges[ci].get((b, a), 0) != 1:
                closed_ok = False
                break
    detL = (Li[0][0] * (Li[1][1] * Li[2][2] - Li[1][2] * Li[2][1])
            - Li[0][1] * (Li[1][0] * Li[2][2] - Li[1][2] * Li[2][0])
            + Li[0][2] * (Li[1][0] * Li[2][1] - Li[1][1] * Li[2][0]))
    vol_ok = (sum(V6) == 6 * abs(detL))
    cones_ok = min_cone > 0
    # exact rationals (in units of the scaled grid)
    sumI = Fraction(0)
    for i in range(n):
        Vi = Fraction(V6[i], 6)
        Mi = [Fraction(m, 24) for m in M24[i]]
        Si = Fraction(S60[i], 60)
        sumI += Si - (Mi[0] ** 2 + Mi[1] ** 2 + Mi[2] ** 2) / Vi
    Atot_up = Fraction(area4_up, 2)
    Vbox = Fraction(abs(detL))
    F3_up = Fraction(20) * Atot_up ** 2 * sumI / (27 * Vbox ** 3)   # scale-invariant
    # decimal upper bound of the cube root
    fl = float(F3_up) ** (1.0 / 3.0)
    q = 10 ** digits
    Fu = Fraction(math.ceil(fl * q), q)
    while Fu ** 3 < F3_up:
        Fu += Fraction(1, q)
    # also a lower bound of F for THIS polyhedral partition is not claimed
    res = {"ntri": len(tri), "cones_positive": cones_ok, "volumes_sum_ok": vol_ok,
           "F_upper": Fu, "F_float": float(F3_up) ** (1 / 3), "time": time.time() - t0}
    ok = cones_ok and vol_ok and closed_ok
    res["closed_surfaces"] = closed_ok
    if verbose:
        print("triangles: %d   exact checks: cell boundaries closed %s, cones>0 %s, sum V_i = det L %s   (%.1f s)"
              % (len(tri), closed_ok, cones_ok, vol_ok, res["time"]))
        print("certified: inf F <= F(Pi) <= %s" % (str(float(Fu)) if False else ("%.*f" % (digits, Fu))))
        print("PASS" if ok else "FAIL", "partition embedded and upper bound certified")
        if threshold is not None:
            thr = Fraction(str(threshold))
            print(("PASS" if (ok and Fu < thr) else "FAIL"),
                  "certified upper bound %.*f < %s" % (digits, Fu, threshold))
    return ok, res


def load(path):
    if path.endswith(".npz"):
        z = np.load(path)
        return {k: z[k] for k in z.files}
    return pickle.load(open(path, "rb"))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("state", help="state.pkl of a run, or an integer certificate .npz")
    ap.add_argument("--threshold", type=str, default=None)
    ap.add_argument("--digits", type=int, default=9)
    ap.add_argument("--export", default=None, help="write the integer certificate (.npz) and exit after checking")
    a = ap.parse_args()
    M = load(a.state)
    if a.export:
        cert = mesh_to_cert(M)
        np.savez_compressed(a.export, **cert)
        M = load(a.export)
    ok, res = certify(M, threshold=a.threshold, digits=a.digits)
    thr_ok = True
    if a.threshold is not None:
        thr_ok = res["F_upper"] < Fraction(a.threshold)
    sys.exit(0 if (ok and thr_ok) else 1)
