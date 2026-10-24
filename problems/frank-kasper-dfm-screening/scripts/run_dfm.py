#!/usr/bin/env python3
"""Resumable driver: relax the diblock-foam-model energy F = (A^2 I)^(1/3) of a
periodic Frank-Kasper partition (fixed Voronoi/contact topology, curved
triangulated faces, free cell volumes) on a sequence of refined meshes.

usage: python3 run_dfm.py NAME [NAME ...] --box fixed|relax --maxlev K [--out DIR]

NAME: BCC A15 C15 sigma Z C14 P delta psigma M mu H (paper boxes), or
      cat:<ID> (84 structures of Dutour Sikiric et al. 2010), cat2:<ID> (part II).
Each level k: refine (k>0), then `cycles` x [restricted L-BFGS (+box strain if
--box relax), tangential smoothing], final restricted L-BFGS.  After every level
a JSON line is appended to DIR/<name>__<box>/results.jsonl and the mesh is
pickled to state.pkl (only the last level is kept), so a killed run resumes.

Environment: Python 3.13, numpy, scipy, numba, gemmi; single thread.
Runtime: seconds (A15, level 4) to ~1 h (56-cell structures, level 5).
"""
import os
import sys
import json
import time
import pickle
import argparse
import platform

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import dfm  # noqa: E402  (pins thread pools to 1)
import numpy as np  # noqa: E402
from dfm import registry, topology, mesh, energy, relax  # noqa: E402

FOURPI = 4.0 * np.pi


def cell_data(M):
    Atot, V, Mo, S = energy.geometry(M)
    I = S - np.sum(Mo ** 2, axis=1) / V
    P = M.X[M.tri] + np.einsum("ij,tkj->tki", M.L, M.tsh)
    ar = 0.5 * np.linalg.norm(np.cross(P[:, 1] - P[:, 0], P[:, 2] - P[:, 0]), axis=1)
    fa = np.bincount(M.tface, weights=ar, minlength=len(M.fi))
    Ai = np.bincount(M.fi, weights=fa, minlength=M.ncell) + np.bincount(M.fj, weights=fa, minlength=M.ncell)
    nf = np.bincount(M.fi, minlength=M.ncell) + np.bincount(M.fj, minlength=M.ncell)
    return V, Ai, I, nf


def initial_mesh(st):
    L, P, tets = st["L"], st["P"], st["tets"]
    faces = topology.build_foam(L, P, tets)
    vp = None
    if st["init"] == "centroid":
        vp = []
        for key in tets:
            X = np.array([L @ (P[k] + np.array(s)) for (k, s) in key])
            vp.append(X.mean(axis=0))
    M = mesh.build_initial_mesh(L, P, tets, faces, vertex_pos=vp)
    return M


def run_one(name, box, maxlev, outdir, cycles=3, maxiter=100000, keep_state=True):
    safe = name.replace(":", "_")
    d = os.path.join(outdir, safe + "__" + box)
    os.makedirs(d, exist_ok=True)
    resf = os.path.join(d, "results.jsonl")
    stf = os.path.join(d, "state.pkl")
    if os.path.exists(stf):
        with open(stf, "rb") as fh:
            M = pickle.load(fh)
        start = M.level + 1
        if start > maxlev:
            print(name, box, "already at level", M.level, flush=True)
            return
    else:
        st = registry.get(name)
        M = initial_mesh(st)
        start = 0
        with open(os.path.join(d, "meta.json"), "w") as fh:
            json.dump({"name": name, "box": box, "ncell": int(M.ncell), "init": st["init"],
                       "info": {k: (v if isinstance(v, (int, float, str, bool, list, tuple)) else str(v))
                                for k, v in st["info"].items()},
                       "L0": st["L"].tolist(), "host": platform.platform(),
                       "python": sys.version.split()[0], "numpy": np.__version__}, fh, indent=1)
    for lev in range(start, maxlev + 1):
        t0 = time.time()
        if lev > 0:
            M = mesh.refine(M)
            M.tsh = M.tsh.astype(np.int8)
        infos = relax.relax_cycles(M, cycles=cycles, maxiter=maxiter, box=(box == "relax"))
        dv = energy.dfm_values(M)
        V, Ai, Ii, nf = cell_data(M)
        q = relax.mesh_quality(M)
        fe = relax.foam_edge_lengths(M)
        cone = relax.cone_check(M)
        V0 = dv["Vbox"] / M.ncell
        Lr = M.L / np.cbrt(abs(np.linalg.det(M.L)))
        rec = {"name": name, "box": box, "level": lev, "ntri": int(len(M.tri)), "nvert": int(len(M.X)),
               "F": dv["F"], "A": dv["A"], "I": dv["I"],
               "Fcycles": [i["F"] for i in infos], "nit": [i["nit"] for i in infos],
               "nfev": [i["nfev"] for i in infos], "msg": [i["msg"] for i in infos],
               "gmax": infos[-1]["gmax"], "cone_min": cone, "min_angle": q["min_angle"],
               "min_foam_edge_over_R0": float(fe.min() / dv["R0"]),
               "Vmin_over_V0": float(V.min() / V0), "Vmax_over_V0": float(V.max() / V0),
               "Vrms_over_V0": float(np.sqrt(np.mean((V / V0 - 1) ** 2))),
               "L_unitvol": Lr.tolist(), "time_s": time.time() - t0}
        if lev == maxlev:
            rec["cells"] = {"V_over_V0": (V / V0).tolist(), "nfaces": nf.tolist(),
                            "A_dimless": (Ai / (FOURPI * (3 * V / FOURPI) ** (2 / 3))).tolist(),
                            "I_dimless": (Ii / (FOURPI * (3 * V / FOURPI) ** (5 / 3) / 5)).tolist()}
        with open(resf, "a") as fh:
            fh.write(json.dumps(rec) + "\n")
        if keep_state:
            tmp = stf + ".tmp"
            M._nbcache = None
            with open(tmp, "wb") as fh:
                pickle.dump(M, fh, protocol=pickle.HIGHEST_PROTOCOL)
            os.replace(tmp, stf)
        print("%-12s %-5s L%d ntri=%8d F=%.10f A=%.8f I=%.8f cone=%.2e minedge/R0=%.3f V/V0=[%.4f,%.4f] t=%.1fs"
              % (name, box, lev, len(M.tri), dv["F"], dv["A"], dv["I"], cone, fe.min() / dv["R0"],
                 V.min() / V0, V.max() / V0, time.time() - t0), flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("names", nargs="+")
    ap.add_argument("--box", default="fixed", choices=["fixed", "relax"])
    ap.add_argument("--maxlev", type=int, default=3)
    ap.add_argument("--cycles", type=int, default=3)
    ap.add_argument("--out", default=os.path.join(HERE, "..", "work", "runs"))
    a = ap.parse_args()
    for nm in a.names:
        try:
            run_one(nm, a.box, a.maxlev, a.out, cycles=a.cycles)
        except Exception as e:  # keep going with the other structures
            import traceback
            traceback.print_exc()
            print("FAILED", nm, repr(e), flush=True)


if __name__ == "__main__":
    main()
