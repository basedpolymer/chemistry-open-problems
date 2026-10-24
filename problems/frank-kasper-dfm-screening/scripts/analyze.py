#!/usr/bin/env python3
"""Collect run results and extrapolate to the continuum (mesh size h -> 0).

For a run with mesh levels 0..K (each refinement halves h), F_L is an exact
upper bound of the infimum for that topology/box (piecewise-flat partition).
Assuming F_L = F* + c h^2 + O(h^3), the Richardson estimate is
    R_L = F_L - (F_{L-1} - F_L)/3,
and the uncertainty is estimated as max(|R_K - R_{K-1}|, |R_K - Aitken_K|),
Aitken_K = F_K - d_K/(rho_K - 1), d_L = F_{L-1}-F_L, rho_K = d_{K-1}/d_K.
This error estimate is heuristic (not rigorous).

usage: analyze.py [--runs DIR] [--box fixed|relax|all] [--out CSV]
"""
import os
import sys
import json
import glob
import argparse

HERE = os.path.dirname(os.path.abspath(__file__))


def load_runs(rundir):
    runs = {}
    for f in sorted(glob.glob(os.path.join(rundir, "*", "results.jsonl"))):
        recs = [json.loads(l) for l in open(f) if l.strip()]
        if not recs:
            continue
        bylev = {}
        for r in recs:
            bylev[r["level"]] = r          # last record per level wins (resumed runs)
        meta = {}
        mf = os.path.join(os.path.dirname(f), "meta.json")
        if os.path.exists(mf):
            meta = json.load(open(mf))
        key = (recs[0]["name"], recs[0]["box"])
        runs[key] = {"levels": bylev, "meta": meta}
    return runs


def extrapolate(F):
    """F: dict level->value (consecutive levels from 0). Returns (R, err, rho, K)."""
    levs = sorted(F)
    K = levs[-1]
    if K < 2:
        return None, None, None, K
    d = {L: F[L - 1] - F[L] for L in levs if L - 1 in F}
    R = {L: F[L] - d[L] / 3.0 for L in d}
    err = None
    rho = None
    if K - 1 in d and d[K] != 0:
        rho = d[K - 1] / d[K]
        cands = [abs(R[K] - R[K - 1])]
        if rho > 1.5:
            aitken = F[K] - d[K] / (rho - 1.0)
            cands.append(abs(R[K] - aitken))
        err = max(cands)
    return R[K], err, rho, K


def summarize(runs, box=None):
    rows = []
    for (name, bx), r in runs.items():
        if box not in (None, "all") and bx != box:
            continue
        F = {L: rec["F"] for L, rec in r["levels"].items()}
        R, err, rho, K = extrapolate(F)
        last = r["levels"][K]
        comp = None
        cells = last.get("cells")
        if cells:
            from collections import Counter
            c = Counter(cells["nfaces"])
            comp = "%d/%d/%d/%d" % (c.get(12, 0), c.get(14, 0), c.get(15, 0), c.get(16, 0))
        rows.append({"name": name, "box": bx, "ncell": r["meta"].get("ncell"), "Kmax": K,
                     "F_K": F[K], "F_extrap": R, "err": err, "rho": rho, "A_K": last["A"], "I_K": last["I"],
                     "Vmin": last["Vmin_over_V0"], "Vmax": last["Vmax_over_V0"], "Vrms": last["Vrms_over_V0"],
                     "cone_min": last["cone_min"], "min_edge": last["min_foam_edge_over_R0"],
                     "comp_Z12_14_15_16": comp, "L_unitvol": last["L_unitvol"], "F_levels": [F[L] for L in sorted(F)]})
    rows.sort(key=lambda x: (x["F_extrap"] if x["F_extrap"] is not None else x["F_K"]))
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--runs", default=os.path.join(HERE, "..", "work", "runs"))
    ap.add_argument("--box", default="all")
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    runs = load_runs(a.runs)
    rows = summarize(runs, a.box)
    hdr = ["name", "box", "ncell", "Kmax", "F_K", "F_extrap", "err", "rho", "A_K", "I_K", "Vrms", "Vmin", "Vmax",
           "cone_min", "min_edge", "comp_Z12_14_15_16", "F_levels"]
    lines = [",".join(hdr)]
    for r in rows:
        vals = []
        for h in hdr:
            v = r[h]
            if h == "F_levels":
                vals.append(";".join("%.10f" % x for x in v))
                continue
            if isinstance(v, float):
                vals.append("%.10f" % v if h in ("F_K", "F_extrap", "A_K", "I_K") else "%.3g" % v)
            else:
                vals.append("" if v is None else str(v))
        lines.append(",".join(vals))
    txt = "\n".join(lines) + "\n"
    if a.out:
        open(a.out, "w").write(txt)
    print(txt)


if __name__ == "__main__":
    main()
