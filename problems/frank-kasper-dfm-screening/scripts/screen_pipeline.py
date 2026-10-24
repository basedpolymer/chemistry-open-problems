#!/usr/bin/env python3
"""Adaptive screen (box relaxed): level 2 for every structure, then level 3 for
those whose extrapolated value is within WIN3 of sigma's, then level 4 within
WIN4.  Resumable (run_dfm checkpoints).
usage: screen_pipeline.py LISTFILE [--win3 1e-3] [--win4 3e-4] [--lev4 4]
Runtime: hours on a loaded 4-core machine (single thread)."""
import os
import sys
import argparse

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import dfm  # noqa: E402
from run_dfm import run_one  # noqa: E402
from analyze import load_runs, extrapolate  # noqa: E402

OUT = os.path.join(HERE, "..", "work", "runs")


def estimates(names):
    runs = load_runs(OUT)
    est = {}
    for nm in names:
        r = runs.get((nm, "relax"))
        if not r:
            continue
        F = {L: rec["F"] for L, rec in r["levels"].items()}
        R, err, rho, K = extrapolate(F)
        est[nm] = R if R is not None else F[max(F)]
    return est


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("listfile")
    ap.add_argument("--win3", type=float, default=1e-3)
    ap.add_argument("--win4", type=float, default=3e-4)
    ap.add_argument("--lev4", type=int, default=4)
    a = ap.parse_args()
    names = [l.strip() for l in open(a.listfile) if l.strip()]
    for nm in names:
        try:
            run_one(nm, "relax", 2, OUT)
        except Exception as e:
            print("FAILED", nm, repr(e), flush=True)
    est = estimates(names)
    ref = est["sigma"]
    sel3 = [nm for nm in names if nm in est and est[nm] <= ref + a.win3]
    print("level-3 selection (%d): %s" % (len(sel3), " ".join(sel3)), flush=True)
    for nm in sel3:
        try:
            run_one(nm, "relax", 3, OUT)
        except Exception as e:
            print("FAILED", nm, repr(e), flush=True)
    est = estimates(names)
    ref = est["sigma"]
    sel4 = [nm for nm in names if nm in est and est[nm] <= ref + a.win4]
    print("level-%d selection (%d): %s" % (a.lev4, len(sel4), " ".join(sel4)), flush=True)
    for nm in sel4:
        try:
            run_one(nm, "relax", a.lev4, OUT)
        except Exception as e:
            print("FAILED", nm, repr(e), flush=True)
    print("done", flush=True)


if __name__ == "__main__":
    main()
