#!/usr/bin/env python3
"""Certify exact coverage, branch labels, and all adjacent-box junctions.

Python 3.12 or later / python-flint 0.9. One process, about a minute. At each common
epsilon endpoint a new Krawczyk enclosure must lie inside BOTH stored boxes;
overlap of boxes alone would not prove that their zeros are the same.
"""
import csv
from fractions import Fraction
from pathlib import Path
import sys
from flint import arb, ctx
import allm_core as C

ctx.prec = 192
CERT = Path(__file__).resolve().parent.parent / "certificates"

def box(r):
    return [arb(r[k + "_mid"], r["rad_" + k]) for k in ("kappa", "x", "y")]

def main():
    failures = 0
    for j in (1, 2, 3):
        path = CERT / ("l1h2_j1.csv" if j == 1 else f"l1_j{j}.csv")
        rows = list(csv.DictReader(path.open()))
        cover = (bool(rows) and Fraction(rows[0]["eps_lo"]) == 0
                 and Fraction(rows[-1]["eps_hi"]) ** 2 * 161 >= 1
                 and all(Fraction(r["eps_lo"]) < Fraction(r["eps_hi"]) for r in rows)
                 and all(a["eps_hi"] == b["eps_lo"] for a, b in zip(rows, rows[1:])))
        good = cover
        for n, r in enumerate(rows):
            e = arb(r["eps_lo"])
            Xm = [arb(r[k + "_mid"]) for k in ("kappa", "x", "y")]
            enc, _ = C.certify_hopf(Xm, e)
            ok = enc is not None and all(b.contains(v) for b, v in zip(box(r), enc))
            if n and ok:
                ok = all(b.contains(v) for b, v in zip(box(rows[n - 1]), enc))
            if not n and ok:
                # On Re Lambda = 0 at eps=0, x=0. The imaginary-root equations
                # imply tan(y)=y, cos(y)<0 and kappa^2=sqrt(1+y^2)-1.
                ok = bool(enc[0] > 0 and enc[2] > (2*j - 1)*arb.pi()
                          and enc[2] < (arb(2*j) - arb("0.5"))*arb.pi())
            good &= ok
            if not ok:
                print("FAIL", path.name, "junction", n, r["eps_lo"], flush=True)
        failures += not good
        print("PASS" if good else "FAIL", f"j={j}: exact coverage, limit label, {len(rows)-1} junctions", flush=True)
    print("SUMMARY", "PASS" if failures == 0 else "FAIL", flush=True)
    return failures

if __name__ == "__main__":
    sys.exit(bool(main()))
