#!/usr/bin/env python3
"""Close the onset and continuation gaps for branches 1--3.

Python 3.12 or later, python-flint 0.9, mpmath; a few minutes per branch. Reuses the
K=infinity chart from higher-branches (no higher-branch result is assumed).
Generates its boxes, re-proves every stored line and every adjacent junction,
and links the chart to the all-M certificate. Usage: check_onsets.py J OUTDIR.
"""
import csv
from fractions import Fraction
from pathlib import Path
import sys
from flint import arb, ctx
import allm_core as C

BASE = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(BASE.parent / "higher-branches/scripts"))
import onset_sweep as O

def main():
    j = int(sys.argv[1]); out = Path(sys.argv[2]).resolve()
    out.mkdir(parents=True, exist_ok=True)
    cert = BASE / "certificates" / ("l1h2_j1.csv" if j == 1 else f"l1_j{j}.csv")
    file = out / f"onset_j{j}.csv"
    if "--reuse" not in sys.argv:
        sys.argv = ["onset_sweep.py", str(j), "162", str(file), str(cert)]
        try:
            O.main()
        except SystemExit as exc:
            if exc.code != 0:
                return 1
    ctx.prec = 192
    rows = list(csv.DictReader(file.open()))
    proved = []
    for n, r in enumerate(rows):
        chart = r["chart"]; lo, hi = r["p_lo"], r["p_hi"]
        assert Fraction(lo) < Fraction(hi)
        P = C.eps_interval(lo, hi)
        Xm = [arb(r[k]) for k in ("X0_mid", "u_mid", "v_mid")]
        rad = [arb(r[k]) for k in ("rad_X0", "rad_u", "rad_v")]
        img, X = O.krawczyk(chart, Xm, rad, P)
        assert img is not None and X[2] > 0
        if chart == "s":
            assert O.slope(chart, X, P) > 0
        else:
            assert X[0] > 0 and X[0] < 1
        proved.append((chart, lo, hi, P, X))
        if n == 0:
            assert chart == "s" and Fraction(lo) == 0
            Mj = (15, 44, 87)[j-1]
            assert X[0] > Mj-1 and X[0] < Mj
            print(f"PASS j={j}: M* in {X[0].str(15)}, integer onset {Mj}", flush=True)
        elif chart == rows[n-1]["chart"]:
            assert lo == rows[n-1]["p_hi"]
            _, _, _, endpoint, _ = O.box_for(chart, lo, lo, Xm)
            assert endpoint is not None
            assert all(a.contains(b) for a, b in zip(X, endpoint))
            assert all(a.contains(b) for a, b in zip(proved[n-1][4], endpoint))
        else:
            # A separately certified point at the first M endpoint is in both charts.
            ch, _, _, Ps, Xs = proved[n-1]
            assert ch == "s" and chart == "M"
            _, _, _, ep, _ = O.box_for(chart, lo, lo, Xm)
            assert ep is not None and all(a.contains(b) for a, b in zip(X, ep))
            assert Ps.contains(ep[0]) and Xs[0].contains(arb(lo))
            assert Xs[1].contains(ep[1]) and Xs[2].contains(ep[2])
    assert Fraction(rows[-1]["p_hi"]) == 162
    assert O.junction(str(cert), proved)
    # Link every finite integer certificate used in the theorem to this same arc.
    import gzip
    count = 0
    for r in csv.DictReader(gzip.open(BASE.parent / f"certificates/sweep_j{j}.csv.gz", "rt")):
        M = int(r["M"])
        if not ( (15, 44, 87)[j-1] <= M <= 161 ):
            continue
        X0 = [arb(r[k+"_mid"]) for k in ("kappa", "x", "y")]
        eps = 1 / arb(M).sqrt()
        import bdlib as bd
        # Load result 8's formula only for a new, rigorous point enclosure.
        res = bd.certify_point(X0, M, prec=192)
        assert res is not None
        X = res["box"]
        sx = 1/(1 + eps*X[0]); ux = X[1]/M; vx = X[2]/M
        assert any(ch == "M" and Pb.contains(arb(M)) and B[0].contains(sx)
                   and B[1].contains(ux) and B[2].contains(vx) for ch, _, _, Pb, B in proved)
        count += 1
    assert count == 162 - (15, 44, 87)[j-1]
    print(f"SUMMARY PASS j={j}: {len(rows)} boxes re-proved, all junctions, {count} integer points linked", flush=True)
    return 0

if __name__ == "__main__":
    sys.exit(main())
