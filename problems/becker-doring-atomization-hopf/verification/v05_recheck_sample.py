#!/usr/bin/env python3
"""v05_recheck_sample.py -- re-proves a SAMPLE of the solver's certificate lines with the solver's own library
(bdlib.certify_point, i.e. the code path of 04_verify_certificates.py), plus checks that 04 does not make.

Selection (per branch j = 1, 2, 3): every M <= 300; NRAND random M in (300, 1e5] (seeded); the boundary values
M = 15, 16 (j=1), 43/44 (j=2), 86/87 (j=3) and 1e5.
For each line: Krawczyk re-proof from the stored midpoints, same kappa (1e-12), sign of l1 in BOTH normalisations
(positive only at (15,1)), simplicity, transversality (as in 04), and in addition
  omega > 0 (Kuznetsov's formula changes sign with omega), Im z > 0 on the whole box (phi non-real, so the root
  of F is not one of the spurious roots 1, 1/A, +-A^{-1/2}), and f(phi) != 0 or f(1/(A phi)) != 0 (the matrix D(phi)
  of the paper's Lemma 5.4 is nonzero; the paper proves this only for real phi).
Usage: python v05_recheck_sample.py SEED NRAND
Environment: WSL Fedora, python-flint 0.9. Runtime: ~20 ms per line.
"""
import sys, os, csv, gzip, random, time
from flint import arb, ctx
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "scripts"))
import bdlib as bd

def load(j):
    out = {}
    for fn in ("sweep_j%d.csv.gz" % j, "sweep_j%d_b.csv.gz" % j):
        for r in csv.DictReader(gzip.open(os.path.join(HERE, "..", "certificates", fn), "rt")):
            out[int(r["M"])] = r
    return out

def main():
    seed, nrand = int(sys.argv[1]), int(sys.argv[2])
    rnd = random.Random(seed)
    nfail = 0; t0 = time.time(); count = 0
    for j, m0 in ((1, 14), (2, 43), (3, 86)):
        rows = load(j)
        sel = set(range(m0, 301)) | set(rnd.sample(range(301, 100001), nrand)) | {m0, m0 + 1, 100000}
        nc = [M for M in sel if rows[M]["status"] == "NO_CROSSING"]
        ok_nc = nc == [m0]
        nfail += not ok_nc
        fails = []
        for M in sorted(sel):
            r = rows[M]
            if r["status"] == "NO_CROSSING":
                continue
            count += 1
            X0 = [arb(r["kappa_mid"]), arb(r["x_mid"]), arb(r["y_mid"])]
            res = bd.certify_point(X0, M, prec=192)
            if res is None:
                fails.append((M, "krawczyk")); continue
            same = abs(float(res["kappa"].mid()) - float(r["kappa_mid"])) < 1e-12*abs(float(r["kappa_mid"]))
            sub = (M, j) == (15, 1)
            sgn = ((res["l1"] > 0) and (res["l1q1"] > 0)) if sub else ((res["l1"] < 0) and (res["l1q1"] < 0))
            extra = bool(res["omega"] > 0) and bool(res["box"][2] > 0)
            # Lemma 5.4 of the paper argues D(phi) != 0 by a real-variable convexity argument; for complex phi we
            # check directly that f(phi) or f(psi) is nonzero (then the paper's proof goes through verbatim)
            ctx.prec = 192
            from flint import acb
            Kb = res["kappa"]/arb(M).sqrt(); Ab = 1 + Kb
            phi = 1 + acb(res["x"], res["y"])/M; psi = 1/(Ab*phi)
            fphi = Ab*phi**(M - 1) + 1 - 1/phi; fpsi = Ab*psi**(M - 1) + 1 - 1/psi
            extra = extra and not (fphi.contains(0) and fpsi.contains(0))
            ok = same and sgn and res["simple"] and bool(res["transv"] > 0) and extra
            if not ok:
                fails.append((M, same, sgn, res["simple"], extra))
        nfail += len(fails)
        print("j=%d: %d lines re-proved (all M<=300, %d random in (300,1e5], boundaries), NO_CROSSING only at M=%s: %s"
              " ; failures: %s" % (j, len([M for M in sel if M not in nc]), nrand, nc, ok_nc, fails or "none"),
              flush=True)
    print("SUMMARY: %d lines, %s, %.0fs" % (count, "PASS" if nfail == 0 else "FAIL", time.time() - t0))
    sys.exit(1 if nfail else 0)

if __name__ == "__main__":
    main()
