#!/usr/bin/env python3
"""07_spectrum_H2.py -- rigorous check of the spectral hypotheses of the Hopf theorem for small M.

Usage: python 07_spectrum_H2.py MMAX sweep_j1.csv [sweep_j2.csv ...]
Environment: WSL Fedora, python-flint 0.9.  Runtime: ~ (M/120)^3 * 4 s per point (about 15 min for MMAX = 120).

For every certified Hopf point (M, j) with M <= MMAX read from the sweep files, the matrix B(K) is built in ball
arithmetic with K = kappa/sqrt(M), kappa the Krawczyk box of the certificate (so B(K) contains the exact Hopf
matrix), and all M eigenvalues are enclosed with Arb's certified eigensolver (acb_mat.eig, which also proves
that the M eigenvalues are simple and isolated).  Claims (PASS/FAIL):
  (a) exactly one enclosure meets the imaginary axis near 0 and it contains 0 (the simple zero eigenvalue
      removed by mass conservation, paper Lemma 3.1);
  (b) exactly two further enclosures meet the imaginary axis, and they contain +-i*omega, omega from the
      certificate (the Hopf pair);
  (c) exactly j-1 conjugate pairs have Re > 0, all other eigenvalues have Re < 0.
(a)-(c) imply: at kappa_j(M) the restriction to the mass hyperplane has exactly one pair of eigenvalues on the
imaginary axis, which is simple (hypotheses of the Hopf theorem besides transversality and l1 != 0), and the
point is the j-th crossing.
"""
import sys, csv, time, gzip
from flint import arb, acb, acb_mat, ctx

def check(M, j, kap, kap_rad, om_lo, om_hi):
    ctx.prec = 128
    kappa = arb(kap, kap_rad)
    K = kappa/arb(M).sqrt(); A = 1 + K; N = M - 1
    rows = [[acb(0)]*M for _ in range(M)]
    rows[0][0] = -A*(M + 2); rows[0][1] = 1 - K
    for c in range(2, N):
        rows[0][c] = -K
    rows[0][M-1] = M*K + 1
    for l in range(1, N):
        rows[l][l-1] += A; rows[l][l] += -A - 1; rows[l][l+1] += 1
    rows[M-1][0] += A; rows[M-1][N-1] += A; rows[M-1][M-1] += -A
    B = acb_mat(rows)
    E = B.eig()                 # certified, raises ValueError if it cannot isolate simple eigenvalues
    on_axis = [e for e in E if e.real.contains(0)]
    pos = [e for e in E if e.real > 0]
    neg = [e for e in E if e.real < 0]
    zero = [e for e in on_axis if e.contains(0)]
    hopf = [e for e in on_axis if not e.contains(0)]
    om = arb(om_lo).union(arb(om_hi))
    om = arb(om.mid(), om.rad() + 1e-12*abs(float(om.mid())))   # csv stores omega in float64
    ok_a = len(zero) == 1
    ok_b = len(hopf) == 2 and all(abs(e.imag).overlaps(om) for e in hopf)
    ok_c = len(pos) == 2*(j - 1) and len(neg) + len(pos) + len(on_axis) == M
    return ok_a and ok_b and ok_c, (len(zero), len(hopf), len(pos), len(neg))

def main():
    Mmax = int(sys.argv[1]); files = sys.argv[2:]
    nfail = 0; n = 0; t0 = time.time()
    for fn in files:
        for r in csv.DictReader((gzip.open if fn.endswith(".gz") else open)(fn, "rt")):
            if r["status"] not in ("PASS", "FAIL") or int(r["M"]) > Mmax:
                continue
            M, j = int(r["M"]), int(r["j"])
            try:
                ok, info = check(M, j, r["kappa_mid"], float(r["box_rad_kappa"]), r["omega_lo"], r["omega_hi"])
            except ValueError as ex:
                ok, info = False, str(ex)
            n += 1; nfail += (not ok)
            print("M=%3d j=%d  (zero, hopf, Re>0, Re<0) = %s  %s" % (M, j, info, "PASS" if ok else "FAIL"), flush=True)
    print("SUMMARY: %d points, %s (%d failures), %.0fs" % (n, "PASS" if nfail == 0 else "FAIL", nfail, time.time() - t0))
    sys.exit(0 if nfail == 0 else 1)

if __name__ == "__main__":
    main()
