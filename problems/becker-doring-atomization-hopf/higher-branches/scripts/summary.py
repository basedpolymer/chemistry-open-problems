#!/usr/bin/env python3
"""summary.py -- the results table, read from the certificate files of every branch present in ../certificates/.

Usage:  python summary.py
Prints, for each branch j: the onset enclosure (first box of onset_j.csv), M_j, the eps-chart range, numbers of
boxes / groups / integer points, the ranges of kappa and of Lhat on the eps-chart, the extreme l1 values at the
integer points, and every integer M where l1 > 0.  No claim is made here (the claims are in the generating scripts
and their logs); the script only reports and checks that all status columns are PASS / SUBCRITICAL.
"""
import os, csv, glob, gzip, math, sys
from fractions import Fraction

D = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "certificates")

def rd(fn):
    op = gzip.open if fn.endswith(".gz") else open
    with op(fn, "rt") as fh:
        return list(csv.DictReader(fh))

def find(name):
    for ext in ("", ".gz"):
        p = os.path.join(D, name + ext)
        if os.path.exists(p):
            return p
    return None

def main():
    js = sorted({int(os.path.basename(p).split("_j")[1].split(".")[0]) for p in glob.glob(os.path.join(D, "l1_j*.csv*"))})
    allok = bool(js)
    print("| j | M*_j (enclosure) | M_j | M_j - M*_j | eps-chart [0, eps_c] | boxes / H2 groups | onset boxes (s/M) |"
          " integer M | kappa on eps-chart | Lhat on eps-chart | l1 > 0 at |")
    print("|---|---|---|---|---|---|---|---|---|---|---|")
    for j in js:
        L1 = rd(find("l1_j%d.csv" % j)); H2 = rd(find("h2_j%d.csv" % j)) if find("h2_j%d.csv" % j) else []
        ON = rd(find("onset_j%d.csv" % j)) if find("onset_j%d.csv" % j) else []
        PW = rd(find("pw_j%d.csv" % j)) if find("pw_j%d.csv" % j) else []
        ok = all(r["status"] == "PASS" for r in L1) and all(r["status"] == "PASS" for r in H2) and \
            all(r["status"] == "PASS" for r in ON) and all(r["status"] in ("PASS", "SUBCRITICAL") for r in PW)
        # Empty/missing files must not pass by vacuous truth. Check that the
        # pointwise and epsilon ranges together cover every integer M.
        def contiguous(rows, lo, hi):
            return bool(rows) and all(Fraction(a[hi]) == Fraction(b[lo])
                                      for a, b in zip(rows, rows[1:]))
        complete = all((L1, H2, ON, PW))
        if complete:
            ms_int = [int(r['M']) for r in PW]
            eps_c = Fraction(L1[-1]['eps_hi'])
            complete = (contiguous(L1, 'eps_lo', 'eps_hi')
                        and contiguous(H2, 'eps_lo', 'eps_hi')
                        and Fraction(L1[0]['eps_lo']) == Fraction(H2[0]['eps_lo']) == 0
                        and Fraction(H2[-1]['eps_hi']) == eps_c
                        and ms_int == list(range(ms_int[0], ms_int[-1] + 1))
                        and (ms_int[-1] + 1) * eps_c**2 >= 1
                        and all(int(r['NP']) == 2*j+2 and int(r['NPm']) == 4 for r in H2))
        ok &= complete
        allok &= ok
        if not ok:
            print('FAIL branch %d: missing, failed or incomplete certificate data' % j)
        kap = [float(r["kappa_mid"]) for r in L1]
        Llo = min(float(r["Lhat_lo"]) for r in L1); Lhi = max(float(r["Lhat_hi"]) for r in L1)
        ms = "%.10f, %.10f" % (float(ON[0]["M_lo"]), float(ON[0]["M_hi"])) if ON else "-"
        Mj = int(PW[0]["M"]) if PW else None
        d = Mj - 0.5 * (float(ON[0]["M_lo"]) + float(ON[0]["M_hi"])) if (ON and PW) else float("nan")
        nA = sum(1 for r in ON if r["chart"] == "s"); nB = len(ON) - nA
        sub = [r["M"] for r in PW if r["status"] == "SUBCRITICAL"]
        print("| %d | [%s] | %s | %.4f | [0, %s] | %d / %d | %d / %d | %s..%s (%d) | [%.4f, %.3f] | [%.3e, %.3e] | %s |" % (
            j, ms, Mj, d, L1[-1]["eps_hi"], len(L1), len(H2), nA, nB, PW[0]["M"] if PW else "-",
            PW[-1]["M"] if PW else "-", len(PW), min(kap), max(kap), Llo, Lhi, ", ".join(sub) if sub else "none"))
        if PW:
            r = PW[0]
            print("|   | l1 at M_j = %s: unit [%s, %s], q1 [%s, %s]; kappa_j(M_j) = %s | | | | | | | | | |" % (
                r["M"], r["l1_unit_lo"], r["l1_unit_hi"], r["l1_q1_lo"], r["l1_q1_hi"], r["kappa_mid"][:14]))
    print("PASS: complete certificate data; all status columns PASS/SUBCRITICAL" if allok
          else "FAIL: certificate data incomplete or failed")
    return 0 if allok else 1

if __name__ == "__main__":
    sys.exit(main())
