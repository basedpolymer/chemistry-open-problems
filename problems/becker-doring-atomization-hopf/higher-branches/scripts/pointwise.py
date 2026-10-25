#!/usr/bin/env python3
"""pointwise.py -- certificate at every integer M of a range, branch j (the range between the onset of branch j
and the start of the eps-box sweep, where kappa_j(M) is large).

Usage:
  python pointwise.py J MMIN MMAX OUT.csv R0 Y cL cR [PREC [cL2 cR2]]   generate (seeds by continuation, then certify)
  python pointwise.py --verify IN.csv [PREC]                     re-prove every line from the stored midpoints
Environment: Python 3.13, python-flint 0.9 (Arb), mpmath.  One thread.  Runtime about 5-10 s per M.

Claims per integer M (all in ball arithmetic):
  (a) result 8's Krawczyk test (bdlib.krawczyk, integer M, unknowns kappa = K sqrt(M), z = M(phi - 1)) proves a
      unique zero of the Hopf system  K^-3 F(phi) = 0, Re(M lambda/K) = 0  in a box; on it 0 is not in dQ/dz
      (simple zero of F), Re dlambda/dkappa > 0 (transversality);
  (b) the first Lyapunov coefficient l1 (restricted to the mass hyperplane; result 8's closed form bdlib.lyapunov_cf,
      normalizations q^H q = 1 and q_1 = 1) is enclosed; its sign is certified (negative: supercritical;
      positive: subcritical);
  (c) the eps-regular Krawczyk test of the all-M extension (allm_core.krawczyk) at the ball eps = M^(-1/2) proves a
      unique zero in a box that contains the box of (a) (same zero); on it Im z > 0 and |A E + eps^2 z| > 0
      (f(phi) != 0, needed in Pego-Velazquez Lemma 5.4);
  (d) the spectral certificate of h2j_lib (all-M Theorem 4.5) at kappa_j(M): the only eigenvalues of B on iR are
      0 and +-i omega, and exactly 2j - 2 zeros of F (with multiplicity) give eigenvalues with Re lambda > 0
      (counts N(P) = 2j + 2, N(P-) = 4).
Status per line: PASS (all checks, l1 < 0), SUBCRITICAL (all checks, l1 > 0 certified), FAIL otherwise.
Exit code 0 iff no line is FAIL and every integer M of [MMIN, MMAX] has a line.
"""
import sys, os, csv, time, math
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from flint import arb, acb, ctx, arb_mat
import mpmath as mp
import allm_core as C
import h2j_lib as HJ
import bdlib as bd

FIELDS = ["M", "j", "kappa_mid", "x_mid", "y_mid", "rad_kappa", "rad_x", "rad_y", "omega_lo", "omega_hi",
          "l1_unit_lo", "l1_unit_hi", "l1_q1_lo", "l1_q1_hi", "simple", "transv_lo", "allm_box_ok", "imz_pos",
          "fphi_lo", "R0", "Y", "cL", "cR", "loc_mode", "mono_ub", "NP", "NPm", "status"]

ALT = []          # optional second pair of level values (cL2, cR2)

def limit_point(j):
    mp.mp.dps = 30
    c = (2 * j - mp.mpf(1) / 2) * mp.pi
    t = mp.findroot(lambda u: mp.tan(u) - u, c - 1 / c)
    return [arb(str(mp.sqrt(mp.sqrt(1 + t * t) - 1))), arb(0), arb(str(t))]

# ------------------------------------------------------------------------------------------------ seeds (float)
def _F3(kap, Y):
    eps, x, y = Y
    Q, Lam = C.QL(kap, acb(x, y), eps)
    f = [Q.v.real, Q.v.imag, Lam.v.real]
    Jm = [[Q.de.real, Q.dz.real, -Q.dz.imag], [Q.de.imag, Q.dz.imag, Q.dz.real],
          [Lam.de.real, Lam.dz.real, -Lam.dz.imag]]
    return f, Jm

def _newton3(kap, Y):
    """Newton in (eps, x, y) at fixed kappa (float midpoints)"""
    Y = [arb(arb(v).mid()) for v in Y]
    for _ in range(60):
        f, Jm = _F3(kap, Y)
        A = arb_mat([[c.mid() for c in r] for r in Jm]); d = A.solve(arb_mat([[c.mid()] for c in f]))
        Y = [arb((Y[i] - d[i, 0]).mid()) for i in range(3)]
        if max(abs(float(d[i, 0].mid())) for i in range(3)) < 1e-40:
            break
    return Y

def seeds(j, Mmin, Mmax):
    """starting guesses for the j-th Hopf point at every integer M in [Mmin, Mmax]: continuation of the Hopf curve
    from eps = 0 (the limit root (kappa_j^0, i t_j)) in eps up to kappa ~ 30, then in kappa (unknowns eps, x, y) up
    to kappa = 1e9; the integer M are read off by interpolation in log kappa.  Exploration only (floating point)."""
    X = C.newton(limit_point(j), arb(0)); e = 0.0; h = 0.002; prev = None; hprev = h
    eps_stop = 1 / math.sqrt(Mmax + 0.5)
    while True:
        en = e + h
        Xp = X if prev is None else [arb(float(X[i].mid()) + (float(X[i].mid()) - float(prev[i].mid())) * h / hprev)
                                     for i in range(3)]
        Xn = C.newton(Xp, arb(en))
        f, _ = C.fvec_jac(Xn, arb(en))
        if not (max(abs(float(v.mid())) for v in f) < 1e-20) or abs(float(Xn[0].mid()) / float(X[0].mid()) - 1) > 0.05:
            h /= 2
            continue
        prev = X; hprev = h; X = Xn; e = en; h = min(h * 1.5, 0.004)
        if e > eps_stop or float(X[0].mid()) > 30:
            break
    hist = []
    kap = float(X[0].mid()); Yv = [arb(e), X[1], X[2]]
    while kap < 1e9:
        Yv = _newton3(arb(kap), Yv)
        hist.append((kap, 1 / float(Yv[0].mid()) ** 2, Yv))
        if hist[-1][1] < Mmin - 0.5 and len(hist) > 3:
            break
        kap *= 1.03
        if len(hist) > 2:
            (k1, _, Y1), (k2, _, Y2) = hist[-2], hist[-1]
            r = math.log(kap / k2) / math.log(k2 / k1)
            Yv = [arb(float(Y2[i].mid()) + (float(Y2[i].mid()) - float(Y1[i].mid())) * r) for i in range(3)]
    out = {}
    for M in range(Mmin, Mmax + 1):
        # bracketing points in M (M decreases along hist), linear interpolation in log kappa
        for (k1, M1, Y1), (k2, M2, Y2) in zip(hist[:-1], hist[1:]):
            if (M1 - M) * (M2 - M) <= 0:
                w = (M1 - M) / (M1 - M2) if M1 != M2 else 0.0
                lk = math.log(k1) + w * (math.log(k2) - math.log(k1))
                out[M] = [arb(math.exp(lk)), arb(float(Y1[1].mid()) * (1 - w) + float(Y2[1].mid()) * w),
                          arb(float(Y1[2].mid()) * (1 - w) + float(Y2[2].mid()) * w)]
                break
    return out

# ------------------------------------------------------------------------------------------------ certificate
def certify_M(j, M, X0, R0, Y, cL, cR, prec, rad0=None):
    """all checks (a)-(d) at integer M from the starting point X0; returns (row dict, status)"""
    ctx.prec = prec
    row = dict(M=M, j=j, R0=R0, Y=Y, cL=cL, cR=cR)
    if rad0 is None:
        res = bd.certify_point(X0, M, prec=prec)
    else:
        # verification mode: Krawczyk on the stored box (midpoint X0, radii rad0), no Newton
        enc, box = bd.krawczyk([arb(x) for x in X0], M, [arb(r) for r in rad0])
        res = None
        if enc is not None:
            kap, x, y = enc
            K = kap / arb(M).sqrt(); phi = 1 + acb(x, y) / M
            l1 = bd.lyapunov_cf(M, K, phi, normalize="unit")
            l1q = bd.lyapunov_cf(M, K, phi, normalize="q1")
            simple, tr, om = bd.hopf_checks(enc, M)
            res = dict(M=M, kappa=kap, x=x, y=y, l1=l1, l1q1=l1q, simple=simple, transv=tr, omega=om, box=box)
    if res is None:
        row["status"] = "FAIL"; row["fail"] = "krawczyk(bdlib)"; return row
    box8 = res["box"]
    Xm = [arb(res["kappa"].mid()), arb(res["x"].mid()), arb(res["y"].mid())] if rad0 is None else [arb(x) for x in X0]
    row.update(kappa_mid=Xm[0].str(40, radius=False), x_mid=Xm[1].str(40, radius=False), y_mid=Xm[2].str(40, radius=False),
               rad_kappa="%.6e" % float(box8[0].rad()), rad_x="%.6e" % float(box8[1].rad()),
               rad_y="%.6e" % float(box8[2].rad()),
               omega_lo="%.17g" % float(res["omega"].lower()), omega_hi="%.17g" % float(res["omega"].upper()),
               l1_unit_lo="%.17g" % float(res["l1"].lower()), l1_unit_hi="%.17g" % float(res["l1"].upper()),
               l1_q1_lo="%.17g" % float(res["l1q1"].lower()), l1_q1_hi="%.17g" % float(res["l1q1"].upper()),
               simple=int(res["simple"]), transv_lo="%.6g" % float(res["transv"].lower()))
    # (c) eps-regular Krawczyk at eps = M^-1/2 around the same midpoint; its box must contain box8
    eb = arb(1) / arb(M).sqrt()
    ok_c = False
    box = None
    for r in (1e-20, 1e-14, 1e-10, 1e-7):
        rad = [arb(max(abs(float(Xm[i].mid())) * r, float(box8[i].rad()) * 4) + 1e-30) for i in range(3)]
        enc, box = C.krawczyk(Xm, rad, eb)
        if enc is not None and all(box[i].contains(box8[i]) for i in range(3)):
            ok_c = True
            break
    row["allm_box_ok"] = int(ok_c)
    if not ok_c:
        row["status"] = "FAIL"; row["fail"] = "krawczyk(allm)"; return row
    simple2, tr2, fphi, laim = C.hopf_checks(box, eb)
    imz = bool(arb(box[2].lower()) > 0)
    row.update(imz_pos=int(imz), fphi_lo="%.6g" % float(fphi.abs_lower()))
    # (d) spectral certificate
    for (cl, cr) in [(cL, cR)] + list(ALT):
        h = HJ.h2j_checks(j, eb, box, R0, Y, cl, cr)
        if h.get("ok", False):
            break
    row.update(cL=h.get("cL", cL), cR=h.get("cR", cR))
    row.update(loc_mode=h.get("loc_mode", ""), mono_ub="%.4g" % h["mono_ub"] if "mono_ub" in h else "",
               NP=h.get("NP", ""), NPm=h.get("NPm", ""))
    checks = (res["simple"] and bool(res["transv"] > 0) and simple2 and bool(tr2 > 0) and imz
              and bool(fphi.abs_lower() > 0) and bool(res["omega"] > 0) and h.get("ok", False))
    if not checks:
        row["status"] = "FAIL"; row["fail"] = h.get("fail", "checks"); return row
    neg = bool(res["l1"] < 0) and bool(res["l1q1"] < 0)
    pos = bool(res["l1"] > 0) and bool(res["l1q1"] > 0)
    row["status"] = "PASS" if neg else ("SUBCRITICAL" if pos else "FAIL")
    if row["status"] == "FAIL":
        row["fail"] = "l1 sign"
    return row

def main_generate():
    j = int(sys.argv[1]); Mmin = int(sys.argv[2]); Mmax = int(sys.argv[3]); out = sys.argv[4]
    R0, Y, cL, cR = (float(a) for a in sys.argv[5:9])
    prec = int(sys.argv[9]) if len(sys.argv) > 9 else 256
    if len(sys.argv) > 11:
        ALT.append((float(sys.argv[10]), float(sys.argv[11])))
    ctx.prec = 192
    t0 = time.time()
    S = seeds(j, Mmin, Mmax)
    print("seeds for %d values of M (%.0fs)" % (len(S), time.time() - t0), flush=True)
    nf = 0; nsub = 0
    with open(out, "w", newline="") as fh:
        wr = csv.DictWriter(fh, fieldnames=FIELDS, extrasaction="ignore"); wr.writeheader()
        for M in range(Mmin, Mmax + 1):
            if M not in S:
                nf += 1; print("M=%d: no seed" % M, flush=True); continue
            row = certify_M(j, M, S[M], R0, Y, cL, cR, prec)
            wr.writerow(row); fh.flush()
            nf += row["status"] == "FAIL"; nsub += row["status"] == "SUBCRITICAL"
            if row["status"] != "PASS" or M % 10 == 0 or M == Mmin:
                print("M=%d j=%d kappa=%s l1(unit) in [%s, %s] NP=%s NPm=%s %s %s (%.0fs)" % (
                    M, j, row.get("kappa_mid", "")[:12], row.get("l1_unit_lo", ""), row.get("l1_unit_hi", ""),
                    row.get("NP", ""), row.get("NPm", ""), row["status"], row.get("fail", ""), time.time() - t0), flush=True)
    print("SUMMARY j=%d M in [%d, %d]: %d lines, %d SUBCRITICAL, %d FAIL -> %s (%.0fs)" % (
        j, Mmin, Mmax, Mmax - Mmin + 1, nsub, nf, "PASS" if nf == 0 else "FAIL", time.time() - t0))
    sys.exit(0 if nf == 0 else 1)

def main_verify():
    fin = sys.argv[2]
    prec = int(sys.argv[3]) if len(sys.argv) > 3 else 256
    rows = list(csv.DictReader(open(fin)))
    t0 = time.time(); nf = 0; nsub = 0
    Ms = [int(r["M"]) for r in rows]
    cover = Ms == list(range(Ms[0], Ms[-1] + 1))
    for r in rows:
        M = int(r["M"]); j = int(r["j"])
        X0 = [r["kappa_mid"], r["x_mid"], r["y_mid"]]
        rad0 = [r["rad_kappa"], r["rad_x"], r["rad_y"]]
        row = certify_M(j, M, X0, float(r["R0"]), float(r["Y"]), float(r["cL"]), float(r["cR"]), prec, rad0=rad0)
        same = row["status"] == r["status"]
        bad = row["status"] == "FAIL" or not same
        nf += bad; nsub += row["status"] == "SUBCRITICAL"
        if bad or row["status"] != "PASS" or M % 25 == 0:
            print("M=%d j=%d stored %s, re-proved %s %s (%.0fs)" % (M, j, r["status"], row["status"], row.get("fail", ""),
                                                               time.time() - t0), flush=True)
    okall = nf == 0 and cover
    print("SUMMARY %s: %d lines re-proved (M = %d..%d contiguous: %s), %d SUBCRITICAL, %d FAIL/mismatch -> %s (%.0fs)" % (
        os.path.basename(fin), len(rows), Ms[0], Ms[-1], cover, nsub, nf, "PASS" if okall else "FAIL", time.time() - t0))
    sys.exit(0 if okall else 1)

if __name__ == "__main__":
    if sys.argv[1] == "--verify":
        main_verify()
    else:
        main_generate()
