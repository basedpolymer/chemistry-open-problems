"""h2_lib -- certificate of the spectral hypothesis H2 (and of the stability count) on an eps-box.

For a box E = [eps_lo, eps_hi] (eps = M^-1/2) the routine `certify_box` proves, for every eps in E:
  (H)  the scaled Hopf system  Q^eps(z; kappa) = 0, Re Lambda(z; kappa, eps) = 0  has a unique zero
       (kappa, z) = (kappa_1(eps), z_H(eps)) in an explicit box X (parametric Krawczyk, allm_core.krawczyk);
       on X: dQ/dz != 0 (simple root), d Re Lambda/d kappa > 0 along the root curve (transversality),
       Im z > 0, |phi f(phi)| > 0 (needed in Pego-Velazquez Lemma 5.4);
  (L)  localization: F has no zero phi with |phi| >= 1, Re lambda(phi) >= 0 and |phi - 1| >= R0/M
       (Lemma L of result.md: part I = ratio bounds on [R0/M, a1], part II = growth of |phi|^M on [a1, a_max];
        or part II alone on [R0 eps_lo^2, a_max] when eps_lo > 0);
  (R)  region: |z| <= R0 and Re z <= -sigma imply Re Lambda < 0;
  (W)  the argument principle on a convex polygon P containing D = {|z| <= R0, Re z >= -sigma} gives exactly
       4 zeros of Q^eps(.; kappa) in P (z = 0 double, z_H, conj(z_H)), for all kappa in the kappa-range of X.
Consequence (result.md, Theorem H): for every integer M with M^-1/2 in E, at kappa = kappa_1(M) the matrix B has
exactly the eigenvalues 0 (simple), +-i omega (simple) in the closed right half-plane, all others Re < 0.
"""
import math
from flint import arb, acb, ctx
import allm_core as C

def hopf_box(eps_lo, eps_hi, X0):
    """parametric Krawczyk over [eps_lo, eps_hi]; radii from the predicted variation of the solution"""
    eb = C.eps_interval(eps_lo, eps_hi)
    sols = []
    for e in (arb(eps_lo), arb(eb.mid()), arb(eps_hi)):
        sols.append(C.newton(X0, e))
    Xm = sols[1]
    dev = [max(abs(float(s[i].mid()) - float(Xm[i].mid())) for s in sols) for i in range(3)]
    for fac in (1.5, 4, 16, 64):
        rad = [arb(fac * dev[i] + 1e-25 * (1 + abs(float(Xm[i].mid())))) for i in range(3)]
        enc, box = C.krawczyk(Xm, rad, eb)
        if enc is not None:
            return eb, Xm, rad, enc, box
    return eb, Xm, None, None, None

def certify_box(eps_lo, eps_hi, X0, R0=8, sigma="0.3", a1="0.1", n_arc=96, n_seg=24):
    """returns dict with all checks (see module docstring); 'ok' is True iff all hold"""
    out = dict(eps_lo=str(eps_lo), eps_hi=str(eps_hi), R0=R0, sigma=sigma, a1=a1)
    eb, Xm, rad, enc, box = hopf_box(eps_lo, eps_hi, X0)
    out["Xm"] = Xm; out["rad"] = rad; out["box"] = box
    if enc is None:
        out["ok"] = False; out["fail"] = "krawczyk"; return out
    simple, tr, fphi, laim = C.hopf_checks(box, eb)
    imz_pos = arb(box[2].lower()) > 0
    out.update(simple=simple, transv_lo=float(arb(tr.lower())), fphi_lo=float(fphi.abs_lower()), imz_pos=imz_pos)
    kap = box[0]
    # localization
    el = arb(eb.lower())
    sI = C.loc_check_I(eb, kap, R0, a1)
    okI = bool(sI < 1)
    okII, nII = C.loc_check_II_range(eb, kap, a1)
    mode = "I+II"
    okL = okI and okII
    if not okL and el > 0:
        okL, nII = C.loc_check_II_range(eb, kap, arb(R0) * el * el, n0=256, depth=18)
        mode = "II"
    out.update(loc_mode=mode, rhoI=float(sI), loc_ok=okL)
    # region
    reg = C.region_check(eb, kap, R0, sigma)
    out["region_ub"] = float(reg); okR = bool(reg < 0)
    # containment of the Hopf zero in D (strictly), with margin
    zb = acb(box[1], box[2])
    inD = bool(acb_abs_up(zb) < R0) and bool(arb(box[1].lower()) > -arb(sigma))
    # winding
    pts, R1 = C.contour_D(R0, float(sigma), n_arc, n_seg)
    n, nseg = C.winding_count(pts, kap, eb)
    out.update(count=n, nseg=nseg if n is not None else str(nseg), inD=inD)
    out["ok"] = bool(simple and tr > 0 and fphi.abs_lower() > 0 and imz_pos and okL and okR and inD and n == 4)
    if not out["ok"]:
        out["fail"] = "checks"
    return out

def acb_abs_up(z):
    return arb(acb(z).abs_upper())

CSV_FIELDS = ["eps_lo", "eps_hi", "kappa_mid", "x_mid", "y_mid", "rad_kappa", "rad_x", "rad_y",
              "kappa_lo", "kappa_hi", "simple", "transv_lo", "fphi_lo", "imz_pos", "loc_mode", "rhoI",
              "loc_ok", "region_ub", "R0", "sigma", "a1", "inD", "count", "nseg", "status"]

def row(out):
    Xm, rad, box = out["Xm"], out["rad"], out["box"]
    return [out["eps_lo"], out["eps_hi"], Xm[0].mid().str(40, radius=False), Xm[1].mid().str(40, radius=False),
            Xm[2].mid().str(40, radius=False), "%.6e" % float(rad[0]), "%.6e" % float(rad[1]), "%.6e" % float(rad[2]),
            "%.17g" % float(box[0].lower()), "%.17g" % float(box[0].upper()), int(out["simple"]),
            "%.6g" % out["transv_lo"], "%.6g" % out["fphi_lo"], int(out["imz_pos"]), out["loc_mode"],
            "%.6g" % out["rhoI"], int(out["loc_ok"]), "%.6g" % out["region_ub"], out["R0"], out["sigma"], out["a1"],
            int(out["inD"]), out["count"], out["nseg"], "PASS" if out["ok"] else "FAIL"]
