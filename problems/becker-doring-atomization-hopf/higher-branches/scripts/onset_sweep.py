#!/usr/bin/env python3
"""onset_sweep.py -- the j-th Hopf curve from K = infinity to the eps-chart, in a chart regular at K = infinity.

Usage:  python onset_sweep.py J M_HI OUT.csv [L1.csv]
  J      branch
  M_HI   the sweep goes from the onset (K = infinity, M = M*_j) up to M = M_HI (real M)
  L1.csv l1_sweep certificate of the same branch (eps-boxes); if given, the junction with its last box is checked
Environment: Python 3.13, python-flint 0.9 (Arb), mpmath.  One thread.  Runtime: about 1-10 min.

Chart.  Real M (phi^M = exp(M log phi), principal branch), s = 1/A = 1/(1+K) in [0, 1), phi = 1 + u + i v.
Exactly (result.md, Section 3), with F = PV's characteristic function (71)-(76),
   F/A^3 = f1 (H1 - s^(M-1) phi^-M H2) - (G1 - phi^M G2) f2,
   f1 = phi^M + s(phi-1),  f2 = s - phi + s^(M-1) phi^(1-M),
   G1 = phi (s-phi)((1-s)phi + (M(1-s)+s)(phi-1)),  G2 = (s-phi)(M(phi-1) + (1-s)phi + phi(phi-1)),
   H1 = (phi-1)(s(1-s) + (M(1-s)+s)(s-phi)),        H2 = (phi-1)((s + M phi)(s-phi) + s(1-s)phi),
and lambda/A = (1 - s/phi)(phi - 1).  The Hopf system  H(s, M, u, v) = (Re F/A^3, Im F/A^3, Re lambda/A) = 0 is
C^1 on [0, 1) x box (the terms s^(M-1) are enclosed for s-balls containing 0); at s = 0 it is G(1 + iv; M) = 0,
u = 0 (kinf.py, claim K4: the onset (M*_j, eta_j)).
Part A (s-parametrized, s in [0, S_A]): for every s in each s-box, exactly one zero (M, u, v) in the box
  (parametric Krawczyk), dM/ds > 0 and v > 0 on the box.
Part B (M-parametrized, M from inside the M-range of the last A box up to M_HI): for every M in each M-box,
  exactly one zero (s, u, v) in the box (so the curve is a graph over M), s > 0, v > 0 on the box.
Junction A/B: the first B box (for its M-range) lies inside the last A box, so the two pieces are one curve.
Junction with L1.csv: one B sub-box whose image in the eps-chart (eps = M^-1/2, kappa = (1/s - 1) M^1/2,
  z = M(phi - 1)) lies inside the last Krawczyk box of L1.csv (for its eps-range); so the curve continues the
  certified eps-chart arc.  Exit 0 iff all claims hold and the boxes cover [0, S_A] and [M_B0, M_HI].
"""
import sys, os, csv, time, math
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from flint import arb, acb, ctx, arb_mat
import mpmath as mp
import allm_core as C

class D3:
    """value + holomorphic derivative in phi + derivatives in M and s (acb)"""
    __slots__ = ("v", "dp", "dM", "ds")
    def __init__(s, v, dp=0, dM=0, ds=0):
        s.v = acb(v); s.dp = acb(dp); s.dM = acb(dM); s.ds = acb(ds)
    @staticmethod
    def c(x): return x if isinstance(x, D3) else D3(x)
    def __add__(a, b): b = D3.c(b); return D3(a.v + b.v, a.dp + b.dp, a.dM + b.dM, a.ds + b.ds)
    __radd__ = __add__
    def __neg__(a): return D3(-a.v, -a.dp, -a.dM, -a.ds)
    def __sub__(a, b): return a + (-D3.c(b))
    def __rsub__(a, b): return D3.c(b) - a
    def __mul__(a, b):
        b = D3.c(b)
        return D3(a.v * b.v, a.dp * b.v + a.v * b.dp, a.dM * b.v + a.v * b.dM, a.ds * b.v + a.v * b.ds)
    __rmul__ = __mul__

def phi_pow(phi, M, sgn=1, shift=0):
    """phi^(sgn M + shift) as D3 (phi, M balls)"""
    lp = acb(phi).log()
    ex = sgn * acb(M) + shift
    P = (ex * lp).exp()
    return D3(P, ex * P / phi, sgn * P * lp, 0)

def s_pow(s, M, shift):
    """s^(M + shift) (exponent p = M + shift >= 2) as D3 in (M, s); valid for s-balls containing 0 (s >= 0):
    for t in [0, b], p in [p_lo, p_hi]: t^p in [0, b^p_lo]; d/ds = p t^(p-1) in [0, p_hi b^(p_lo-1)];
    d/dM = t^p log t in [-b^p_lo |log b|, 0] (b <= exp(-1/p_lo))."""
    s = arb(s); M = arb(M)
    p = M + shift
    if arb(s.lower()) > 0:
        ls = s.log(); P = (p * ls).exp()
        return D3(P, 0, P * ls, p * P / s)
    b = C.arb_up(s)
    plo = arb(p.lower()); phi_ = C.arb_up(p)
    if b == 0:
        return D3(0, 0, 0, 0)                      # s = 0 exactly (p >= 2: value and derivatives vanish)
    assert b < 1 and b <= (-1 / plo).exp()
    v = b ** plo
    val = arb(0).union(v)
    dsv = arb(0).union(phi_ * b ** (plo - 1))
    dMv = arb(0).union(-v * (-b.log()))
    return D3(acb(val), 0, acb(dMv), acb(dsv))

def system(s, M, phi):
    """(F/A^3, lambda/A) as D3 at s (real ball >= 0), M (real ball), phi (complex ball)"""
    s = arb(s)
    sd = D3(acb(s), 0, 0, 1); Md = D3(acb(M), 0, 1, 0); pd = D3(acb(phi), 1, 0, 0)
    one = D3(1)
    PM = phi_pow(phi, M, 1, 0)                     # phi^M
    PmM = phi_pow(phi, M, -1, 0)                   # phi^-M
    P1mM = phi_pow(phi, M, -1, 1)                  # phi^(1-M)
    SM1 = s_pow(s, M, -1)                          # s^(M-1)
    f1 = PM + sd * (pd - 1)
    f2 = sd - pd + SM1 * P1mM
    a = Md * (one - sd) + sd                       # M(1-s) + s
    G1 = pd * (sd - pd) * ((one - sd) * pd + a * (pd - 1))
    G2 = (sd - pd) * (Md * (pd - 1) + (one - sd) * pd + pd * (pd - 1))
    H1 = (pd - 1) * (sd * (one - sd) + a * (sd - pd))
    H2 = (pd - 1) * ((sd + Md * pd) * (sd - pd) + sd * (one - sd) * pd)
    F = f1 * (H1 - SM1 * PmM * H2) - (G1 - PM * G2) * f2
    # lambda/A = (1 - s/phi)(phi - 1) = phi - 1 - s + s/phi
    iphi = 1 / acb(phi)
    lam = D3(acb(phi) - 1 - acb(s) + acb(s) * iphi, 1 - acb(s) * iphi * iphi, 0, -1 + iphi)
    return F, lam

# ------------------------------------------------------------------------------------------------ the two charts
# chart "s": parameter s, unknowns X = (M, u, v);  chart "M": parameter M, unknowns X = (s, u, v)
def fjac(chart, X, p):
    """f = (Re F, Im F, Re lam), Jacobian in the unknowns X, and d f / d(parameter)"""
    if chart == "s":
        M, u, v = X; s = p
    else:
        s, u, v = X; M = p
    F, L = system(s, M, acb(1 + arb(u), v))
    f = [F.v.real, F.v.imag, L.v.real]
    Ju = [F.dp.real, F.dp.imag, L.dp.real]           # d/du (holomorphic in phi)
    Jv = [-F.dp.imag, F.dp.real, -L.dp.imag]         # d/dv = i d/dphi
    Js = [F.ds.real, F.ds.imag, L.ds.real]
    JM = [F.dM.real, F.dM.imag, L.dM.real]
    first, par = (JM, Js) if chart == "s" else (Js, JM)
    J = [[first[i], Ju[i], Jv[i]] for i in range(3)]
    return f, J, par

def newton(chart, X0, p, it=80):
    X = [arb(arb(x).mid()) for x in X0]
    pm = arb(arb(p).mid())
    for _ in range(it):
        f, J, _ = fjac(chart, X, pm)
        Jm = arb_mat([[c.mid() for c in r] for r in J]); d = Jm.solve(arb_mat([[c.mid()] for c in f]))
        X = [arb((X[i] - d[i, 0]).mid()) for i in range(3)]
        if max(abs(float(d[i, 0].mid())) / (1 + abs(float(X[i].mid()))) for i in range(3)) < 2.0 ** (-ctx.prec + 40):
            break
    return X

def krawczyk_img(chart, Xm, rad, P):
    """parametric Krawczyk image over the parameter ball P:  m - Y f(m; p_m) - (Y f_p(m; P))(P - p_m)
    + (I - Y J(X; P))(X - m)  (Y f_p formed before multiplying by the scalar interval)"""
    P = arb(P); pm = arb(P.mid())
    box = [arb(Xm[i], rad[i]) for i in range(3)]
    f0, Jm, _ = fjac(chart, Xm, pm)
    _, _, fpP = fjac(chart, Xm, P)
    _, JX, _ = fjac(chart, box, P)
    Y = arb_mat([[c.mid() for c in r] for r in Jm]).inv()
    Y = arb_mat([[Y[i, k].mid() for k in range(3)] for i in range(3)])
    I3 = arb_mat([[1 if i == k else 0 for k in range(3)] for i in range(3)])
    dev = arb_mat([[box[i] - Xm[i]] for i in range(3)])
    Yfp = Y * arb_mat([[c] for c in fpP])
    Kx = arb_mat([[Xm[i]] for i in range(3)]) - Y * arb_mat([[c] for c in f0]) - Yfp * (P - pm) + (I3 - Y * arb_mat(JX)) * dev
    return [Kx[i, 0] for i in range(3)], box

def krawczyk(chart, Xm, rad, P):
    img, box = krawczyk_img(chart, Xm, rad, P)
    if all(box[i].contains_interior(img[i]) for i in range(3)):
        return img, box
    return None, box

def box_for(chart, plo, phi_, X0):
    """Newton at the ends and the middle of [plo, phi_], then epsilon-inflation of the Krawczyk box"""
    P = C.eps_interval(plo, phi_)
    sols = [newton(chart, X0, arb(e)) for e in (arb(plo), arb(P.mid()), arb(phi_))]
    Xm = sols[1]
    dev = [max(abs(float(q[i].mid()) - float(Xm[i].mid())) for q in sols) for i in range(3)]
    rad = [arb(1.5 * dev[i] + 1e-25 * (1 + abs(float(Xm[i].mid())))) for i in range(3)]
    for _ in range(30):
        img, box = krawczyk_img(chart, Xm, rad, P)
        if all(box[i].contains_interior(img[i]) for i in range(3)):
            return P, Xm, rad, img, box
        new = []
        for i in range(3):
            d = abs(img[i] - Xm[i])
            dup = float(arb(d.upper())) if d.is_finite() else float("inf")
            new.append(arb(max(float(rad[i]) * 1.1, 1.25 * dup + 1e-30)))
        if not all(abs(float(r)) < 1e2 for r in new):
            break
        rad = new
    return P, Xm, None, None, None

def slope(chart, box, P):
    """enclosure of d(first unknown)/d(parameter) along the curve, -[J^-1 f_p]_0 with J, f_p over the box;
    the linear system is preconditioned by the inverse of the midpoint Jacobian"""
    _, JX, fp = fjac(chart, box, P)
    Jm = arb_mat([[c.mid() for c in r] for r in JX])
    try:
        Y = Jm.inv()
        Y = arb_mat([[Y[i, k].mid() for k in range(3)] for i in range(3)])
        sol = (Y * arb_mat(JX)).solve(Y * arb_mat([[c] for c in fp]))
    except ZeroDivisionError:
        return None
    return -sol[0, 0]

def kinf_point(j):
    """float onset (M*_j, eta_j) from G(1 + i eta; M) = 0 (kinf.py): eta from a scan of |G| phi^-M on the line at a
    guess of M (table, or 5.31 j^2), then a 2-d Newton iteration"""
    mp.mp.dps = 80
    def Gm(M, e):
        p = 1 + 1j * e
        return -mp.exp((M - 2) * mp.log(p)) * (p ** 3 + M * p ** 2 - 2 * M * p + M) + (M + 1) * p - M
    guess = {1: 14.893, 2: 43.664, 3: 86.387, 4: 141.98, 5: 209.63, 6: 288.74, 7: 378.87, 8: 479.66, 9: 590.83,
             10: 712.15, 11: 843.40, 12: 984.40, 13: 1134.995}
    M0 = mp.mpf(guess.get(j, 5.31 * j * j))
    best = None
    for i in range(4001):
        zz = 2 * mp.pi * j - 3 + 3.5 * i / 4000
        e = zz / M0
        v = abs(Gm(M0, e)) / abs(mp.exp(M0 * mp.log(1 + 1j * e)))
        if best is None or v < best[0]:
            best = (v, e)
    r = mp.findroot(lambda M, e: [mp.re(Gm(M, e)), mp.im(Gm(M, e))], (M0, best[1]), verify=False)
    r = mp.findroot(lambda M, e: [mp.re(Gm(M, e)), mp.im(Gm(M, e))], (r[0], r[1]), verify=False)
    return [arb(mp.nstr(r[0], 50)), arb(0), arb(mp.nstr(r[1], 50))]

FIELDS = ["j", "chart", "p_lo", "p_hi", "X0_mid", "u_mid", "v_mid", "rad_X0", "rad_u", "rad_v", "M_lo", "M_hi",
          "s_lo", "s_hi", "slope_lo", "status"]

def sweep(chart, p0, p1, X0, w0, wmax, wr, rows, t0):
    lo = p0; w = w0; prev = None; nf = 0
    while float(lo) < float(p1):
        hi = "%.17g" % min(float(lo) + w, float(p1))
        try:
            P, Xm, rad, img, box = box_for(chart, lo, hi, X0)
        except (AssertionError, ZeroDivisionError, ValueError):
            img = None                      # Newton left the domain (e.g. s < 0): treat as a failed box
        ok = img is not None
        if ok:
            # chart s: M must increase along the curve (dM/ds > 0), so that every M is met once;
            # chart M: the curve is a graph over M by construction; s > 0 (K finite) is required
            sl = slope(chart, box, P) if chart == "s" else arb(0)
            vpos = bool(arb(box[2].lower()) > 0)
            ok = vpos and ((chart == "s" and sl is not None and bool(sl > 0)) or
                           (chart == "M" and bool(arb(box[0].lower()) > 0)))
        if not ok:
            if w > 1e-14 * max(1.0, abs(float(lo))):
                w /= 2; continue
            nf += 1
            print("%s-box [%s, %s] FAIL" % (chart, lo, hi), flush=True)
            break
        if chart == "s":
            Mr = box[0]; Sr = P
        else:
            Mr = P; Sr = box[0]
        wr.writerow([J_GLOBAL, chart, lo, hi, Xm[0].mid().str(30, radius=False), Xm[1].mid().str(30, radius=False),
                     Xm[2].mid().str(30, radius=False), "%.6e" % float(rad[0]), "%.6e" % float(rad[1]),
                     "%.6e" % float(rad[2]), "%.17g" % float(Mr.lower()), "%.17g" % float(Mr.upper()),
                     "%.17g" % float(Sr.lower()), "%.17g" % float(Sr.upper()), "%.6g" % float(sl.lower()), "PASS"])
        rows.append((chart, lo, hi, P, box))
        if len(rows) % 100 == 0:
            print("%d boxes; %s-box [%s, %s] M in [%.6f, %.6f] s in [%.6g, %.6g] (%.0fs)" % (
                len(rows), chart, lo, hi, float(Mr.lower()), float(Mr.upper()), float(Sr.lower()), float(Sr.upper()),
                time.time() - t0), flush=True)
        X0 = [arb(2 * float(Xm[i].mid()) - float(prev[i].mid())) for i in range(3)] if prev is not None else Xm
        prev = Xm
        lo = hi
        w = min(w * 1.4, wmax)
    return nf == 0 and float(lo) >= float(p1)

J_GLOBAL = 0

def main():
    global J_GLOBAL
    j = int(sys.argv[1]); Mhi = float(sys.argv[2]); out = sys.argv[3]
    J_GLOBAL = j
    l1file = sys.argv[4] if len(sys.argv) > 4 else None
    ctx.prec = 192
    t0 = time.time()
    X0 = kinf_point(j)
    S_A = 2e-4
    rows = []
    with open(out, "w", newline="") as fh:
        wr = csv.writer(fh); wr.writerow(FIELDS)
        okA = sweep("s", "0", "%.17g" % S_A, X0, 1e-6, 2e-5, wr, rows, t0)
        nA = len(rows)
        okAB = okB = False
        if okA:
            # first M-box: a small interval around the middle of the last A box's M-range, inside it
            chart, lo, hi, P, boxA = rows[-1]
            Mmid = float(boxA[0].mid()); h = 0.05 * float(boxA[0].rad())
            Xs = [arb((1 / P).mid()) if False else arb(P.mid()), boxA[1], boxA[2]]
            Xs = [arb(P.mid()), arb(boxA[1].mid()), arb(boxA[2].mid())]
            PB, XmB, radB, imgB, boxB = box_for("M", "%.17g" % (Mmid - h), "%.17g" % (Mmid + h), Xs)
            okAB = (imgB is not None and P.contains(boxB[0]) and boxA[0].contains(PB) and boxA[1].contains(boxB[1])
                    and boxA[2].contains(boxB[2]))
            print("junction A/B: M-box %s inside the last s-box: %s" % (PB.str(12), okAB), flush=True)
            if okAB and Mhi > Mmid + h:
                okB = sweep("M", "%.17g" % (Mmid - h), "%.17g" % Mhi, Xs, 1e-3, 0.5, wr, rows, t0)
    cover = okA and okAB and okB
    Mlo = float(rows[0][4][0].lower()) if rows else float("nan")
    print("onset sweep j=%d: %d s-boxes (s in [0, %g]) and %d M-boxes up to M = %s; covered: %s; M*_j in "
          "the first box M-range [%.12f, %.12f]" % (j, nA, S_A, len(rows) - nA, Mhi, cover, Mlo,
                                                   float(rows[0][4][0].upper()) if rows else float("nan")), flush=True)
    okj = True
    if l1file is not None and cover:
        okj = junction(l1file, rows)
    okall = cover and okj
    jtxt = ("PASS" if okj else "FAIL") if l1file is not None else "not requested"
    print("SUMMARY j=%d: onset sweep %s, junction with the eps-chart %s -> %s (%.0fs)" % (
        j, "PASS" if cover else "FAIL", jtxt, "PASS" if okall else "FAIL", time.time() - t0))
    sys.exit(0 if okall else 1)

def junction(l1file, rows):
    """an M-sub-box whose image lies in the last Krawczyk box of the l1 file (eps in its eps-range)"""
    r = list(csv.DictReader(open(l1file)))[-1]
    eb = C.eps_interval(r["eps_lo"], r["eps_hi"])
    Xm = [arb(r["kappa_mid"]), arb(r["x_mid"]), arb(r["y_mid"])]
    rad = [arb(r["rad_kappa"]), arb(r["rad_x"]), arb(r["rad_y"])]
    enc, X = C.krawczyk(Xm, rad, eb)
    if enc is None:
        print("junction: l1 box not re-proved"); return False
    Mt = 1 / float(eb.mid()) ** 2
    for (chart, lo, hi, P, box) in rows:
        if chart == "M" and float(lo) <= Mt <= float(hi):
            Xs = [arb(box[i].mid()) for i in range(3)]
            break
    else:
        print("junction: M = %.6f not in the onset sweep" % Mt); return False
    for wloc in (1e-6, 1e-7, 1e-8, 1e-9):
        h = wloc * Mt
        PB, XmB, radB, img, bx = box_for("M", "%.17g" % (Mt - h), "%.17g" % (Mt + h), Xs)
        if img is None:
            continue
        M = PB; s = bx[0]; phi = acb(1 + bx[1], bx[2])
        eps_img = 1 / M.sqrt()
        kap_img = (1 / s - 1) * M.sqrt()
        z_img = M * (phi - 1)
        inside = (eb.contains(eps_img) and X[0].contains(kap_img) and X[1].contains(z_img.real)
                  and X[2].contains(z_img.imag))
        insweep = any(ch == "M" and arb(P2.lower()) <= arb(PB.lower()) and arb(PB.upper()) <= arb(P2.upper()) and
                      all(b2[i].contains(bx[i]) for i in range(3)) for (ch, lo2, hi2, P2, b2) in rows)
        print("junction: M in %s -> eps in %s, kappa in %s; inside the last l1 box: %s; inside a sweep box: %s" % (
            PB.str(12), eps_img.str(10), kap_img.str(10), inside, insweep), flush=True)
        if inside and insweep:
            return True
    return False

if __name__ == "__main__":
    main()
