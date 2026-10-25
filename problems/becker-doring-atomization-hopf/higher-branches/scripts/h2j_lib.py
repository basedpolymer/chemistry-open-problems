"""h2j_lib -- spectral certificate at the j-th Hopf point (any j >= 1), with contours following Re Lambda = const.

For an eps-box E and the Krawczyk box X of the Hopf point (kappa_j, z_H) (all eps in E), we prove:
  (L)  localization (as in h2_lib): every zero phi of F with |phi| >= 1, Re lambda >= 0 has |z| < R0, z = M(phi-1);
  (Mo) monotonicity: Re Lambda'(z) > 0 on the rectangle Rect = [x_min, R0] x [-Y, Y] (R0 <= Y), where
       Lambda = z + eps z^2/(kappa phi), so x -> Re Lambda(x + i y) is strictly increasing there;
  (R)  region: |z| <= R0 and Re z <= x_min imply Re Lambda < 0 (allm_core.region_check_sharp, sigma = -x_min);
  (P)  a polygon P = {l(y) < x < R0, |y| < Y} whose left chain l is a piecewise linear graph over y with
       Re Lambda < c_L < 0 on it (interval check on every edge);  winding number count N(P) = 2j + 2;
  (P-) a polygon P- = {l(y) < x < r(y), |y| < Y} whose right chain r has Re Lambda > c_R > 0 on it;
       N(P-) = 4  (z = 0 double, z_H, conj z_H).
Consequences (result.md, Theorem H_j): every zero with Re Lambda >= 0 lies in P; every zero with Re Lambda = 0 lies
in P- (monotonicity), hence is 0, z_H or conj z_H; exactly 2j - 2 zeros (with multiplicity) have Re Lambda > 0.
So at kappa_j(M): +- i omega are the only nonzero eigenvalues of B on the imaginary axis (H2), and the zeros of F
corresponding to eigenvalues in Re lambda > 0 number 2(j - 1).
"""
import math
from flint import arb, acb, ctx
import allm_core as C

def ReLam_box(zbox, kap, eps):
    """enclosure of Re Lambda over a complex box (all kappa in kap, eps in eps)"""
    e = acb(arb(eps)); z = acb(zbox); k = acb(arb(kap))
    phi = 1 + e * e * z
    return (z + e * z * z / (k * phi)).real

def mono_check(eps, kap, xmin, R0, Y):
    """upper bound of |Lambda'(z) - 1| on Rect (analytic bound); must be < 1"""
    e = C.arb_up(arb(eps)); kl = arb(arb(kap).lower())
    rho = C.arb_up((arb(max(abs(float(xmin)), float(R0))) ** 2 + arb(Y) ** 2).sqrt())
    q = e * e * rho
    assert q < 1
    return C.arb_up((e / kl) * rho * (2 + q) / (1 - q) ** 2)

def _solve_x(c, y, kap, eps, x0):
    """float Newton for Re Lambda(x + i y) = c at mid parameters"""
    k = float(arb(kap).mid()); e = float(arb(eps).mid())
    x = x0
    for _ in range(60):
        z = complex(x, y)
        L = z + e * z * z / (k * (1 + e * e * z))
        dL = 1 + (e / k) * z * (2 + e * e * z) / (1 + e * e * z) ** 2
        dx = (L.real - c) / dL.real
        x -= dx
        if abs(dx) < 1e-14:
            break
    return x

def level_chain(c, Y, n, kap, eps, side, margin):
    """vertices (x_k, y_k), y_k = -Y..Y increasing, on Re Lambda = c, shifted by margin (side=-1: left, +1: right)"""
    pts = []
    for i in range(n + 1):
        y = -Y + 2 * Y * i / n
        x = _solve_x(c, y, kap, eps, c)
        pts.append(complex(x + side * margin, y))
    return pts

def chain_check(pts, kap, eps, sign, thr, depth=12):
    """verify sign*(Re Lambda - thr) > 0 on every edge (refining edges by bisection, inserting the bisection point
    without moving it); returns (ok, refined vertex list)"""
    out = [pts[0]]
    for a, b in zip(pts[:-1], pts[1:]):
        stack = [(a, b, 0)]
        seg = []
        while stack:
            p, q, d = stack.pop()
            box = acb(p).union(acb(q))
            v = ReLam_box(box, kap, eps) - arb(thr)
            if (sign > 0 and v > 0) or (sign < 0 and v < 0):
                seg.append((p, q))
                continue
            if d >= depth:
                return False, None
            m = (p + q) / 2
            stack.append((m, q, d + 1)); stack.append((p, m, d + 1))
        # seg pieces are produced in order p..q
        for p, q in seg:
            out.append(q)
    return True, out

def polygon_P(left, R0, Y):
    """counterclockwise polygon: right edge x = R0 upward, top edge, left chain downward, bottom edge"""
    pts = [complex(R0, -Y), complex(R0, Y)]
    pts += list(reversed(left))            # left chain from (l(Y), Y) down to (l(-Y), -Y)
    pts.append(complex(R0, -Y))
    return pts

def polygon_Pm(left, right):
    pts = list(right)                       # right chain upward (-Y .. Y)
    pts += list(reversed(left))             # left chain downward
    pts.append(right[0])
    return pts

def subdivide(pts, nper):
    """insert nper-1 equally spaced points on every edge (keeps the same polygon)"""
    out = []
    for a, b in zip(pts[:-1], pts[1:]):
        for k in range(nper):
            out.append(a + (b - a) * k / nper)
    out.append(pts[-1])
    return out

def h2j_checks(j, eb, box, R0, Y, cL, cR, a1="0.1", nchain=48):
    out = dict(j=j, R0=R0, Y=Y, cL=cL, cR=cR)
    kap = box[0]
    el = arb(eb.lower())
    # localization (Lemma 3.3 with part II replaced by loc_check_II2_range: no Gershgorin, any K)
    sI = C.loc_check_I(eb, kap, R0, a1)
    okL = bool(sI < 1) and C.loc_check_II2_range(eb, kap, a1)[0]
    mode = "I+II2"
    if not okL and el > 0:
        okL = C.loc_check_II2_range(eb, kap, arb(R0) * el * el, n0=128, depth=18)[0]; mode = "II2"
    out.update(loc_ok=okL, loc_mode=mode, rhoI=float(sI))
    if not okL:
        out.update(ok=False, fail="localization"); return out
    margin = 0.02
    left = level_chain(cL, Y, nchain, kap, eb, -1, margin)
    okc, left = chain_check(left, kap, eb, -1, 0)
    if not okc:
        out.update(ok=False, fail="left chain"); return out
    right = level_chain(cR, Y, nchain, kap, eb, +1, margin)
    okc, right = chain_check(right, kap, eb, +1, 0)
    if not okc:
        out.update(ok=False, fail="right chain"); return out
    xmin = min(p.real for p in left) - 1e-9
    mono = mono_check(eb, kap, xmin, R0, Y)
    out["mono_ub"] = float(mono)
    if not (mono < 1) or not (R0 <= Y):
        out.update(ok=False, fail="monotonicity"); return out
    # zeros with Re Lambda >= 0 and |z| <= R0 have Re z > xmin (Lemma 4.2' with sigma = -xmin), so they lie in Rect
    if not (xmin < 0):
        out.update(ok=False, fail="xmin"); return out
    reg = C.region_check_sharp(eb, kap, R0, arb(-xmin))
    out["region_ub"] = float(reg)
    if not (reg < 0):
        out.update(ok=False, fail="region"); return out
    # the Hopf zero must be strictly inside Rect (|y| < Y) -- it then lies in P- by monotonicity
    if not (C.arb_up(abs(box[2])) < Y):
        out.update(ok=False, fail="Hopf zero outside"); return out
    P = subdivide(polygon_P(left, R0, Y), 1)
    Pm = polygon_Pm(left, right)
    nP, sP = C.winding_count(P, kap, eb)
    nPm, sPm = C.winding_count(Pm, kap, eb)
    out.update(NP=nP, NPm=nPm, segP=sP, segPm=sPm)
    out["ok"] = bool(nP == 2 * j + 2 and nPm == 4)
    if not out["ok"]:
        out["fail"] = "counts"
    return out
