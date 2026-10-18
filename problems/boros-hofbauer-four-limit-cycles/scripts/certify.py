"""
certify.py -- computer-assisted certificate for the Boros-Hofbauer four-limit-cycle conjecture.

Environment : WSL Fedora, Python 3.12, python-flint 0.9.0 (fmpq_mpoly, arb), mpmath 1.3.0.
Runtime     : ~2 s.   Deterministic.
Usage       : python certify.py [--write]
              (--write stores the polynomials and the box in ../certificates/; without it, the
               stored files are re-checked against the recomputation.)

What is certified (all claims below are checked with exact rational arithmetic or with Arb ball
arithmetic, i.e. rigorous interval arithmetic with outward rounding):

  Notation.  For the system  x' = 1/y - x^a y^b,  y' = K(-(d-1) - 1/y + d x^a y^b)  with K = a/(1+bd)
  (trace of the Jacobian at (1,1) equal to 0), put
        c = 1 + b,   s = (d-1)/(1+bd) = omega^2/a^2   (so d = (1+s)/(1-(c-1)s)).
  focal_exact.py produces V~_k in Q[a,c,w] (w = sqrt(s)); V~_k = w * G_k(a,c,w^2) with G_k in Q[a,c,s]
  (checked).  The focal values in the normalisation of the paper are
        L_k = pi * G_k(a,c,s) / (a^(2k) c^(4k) s^(k-1/2)),
  so on the parameter domain (a, c, s > 0) sign L_k = sign G_k.

  C0  the polynomials G_1..G_4 recomputed now coincide with certificates/G_k.json;
  C1  domain: on the box X, a > 0, c > 0 (b > -1), s > 0, 1-(c-1)s > 0, hence d > 1, 1+bd > 0, K > 0;
  C2  Krawczyk: K(X) subset of int(X) for G = (G_1,G_2,G_3) => unique zero mu* of G in X;
  C3  || I - C J(X) ||_inf < 1  => every matrix in the interval Jacobian J(X) = dG/d(a,c,s)(X) is
      nonsingular (in particular dG(mu*) is invertible);
  C4  G_4 < 0 on X (hence L_4(mu*) < 0);
  C5  the Jacobian of (L_1,L_2,L_3) w.r.t. (a,b,d) at mu* is nonsingular (interval determinant
      excludes 0, computed from J(X) via the chain rule and the row scalings).
"""
import sys
import os
import json
import mpmath as mp
from flint import fmpq, arb, ctx, fmpq_mpoly_ctx, Ordering

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from focal_exact import lyapunov_exact  # noqa: E402

CERT = os.path.join(HERE, "..", "certificates")
SCTX = fmpq_mpoly_ctx.get(('a', 'c', 's'), Ordering.degrevlex)
results = []


def report(name, ok, info=""):
    results.append(bool(ok))
    print(("PASS" if ok else "FAIL") + "  " + name + (("  [" + info + "]") if info else ""), flush=True)


def to_G(V):
    """V~ = w * G(a,c,w^2)  ->  G in Q[a,c,s]."""
    d = {}
    for e, cf in V.to_dict().items():
        ea, ec, ew = int(e[0]), int(e[1]), int(e[2])
        if ew % 2 != 1:
            raise ValueError("V~ is not w * polynomial(a, c, w^2)")
        d[(ea, ec, (ew - 1) // 2)] = cf
    return SCTX.from_dict(d)


def poly_to_json(P):
    terms = [[int(e[0]), int(e[1]), int(e[2]), str(cf)] for e, cf in P.to_dict().items()]
    terms.sort()
    return {"vars": ["a", "c", "s"], "terms": terms}


def poly_from_json(data):
    d = {}
    for t in data["terms"]:
        num, _, den = t[3].partition('/')
        d[(t[0], t[1], t[2])] = fmpq(int(num), int(den) if den else 1)
    return SCTX.from_dict(d)


def ev_arb(P, x):
    """Evaluate P in Q[a,c,s] on arb balls x = (a, c, s) (naive term-by-term, rigorous)."""
    degs = P.degrees()
    pw = []
    for i in range(3):
        L = [arb(1)]
        for _ in range(int(degs[i])):
            L.append(L[-1] * x[i])
        pw.append(L)
    tot = arb(0)
    for e, cf in P.to_dict().items():
        tot += arb(cf) * pw[0][int(e[0])] * pw[1][int(e[1])] * pw[2][int(e[2])]
    return tot


def ev_mp(P, x):
    tot = mp.mpf(0)
    for e, cf in P.to_dict().items():
        tot += mp.mpf(int(cf.p)) / int(cf.q) * x[0]**int(e[0]) * x[1]**int(e[1]) * x[2]**int(e[2])
    return tot


def main(write):
    ctx.prec = 512
    mp.mp.dps = 150
    # ---------------- C0: polynomials ----------------
    Vt = lyapunov_exact(4, verbose=False)
    G = [to_G(V) for V in Vt]
    names = ["G_%d.json" % k for k in range(1, 5)]
    if write:
        os.makedirs(CERT, exist_ok=True)
        for P, nm in zip(G, names):
            with open(os.path.join(CERT, nm), "w") as fh:
                json.dump(poly_to_json(P), fh)
    same = True
    for P, nm in zip(G, names):
        with open(os.path.join(CERT, nm)) as fh:
            same = same and (poly_from_json(json.load(fh)) == P)
    report("C0 recomputed G_1..G_4 == certificates/G_k.json", same,
           "terms: " + ", ".join(str(len(P)) for P in G))
    # derivative polynomials
    dG = [[G[k].derivative(i) for i in range(3)] for k in range(3)]

    # ---------------- approximate zero (Newton, 150 digits; not part of the proof) ----------------
    if write:
        a0, b0, d0 = mp.mpf('1.01282'), mp.mpf('0.65463'), mp.mpf('3.28862')
        x = [a0, 1 + b0, (d0 - 1) / (1 + b0 * d0)]
        for _ in range(60):
            F = mp.matrix([ev_mp(G[k], x) for k in range(3)])
            J = mp.matrix([[ev_mp(dG[k][i], x) for i in range(3)] for k in range(3)])
            dx = mp.lu_solve(J, F)
            x = [x[i] - dx[i] for i in range(3)]
            if mp.norm(dx) < mp.mpf(10)**-140:
                break
        # centre: rationals with denominator 10^60; radius 10^-40
        den = 10**60
        m = [fmpq(int(mp.nint(xi * den)), den) for xi in x]
        r = fmpq(1, 10**40)
        # preconditioner: approximate inverse of J(m), rounded to rationals (any real matrix is allowed)
        J = mp.matrix([[ev_mp(dG[k][i], x) for i in range(3)] for k in range(3)])
        Ci = J**-1
        Cq = [[str(fmpq(int(mp.nint(Ci[i, j] * 10**50)), 10**50)) for j in range(3)] for i in range(3)]
        box = {"variables": ["a", "c", "s"],
               "centre": [str(v) for v in m], "radius": str(r), "preconditioner_C": Cq,
               "comment": "box X = prod [centre_i - radius, centre_i + radius]; c = 1+b, s = (d-1)/(1+bd)"}
        with open(os.path.join(CERT, "box.json"), "w") as fh:
            json.dump(box, fh, indent=1)
    with open(os.path.join(CERT, "box.json")) as fh:
        box = json.load(fh)

    def q(sv):
        num, _, den = sv.partition('/')
        return fmpq(int(num), int(den) if den else 1)
    m = [q(v) for v in box["centre"]]
    r = q(box["radius"])
    Cq = [[q(v) for v in row] for row in box["preconditioner_C"]]
    C = [[arb(v) for v in row] for row in Cq]
    X = [arb(mi, r) for mi in m]          # ball with exact rational centre and radius (rounded outward)
    Xlo = [mi - r for mi in m]
    Xhi = [mi + r for mi in m]
    M = [arb(mi) for mi in m]
    print("box centre (a, c, s) =", [str(arb(v).str(45)) for v in m], " radius =", str(r))
    report("C0' the Arb balls X_i contain the exact box [m_i - r, m_i + r]",
           all(X[i].contains(arb(Xlo[i])) and X[i].contains(arb(Xhi[i])) for i in range(3)))

    # ---------------- C1: domain ----------------
    a, c, s = X
    rho = 1 - (c - 1) * s
    dom = (a > 0) and (c > 0) and (s > 0) and (rho > 0)
    b = c - 1
    d = (1 + s) / rho
    Kc = a / (1 + b * d)
    dom = dom and (d > 1) and (1 + b * d > 0) and (Kc > 0)
    report("C1 domain a>0, b>-1, s>0, 1-(c-1)s>0, d>1, 1+bd>0, K>0 on X", dom,
           "a in %s, b in %s, d in %s, K in %s" % (a.str(20), b.str(20), d.str(20), Kc.str(20)))

    # ---------------- C2: Krawczyk ----------------
    Gm = [ev_arb(G[k], M) for k in range(3)]
    JX = [[ev_arb(dG[k][i], X) for i in range(3)] for k in range(3)]
    CG = [sum((C[i][j] * Gm[j] for j in range(3)), arb(0)) for i in range(3)]
    Mmat = [[(arb(1) if i == j else arb(0)) - sum((C[i][k] * JX[k][j] for k in range(3)), arb(0))
             for j in range(3)] for i in range(3)]
    D = [X[i] - M[i] for i in range(3)]
    Kx = [M[i] - CG[i] + sum((Mmat[i][j] * D[j] for j in range(3)), arb(0)) for i in range(3)]
    inside = all((Kx[i] > arb(Xlo[i])) and (Kx[i] < arb(Xhi[i])) for i in range(3))
    report("C2 Krawczyk K(X) in int(X)  (unique zero of (G1,G2,G3) in X)", inside,
           "max |K_i - m_i| <= %s, radius %s" % (max((Kx[i] - M[i]).abs_upper() for i in range(3)).str(5),
                                                 str(r)))
    # ---------------- C3: nonsingularity ----------------
    rown = [sum((Mmat[i][j].abs_upper() for j in range(3)), arb(0)) for i in range(3)]
    nrm = max(rown[0], rown[1], rown[2])
    report("C3 ||I - C J(X)||_inf < 1  (all Jacobians in J(X) nonsingular)", nrm < 1, "norm <= %s" % nrm.str(5))
    detJ = (JX[0][0] * (JX[1][1] * JX[2][2] - JX[1][2] * JX[2][1])
            - JX[0][1] * (JX[1][0] * JX[2][2] - JX[1][2] * JX[2][0])
            + JX[0][2] * (JX[1][0] * JX[2][1] - JX[1][1] * JX[2][0]))
    report("C3' det dG/d(a,c,s) on X excludes 0", not detJ.contains(0), "det in %s" % detJ.str(15))

    # ---------------- C4: sign of the fourth focal value ----------------
    G4 = ev_arb(G[3], X)
    L4 = arb.pi() * G4 / (a**8 * c**16 * s**3 * s.sqrt())
    report("C4 G_4 < 0 on X  (L_4 < 0 at the zero)", G4 < 0, "G_4 in %s ; L_4 in %s" % (G4.str(15), L4.str(15)))
    # also the 5th one, for information
    # ---------------- C5: Jacobian w.r.t. (a, b, d) ----------------
    # L_k = pi G_k / (a^(2k) c^(4k) s^(k-1/2)).  At the zero (G_k = 0):  dL_k = pi dG_k / (...)  (positive
    # row scaling), so it suffices that dG/d(a,b,d) = dG/d(a,c,s) * d(a,c,s)/d(a,b,d) is nonsingular;
    # d(a,c,s)/d(a,b,d) = [[1,0,0],[0,1,0],[0, ds/db, ds/dd]] with ds/dd = (1+b)/(1+bd)^2 > 0.
    dsdd = (1 + b) / (1 + b * d)**2
    det_abd = detJ * dsdd
    report("C5 det d(L1,L2,L3)/d(a,b,d) at the zero is nonzero", not det_abd.contains(0) and not dsdd.contains(0),
           "det dG/d(a,b,d) in %s" % det_abd.str(10))
    for k in range(3):
        print("     G_%d(centre) in %s" % (k + 1, Gm[k].str(5)))
    out = {"box": box, "enclosures": {
        "a": a.str(40), "b": b.str(40), "d": d.str(40), "K": Kc.str(40),
        "omega": (a * s.sqrt()).str(30),
        "G4_on_box": G4.str(30), "L4_on_box": L4.str(30),
        "det_dG_dacs_on_box": detJ.str(20), "krawczyk_norm_bound": nrm.str(10)},
        "all_pass": all(results)}
    with open(os.path.join(CERT, "certificate_output.json"), "w") as fh:
        json.dump(out, fh, indent=1)
    return all(results)


if __name__ == "__main__":
    ok = main("--write" in sys.argv)
    print("ALL PASS" if ok else "SOME CHECKS FAILED")
    sys.exit(0 if ok else 1)
