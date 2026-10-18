"""
uniqueness.py -- the certified zero mu* is the ONLY point of the whole admissible parameter domain
                 where tr J = L1 = L2 = L3 = 0  (bonus result; exact elimination + Arb balls).

Environment: Python 3.12, python-flint 0.9.0.  Runtime: a few seconds.  Requires certificates/G_k.json
and certificates/box.json (written by certify.py --write).

Admissible domain (Boros-Hofbauer Sec. 5.1 with a positive equilibrium): a > 0, b > -1, d > 1, 1+bd > 0.
In the variables (a, c, s), c = 1+b, s = (d-1)/(1+bd), it is  a > 0, c > 0, s > 0, rho := 1-(c-1)s > 0.

Steps (each printed with PASS/FAIL):
 1. On the domain, G_1 = 0  <=>  rho = 1/(t (2t(1+s)+1)) with t := a/c > 0.  Hence
        (t, s) -> (a, c, s) = (t c, c, s),   c = (1+s-rho)/s
    is a bijection from {t > 0, s > 0, 2t(1+s) > 1} onto {G_1 = 0} (in the domain).
    [check: G_1(t c, c, s) == 0 identically after the substitution.]
 2. H_k(t,s) := G_k(tc, c, s) * (s t Q)^m  (Q = 2t(1+s)+1) is a polynomial; we factor it as
    (product of factors in {s, s+1, t, ts+t+1, 2ts+2t-1}, all nonzero on the domain) * h_k(t,s).
 3. Any common zero (t0, s0) of h_2, h_3 satisfies R_t(t0) = 0 and R_s(s0) = 0 where
    R_t = Res_s(h_2, h_3), R_s = Res_t(h_2, h_3) (computed exactly).
 4. All complex roots of R_t, R_s are enclosed in Arb balls (fmpz_poly.complex_roots).  For every pair of
    roots with a real part that may be positive and an imaginary part that may be 0, evaluate h_2, h_3 on the
    pair of balls; pairs where an enclosure excludes 0 are discarded.
 5. Every remaining pair is shown to lie inside the Krawczyk box of certify.py (mapped to (a,c,s)), where the
    zero is unique.  Hence {tr J = L1 = L2 = L3 = 0} in the domain is exactly {mu*}.
"""
import sys
import os
import json
from flint import fmpq_mpoly_ctx, Ordering, fmpq, fmpz_poly, fmpq_poly, acb, arb, ctx

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from certify import poly_from_json  # noqa: E402

CERT = os.path.join(HERE, "..", "certificates")
T = fmpq_mpoly_ctx.get(('t', 's'), Ordering.degrevlex)
t, s = T.gens()
results = []


def report(name, ok, info=""):
    results.append(bool(ok))
    print(("PASS" if ok else "FAIL") + "  " + name + (("  [" + info + "]") if info else ""), flush=True)


def subst(P):
    """G(a,c,s) -> G(tc, c, s) * (s t Q)^dmax with c = Nc/(s t Q)."""
    Q = 2 * t * (1 + s) + 1
    Nc = t * Q * (1 + s) - 1
    Dc = s * t * Q
    terms = P.to_dict()
    dmax = max(int(e[0]) + int(e[1]) for e in terms)
    Np = [T.from_dict({(0, 0): 1})]
    Dp = [T.from_dict({(0, 0): 1})]
    for _ in range(dmax):
        Np.append(Np[-1] * Nc)
        Dp.append(Dp[-1] * Dc)
    tot = T.from_dict({})
    for e, cf in terms.items():
        i, j, k = int(e[0]), int(e[1]), int(e[2])
        tot += cf * t**i * Np[i + j] * Dp[dmax - i - j] * s**k
    return tot


def univariate(P, var):
    """P in Q[t,s] depending on one variable only -> fmpq_poly."""
    idx = 0 if var == 't' else 1
    coeffs = {}
    for e, cf in P.to_dict().items():
        if int(e[1 - idx]) != 0:
            raise ValueError("not univariate")
        coeffs[int(e[idx])] = cf
    deg = max(coeffs)
    return fmpq_poly([coeffs.get(i, 0) for i in range(deg + 1)])


def ev2(P, tb, sb):
    tot = acb(0)
    for e, cf in P.to_dict().items():
        tot += acb(arb(cf)) * tb**int(e[0]) * sb**int(e[1])
    return tot


def roots_of(P, var):
    up = univariate(P, var)
    num = fmpz_poly([int(c) for c in (up * up.denom()).coeffs()])
    out = []
    for fac, m in num.factor()[1]:
        for r, mult in fac.complex_roots():
            out.append((r, str(fac)[:60]))
    return out


def main():
    ctx.prec = 400
    G = [poly_from_json(json.load(open(os.path.join(CERT, "G_%d.json" % k)))) for k in range(1, 5)]
    H = [subst(P) for P in G]
    report("step 1: G_1 vanishes identically on the parametrisation of {G_1 = 0}", H[0] == 0)
    allowed = {str(s), str(s + 1), str(t), str(t * s + t + 1), str(2 * t * s + 2 * t - 1)}
    h = {}
    for k in (1, 2):
        const, facs = H[k].factor()
        big = [(q, m) for q, m in facs if str(q) not in allowed]
        ok = len(big) == 1 and big[0][1] == 1
        h[k] = big[0][0]
        report("step 2: H_%d = const * (domain-nonvanishing factors) * h_%d" % (k + 1, k + 1), ok,
               "h_%d degrees (t,s) = %s, %d terms" % (k + 1, h[k].degrees(), len(h[k])))
    h2, h3 = h[1], h[2]
    Rt = h2.resultant(h3, 's')
    Rs = h2.resultant(h3, 't')
    print("     deg R_t = %s, deg R_s = %s" % (Rt.degrees(), Rs.degrees()))
    rt = roots_of(Rt, 't')
    rs = roots_of(Rs, 's')
    zero = arb(0)
    cand_t = [(r, f) for r, f in rt if not (r.real < zero) and r.imag.contains(0) and not (r.real == zero)]
    cand_s = [(r, f) for r, f in rs if not (r.real < zero) and r.imag.contains(0) and not (r.real == zero)]
    # exact rational roots t = 0 / s = 0 are on the boundary (t > 0, s > 0 strictly required)
    print("     candidate t-roots (Re possibly > 0, Im possibly 0):", [r.real.str(8) for r, f in cand_t])
    print("     candidate s-roots (Re possibly > 0, Im possibly 0):", [r.real.str(8) for r, f in cand_s])
    surviving = []
    for tr_, ft in cand_t:
        for sr_, fs in cand_s:
            v2 = ev2(h2, tr_, sr_)
            v3 = ev2(h3, tr_, sr_)
            if v2.contains(0) and v3.contains(0):
                surviving.append((tr_, sr_))
    report("step 4: exactly one candidate pair (t, s) survives the interval exclusion", len(surviving) == 1,
           "surviving: %s" % [(x.real.str(15), y.real.str(15)) for x, y in surviving])
    # step 5: the survivor lies in the Krawczyk box
    box = json.load(open(os.path.join(CERT, "box.json")))

    def q(sv):
        num, _, den = sv.partition('/')
        return fmpq(int(num), int(den) if den else 1)
    m = [q(v) for v in box["centre"]]
    r = q(box["radius"])
    inside = False
    if len(surviving) == 1:
        tb, sb = surviving[0]
        tb, sb = tb.real, sb.real
        Qb = 2 * tb * (1 + sb) + 1
        rho = 1 / (tb * Qb)
        cb = (1 + sb - rho) / sb
        ab = tb * cb
        X = [arb(mi, r) for mi in m]
        inside = X[0].contains(ab) and X[1].contains(cb) and X[2].contains(sb)
        print("     survivor in (a,c,s): a in %s, c in %s, s in %s" % (ab.str(20), cb.str(20), sb.str(20)))
    report("step 5: the surviving pair lies in the Krawczyk box (zero unique there)", inside)
    return all(results)


if __name__ == "__main__":
    ok = main()
    print("ALL PASS" if ok else "SOME CHECKS FAILED")
    sys.exit(0 if ok else 1)
