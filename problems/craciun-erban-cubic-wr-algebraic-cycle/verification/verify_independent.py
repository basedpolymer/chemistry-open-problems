#!/usr/bin/env python3
"""Independent verification (written by the verifier, without the solver's code or reasoning).

Claim: the 10-reaction network N10 below lies in W_3 of Craciun & Erban (J. Math. Biol. 2025, arXiv:2406.05057,
Definitions 1-3) and its mass-action system has the circle (x-2)^2 + (y-2)^2 = 1 as a hyperbolic stable limit cycle.
The same checks are run on the 17-reaction network N17 quoted in the manuscript.

Methods (deliberately different from scripts/certify_cycle.py, which uses a resultant and Arb integration):
  A. Network membership in W_3, checked against the source definitions: every reactant AND product complex has
     molecularity <= 3; no reaction is trivial (source != product) or repeated; both species occur; rate constants
     positive; weak reversibility checked reaction by reaction (edge u -> v lies on a directed cycle iff u is
     reachable from v).
  B. Mass-action vector field rebuilt from the reaction list; compared with the polynomials printed in the paper.
  C. Invariance of the circle: (i) multivariate division of f h_x + g h_y by h, (ii) exact substitution of the rational
     parametrisation x = 2 + (1-t^2)/(1+t^2), y = 2 + 2t/(1+t^2) into f h_x + g h_y.
  D. No equilibrium on the circle: Q(t) = (1+t^2)^4 P(x(t), y(t)), P = -(y-2) f + (x-2) g the tangential speed, has
     no real root (exact Sturm count over all of R) and P(1,2) > 0 (the point t = infinity); also an exact lower
     bound P >= c on the circle by the same method applied to P - c.
  E. Characteristic exponent Lambda = oint div F dt = int_{-pi}^{pi} divF/P d theta:
     E1. rigorous enclosure with mpmath interval arithmetic (library independent of Arb): composite midpoint rule
         with the remainder term h^3/24 * f''(I_k) enclosed by interval evaluation of the symbolic second derivative;
     E2. closed form by residues: Lambda = int_R 2 D(t)(1+t^2)/Q(t) dt = 2 pi i * sum of residues at the poles
         in the upper half plane after cancellation (Q = 4 (1+t^2)^2 q4(t) for N10), poles isolated exactly by sympy
         CRootOf and evaluated to 40 digits.
     The period T = oint dt is computed the same two ways.
  F. Non-rigorous sanity check: scipy DOP853 integration of the ODE and of its variational equation over one period
     from (3, 2): return to the start point, and det(monodromy) = exp(Lambda) (Liouville).
  G. Side claim of the manuscript (Remark "construction"): on the 10-reaction support, the solutions of the linear
     system "f h_x + g h_y = s h" in the unknown rate constants and cofactor coefficients form a one-dimensional space.

Environment: Python 3.12, sympy 1.14, mpmath 1.3, scipy 1.18 (WSL Fedora). Single thread. Runtime ~ 2-4 min.
Prints PASS/FAIL per claim; exits 1 on any FAIL.
"""
import sys
import sympy as sp
import mpmath
from mpmath import iv

x, y, t, th = sp.symbols("x y t theta", real=True)
Rat = sp.Rational
fails = 0


def check(cond, msg):
    global fails
    print(("PASS " if cond else "FAIL ") + msg, flush=True)
    fails += (not cond)


def parse_complex(c):
    """'2X+Y' -> (2, 1); '0' -> (0, 0)."""
    c = c.replace(" ", "")
    if c == "0":
        return (0, 0)
    a = b = 0
    for term in c.split("+"):
        coef = term[:-1]
        k = int(coef) if coef else 1
        if term[-1] == "X":
            a += k
        elif term[-1] == "Y":
            b += k
        else:
            raise ValueError(term)
    return (a, b)


# Reaction lists typed in from the task statement / manuscript (not copied from the solver's script).
N10_TXT = """
0 -> 2Y : 434
0 -> X+Y : 448
Y -> 0 : 3374
2Y -> 3Y : 1728
3Y -> Y : 194
3Y -> 2X+Y : 11
X+Y -> X+2Y : 328
X+2Y -> 3Y : 250
2X+Y -> 3X : 216
3X -> X+2Y : 28
"""
N17_TXT = """
0 -> X : 801/1000000
Y -> 0 : 104261883/56000000
2Y -> 3Y : 14893029/14000000
3Y -> 0 : 14847677/168000000
3Y -> Y : 211/1000000
3Y -> 2Y : 213/500000
X -> 2X : 1/560000
X -> 3X : 115777/250000
X+Y -> Y : 37358609/28000000
X+Y -> 2X : 748283/560000
X+2Y -> 3Y : 213551/400000
2X -> 3X : 7497783/5600000
2X+Y -> 3X : 7477489/56000000
3X -> 3Y : 11166451/84000000
3X -> X+Y : 59813/56000000
3X -> X+2Y : 107/500000
3X -> 2X+Y : 53/1000000
"""
# Polynomials printed in the manuscript (Theorem 1) for N10.
PAPER_F = -56*x**3 + 216*x**2*y - 250*x*y**2 + 22*y**3 + 448
PAPER_G = 56*x**3 - 216*x**2*y + 250*x*y**2 + 328*x*y - 410*y**3 + 1728*y**2 - 3374*y + 1316
PAPER_S = -112*x**2 + 544*x*y - 448*x - 820*y**2 + 1728*y - 1008


def parse(txt):
    out = []
    for line in txt.strip().splitlines():
        lhs, k = line.split(":")
        src, prod = lhs.split("->")
        out.append((parse_complex(src), parse_complex(prod), sp.nsimplify(k.strip(), rational=True)))
    return out


# ---------------------------------------------------------------------------------------------------- interval tools
def iv_eval(expr, env):
    """Evaluate a sympy expression tree with mpmath interval arithmetic (env maps symbols to iv.mpf)."""
    if expr.is_Symbol:
        return env[expr]
    if expr.is_Integer:
        return iv.mpf(int(expr))
    if expr.is_Rational:
        return iv.mpf(int(expr.p)) / int(expr.q)
    if expr is sp.pi:
        return iv.pi
    if expr.is_Add:
        r = iv.mpf(0)
        for a in expr.args:
            r = r + iv_eval(a, env)
        return r
    if expr.is_Mul:
        r = iv.mpf(1)
        for a in expr.args:
            r = r * iv_eval(a, env)
        return r
    if expr.is_Pow:
        b, e = expr.args
        if e.is_Integer:
            bb = iv_eval(b, env)
            n = int(e)
            if n >= 0:
                return bb ** n
            return iv.mpf(1) / (bb ** (-n))
        raise ValueError(expr)
    if isinstance(expr, sp.cos):
        return iv.cos(iv_eval(expr.args[0], env))
    if isinstance(expr, sp.sin):
        return iv.sin(iv_eval(expr.args[0], env))
    raise ValueError(f"unsupported node {expr.func}")


def rigorous_integral(fexpr, nsub):
    """Enclosure of int_{-pi}^{pi} fexpr(theta) d theta by the composite midpoint rule with interval remainder."""
    f2 = sp.diff(fexpr, th, 2)
    two_pi = 2 * iv.pi
    hk = two_pi / nsub
    total = iv.mpf(0)
    for k in range(nsub):
        lo = -iv.pi + hk * k
        mid = lo + hk / 2
        Ik = lo + hk * iv.mpf([0, 1])
        # exact identity on each cell: int f = h f(m) + h^3/24 f''(xi), xi in cell
        total = total + hk * iv_eval(fexpr, {th: mid}) + hk ** 3 / 24 * iv_eval(f2, {th: Ik})
    return total


# ---------------------------------------------------------------------------------------------------- main checks
def verify(label, txt, paper=None, nsub=2048, ode=True):
    print(f"========== {label} ==========", flush=True)
    net = parse(txt)
    print(f"{len(net)} reactions")
    # ---- A. membership in W_3
    cx = sorted({c for s_, p_, _ in net for c in (s_, p_)})
    check(all(k > 0 for *_, k in net), "A1 all rate constants > 0")
    check(all(sum(s_) <= 3 and sum(p_) <= 3 for s_, p_, _ in net),
          "A2 every reactant and every product complex has molecularity <= 3 (Definition 2, M_3)")
    check(max(sum(s_) for s_, _, _ in net) == 3, "A3 order (2.2) = max reactant molecularity = 3 (cubic)")
    check(all(s_ != p_ for s_, p_, _ in net), "A4 no trivial reaction (source != product)")
    check(len({(s_, p_) for s_, p_, _ in net}) == len(net), "A5 no repeated reaction")
    check(any(s_[0] != p_[0] for s_, p_, _ in net) and any(s_[1] != p_[1] for s_, p_, _ in net),
          "A6 both species change in some reaction")
    succ = {c: set() for c in cx}
    for s_, p_, _ in net:
        succ[s_].add(p_)

    def reach(u):
        seen, st = {u}, [u]
        while st:
            a = st.pop()
            for b in succ[a]:
                if b not in seen:
                    seen.add(b)
                    st.append(b)
        return seen

    on_cycle = [s_ in reach(p_) for s_, p_, _ in net]
    check(all(on_cycle), "A7 weak reversibility: each reaction u->v lies on a directed cycle (u reachable from v)")
    # linkage classes and deficiency
    und = {c: set() for c in cx}
    for s_, p_, _ in net:
        und[s_].add(p_)
        und[p_].add(s_)
    seen, ncl = set(), 0
    for c in cx:
        if c not in seen:
            ncl += 1
            st = [c]
            seen.add(c)
            while st:
                a = st.pop()
                for b in und[a]:
                    if b not in seen:
                        seen.add(b)
                        st.append(b)
    rank = sp.Matrix([[p_[0] - s_[0], p_[1] - s_[1]] for s_, p_, _ in net]).rank()
    print(f"   complexes = {len(cx)}, linkage classes = {ncl}, rank = {rank}, deficiency = {len(cx) - ncl - rank}")

    # ---- B. mass-action field
    f = sp.expand(sum(k * (p_[0] - s_[0]) * x**s_[0] * y**s_[1] for s_, p_, k in net))
    g = sp.expand(sum(k * (p_[1] - s_[1]) * x**s_[0] * y**s_[1] for s_, p_, k in net))
    print("   f =", f)
    print("   g =", g)
    if paper is not None:
        check(sp.expand(f - paper[0]) == 0 and sp.expand(g - paper[1]) == 0,
              "B1 mass-action field equals the polynomials printed in the manuscript")
    check(sp.Poly(f, x, y).total_degree() <= 3 and sp.Poly(g, x, y).total_degree() <= 3, "B2 field has degree <= 3")

    # ---- C. invariance
    h = (x - 2)**2 + (y - 2)**2 - 1
    L = sp.expand(f * sp.diff(h, x) + g * sp.diff(h, y))
    q, r = sp.reduced(L, [h], x, y)
    s = sp.expand(q[0])
    check(r == 0 and sp.expand(L - s * h) == 0, "C1 f h_x + g h_y = s h exactly")
    print("   cofactor s =", s)
    if paper is not None:
        check(sp.expand(s - paper[2]) == 0, "C2 cofactor equals the one printed in the manuscript")
    xt = 2 + (1 - t**2) / (1 + t**2)
    yt = 2 + 2 * t / (1 + t**2)
    check(sp.simplify(sp.together(L.subs({x: xt, y: yt}))) == 0,
          "C3 f h_x + g h_y vanishes identically on the rational parametrisation of the circle")

    # ---- D. tangential speed
    P = sp.expand(-(y - 2) * f + (x - 2) * g)
    Q = sp.Poly(sp.cancel(sp.together(P.subs({x: xt, y: yt}) * (1 + t**2)**4)), t)
    P12 = P.subs({x: 1, y: 2})
    check(Q.count_roots() == 0, f"D1 Q(t) = (1+t^2)^4 P(circle(t)) (degree {Q.degree()}) has no real root (Sturm)")
    check(P12 > 0 and Q.LC() == P12, f"D2 P(1,2) = {P12} > 0 (point t = oo); hence P > 0 on the whole circle")
    # exact lower bound
    num = sp.lambdify(th, P.subs({x: 2 + sp.cos(th), y: 2 + sp.sin(th)}), "mpmath")
    approx_min = min(num(mpmath.mpf(2) * mpmath.pi * i / 20000) for i in range(20000))
    c = sp.Rational(int(mpmath.floor(approx_min * 100)), 100) - Rat(1, 100)
    Qc = sp.Poly(sp.cancel(sp.together((P - c).subs({x: xt, y: yt}) * (1 + t**2)**4)), t)
    check(Qc.count_roots() == 0 and P12 - c > 0, f"D3 exact bound: P > {c} on the circle (sampled min {float(approx_min):.6f})")
    print("   factorisation of Q:", sp.factor_list(Q.as_expr()))

    # ---- E. integrals
    divF = sp.expand(sp.diff(f, x) + sp.diff(g, y))
    print("   div F =", divF)
    sub_th = {x: 2 + sp.cos(th), y: 2 + sp.sin(th)}
    lam_expr = divF.subs(sub_th) / P.subs(sub_th)
    per_expr = 1 / P.subs(sub_th)
    iv.prec = 100
    lam_iv = rigorous_integral(lam_expr, nsub)
    per_iv = rigorous_integral(per_expr, nsub)
    print(f"   E1 mpmath-iv  Lambda in {lam_iv}")
    print(f"   E1 mpmath-iv  T      in {per_iv}")
    check(lam_iv.b < 0, "E1 rigorous interval enclosure: oint div F dt < 0")
    check(per_iv.a > 0, "E1 rigorous interval enclosure: period T > 0")
    # residues
    D = sp.Poly(sp.cancel(sp.together(divF.subs({x: xt, y: yt}) * (1 + t**2)**2)), t)
    A_lam = sp.Poly(2 * D.as_expr() * (1 + t**2), t)      # Lambda = int_R A_lam / Q dt
    A_per = sp.Poly(2 * (1 + t**2)**3, t)                 # T      = int_R A_per / Q dt
    mpmath.mp.dps = 45

    def coeffs_mp(poly):
        return [mpmath.mpf(int(sp.Rational(cc).p)) / int(sp.Rational(cc).q) for cc in poly.all_coeffs()]

    def residue_integral(A, tag):
        # int_R A/Q dt with A/Q reduced to num/den in lowest terms; if den is squarefree with no real root and
        # deg num <= deg den - 2:  integral = 2 pi i sum_{Im z > 0} num(z)/den'(z).
        num, den = sp.fraction(sp.cancel(A.as_expr() / Q.as_expr()))
        num, den = sp.Poly(num, t), sp.Poly(den, t)
        ok = (num.degree() <= den.degree() - 2 and sp.gcd(den, den.diff(t)).degree() == 0
              and den.count_roots() == 0)
        roots = den.all_roots()
        upper = [z for z in roots if sp.im(z.evalf(30)) > 0]
        ok = ok and 2 * len(upper) == den.degree()
        check(ok, f"E2 [{tag}] reduced integrand num/den: deg num = {num.degree()} <= deg den - 2 = "
                  f"{den.degree() - 2}, den squarefree, no real pole, {len(upper)} simple poles with Im > 0")
        dden = den.diff(t)
        tot = mpmath.mpc(0)
        for z in upper:
            zv = z.evalf(50)
            zz = mpmath.mpc(mpmath.mpf(str(sp.re(zv))), mpmath.mpf(str(sp.im(zv))))
            tot += mpmath.polyval(coeffs_mp(num), zz) / mpmath.polyval(coeffs_mp(dden), zz)
        return 2j * mpmath.pi * tot

    lam_res = residue_integral(A_lam, "Lambda")
    per_res = residue_integral(A_per, "period")
    print(f"   E2 residues   Lambda = {mpmath.nstr(lam_res.real, 30)}  (imag {mpmath.nstr(lam_res.imag, 3)})")
    print(f"   E2 residues   T      = {mpmath.nstr(per_res.real, 30)}  (imag {mpmath.nstr(per_res.imag, 3)})")
    check(lam_iv.a <= lam_res.real <= lam_iv.b and per_iv.a <= per_res.real <= per_iv.b,
          "E3 residue values lie inside the rigorous interval enclosures")
    print(f"   multiplier exp(Lambda) = {mpmath.nstr(mpmath.exp(lam_res.real), 20)}")

    # ---- F. ODE sanity check (floating point, not a proof)
    if ode:
        import numpy as np
        from scipy.integrate import solve_ivp
        fx, fy = sp.diff(f, x), sp.diff(f, y)
        gx, gy = sp.diff(g, x), sp.diff(g, y)
        Fn = sp.lambdify((x, y), [f, g, fx, fy, gx, gy], "math")

        def rhs(_, u):
            a, b, c_, d, e_, ff = Fn(u[0], u[1])
            M = np.array([[c_, d], [e_, ff]])
            Phi = u[2:].reshape(2, 2)
            return np.concatenate([[a, b], (M @ Phi).ravel()])

        T = float(per_res.real)
        sol = solve_ivp(rhs, (0, T), [3.0, 2.0, 1, 0, 0, 1], method="DOP853", rtol=1e-13, atol=1e-13)
        end = sol.y[:, -1]
        Mono = end[2:].reshape(2, 2)
        ev = np.linalg.eigvals(Mono)
        print(f"   F  return error |u(T)-u(0)| = {np.hypot(end[0]-3, end[1]-2):.2e}; det(monodromy) = "
              f"{np.linalg.det(Mono):.12f}; eigenvalues = {np.sort(ev.real)}")
        check(np.hypot(end[0] - 3, end[1] - 2) < 1e-8 and abs(np.linalg.det(Mono) - float(mpmath.exp(lam_res.real))) < 1e-7,
              "F  (numerical) orbit closes after T and det(monodromy) = exp(Lambda)")
        for x0 in (2.8, 3.2):
            s2 = solve_ivp(lambda _, u: Fn(u[0], u[1])[:2], (0, 40 * T), [x0, 2.0], method="DOP853", rtol=1e-11,
                           atol=1e-12)
            hv = (s2.y[0, -1] - 2) ** 2 + (s2.y[1, -1] - 2) ** 2 - 1
            print(f"   F  start ({x0}, 2): h after 40 periods = {hv:.3e}")
    return f, g, s


fN10, gN10, sN10 = verify("N10", N10_TXT, paper=(PAPER_F, PAPER_G, PAPER_S))

# ---- G. dimension of the solution space on the N10 support
print("========== G: uniqueness on the N10 support ==========")
net = parse(N10_TXT)
ks = sp.symbols("k0:%d" % len(net))
cs = sp.symbols("c0:6")
h = (x - 2)**2 + (y - 2)**2 - 1
f = sum(ks[i] * (p_[0] - s_[0]) * x**s_[0] * y**s_[1] for i, (s_, p_, _) in enumerate(net))
g = sum(ks[i] * (p_[1] - s_[1]) * x**s_[0] * y**s_[1] for i, (s_, p_, _) in enumerate(net))
sgen = cs[0] + cs[1]*x + cs[2]*y + cs[3]*x**2 + cs[4]*x*y + cs[5]*y**2
E = sp.Poly(sp.expand(f * sp.diff(h, x) + g * sp.diff(h, y) - sgen * h), x, y)
M = sp.Matrix([[sp.diff(cf, v) for v in ks + cs] for cf in E.coeffs()])
null = M.nullspace()
print(f"   linear system: {M.shape[0]} equations, {M.shape[1]} unknowns, nullity {len(null)}")
check(len(null) == 1, "G1 solution space (rates, cofactor) on the 10-reaction support is one-dimensional")
if len(null) == 1:
    v = null[0]
    v = v / v[0] * 434
    check(all(v[i] == net[i][2] for i in range(len(net))), "G2 its positive generator is the stated integer rate vector")

verify("N17", N17_TXT, paper=None, nsub=2048, ode=True)

print("ALL PASS" if not fails else f"{fails} FAILURES")
sys.exit(1 if fails else 0)
