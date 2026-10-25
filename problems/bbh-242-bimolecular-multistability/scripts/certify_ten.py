#!/usr/bin/env python3
"""Exact certificates that each of the 10 networks of BBH (2024, Table 1, orange bullets) has, for explicit
rational rate constants, two asymptotically stable equilibria in the closed quadrant: the origin and a positive
equilibrium (the positive equilibrium is placed at z* = (1,1) by construction).

For each network the script checks in exact rational arithmetic (sympy):
  E1  kappa > 0 and F(1,1) = 0 (z* = (1,1) is a positive equilibrium);
  E2  J(z*) has det > 0 and trace < 0 (both eigenvalues have negative real part: z* is asymptotically stable);
  E3  F(0,0) = 0 and the origin is asymptotically stable relative to the closed quadrant:
        - either J(0) has det > 0 and trace < 0 (hyperbolic sink), or
        - (networks with 2X -> 0, where J(0) is singular) the linear function V = x + 3y satisfies
          V'(z) <= y*(l2 + |E| x + max(F,0) y) + D x^2 with l2 < 0, D < 0, so V' < 0 on
          {z >= 0, z != 0, |E| x + max(F,0) y < -l2}: a strict Lyapunov function on a neighbourhood of 0 in
          the quadrant (sublevel sets {V <= r} are compact triangles);
  E4  (illustration) the exact number of positive equilibria for these rate constants is 2 (resultant +
      Sturm count) and the other one is a saddle (det J < 0, checked on an isolating interval).

Runtime: ~20 s.  Environment: Python 3.12, sympy 1.14 (WSL Fedora 44).  Exit status 1 on any FAIL.
The rate constants are written to ../certificates/ten_networks_rates.txt.
"""
import sys
import os
import sympy as sp
from networks import *

x, y, t = sp.symbols('x y t')
FAILS = []


def check(cond, msg):
    print(("PASS " if cond else "FAIL ") + msg)
    if not cond:
        FAILS.append(msg)


def field(net, k):
    F = sp.Matrix([0, 0])
    for kj, r in zip(k, net):
        s = r[0]
        v = vec(r)
        F += kj * x**s[0] * y**s[1] * sp.Matrix(v)
    return sp.expand(F)


def kernel_rays(net):
    from math import gcd
    G = sp.Matrix([[vec(r)[0] for r in net], [vec(r)[1] for r in net]])
    b1, b2 = G.nullspace()
    cands = set()
    for i in range(4):
        a, b = b1[i], b2[i]
        if a == 0 and b == 0:
            continue
        w = b * b1 - a * b2
        w = w * sp.lcm([sp.fraction(q)[1] for q in w])
        wi = [int(q) for q in w]
        g = 0
        for q in wi:
            g = gcd(g, abs(q))
        wi = [q // g for q in wi]
        for s in (wi, [-q for q in wi]):
            if all(q >= 0 for q in s):
                cands.add(tuple(s))
    return sorted(cands)


def choose_rates(net):
    """kappa = h(t) for a rational t with det J(1,1) > 0 (rates on the stable side of the fold)."""
    u, v = kernel_rays(net)
    for tt in [sp.Rational(i, 64) for i in range(1, 64)]:
        k = [(1 - tt) * u[i] + tt * v[i] for i in range(4)]
        if any(q <= 0 for q in k):
            continue
        F = field(net, k)
        J = F.jacobian([x, y]).subs({x: 1, y: 1})
        if J.det() > 0 and J.trace() < 0:
            J0 = F.jacobian([x, y]).subs({x: 0, y: 0})
            return k, tt
    return None, None


def main():
    lines = []
    for s in ORANGE10:
        net = parse_net(s)
        print("== " + s)
        k, tt = choose_rates(net)
        check(k is not None, "rational rate constants found (t = %s)" % tt)
        if k is None:
            continue
        # make integer rates (scaling kappa by a positive constant only rescales time)
        L = sp.lcm([sp.fraction(q)[1] for q in k])
        k = [q * L for q in k]
        kd = {r: q for r, q in zip(net, k)}
        listed = [parse_net(part)[0] for part in s.split(",")]
        lines.append(", ".join(rname(r) for r in listed) + "  |  kappa = " + ", ".join(str(kd[r]) for r in listed))
        F = field(net, k)
        print("   xdot =", F[0], "   ydot =", F[1])
        J = F.jacobian([x, y])
        check(all(q > 0 for q in k) and F.subs({x: 1, y: 1}) == sp.Matrix([0, 0]), "E1 (1,1) is an equilibrium, kappa > 0")
        J1 = J.subs({x: 1, y: 1})
        check(J1.det() > 0 and J1.trace() < 0, "E2 det J(1,1) = %s > 0, tr J(1,1) = %s < 0" % (J1.det(), J1.trace()))
        J0 = J.subs({x: 0, y: 0})
        check(F.subs({x: 0, y: 0}) == sp.Matrix([0, 0]), "E3a origin is an equilibrium")
        if J0.det() > 0 and J0.trace() < 0:
            check(True, "E3b origin hyperbolic sink: det J(0) = %s, tr J(0) = %s" % (J0.det(), J0.trace()))
        else:
            V = F[0] + 3 * F[1]          # derivative of V = x + 3y along the flow
            P = sp.Poly(V, x, y)
            l1 = P.coeff_monomial(x); l2 = P.coeff_monomial(y)
            D = P.coeff_monomial(x**2); E = P.coeff_monomial(x * y); Fc = P.coeff_monomial(y**2)
            ok = (P.coeff_monomial(1) == 0 and l1 == 0 and l2 < 0 and D < 0)
            check(ok, "E3c V = x+3y: V' = %s*y + %s*x^2 + %s*x*y + %s*y^2 with l2 < 0, D < 0, no x term" % (l2, D, E, Fc))
        # E4: positive equilibria <-> roots t in (0,1) of prod h_i(t)^c_i = K := prod kappa_i^c_i (BBH Lemma 11)
        u, v = kernel_rays(net)
        A1 = sp.Matrix([[r[0][0], r[0][1], 1] for r in net])
        c = A1.T.nullspace()[0]
        c = c * sp.lcm([sp.fraction(q)[1] for q in c])
        h = [(1 - t) * u[i] + t * v[i] for i in range(4)]
        K = sp.prod([sp.Rational(k[i]) ** c[i] for i in range(4)])
        lhs = sp.prod([h[i] ** c[i] for i in range(4) if c[i] > 0])
        rhs = sp.prod([h[i] ** (-c[i]) for i in range(4) if c[i] < 0])
        P = sp.Poly(sp.expand(lhs - K * rhs), t)
        Q = sp.Poly(sp.expand((sp.Matrix([[vec(r)[0] for r in net], [vec(r)[1] for r in net]]) * sp.diag(*h)
                               * sp.Matrix([[r[0][0], r[0][1]] for r in net])).det()), t)
        ivs = [(a, b) for (a, b), m in P.intervals(eps=sp.Rational(1, 10**9)) if a >= 0 and b <= 1 and not (a == b == 0) and not (a == b == 1)]
        signs = []
        for a, b in ivs:
            if a == b:
                signs.append(sp.sign(Q.eval(a)))
            else:
                assert Q.count_roots(a, b) == 0
                signs.append(sp.sign(Q.eval((a + b) / 2)))
        check(len(ivs) == 2 and sorted(signs) == [-1, 1],
              "E4 exactly two positive equilibria (roots of a degree-%d polynomial in t on (0,1)), det J signs %s" % (P.degree(), signs))
    out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "certificates", "ten_networks_rates.txt")
    with open(out, "w") as f:
        f.write("# Integer rate constants (in the order of the reactions as listed) for which the origin and (1,1)\n")
        f.write("# are both asymptotically stable (checked exactly by scripts/certify_ten.py)\n")
        for l in lines:
            f.write(l + "\n")
    print("SUMMARY:", "ALL PASS" if not FAILS else "FAILURES: %d" % len(FAILS))
    sys.exit(1 if FAILS else 0)


if __name__ == "__main__":
    main()
