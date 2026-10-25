"""check_rook_hubs.py -- Proposition 6.1: the rook graph K_a x K_b with one hub per row.

R(a,b): vertices (i,j) (i in [a], j in [b]) and hubs h_i; (i,j) ~ (i,j') (rows are cliques K_b), (i,j) ~ (i',j)
(columns are cliques K_a), h_i ~ (i,j) for all j.  n = a(b+1).  The ab rook vertices form one orbit of
Aut(R(a,b)) >= S_a x S_b.  We derive Delta(a,b) = Kf(R) - Kf(R - (1,1)) symbolically from the isotypic
decomposition of R^{V} under S_a x S_b (resp. of R^{V - v} under S_{a-1} x S_{b-1}):
  tr L(R)^+    = 1/(b+1) [trivial, nonzero part] + (a-1) tr M_A^{-1} + (b-1)/(b+1) + (a-1)(b-1)/(a+b+1),
                 M_A = [[a+1, -1], [-b, b]];
  tr L(R-v)^+  = (sum of 1/lambda over the nonzero eigenvalues of the 5x5 quotient Q) + (a-2) tr M_A'^{-1}
                 + (b-2) tr M_B'^{-1} + (a-2)(b-2)/(a+b+1),
  with Q, M_A' (3x3), M_B' (2x2) as written below (cells: row-mates, column-mates, others, own hub, other hubs).
Claims (exact):
  (S1) the symbolic Kf(R) and Kf(R - v) agree with brute-force exact Kf for 12 pairs (a,b) with n <= 30;
  (S2) Delta = -P(a,b)/Den(a,b) with Den > 0 for a, b >= 2 and P the printed integer polynomial (irreducible over Q);
  (S3) for every 3 <= b <= 400, P(., b) has exactly one real root a*(b) > 1.5, and Delta changes sign between
       floor(a*) and ceil(a*) (exact evaluation): the rook class is good on a real curve;
  (S4) P(a, b) != 0 for all integers a >= 2 and 2 <= b <= 400 (no linear factor of P(., b) over Q has an integer
       root >= 2): no member of the family with b <= 400 has its rook vertices good;
  (S5) a*(b)/b^3 -> 1/2 (printed for b = 10, 50, 100, 200, 400): along the branch the proportion b/(b+1) -> 1.
Runtime: about 3 min (Fedora 44 / WSL, sympy 1.14, python-flint 0.6).
"""
import sympy as sp
import networkx as nx
from fractions import Fraction
from flint import fmpq
from kfexact import kf, minus, Checker

C = Checker()
a, b, x = sp.symbols('a b x')


def sum_inv_nonzero(Q):
    """sum of 1/lambda over the nonzero eigenvalues of a matrix with a simple zero eigenvalue."""
    p = sp.Poly(sp.expand((x * sp.eye(Q.shape[0]) - Q).det()), x)
    c = p.all_coeffs()
    return sp.cancel(-c[-3] / c[-2])


n = a * (b + 1)
MA = sp.Matrix([[a + 1, -1], [-b, b]])
trG = sp.Rational(1) / (b + 1) + (a - 1) * sp.cancel(MA.inv().trace()) + (b - 1) / (b + 1) + (a - 1) * (b - 1) / (a + b + 1)
KfG = sp.cancel(n * trG)
Q = sp.Matrix([[a, 0, -(a - 1), -1, 0], [0, b, -(b - 1), 0, -1], [-1, -1, 3, 0, -1], [-(b - 1), 0, 0, b - 1, 0],
               [0, -1, -(b - 1), 0, b]])
MA2 = sp.Matrix([[a + b - 1, -(b - 1), -1], [-1, a + 2, -1], [-1, -(b - 1), b]])
MB2 = sp.Matrix([[a + b - 1, -(a - 1)], [-1, b + 2]])
trGv = sum_inv_nonzero(Q) + (a - 2) * sp.cancel(MA2.inv().trace()) + (b - 2) * sp.cancel(MB2.inv().trace()) \
    + (a - 2) * (b - 2) / (a + b + 1)
KfGv = sp.cancel((n - 1) * trGv)
Delta = sp.factor(sp.cancel(KfG - KfGv))


def rook(aa, bb):
    G = nx.Graph()
    for i in range(aa):
        for j in range(bb):
            G.add_edge(('r', i, j), ('h', i))
            for jj in range(j + 1, bb):
                G.add_edge(('r', i, j), ('r', i, jj))
            for ii in range(i + 1, aa):
                G.add_edge(('r', i, j), ('r', ii, j))
    return G


def to_fmpq(e):
    r = sp.Rational(e)
    return fmpq(int(r.p), int(r.q))


ok = True
pairs = [(2, 2), (2, 3), (3, 2), (3, 3), (2, 4), (4, 2), (3, 4), (4, 3), (2, 5), (5, 2), (2, 6), (2, 9)]
for aa, bb in pairs:
    G = rook(aa, bb)
    k1, k2 = kf(G), kf(minus(G, ('r', 0, 0)))
    s1, s2 = to_fmpq(KfG.subs({a: aa, b: bb})), to_fmpq(KfGv.subs({a: aa, b: bb}))
    if not (k1 == s1 and k2 == s2):
        ok = False
        print('   mismatch', aa, bb, k1, s1, k2, s2)
C.check(ok, '(S1) symbolic Kf(R) and Kf(R - v) = brute-force exact values for %d pairs (a,b)' % len(pairs))

num, den = sp.fraction(sp.together(Delta))
P = sp.Poly(sp.expand(-num), a, b)
Den = sp.factor(den)
print('   Delta = -P(a,b) / Den(a,b), Den =', Den)
print('   P(a,b) =', P.as_expr())
den_factors = [f for f, _ in sp.factor_list(den)[1]]
pos = all(sp.Poly(f, a, b).eval({a: 2, b: 2}) > 0 and all(c > 0 for c in sp.Poly(sp.expand(f.subs({a: a + 2, b: b + 2})), a, b).coeffs())
          for f in den_factors)
irr = len(sp.factor_list(P.as_expr())[1]) == 1 and sp.factor_list(P.as_expr())[1][0][1] == 1
C.check(pos and irr, '(S2) every factor of Den has positive coefficients after a,b -> a+2,b+2 (Den > 0); P irreducible over Q')

ok3 = True
ok4 = True
ratios = {}
for bb in range(2, 401):
    pa = sp.Poly(P.as_expr().subs(b, bb), a)
    # integer roots >= 2: linear factors over Q
    for f, _ in sp.factor_list(pa.as_expr())[1]:
        fp = sp.Poly(f, a)
        if fp.degree() == 1:
            r = -fp.all_coeffs()[1] / fp.all_coeffs()[0]
            if r.is_integer and r >= 2:
                ok4 = False
                print('   integer root', r, bb)
    if bb >= 3:
        rts = [r for r in pa.real_roots() if r > sp.Rational(3, 2)]
        if len(rts) != 1:
            ok3 = False
            print('   b =', bb, 'real roots > 1.5:', rts)
            continue
        r0 = rts[0]
        lo = int(sp.floor(r0))
        hi = lo + 1
        dlo = Delta.subs({a: lo, b: bb})
        dhi = Delta.subs({a: hi, b: bb})
        if not (sp.sign(dlo) * sp.sign(dhi) < 0):
            ok3 = False
        ratios[bb] = float(r0) / bb ** 3
C.check(ok3, '(S3) for 3 <= b <= 400: one real root a*(b) > 1.5 of P(., b), and Delta changes sign between floor and ceil')
C.check(ok4, '(S4) no integer a >= 2 with P(a, b) = 0 for 2 <= b <= 400')
print('   a*(b)/b^3:', ', '.join('b=%d: %.6f' % (bb, ratios[bb]) for bb in (3, 10, 50, 100, 200, 400)))
# With u=1/b and z=a/b^3, the polynomial below is u^18 P(z/u^3,1/u).
# Its simple zero at (z,u)=(1/2,0) yields a real analytic zero branch
# z(u) -> 1/2 by the implicit function theorem. Finite samples do not prove a limit.
z, u = sp.symbols('z u')
weight = max(3*ia+ib for (ia, ib), coeff in P.terms())
scaled = sp.Poly(sum(coeff*z**ia*u**(weight-3*ia-ib)
                     for (ia, ib), coeff in P.terms()), z, u).as_expr()
limit = sp.expand(scaled.subs(u, 0))
simple = sp.diff(limit, z).subs(z, sp.Rational(1, 2))
C.check(weight == 18 and sp.expand(limit-z**4*(2*z-1)) == 0 and simple == sp.Rational(1, 8),
        '(S5) u^18 P(z/u^3,1/u) at u=0 equals z^4(2z-1), derivative 1/8 at z=1/2; a real zero branch has a/b^3 -> 1/2')
C.finish()
