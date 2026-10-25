"""check_rook_hubs_c.py -- Proposition 13(5): the rook graph K_a x K_b with c hubs per row (R_c(a,b)).

R_c(a,b): rook vertices (i,j) with row and column cliques, and c independent hubs per row, each adjacent to the b
vertices of its row.  n = a(b+c); the ab rook vertices form one orbit (proportion b/(b+c) if they are good).
Isotypic decomposition as in check_rook_hubs.py, plus the hub-twin eigenvalues (b with multiplicity a(c-1) in R;
b-1 with multiplicity c-1 and b with multiplicity (a-1)(c-1) in R - v):
  tr L(R)^+   = 1/(b+c) + a(c-1)/b + (a-1)(a+b+c)/(ab) + (b-1)/(b+c) + (a-1)(b-1)/(a+b+c);
  tr L(R-v)^+ = sum_{nonzero} 1/lambda(Q_c) + (c-1)/(b-1) + (a-1)(c-1)/b + (a-2) tr M_A^{-1} + (b-2) tr M_B^{-1}
                + (a-2)(b-2)/(a+b+c).
Claims (exact):
  (C1) the symbolic Kf(R_c) and Kf(R_c - v) agree with brute-force exact values for 8 triples (a,b,c), n <= 33;
  (C2) for c = 1 the formula reduces to the one of check_rook_hubs.py;
  (C3) for 3 <= b <= BMAX and 1 <= c <= b there is no integer a >= 2 with Kf(R_c) = Kf(R_c - v) (factorisation over
       Q of the numerator, degree 5 in a); in particular no R_c(a,b) with c < b/2 (proportion > 2/3) and b <= BMAX
       has good rook vertices.
  (C4) printed: the real roots a*(b,c) > 1.5 for a few (b,c), showing the real zero set for c < b/2.
Usage: python3 check_rook_hubs_c.py [BMAX]  (default 60, about 2 min; the log for BMAX = 200 took about 15 min).
"""
import sys
import sympy as sp
import networkx as nx
from sympy import Poly, ZZ
from flint import fmpq
from kfexact import kf, minus, Checker

BMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 60
C = Checker()
a, b, c, x = sp.symbols('a b c x')


def sum_inv_nonzero(Q):
    p = sp.Poly(sp.expand((x * sp.eye(Q.shape[0]) - Q).det()), x)
    co = p.all_coeffs()
    return sp.cancel(-co[-3] / co[-2])


n = a * (b + c)
MA = sp.Matrix([[a + c, -c], [-b, b]])
trG = 1 / (b + c) + a * (c - 1) / b + (a - 1) * sp.cancel(MA.inv().trace()) + (b - 1) / (b + c) \
    + (a - 1) * (b - 1) / (a + b + c)
KfG = sp.cancel(n * trG)
Q = sp.Matrix([[a + c - 1, 0, -(a - 1), -c, 0], [0, b + c - 1, -(b - 1), 0, -c], [-1, -1, c + 2, 0, -c],
               [-(b - 1), 0, 0, b - 1, 0], [0, -1, -(b - 1), 0, b]])
MA2 = sp.Matrix([[a + b + c - 2, -(b - 1), -c], [-1, a + c + 1, -c], [-1, -(b - 1), b]])
MB2 = sp.Matrix([[a + b + c - 2, -(a - 1)], [-1, b + c + 1]])
trGv = sum_inv_nonzero(Q) + (c - 1) / (b - 1) + (a - 1) * (c - 1) / b + (a - 2) * sp.cancel(MA2.inv().trace()) \
    + (b - 2) * sp.cancel(MB2.inv().trace()) + (a - 2) * (b - 2) / (a + b + c)
KfGv = sp.cancel((n - 1) * trGv)
Delta = sp.factor(sp.cancel(KfG - KfGv))


def build(aa, bb, cc):
    G = nx.Graph()
    for i in range(aa):
        for j in range(bb):
            for t in range(cc):
                G.add_edge(('r', i, j), ('h', i, t))
            for jj in range(j + 1, bb):
                G.add_edge(('r', i, j), ('r', i, jj))
            for ii in range(i + 1, aa):
                G.add_edge(('r', i, j), ('r', ii, j))
    return G


def q(e):
    r = sp.Rational(e)
    return fmpq(int(r.p), int(r.q))


ok = True
triples = [(2, 2, 1), (2, 3, 2), (3, 3, 2), (3, 4, 3), (4, 3, 2), (2, 5, 3), (3, 2, 3), (3, 6, 5)]
for aa, bb, cc in triples:
    G = build(aa, bb, cc)
    sub = {a: aa, b: bb, c: cc}
    if not (kf(G) == q(KfG.subs(sub)) and kf(minus(G, ('r', 0, 0))) == q(KfGv.subs(sub))):
        ok = False
        print('   mismatch', aa, bb, cc)
C.check(ok, '(C1) symbolic Kf(R_c) and Kf(R_c - v) = brute-force exact values for %d triples' % len(triples))

# reference from check_rook_hubs.py (c = 1)
nC = a * (b + 1)
MAr = sp.Matrix([[a + 1, -1], [-b, b]])
KfR = sp.cancel(nC * (sp.Rational(1) / (b + 1) + (a - 1) * sp.cancel(MAr.inv().trace()) + (b - 1) / (b + 1)
                       + (a - 1) * (b - 1) / (a + b + 1)))
C.check(sp.simplify(KfG.subs(c, 1) - KfR) == 0, '(C2) c = 1 reproduces Kf(R(a,b)) of check_rook_hubs.py')

num, den = sp.fraction(sp.together(Delta))
P = Poly(sp.expand(num), a, b, c)
print('   Delta denominator:', sp.factor(den))
hits = []
for bb in range(3, BMAX + 1):
    Pb = Poly(P.as_expr().subs(b, bb), a, c)
    for cc in range(1, bb + 1):
        pa = Poly(Pb.as_expr().subs(c, cc), a, domain=ZZ)
        for f, _ in sp.factor_list(pa.as_expr())[1]:
            fp = Poly(f, a)
            if fp.degree() == 1:
                c1, c0 = fp.all_coeffs()
                if c0 % c1 == 0 and -c0 // c1 >= 2:
                    hits.append((-c0 // c1, bb, cc))
C.check(not hits, '(C3) no integer a >= 2 with Delta = 0 for 3 <= b <= %d, 1 <= c <= b (hits: %s)' % (BMAX, hits))
for bb, cc in [(6, 3), (10, 3), (10, 4), (20, 9), (40, 19)]:
    rts = [float(r) for r in Poly(P.as_expr().subs({b: bb, c: cc}), a).real_roots() if float(r) > 1.5]
    print('   (C4) b=%d c=%d: real roots a* > 1.5: %s  (proportion b/(b+c) = %.4f)' % (bb, cc, ['%.3f' % r for r in rts], bb / (bb + cc)))
C.finish()
