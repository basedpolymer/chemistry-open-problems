"""check_exact_proportions.py -- exact checks for Theorem 3.1 (every rational beta in (0, 2/3) is the exact proportion
of good vertices of infinitely many graphs) and Corollary 3.3.

Building blocks (rooted at x):
  D  = K4 - e rooted at a vertex of degree 2 (the gadget G_1 of Gottwald et al.);  M = 3, T = 9/4;
  F5 = (K4 - e rooted at a vertex of degree 3) + 2 pendant edges at x;              M = 5, T = 15/4;
  F6 = K4 rooted at a vertex + 3 pendant edges at x;                                M = 6, T = 9/2.
  G(k, al, ga) = splice at x of k copies of D, al copies of F5, ga copies of F6;  n = 1 + 3k + 5 al + 6 ga.
Claims:
  (V) the rooted constants (Rt, delta, rho) of D, D', K4, K2 used in the proof;
  (N) the three branch types are neutral: T = (3/4) M;
  (F) for every (k, al, ga) with n <= 46 (and k >= 1): the good vertices of G(k, al, ga) are exactly the 2k
      vertices of degree 3 of the copies of D (brute force over all vertices, exact rationals);
  (R) every integer N >= 20 is 5 al + 6 ga with al, ga >= 0, and 19 is not;
  (B) for a list of rationals beta = p/q < 2/3: the parameters of the proof (k = p j, 5 al + 6 ga = (2q-3p) j - 1)
      give n = 2qj and 2k/n = beta; for the smallest admissible j of several beta the graph is built and its good
      set computed by brute force (n <= 80).
Runtime: about 2 min (Fedora 44 / WSL, python-flint 0.6).
"""
import itertools
import networkx as nx
from fractions import Fraction
from flint import fmpq
from kfexact import kf, rt, delta, minus, splice, good_vertices, resistance, Checker

C = Checker()

D = nx.Graph([('x', 'a'), ('x', 'b'), ('a', 'b'), ('a', 'y'), ('b', 'y')])        # K4 - e, root x of degree 2
Dp = nx.Graph([('r', 'b'), ('r', 'u'), ('r', 'w'), ('b', 'u'), ('b', 'w')])       # K4 - e, root r of degree 3
K4 = nx.complete_graph(4)
E = nx.path_graph(2)


def consts(B, x):
    K, R = kf(B), rt(B, x)
    out = {}
    for v in B.nodes():
        if v != x:
            d = delta(B, v, K)
            out[v] = (d, rt(minus(B, v), x) - R)
    return R, out


R, c = consts(D, 'x')
C.check(R == fmpq(9, 4) and c['a'] == c['b'] == (0, fmpq(3, 4)) and c['y'] == (2, fmpq(-11, 12)),
        '(V) D rooted at degree-2 vertex: Rt = 9/4; a,b: delta 0, rho 3/4; y: delta 2, rho -11/12')
R, c = consts(Dp, 'r')
C.check(R == fmpq(7, 4) and c['b'] == (0, fmpq(1, 4)) and c['u'] == c['w'] == (2, fmpq(-5, 12)),
        "(V) D' rooted at degree-3 vertex: Rt = 7/4; other degree-3: delta 0, rho 1/4; degree-2: delta 2, rho -5/12")
R, c = consts(K4, 0)
C.check(R == fmpq(3, 2) and all(c[v] == (1, fmpq(-1, 6)) for v in (1, 2, 3)),
        '(V) K4 rooted: Rt = 3/2; others: delta 1, rho -1/6')
R, c = consts(E, 0)
C.check(R == 1 and c[1] == (1, -1), '(V) pendant edge: Rt = 1; leaf: delta 1, rho -1')
# minimum resistance from the root inside each branch type (used for the vertex b' of D')
mins = [min(resistance(B, u, x) for u in B.nodes() if u != x) for B, x in [(D, 'x'), (Dp, 'r'), (K4, 0), (E, 0)]]
C.check(min(mins) == fmpq(1, 2), '(V) every resistance from x inside a branch is >= 1/2 (minimum %s)' % min(mins))

F5 = splice([(Dp, 'r'), (E, 0), (E, 0)])
F6 = splice([(K4, 0), (E, 0), (E, 0), (E, 0)])
for name, (B, x) in [('D', (D, 'x')), ('F5', (F5, 'x')), ('F6', (F6, 'x'))]:
    M = B.number_of_nodes() - 1
    C.check(rt(B, x) == fmpq(3, 4) * M, '(N) %s is neutral: Rt = (3/4) M with M = %d' % (name, M))


def family(k, al, ga):
    return splice([(D, 'x')] * k + [(F5, 'x')] * al + [(F6, 'x')] * ga)


def predicted(k):
    return {(i, 'a') for i in range(k)} | {(i, 'b') for i in range(k)}


# (F) brute force on all small members
cnt = 0
ok = True
for k in range(1, 15):
    for al in range(0, 9):
        for ga in range(0, 8):
            n = 1 + 3 * k + 5 * al + 6 * ga
            if n > 46:
                continue
            G = family(k, al, ga)
            assert G.number_of_nodes() == n
            g = set(good_vertices(G))
            cnt += 1
            if g != predicted(k):
                ok = False
                print('   mismatch at', (k, al, ga), sorted(map(str, g)))
C.check(ok, '(F) good set of G(k, al, ga) = the 2k degree-3 vertices of the copies of D, %d graphs with n <= 46' % cnt)

# (R) representability
rep = {5 * s + 6 * t for s in range(81) for t in range(68) if 5 * s + 6 * t < 400}
C.check(all(N in rep for N in range(20, 400)) and 19 not in rep, '(R) 5a + 6b represents every N >= 20, not 19')


def params(p, q, j):
    d = 2 * q - 3 * p
    Nf = d * j - 1
    for ga in range(Nf // 6 + 1):
        if (Nf - 6 * ga) % 5 == 0:
            return p * j, (Nf - 6 * ga) // 5, ga
    return None


betas = [Fraction(1, 2), Fraction(1, 3), Fraction(2, 5), Fraction(3, 5), Fraction(5, 8), Fraction(1, 7),
         Fraction(7, 11), Fraction(13, 20), Fraction(65, 98), Fraction(1, 1000), Fraction(333, 500)]
ok = True
for be in betas:
    p, q = be.numerator, be.denominator
    assert 3 * p < 2 * q
    for j in range(21, 41):
        k, al, ga = params(p, q, j)
        n = 1 + 3 * k + 5 * al + 6 * ga
        if not (n == 2 * q * j and Fraction(2 * k, n) == be and al >= 0 and ga >= 0):
            ok = False
C.check(ok, '(B) parameters of the proof give proportion exactly beta for j = 21..40, %d values of beta' % len(betas))

# smallest admissible j and brute-force check of the graph when n <= 80
for be in [Fraction(1, 2), Fraction(1, 3), Fraction(2, 5), Fraction(3, 5), Fraction(5, 8), Fraction(4, 7)]:
    p, q = be.numerator, be.denominator
    for j in range(1, 60):
        r = params(p, q, j)
        if r is None:
            continue
        k, al, ga = r
        if k >= 1 and (k + al + ga >= 2 or (k == 1 and al == ga == 0)):
            break
    n = 1 + 3 * k + 5 * al + 6 * ga
    if n <= 80:
        G = family(k, al, ga)
        g = set(good_vertices(G))
        C.check(g == predicted(k) and Fraction(len(g), n) == be,
                '(B) beta = %s: j = %d, (k, al, ga) = (%d, %d, %d), n = %d, exactly %d good vertices' % (be, j, k, al, ga, n, len(g)))
C.finish()
