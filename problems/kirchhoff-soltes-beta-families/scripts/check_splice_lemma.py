"""check_splice_lemma.py -- exact checks of Lemma 2.1 (splice formula with an arbitrary rest), Corollary 2.2 (k copies)
and Example 2.5 (C5 copies tuned by a pendant path, the regime a < 0).

Claims checked (exact rational arithmetic, python-flint):
  (L1) Kf(B x W) = Kf(B) + Kf(W) + M Rt_B(x) + (|B|-1) T,   M = |W|-1, T = Rt_W(x);
  (L2) for v in B - x with B - v connected: Kf(G) - Kf(G - v) = delta_B(v) - M rho_B(v) + T;
  (C1) G = k copies of B at x plus W0:  Delta_G(v) = delta + (k-1) a - M0 rho + T0,  a = Rt_B(x) - (|B|-1) rho;
  (E1) C5 rooted at a vertex: the two neighbours of the root have rho = 2, delta = 0, a = -4; the other two
       have rho = 0, delta = 0;
  (E2) k copies of C5 at x plus a pendant path of length L with L(L-3) = 8(k-1): exactly 2k good vertices
       (the neighbours of x in the copies), for (L,k) = (3,1), (8,6), (11,12).
Random instances use a fixed seed (deterministic).  Runtime: about 1 min (Fedora 44 / WSL, python-flint 0.6).
"""
import random
import networkx as nx
from flint import fmpq
from kfexact import kf, rt, delta, minus, splice, good_vertices, Checker

C = Checker()
rng = random.Random(20261025)


def random_connected(nv, p):
    while True:
        G = nx.gnp_random_graph(nv, p, seed=rng.randrange(10**9))
        if nx.is_connected(G):
            return G


# (L1), (L2)
cnt = 0
ok1 = ok2 = True
for trial in range(60):
    B = random_connected(rng.randint(3, 7), 0.5)
    xB = rng.choice(list(B.nodes()))
    W = random_connected(rng.randint(1, 6), 0.5) if rng.random() < 0.9 else nx.path_graph(1)
    xW = rng.choice(list(W.nodes()))
    G = splice([(B, xB), (W, xW)])
    M = W.number_of_nodes() - 1
    T = rt(W, xW)
    KB, KW, KG = kf(B), kf(W), kf(G)
    RB = rt(B, xB)
    if KG != KB + KW + M * RB + (B.number_of_nodes() - 1) * T:
        ok1 = False
    for v in B.nodes():
        if v == xB:
            continue
        dB = delta(B, v, KB)
        if dB is None:
            continue
        rho = rt(minus(B, v), xB) - RB
        dG = delta(G, (0, v), KG)
        cnt += 1
        if dG != dB - M * rho + T:
            ok2 = False
C.check(ok1, '(L1) Kf of a splice, 60 random instances')
C.check(ok2, '(L2) Delta_G(v) = delta_B(v) - M rho_B(v) + T, %d random deletions' % cnt)

# (C1)
okc = True
cnt = 0
for trial in range(25):
    B = random_connected(rng.randint(3, 6), 0.55)
    xB = rng.choice(list(B.nodes()))
    W0 = random_connected(rng.randint(1, 5), 0.5)
    x0 = rng.choice(list(W0.nodes()))
    k = rng.randint(1, 3)
    G = splice([(B, xB)] * k + [(W0, x0)])
    KB, KG = kf(B), kf(G)
    RB = rt(B, xB)
    m = B.number_of_nodes() - 1
    M0 = W0.number_of_nodes() - 1
    T0 = rt(W0, x0)
    for v in B.nodes():
        if v == xB:
            continue
        d = delta(B, v, KB)
        if d is None:
            continue
        rho = rt(minus(B, v), xB) - RB
        a = RB - m * rho
        cnt += 1
        if delta(G, (0, v), KG) != d + (k - 1) * a - M0 * rho + T0:
            okc = False
C.check(okc, '(C1) k-copy formula Delta = delta + (k-1)a - M0 rho + T0, %d deletions' % cnt)

# (E1) C5 rooted at vertex 0
C5 = nx.cycle_graph(5)
K5 = kf(C5)
R5 = rt(C5, 0)
C.check(K5 == 10 and R5 == 4, '(E1) Kf(C5) = 10, Rt_C5(x) = 4')
vals = {}
for v in range(1, 5):
    rho = rt(minus(C5, v), 0) - R5
    vals[v] = (delta(C5, v, K5), rho, R5 - 4 * rho)
C.check(vals[1] == vals[4] == (0, 2, -4), '(E1) neighbours of the root: delta = 0, rho = 2, a = -4')
C.check(vals[2] == vals[3] == (0, 0, 4), '(E1) non-neighbours of the root: delta = 0, rho = 0, a = 4')

# (E2) C5 copies + pendant path
for (L, k) in [(3, 1), (8, 6), (11, 12)]:
    assert L * (L - 3) == 8 * (k - 1)
    P = nx.path_graph(L + 1)
    G = splice([(C5, 0)] * k + [(P, 0)])
    good = set(good_vertices(G))
    pred = {(i, 1) for i in range(k)} | {(i, 4) for i in range(k)}
    n = G.number_of_nodes()
    C.check(good == pred, '(E2) L=%d, k=%d: n=%d, good vertices = the 2k neighbours of x (proportion %s)'
            % (L, k, n, fmpq(2 * k, n)))
C.finish()
