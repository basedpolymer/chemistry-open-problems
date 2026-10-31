"""hubgraph.py -- shared library for the all-beta certificates (exact rational arithmetic only).

The graph H(N; b_1, ..., b_a): a clique K_N (the "base"); the base is split into blocks B_1, ..., B_a of sizes b_i >= 1
and a set C of "hubless" vertices; each block has a hub h_i adjacent to the b_i vertices of B_i and to nothing else.
Write a = number of hubs, a' = a + |C| (blocks plus hubless vertices), Y = sum 1/b_i.

Closed forms proved in paper.tex, Lemma 3:
    Kf(H) = Y (N+1)(N+a-1)/N + (N+a-2) + (1+a')(1-a)/(N+1) + a(a'-2)/N
    Rt_H(c) = 2 + Y + (Y+a+a'-2)/N - (1+a')/(N+1)     (c hubless)

Environment: Python 3.12+, sympy, python-flint, networkx (Fedora 44 under WSL).  Fractions only; no floating point in a
claim.  Functions marked "exact" return fractions.Fraction.
"""
from fractions import Fraction as Fr
from math import gcd, isqrt
import networkx as nx


class Checker:
    """PASS/FAIL bookkeeping; exit status 1 if any claim fails."""
    def __init__(self):
        self.n_pass = 0
        self.n_fail = 0

    def check(self, ok, msg):
        if ok:
            self.n_pass += 1
            print('PASS', msg, flush=True)
        else:
            self.n_fail += 1
            print('FAIL', msg, flush=True)

    def finish(self):
        print('summary: %d PASS, %d FAIL' % (self.n_pass, self.n_fail), flush=True)
        raise SystemExit(1 if self.n_fail else 0)


# ---- closed forms (exact) ----------------------------------------------------------------------------------------

def kf_closed(N, a, ap, Y):
    """Kf of H with N base vertices, a hubs, ap = a + #hubless, Y = sum 1/b_i."""
    return Fr(Y)*Fr((N+1)*(N+a-1), N) + (N+a-2) + Fr((1+ap)*(1-a), N+1) + Fr(a*(ap-2), N)


def rt_hubless_closed(N, a, ap, Y):
    """Resistance transmission of a hubless base vertex."""
    return 2 + Fr(Y) + (Fr(Y)+a+ap-2)/N - Fr(1+ap, N+1)


def delta_closed(a, b):
    """delta = Kf(H) - Kf(H-v), H = H_{a,b} (N = ab+1, one hubless vertex w), v a block vertex; partial fractions."""
    B = b*(b-1)
    return (Fr(-2*a, B) + Fr(b**4-b**3-2*b+1, b**3*(b-1)) + Fr(2*(b-1)*(b+2), b*b*(a*b+2))
            - Fr((b+1)*(2*b*b-2*b-1), b**3*(a*b+1)) + Fr(1, a*b*b*(b-1)))


def rho_closed(a, b):
    """rho = Rt_{H-v}(w) - Rt_H(w); partial fractions."""
    B = b*(b-1)
    return (Fr(1, B) + Fr(2*(b-1), b*(a*b+2)) - Fr(b*b-3*b-1, b*b*(a*b+1)) - Fr(b*b-b-1, a*b*b*(b-1)))


def delta_hub_closed(a, b):
    """delta_h = Kf(H) - Kf(H-h), h a hub of H_{a,b}."""
    N = a*b+1
    return Fr((N+1)*(N+2*a-2), b*N) + 1 + Fr(a*b-2*a-2*b, N+1) + Fr((a-1)*(2-b), N)


def rho_hub_closed(a, b):
    """rho_h = Rt_{H-h}(w) - Rt_H(w)."""
    N = a*b+1
    return Fr(-1, b) + (b-2-Fr(1, b))/N - Fr(b-1, N+1)


# ---- explicit graphs ------------------------------------------------------------------------------------------------

def hub_graph(blocks, hubless=0):
    """H with block sizes `blocks` and `hubless` hubless base vertices; N = sum(blocks)+hubless."""
    G = nx.Graph()
    base = []
    for i, b in enumerate(blocks):
        for j in range(b):
            base.append(('u', i, j))
    for j in range(hubless):
        base.append(('w', j))
    G.add_nodes_from(base)
    for p in range(len(base)):
        for q in range(p+1, len(base)):
            G.add_edge(base[p], base[q])
    for i, b in enumerate(blocks):
        for j in range(b):
            G.add_edge(('h', i), ('u', i, j))
    return G


def tree_depth_data(l, t):
    """Lemma 5: a rooted tree with l non-root vertices and sum of depths t, l <= t <= l(l+1)/2.
    Returns (k, e): a path of k vertices at depths 1..k, the other l-k at depth 1, except one at depth 1+e (0 <= e < k)."""
    if l == 0:
        assert t == 0
        return (0, 0)
    assert l <= t <= l*(l+1)//2
    k = 1
    while k+1 <= l and (k+1)*(k+2)//2 + l-k-1 <= t:
        k += 1
    e = t - (k*(k+1)//2 + l-k)
    assert 0 <= e and (e == 0 or (e < k and k < l))
    return (k, e)


def add_tree(G, root, l, t, tag='t'):
    """Attach the tree of Lemma 5 at root; returns the list of its non-root vertices."""
    if l == 0:
        return []
    k, e = tree_depth_data(l, t)
    path = [root] + [(tag, 'p', i) for i in range(k)]
    for i in range(k):
        G.add_edge(path[i], path[i+1])
    rest = l-k
    new = path[1:]
    for i in range(rest):
        z = (tag, 'q', i)
        G.add_edge(path[e] if (i == 0 and e > 0) else root, z)
        new.append(z)
    return new


def add_clique(G, root, c, tag):
    """Glue K_c at root; returns the c-1 new vertices."""
    vs = [root] + [(tag, i) for i in range(c-1)]
    for p in range(len(vs)):
        for q in range(p+1, len(vs)):
            G.add_edge(vs[p], vs[q])
    return vs[1:]


def add_gadget(G, root, c, j, tag):
    """Gadget Z: K_c at root with j pendant leaves on one fixed non-root vertex y."""
    vs = add_clique(G, root, c, tag)
    y = vs[0]
    for i in range(j):
        G.add_edge(y, (tag, 'leaf', i))
    return vs


# ---- the construction of Theorem 1 ----------------------------------------------------------------------------------

def q_list(b, s):
    """x = b(b-1)s and the pairwise coprime moduli Q_0, ..., Q_4 of Lemma 6."""
    x = b*(b-1)*s
    Q = [b*b*(b-1)*x*(x-2), x*x-3*x+1, x*x-2*x-1, x*x-x-1, x-1]
    return x, Q


def construct(p, r, b, s):
    """Run the proof of Theorem 1 for beta = p/r, block size b, parameter s (multiple of p), exact arithmetic.
    Returns a dict with every number of the proof, or None if s is too small for the size conditions."""
    assert gcd(p, r) == 1 and 0 < p < r and s % p == 0 and b*(r-p) > p
    B = b*(b-1)
    x, Q = q_list(b, s)
    if x < 4:
        return None
    a = (b-1)*s*Q[1]
    assert a*b == x*Q[1]
    N = a*b+1
    n = a*b*r//p
    assert n*p == a*b*r
    m = n-N-a
    delta = delta_closed(a, b)
    rho = rho_closed(a, b)
    E = -delta
    Sig = sum(Q)
    # L = ceiling of 2 sqrt(E+1)
    Lval = E+1
    L = isqrt(4*(Lval.numerator // Lval.denominator + 1)) + 1
    k = (m-3*Sig-L) // (B*(2*B-1))
    if k < 0:
        return None
    m_r = m-B*k*(2*B-1)
    tau = m*rho-delta
    taup = tau-k*(2*B-1)
    D = 1
    for q in Q:
        D *= q
    assert (taup*D).denominator == 1
    A = int(taup*D)
    j = []
    for q in Q:
        c = A*pow(D//q, -1, q) % q
        j.append((c+1) % q)
    T_g = sum((2-Fr(1, q)+ji*(1+Fr(1, q)) for q, ji in zip(Q, j)), Fr(0))
    m_g = sum(2*q-1+ji for q, ji in zip(Q, j))
    l = m_r-m_g
    t = taup-T_g
    if t.denominator != 1:
        raise AssertionError('t not integral')
    t = int(t)
    if not (l <= t <= l*(l+1)//2):
        return None
    return dict(p=p, r=r, b=b, s=s, x=x, Q=Q, a=a, N=N, n=n, m=m, delta=delta, rho=rho, E=E, Sig=Sig, L=L, k=k,
                m_r=m_r, tau=tau, j=j, T_g=T_g, m_g=m_g, l=l, t=t, B=B,
                T=k*(2*B-1)+T_g+t)
