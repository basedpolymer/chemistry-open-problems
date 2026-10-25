"""kfexact.py -- exact Kirchhoff-index utilities shared by the certificate scripts of this folder.

Graphs are networkx graphs.  For a connected graph H on N vertices and a vertex g, let A be the Laplacian of H with
the row and column of g deleted and M = A^{-1} (exact rationals, python-flint fmpq_mat).  Then
    r(u, g) = M_uu,  r(u, w) = M_uu + M_ww - 2 M_uw,  Kf(H) = N tr(M) - 1^T M 1,  Rt_H(g) = tr(M).
(Grounded-Laplacian formula, see e.g. Lemma 1 of ../kirchhoff-soltes-problem/result.md.)

Conventions of the manuscript: a vertex v is good if H - v is connected and Kf(H - v) = Kf(H);
delta_H(v) = Kf(H) - Kf(H - v);  for a root x != v, rho_x(v) = Rt_{H-v}(x) - Rt_H(x).
Environment: Python 3.12+, python-flint 0.6+, networkx 3.x (Fedora 44 under WSL).
"""
from flint import fmpq_mat, fmpq
import networkx as nx


def _grounded_inverse(H, g):
    nodes = [u for u in H.nodes() if u != g]
    idx = {u: i for i, u in enumerate(nodes)}
    N = len(nodes)
    A = fmpq_mat(N, N)
    for u in nodes:
        i = idx[u]
        A[i, i] = H.degree(u)
        for w in H[u]:
            if w != g:
                A[i, idx[w]] = -1
    return nodes, idx, A.inv()


def kf(H):
    """Exact Kirchhoff index of a connected graph (fmpq)."""
    if H.number_of_nodes() <= 1:
        return fmpq(0)
    if not nx.is_connected(H):
        raise ValueError('disconnected graph')
    g = next(iter(H.nodes()))
    nodes, idx, M = _grounded_inverse(H, g)
    N = len(nodes)
    tr = sum((M[i, i] for i in range(N)), fmpq(0))
    s = sum((M[i, j] for i in range(N) for j in range(N)), fmpq(0))
    return H.number_of_nodes() * tr - s


def rt(H, x):
    """Exact resistance transmission Rt_H(x) = sum_u r(x, u)."""
    if H.number_of_nodes() <= 1:
        return fmpq(0)
    nodes, idx, M = _grounded_inverse(H, x)
    return sum((M[i, i] for i in range(len(nodes))), fmpq(0))


def resistance(H, u, w):
    nodes, idx, M = _grounded_inverse(H, w)
    return M[idx[u], idx[u]]


def minus(H, v):
    K = H.copy()
    K.remove_node(v)
    return K


def delta(H, v, K_H=None):
    """Kf(H) - Kf(H - v), or None if H - v is disconnected (v is then never good)."""
    Hv = minus(H, v)
    if Hv.number_of_nodes() > 0 and not nx.is_connected(Hv):
        return None
    if K_H is None:
        K_H = kf(H)
    return K_H - kf(Hv)


def good_vertices(H):
    K = kf(H)
    out = []
    for v in H.nodes():
        d = delta(H, v, K)
        if d is not None and d == 0:
            out.append(v)
    return out


def splice(branches):
    """Splice of rooted graphs at a common root.  branches: list of (graph, root).  The root of every branch is
    renamed 'x'; the other vertices become (branch index, original name)."""
    G = nx.Graph()
    G.add_node('x')
    for i, (B, r) in enumerate(branches):
        ren = {u: ('x' if u == r else (i, u)) for u in B.nodes()}
        for u in B.nodes():
            G.add_node(ren[u])
        for u, w in B.edges():
            G.add_edge(ren[u], ren[w])
    return G


class Checker:
    """PASS/FAIL bookkeeping; exit status 1 if any claim fails."""
    def __init__(self):
        self.n_pass = 0
        self.n_fail = 0

    def check(self, ok, msg):
        if ok:
            self.n_pass += 1
            print('PASS', msg)
        else:
            self.n_fail += 1
            print('FAIL', msg)

    def finish(self):
        print('summary: %d PASS, %d FAIL' % (self.n_pass, self.n_fail))
        raise SystemExit(1 if self.n_fail else 0)
