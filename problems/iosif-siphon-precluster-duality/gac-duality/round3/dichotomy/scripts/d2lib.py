"""Agent D2 (round 3): exact boundary analysis of the rank-2 dual N* of a network N (s = 2).

Environment: Python 3.12, standard library only (fractions.Fraction).  Imports nothing from other agents except
through the callers.

Conventions.  N has species X_1..X_n (index i) and reactions R_j: y_j -> y'_j (index j), Gamma = Yp - Ye (n x r), rank 2.
The dual N* has species k_1..k_r and reactions rho_i: sum_j (y_j)_i k_j -> sum_j (y'_j)_i k_j, rate kappa_i k^{s_i},
s_i = row i of Ye.  Planar coordinates: Gamma = Pi * Ghat with Pi (n x 2, rows pi_i) and Ghat (2 x r, columns g_j), so
Gamma_ij = <pi_i, g_j>.  A class of N* is the polygon P(k0) = {X in R^2 : k_j(X) = k0_j + <g_j, X> >= 0}, and the
dual dynamics reads X' = sum_i kappa_i k(X)^{s_i} pi_i.

Generic boundary pieces (every positive class whose polygon has no tie):
  edge j    : zero set {j}; every point p >= 0 with zero set exactly {j} lies in the relative interior of the j-edge
              of some class (proved in notes.md, Lemma 2.1);
  vertex jj': zero set {j, j'}, g_j and g_j' linearly independent; the other coordinates of the vertex can be any
              positive vector (Lemma 2.1).
"""
from fractions import Fraction
from itertools import combinations


def frac(x):
    return Fraction(x)


class Net:
    def __init__(self, n, rx):
        self.n = n
        self.rx = tuple(rx)
        self.r = r = len(rx)
        self.Y = [[rx[j][0][i] for j in range(r)] for i in range(n)]          # Ye, n x r
        self.Yp = [[rx[j][1][i] for j in range(r)] for i in range(n)]
        self.G = [[self.Yp[i][j] - self.Y[i][j] for j in range(r)] for i in range(n)]
        self.changed = [i for i in range(n) if any(self.G[i])]
        # planar factorization
        cols = [[Fraction(self.G[i][j]) for i in range(n)] for j in range(r)]
        j1 = next(j for j in range(r) if any(cols[j]))
        j2 = None
        for j in range(r):
            if any(cols[j1][a] * cols[j][b] - cols[j1][b] * cols[j][a] for a in range(n) for b in range(a + 1, n)):
                j2 = j
                break
        assert j2 is not None, "rank < 2"
        u, w = cols[j1], cols[j2]
        guu = sum(x * x for x in u)
        guw = sum(x * y for x, y in zip(u, w))
        gww = sum(x * x for x in w)
        det = guu * gww - guw * guw
        self.gh = []
        for j in range(r):
            ru = sum(x * y for x, y in zip(u, cols[j]))
            rw = sum(x * y for x, y in zip(w, cols[j]))
            al = (gww * ru - guw * rw) / det
            be = (-guw * ru + guu * rw) / det
            assert all(al * u[i] + be * w[i] == cols[j][i] for i in range(n)), "rank > 2"
            self.gh.append((al, be))
        self.pi = [(u[i], w[i]) for i in range(n)]

    # ------------------------------------------------------------------ dual
    def dual_rx(self):
        return tuple((tuple(self.Y[i]), tuple(self.Yp[i])) for i in range(self.n))

    def parallel(self, j, l):
        a, b = self.gh[j], self.gh[l]
        return a[0] * b[1] - a[1] * b[0] == 0 and a[0] * b[0] + a[1] * b[1] > 0

    def independent(self, j, l):
        a, b = self.gh[j], self.gh[l]
        return a[0] * b[1] - a[1] * b[0] != 0

    # ------------------------------------------------------------------ generic edge j
    def edge(self, j):
        """Data of the generic edge {j}: siphon flag, transverse type, minimal order mu, the minimal-order terms
        (i, sign, beta_i) and the edge-flow terms (i, c_i) with pi_i = c_i * d, d = (-g_j[1], g_j[0])."""
        n, r = self.n, self.r
        y = [self.Y[i][j] for i in range(n)]
        ch = [i for i in range(n) if self.G[i][j] != 0]
        siphon = all(y[i] > 0 for i in ch)
        mu = min(y[i] for i in ch)
        J = [i for i in ch if y[i] == mu]
        terms = [(i, 1 if self.G[i][j] > 0 else -1, tuple(self.Y[i][l] for l in range(r) if l != j)) for i in J]
        good = [t for t in terms if t[1] > 0]
        bad = [t for t in terms if t[1] < 0]
        if not bad:
            typ = "IN"
        elif not good:
            typ = "OUT"
        else:
            typ = "MIX"  # ESC or DOM, decided by the caller with an exact convex-hull test
        d = (-self.gh[j][1], self.gh[j][0])
        flow = []
        for i in range(n):
            if self.Y[i][j] == 0 and any(self.pi[i]):
                p = self.pi[i]
                # pi_i parallel to d on a siphon edge (Gamma_ij = 0 for order-0 species)
                c = p[0] / d[0] if d[0] != 0 else p[1] / d[1]
                assert (c * d[0], c * d[1]) == p or not siphon
                flow.append((i, c))
        return {"j": j, "siphon": siphon, "type": typ, "mu": mu, "terms": terms, "flow": flow, "d": d}

    # ------------------------------------------------------------------ generic vertex {j, j'}
    def vertex(self, j, jj):
        n, r = self.n, self.r
        eq = all(self.Y[i][j] + self.Y[i][jj] >= 1 for i in self.changed)
        data = []
        for i in self.changed:
            p, q = self.Y[i][j], self.Y[i][jj]
            w = tuple(self.Y[i][l] for l in range(r) if l not in (j, jj))
            data.append((i, p, q, w, self.G[i][j], self.G[i][jj]))
        return {"j": j, "jj": jj, "eq": eq, "data": data}


def describe(n, rx, names=None):
    names = names or [f"X{i+1}" for i in range(n)]

    def cx(v):
        t = []
        for i, c in enumerate(v):
            if c:
                t.append((f"{c}" if c > 1 else "") + names[i])
        return " + ".join(t) if t else "0"
    return ", ".join(f"{cx(a)} -> {cx(b)}" for a, b in rx)
