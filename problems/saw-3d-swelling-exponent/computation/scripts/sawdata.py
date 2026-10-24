"""sawdata.py -- parse the output of src/saw_enum.cpp (plain or .gz) and compute exact derived quantities.

Python 3 standard library only (fractions, gzip).  Imported by the other scripts of this folder.

Data model (all integers exact):
  D.N                      largest n
  D.c[n], D.h[n], D.b[n]   number of n-step SAWs / strict half-space walks / bridges (Hammersley-Welsh convention)
  D.A[n][(a,b,c)] = (Nsum, orb)   a>=b>=c>=0: Nsum = number of n-step SAWs ending in the B_3-orbit of (a,b,c),
                                   orb = orbit size; c_n(x) = Nsum // orb for every x in the orbit
  D.H[n][(x1,b,c)], D.B[n][(x1,b,c)] = (Nsum, orb)   same for half-space walks / bridges, orbits of Stab(e1)
  D.M[(k,n)][(a,b,c)] = (Nsum, orb)   number of n-step SAWs whose k-th vertex lies in the orbit (1 <= k < n)
"""
import gzip
from fractions import Fraction


class SawData:
    def __init__(self, path):
        op = gzip.open if path.endswith(".gz") else open
        self.c, self.h, self.b = {}, {}, {}
        self.A, self.H, self.B, self.M = {}, {}, {}, {}
        self.header = []
        with op(path, "rt") as f:
            for line in f:
                if line.startswith("#"):
                    self.header.append(line.strip())
                    continue
                t = line.split()
                if not t:
                    continue
                tag = t[0]
                if tag == "T":
                    n = int(t[1])
                    self.c[n], self.h[n], self.b[n] = int(t[2]), int(t[3]), int(t[4])
                elif tag in "AHB":
                    n = int(t[1])
                    key = (int(t[2]), int(t[3]), int(t[4]))
                    getattr(self, tag).setdefault(n, {})[key] = (int(t[5]), int(t[6]))
                elif tag == "M":
                    k, n = int(t[1]), int(t[2])
                    key = (int(t[3]), int(t[4]), int(t[5]))
                    self.M.setdefault((k, n), {})[key] = (int(t[6]), int(t[7]))
        self.N = max(self.c)

    # ---- endpoint distribution -------------------------------------------------------------------------------
    def cnx(self, n, x):
        """c_n(x) for x a 3-tuple of integers (exact)."""
        a = sorted((abs(v) for v in x), reverse=True)
        e = self.A.get(n, {}).get(tuple(a))
        return 0 if e is None else e[0] // e[1]

    def S(self, n, power=1):
        """sum over n-step SAWs of |omega_n|^(2*power)."""
        return sum(Ns * (a * a + b * b + c * c) ** power for (a, b, c), (Ns, orb) in self.A[n].items())

    def mean_sq(self, n):
        return Fraction(self.S(n), self.c[n])

    def ending_on_plane(self, n):
        """number of n-step SAWs with x_1(omega_n) = 0.  For an orbit with z zero entries, z/3 of its points have
        first coordinate 0 (the orbit is invariant under coordinate permutations)."""
        tot = 0
        for (a, b, c), (Ns, orb) in self.A[n].items():
            z = (a == 0) + (b == 0) + (c == 0)
            assert (Ns * z) % 3 == 0
            tot += Ns * z // 3
        return tot

    def max_point(self, n):
        """(max_x c_n(x), list of orbit representatives attaining it)."""
        best, arg = -1, []
        for key, (Ns, orb) in self.A[n].items():
            v = Ns // orb
            if v > best:
                best, arg = v, [key]
            elif v == best:
                arg.append(key)
        return best, arg

    # ---- intermediate points ---------------------------------------------------------------------------------
    def S_mid(self, k, n):
        """sum over n-step SAWs of |omega_k|^2 (k <= n)."""
        if k == n:
            return self.S(n)
        if k == 0:
            return 0
        return sum(Ns * (a * a + b * b + c * c) for (a, b, c), (Ns, orb) in self.M[(k, n)].items())

    def max_point_mid(self, k, n):
        """max_x #{n-step SAWs with omega_k = x}."""
        tab = self.A[n] if k == n else self.M[(k, n)]
        return max(Ns // orb for (Ns, orb) in tab.values())

    # ---- bridges / half-space walks --------------------------------------------------------------------------
    def bridge_span_dist(self, n):
        """dict span -> number of n-step bridges with x_1(omega_n) = span."""
        d = {}
        for (x1, b, c), (Ns, orb) in self.B[n].items():
            d[x1] = d.get(x1, 0) + Ns
        return d

    def S_bridge(self, n):
        return sum(Ns * (x1 * x1 + b * b + c * c) for (x1, b, c), (Ns, orb) in self.B[n].items())

    def S_half(self, n):
        return sum(Ns * (x1 * x1 + b * b + c * c) for (x1, b, c), (Ns, orb) in self.H[n].items())


def orbit_size(a, b, c):
    perms = 1 if a == b == c else (3 if (a == b or b == c) else 6)
    return perms << ((a > 0) + (b > 0) + (c > 0))


def read_reference(path):
    """read a whitespace-separated reference table with '#' comments; returns list of rows (strings)."""
    rows = []
    for line in open(path):
        if line.startswith("#") or not line.strip():
            continue
        rows.append(line.split())
    return rows
