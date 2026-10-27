"""Shared helpers of agent E2 (round 2 of the gac-duality work on result 17).

Environment: Python 3.12 (WSL Fedora 44), standard library only.  Exact rational arithmetic (fractions.Fraction);
the network classes and the first-round library are imported unchanged from ../../../computation/scripts/ecore.py
(which imports ../../../../scripts/crn.py).

Contents
  * exact linear feasibility with a witness point: Fourier-Motzkin elimination with back-substitution
    (fm_point, lp_point), used for convex-hull membership, rate-constant LPs and sign-vector realisability;
  * record loading (first-round pickles in ~/gacE-data/, round-2 pickles in ~/gacE2-data/), class labels;
  * small structural helpers (Petri connectivity, EC, WF-EC, describe).
"""
import os
import sys
import gzip
import pickle
from fractions import Fraction
from itertools import combinations

HERE = os.path.dirname(os.path.abspath(__file__))
ESCRIPTS = os.path.normpath(os.path.join(HERE, "..", "..", "..", "computation", "scripts"))
PSCRIPTS = os.path.normpath(os.path.join(HERE, "..", "..", "..", "petri", "scripts"))
sys.path.insert(0, ESCRIPTS)
import ecore  # noqa: E402
import endo  # noqa: E402
from ecore import CNet, crn  # noqa: E402,F401

DATA_E = os.path.expanduser("~/gacE-data")
DATA_E2 = os.path.expanduser("~/gacE2-data")

E_FILES = ["wr_2_3_6", "wr_3_2_6", "wr_3_3_5", "wr_4_2_6", "wr_5_2_5", "wr_5_2_6_homog", "wr_6_2_5_homog",
           "wr_2_2_5_zero", "wr_3_2_5_zero", "wr_4_3_4", "wr_5_3_4", "sample_5_2_4_8_s1", "sample_6_2_5_8_s2",
           "all_2_2_5", "all_2_2_6", "all_2_3_4", "all_3_2_4", "all_3_3_3", "all_4_1_12", "all_4_2_3",
           "fc_3_6", "fc_4_5"]


def rel(path):
    """Path printed in logs: relative to the data directory or the repository, never with a user name."""
    home = os.path.expanduser("~")
    if path.startswith(home):
        return "~" + path[len(home):]
    i = path.find("problems" + os.sep)
    return path[i:] if i >= 0 else os.path.basename(path)


def load(names, base=DATA_E):
    seen = {}
    for nm in names:
        fn = nm if os.path.isabs(nm) else os.path.join(base, nm + (".pkl.gz" if not nm.endswith(".gz") else ""))
        with gzip.open(fn, "rb") as fh:
            for r in pickle.load(fh):
                key = (r["n"], r["rx"])
                if key not in seen:
                    seen[key] = r
    return list(seen.values())


def H(rec):
    N = rec["N"]
    return N["cons"] and N["consist"] and N["rows_distinct"]


def Hd(rec):
    D = rec["D"]
    return D["cons"] and D["consist"] and D["rows_distinct"]


def classes(rec):
    out = ["all"]
    if H(rec):
        out.append("H")
        if rec["N"]["WR"]:
            out.append("WR")
            if rec["N"]["deficiency"] == 0:
                out.append("WR0")
    return out


# ----------------------------------------------------------------------------------------------
# exact linear feasibility with witness
# ----------------------------------------------------------------------------------------------
def _norm(a, b):
    """Scale a constraint a.x <= b so that its first nonzero coefficient has absolute value 1 (for dedup)."""
    for x in a:
        if x != 0:
            s = abs(x)
            return tuple(y / s for y in a), b / s
    return tuple(a), b


def fm_point(ineqs, d, cap=20000):
    """ineqs: list of (a, b) meaning a.x <= b, a of length d.  Returns an exact rational point or None."""
    cur = []
    seen = set()
    for a, b in ineqs:
        a, b = _norm([Fraction(x) for x in a], Fraction(b))
        if (a, b) not in seen:
            seen.add((a, b))
            cur.append((a, b))
    stages = [cur]
    for v in reversed(range(d)):
        pos = [c for c in cur if c[0][v] > 0]
        neg = [c for c in cur if c[0][v] < 0]
        new = [c for c in cur if c[0][v] == 0]
        seen = set(new)
        for ap, bp in pos:
            for an, bn in neg:
                fp, fn = ap[v], -an[v]
                a = [x / fp + y / fn for x, y in zip(ap, an)]
                a[v] = Fraction(0)
                c = _norm(a, bp / fp + bn / fn)
                if c not in seen:
                    seen.add(c)
                    new.append(c)
        if len(new) > cap:
            raise RuntimeError("Fourier-Motzkin blow-up")
        cur = new
        stages.append(cur)
    if any(b < 0 for a, b in cur):
        return None
    x = [Fraction(0)] * d
    for v in range(d):
        sysv = stages[d - 1 - v]
        lo, hi = None, None
        for a, b in sysv:
            if a[v] == 0:
                continue
            rhs = (b - sum(a[u] * x[u] for u in range(v))) / a[v]
            if a[v] > 0:
                hi = rhs if hi is None else min(hi, rhs)
            else:
                lo = rhs if lo is None else max(lo, rhs)
        if lo is not None and hi is not None:
            assert lo <= hi
            x[v] = lo if lo == hi else (lo + hi) / 2
        elif lo is not None:
            x[v] = lo
        elif hi is not None:
            x[v] = hi
    for a, b in ineqs:
        assert sum(Fraction(p) * q for p, q in zip(a, x)) <= b
    return x


def lp_point(eqs, ineqs, d):
    """Exists x in Q^d with a.x = b ((a,b) in eqs) and a.x <= b ((a,b) in ineqs)?  Returns a point or None."""
    if not eqs:
        return fm_point(ineqs, d)
    M = [[Fraction(x) for x in a] + [Fraction(b)] for a, b in eqs]
    piv = []
    row = 0
    for col in range(d):
        p = next((i for i in range(row, len(M)) if M[i][col] != 0), None)
        if p is None:
            continue
        M[row], M[p] = M[p], M[row]
        pv = M[row][col]
        M[row] = [x / pv for x in M[row]]
        for i in range(len(M)):
            if i != row and M[i][col] != 0:
                f = M[i][col]
                M[i] = [x - f * y for x, y in zip(M[i], M[row])]
        piv.append(col)
        row += 1
    if any(all(x == 0 for x in M[i][:d]) and M[i][d] != 0 for i in range(len(M))):
        return None
    free = [c for c in range(d) if c not in piv]
    # x_piv = rhs - sum_free M * x_free
    def expand(t):
        x = [Fraction(0)] * d
        for k, c in enumerate(free):
            x[c] = t[k]
        for i, c in enumerate(piv):
            x[c] = M[i][d] - sum(M[i][f] * x[f] for f in free)
        return x
    # substitute into inequalities: a.x = a0 + sum_k (coef_k) t_k
    base = expand([Fraction(0)] * len(free))
    units = []
    for k in range(len(free)):
        t = [Fraction(0)] * len(free)
        t[k] = Fraction(1)
        e = expand(t)
        units.append([e[i] - base[i] for i in range(d)])
    red = []
    for a, b in ineqs:
        a = [Fraction(x) for x in a]
        a0 = sum(p * q for p, q in zip(a, base))
        red.append(([sum(p * q for p, q in zip(a, u)) for u in units], Fraction(b) - a0))
    t = fm_point(red, len(free)) if free else ([] if all(b >= 0 for a, b in red) else None)
    if t is None:
        return None
    x = expand(t)
    for a, b in eqs:
        assert sum(Fraction(p) * q for p, q in zip(a, x)) == b
    for a, b in ineqs:
        assert sum(Fraction(p) * q for p, q in zip(a, x)) <= b
    return x


def in_convex_hull(e, pts):
    """Exact: is e a convex combination of pts?  Returns the weights or None."""
    k = len(pts)
    d = len(e)
    eqs = [([pts[j][i] for j in range(k)], e[i]) for i in range(d)] + [([1] * k, 1)]
    ineqs = [([-int(i == j) for i in range(k)], 0) for j in range(k)]
    return lp_point(eqs, ineqs, k)


# ----------------------------------------------------------------------------------------------
# structure
# ----------------------------------------------------------------------------------------------
def petri_connected(n, rx):
    parent = list(range(n + len(rx)))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    for j, (a, b) in enumerate(rx):
        for i in range(n):
            if a[i] or b[i]:
                parent[find(i)] = find(n + j)
    return len({find(x) for x in range(n + len(rx))}) == 1


def is_EC(rx):
    for (a1, b1), (a2, b2) in combinations(rx, 2):
        if a1 != a2 and any(x and y for x, y in zip(a1, a2)):
            return False
    return True


def n_sources(rx):
    return len({a for a, b in rx})


def wfec(rec):
    """WF-EC (agent P): EC, connected, (H1), (H2), s = m - 1 with m the number of source complexes."""
    rx, n = rec["rx"], rec["n"]
    N = rec["N"]
    return (N["cons"] and N["consist"] and is_EC(rx) and petri_connected(n, rx)
            and rec["s"] == n_sources(rx) - 1)


def dual_rx(n, rx):
    Ye = [tuple(a[i] for a, b in rx) for i in range(n)]
    Yp = [tuple(b[i] for a, b in rx) for i in range(n)]
    return tuple(zip(Ye, Yp))


def describe(n, rx, names=None):
    return CNet(n, rx, names).describe()
