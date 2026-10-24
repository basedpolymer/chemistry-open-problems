"""Exact (rational) representation of minimum-number-of-matches instances.

Independent re-implementation (does not import or copy the Letsios-Kouyialis-Misener Pyomo code).

Instance format (data/mip_instances/*/*.dat of github.com/dimletsios/min_matches_heuristics):
    n=.., m=.., k=..                      hot streams (incl. hot utilities), cold streams (incl. cold utilities),
                                         temperature intervals T0 (hottest) .. T{k-1} (coldest)
    QH[i]: T<t> <heat> ...                heat supplied by hot stream i in interval t
    QC[j]: T<t> <heat> ...                heat demanded by cold stream j in interval t
    R[t]= ..                              residual heat entering interval t (t = 0..k)
All numbers are read as exact decimal rationals (fractions.Fraction of the printed decimal string).

Heat of hot stream i in interval s may be sent to cold stream j in interval t only if s <= t
(heat flows from hotter to colder intervals; the DTmin shift is already built into the intervals).

Environment: Python >= 3.10, standard library only.
"""
from fractions import Fraction
import re
import math


class Instance:
    def __init__(self, name, n, m, k, QH, QC, R, cost=None):
        self.name = name
        self.n, self.m, self.k = n, m, k
        self.QH = QH  # list (len n) of lists (len k) of Fraction
        self.QC = QC  # list (len m) of lists (len k) of Fraction
        self.R = R    # list (len k+1) of Fraction, as printed in the file
        self.cost = cost

    # ---- aggregated quantities -------------------------------------------------
    def hot_total(self, i):
        return sum(self.QH[i])

    def cold_total(self, j):
        return sum(self.QC[j])

    def hot_cum(self, i, t):
        """heat of hot stream i in intervals 0..t (inclusive)"""
        return sum(self.QH[i][: t + 1])

    def cold_cum(self, j, t):
        return sum(self.QC[j][: t + 1])

    def residuals(self):
        """exact cascade residual R'[t] = sum_{s<t} (sum_i QH[i][s] - sum_j QC[j][s]), t=0..k"""
        res = [Fraction(0)]
        for t in range(self.k):
            res.append(res[-1] + sum(self.QH[i][t] for i in range(self.n)) - sum(self.QC[j][t] for j in range(self.m)))
        return res

    def denominators_lcm(self):
        L = 1
        for row in self.QH + self.QC:
            for v in row:
                L = L * v.denominator // math.gcd(L, v.denominator)
        return L

    # ---- valid upper bounds on exchanged heat -----------------------------------
    def U_pair(self, i, j):
        """Maximum heat that hot stream i can transfer to cold stream j in ANY feasible solution, when only
        the two streams are considered (a relaxation, hence a valid upper bound).  Single-hot single-cold
        transshipment: greedy cascade from the hottest interval is optimal (exchange argument: serving the
        current interval as much as possible never hurts, since unused hot heat stays available below)."""
        avail = Fraction(0)
        tot = Fraction(0)
        for t in range(self.k):
            avail += self.QH[i][t]
            x = min(avail, self.QC[j][t])
            tot += x
            avail -= x
        return tot

    def u_interval(self, i, j, t):
        """valid bound on heat from hot i received by cold j inside interval t"""
        return min(self.QC[j][t], self.hot_cum(i, t))

    def allowed(self, i, j):
        return self.U_pair(i, j) > 0


def parse_dat(path, name=None):
    num = r"[-+]?\d+(?:\.\d*)?(?:[eE][-+]?\d+)?"
    with open(path) as f:
        txt = f.read()
    def grab(key):
        mo = re.search(r"\b" + key + r"\s*=\s*(" + num + ")", txt)
        return mo.group(1) if mo else None
    n, m, k = int(grab("n")), int(grab("m")), int(grab("k"))
    cost = grab("Cost")
    QH = [[Fraction(0)] * k for _ in range(n)]
    QC = [[Fraction(0)] * k for _ in range(m)]
    R = [None] * (k + 1)
    seenH, seenC = set(), set()
    for line in txt.splitlines():
        line = line.strip()
        mo = re.match(r"Q([HC])\[(\d+)\]:(.*)$", line)
        if mo:
            kind, idx, rest = mo.group(1), int(mo.group(2)), mo.group(3).split()
            assert len(rest) % 2 == 0, line
            target = QH if kind == "H" else QC
            (seenH if kind == "H" else seenC).add(idx)
            for a in range(0, len(rest), 2):
                assert rest[a][0] == "T"
                t = int(rest[a][1:])
                v = Fraction(rest[a + 1])
                assert v >= 0
                assert target[idx][t] == 0, "duplicate entry " + line
                target[idx][t] = v
            continue
        mo = re.match(r"R\[(\d+)\]\s*=\s*(" + num + ")", line)
        if mo:
            R[int(mo.group(1))] = Fraction(mo.group(2))
    assert seenH == set(range(n)) and seenC == set(range(m)), (path, seenH, seenC)
    return Instance(name or path.split("/")[-1].replace(".dat", ""), n, m, k, QH, QC, R,
                    Fraction(cost) if cost else None)


def consistency_report(inst):
    """Exact consistency checks of the data file.  Returns dict of facts."""
    res = inst.residuals()
    out = {}
    out["total_hot"] = sum(inst.hot_total(i) for i in range(inst.n))
    out["total_cold"] = sum(inst.cold_total(j) for j in range(inst.m))
    out["balanced"] = out["total_hot"] == out["total_cold"]
    out["cascade_nonneg"] = all(r >= 0 for r in res)
    out["R_matches_file"] = all(inst.R[t] is not None and inst.R[t] == res[t] for t in range(inst.k + 1))
    out["R_maxdev"] = max(abs((inst.R[t] or 0) - res[t]) for t in range(inst.k + 1))
    out["empty_hot"] = [i for i in range(inst.n) if inst.hot_total(i) == 0]
    out["empty_cold"] = [j for j in range(inst.m) if inst.cold_total(j) == 0]
    out["lcm_den"] = inst.denominators_lcm()
    return out
