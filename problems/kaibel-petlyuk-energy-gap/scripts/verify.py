#!/usr/bin/env python3
"""Check the Kaibel/Petlyuk examples with exact rational arithmetic.

Python >= 3.10, standard library only; expected runtime < 1 second.
No numerical root finder is used by this checker. All brackets are checked
by endpoint signs and the strictly positive derivatives proved in paper.tex.
Run from any directory: python scripts/verify.py
"""
from fractions import Fraction as Q
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]


class Interval:
    def __init__(self, lo, hi=None):
        self.lo, self.hi = Q(lo), Q(lo if hi is None else hi)
        if self.lo > self.hi:
            raise ValueError("reversed interval")

    @staticmethod
    def coerce(other):
        return other if isinstance(other, Interval) else Interval(other)

    def __add__(self, other):
        other = self.coerce(other)
        return Interval(self.lo + other.lo, self.hi + other.hi)

    __radd__ = __add__

    def __neg__(self):
        return Interval(-self.hi, -self.lo)

    def __sub__(self, other):
        return self + -self.coerce(other)

    def __rsub__(self, other):
        return self.coerce(other) + -self

    def __mul__(self, other):
        other = self.coerce(other)
        ends = [a * b for a in (self.lo, self.hi)
                for b in (other.lo, other.hi)]
        return Interval(min(ends), max(ends))

    __rmul__ = __mul__

    def __truediv__(self, other):
        other = self.coerce(other)
        if other.lo <= 0 <= other.hi:
            raise ValueError("division by an interval containing zero")
        return self * Interval(1 / other.hi, 1 / other.lo)

    def __rtruediv__(self, other):
        return self.coerce(other) / self


def maximum(*xs):
    return Interval(max(x.lo for x in xs), max(x.hi for x in xs))


def functions(alpha, z):
    def h(x):
        return sum(alpha[i] * z[i] / (alpha[i] - x) for i in (0, 1))

    def g(x):
        return sum(alpha[i] * z[i] / (alpha[i] - x) for i in (2, 3))

    return h, g


def values(alpha, z, q, roots):
    a, b, c, d = alpha
    s = 1 - q
    h, g = functions(alpha, z)
    ta, t, tc, phi, psi = [roots[k] for k in
                          ("theta_A", "theta_B", "theta_C", "phi", "psi")]
    u = h(t)
    v = u - s
    ha = a * z[0] / (a - ta)
    hc = s + d * z[3] / (tc - d)
    ka = a * z[0] / (a - phi)
    kc = s + d * z[3] / (psi - d)
    ka_reduced = u * (a - t) / (a - b)
    kc_reduced = s + v * (t - d) / (c - d)
    petlyuk = maximum(ha, u, hc)
    kaibel = maximum(ka, kc)
    return dict(H_A=ha, H_B=u, H_C=hc, K_A=ka, K_C=kc,
                Petlyuk=petlyuk, Kaibel=kaibel, gap=kaibel - petlyuk,
                K_A_reduced=ka_reduced, K_C_reduced=kc_reduced)


def require(condition, label):
    if not condition:
        raise ValueError(label)
    print("PASS", label)


def verify_case(case):
    name = case["name"]
    alpha, z = [list(map(Q, case[k])) for k in ("alpha", "z")]
    q = Q(case["q"])
    a, b, c, d = alpha
    require(a > b > c > d > 0 and all(x > 0 for x in z)
            and sum(z) == 1, name + ": admissible Underwood input")
    h, g = functions(alpha, z)
    s = 1 - q
    roots = {k: Interval(*pair) for k, pair in case["roots"].items()}
    for key, left, right in (("theta_A", b, a), ("theta_B", c, b),
                             ("theta_C", d, c)):
        r = roots[key]
        require(left < r.lo <= r.hi < right
                and h(r.lo) + g(r.lo) <= s <= h(r.hi) + g(r.hi),
                name + ": unique " + key)
    t, phi, psi = [roots[k] for k in ("theta_B", "phi", "psi")]
    u = h(t)
    v = u - s
    require(u.lo > 0 and v.lo > 0, name + ": positive prefractionator flows")
    require(b < phi.lo <= phi.hi < a and h(phi.lo) < u.lo
            and h(phi.hi) > u.hi, name + ": unique phi")
    require(d < psi.lo <= psi.hi < c and g(psi.lo) < -v.hi
            and g(psi.hi) > -v.lo, name + ": unique psi")
    result = values(alpha, z, q, roots)
    for key in ("H_A", "H_B", "H_C", "K_A", "K_C", "Petlyuk", "Kaibel", "gap"):
        bound = Interval(*case["bounds"][key])
        actual = result[key]
        require(bound.lo <= actual.lo <= actual.hi <= bound.hi
                and bound.hi - bound.lo <= Q(1, 10**10),
                name + ": certified " + key)
    for key in ("K_A", "K_C"):
        direct, reduced = result[key], result[key + "_reduced"]
        bound = Interval(*case["bounds"][key])
        require(bound.lo <= reduced.lo <= reduced.hi <= bound.hi
                and max(direct.lo, reduced.lo) <= min(direct.hi, reduced.hi),
                name + ": direct and reduced " + key)
    require(result["K_A"].lo > max(result["H_A"].hi, result["H_B"].hi)
            and result["K_C"].lo > max(result["H_C"].hi, result["H_B"].hi)
            and result["gap"].lo > 0, name + ": all strict comparisons")
    expected = case["bottleneck"]
    if expected != "either":
        other = "K_C" if expected == "K_A" else "K_A"
        require(result[expected].lo > result[other].hi,
                name + ": bottleneck " + expected)
    if name == "source_table_1":
        require(Q("1.375") < result["Petlyuk"].lo
                and result["Petlyuk"].hi < Q("1.385")
                and Q("1.825") < result["Kaibel"].lo
                and result["Kaibel"].hi < Q("1.835"),
                name + ": published rounded values 1.38 and 1.83")
    return result


def main():
    data = json.loads((ROOT / "certificates" / "examples.json").read_text())
    require(data["schema"] == 1, "certificate schema")
    require(len(data["cases"]) == 7
            and len({c["name"] for c in data["cases"]}) == 7,
            "seven distinct examples")
    for case in data["cases"]:
        verify_case(case)
    print("PASS all certificates (exact Fraction arithmetic)")


if __name__ == "__main__":
    try:
        main()
    except (ValueError, KeyError, TypeError, ZeroDivisionError) as exc:
        print("FAIL", exc)
        sys.exit(1)
