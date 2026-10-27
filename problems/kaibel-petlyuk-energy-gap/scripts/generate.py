#!/usr/bin/env python3
"""Generate exact brackets; Python >= 3.10, standard library, < 2 seconds.

The separate default checker reads the output and checks signs, not this
root-finding procedure. Run only to regenerate certificates/examples.json.
"""
import json
from fractions import Fraction as Q
from verify import ROOT, Interval, functions, values


def bisect(f, left, right, iterations):
    # Endpoints may be poles, so evaluate only strictly interior midpoints.
    left, right = Q(left), Q(right)
    for _ in range(iterations):
        middle = (left + right) / 2
        value = f(middle)
        if value == 0:
            return Interval(middle)
        if value < 0:
            left = middle
        else:
            right = middle
    return Interval(left, right)


def enclosing_strings(x, digits=12):
    scale = 10**digits
    lo = (x.lo * scale).__floor__() - 1
    hi = (x.hi * scale).__ceil__() + 1
    return [str(Q(lo, scale)), str(Q(hi, scale))]


def generate(name, z, q, bottleneck):
    alpha = list(map(Q, (6, 4, 2, 1)))
    z, q = list(map(Q, z)), Q(q)
    a, b, c, d = alpha
    h, g = functions(alpha, z)
    s = 1 - q
    f = lambda x: h(x) + g(x) - s
    roots = {key: bisect(f, left, right, 180)
             for key, left, right in (("theta_A", b, a),
                                      ("theta_B", c, b), ("theta_C", d, c))}
    t = roots["theta_B"]
    u = h((t.lo + t.hi) / 2)
    # Round outwards at 30 decimal places to dominate uncertainty in u.
    roots["phi"] = Interval(*enclosing_strings(bisect(lambda x: h(x) - u, b, a, 140), 30))
    roots["psi"] = Interval(*enclosing_strings(bisect(lambda x: g(x) + u - s, d, c, 140), 30))
    result = values(alpha, z, q, roots)
    return dict(name=name, alpha=list(map(str, alpha)), z=list(map(str, z)),
                q=str(q), bottleneck=bottleneck,
                roots={k: [str(v.lo), str(v.hi)] for k, v in roots.items()},
                bounds={k: enclosing_strings(v) for k, v in result.items()
                        if not k.endswith("_reduced")})


if __name__ == "__main__":
    cases = [generate("source_table_1", ("3/10", "1/5", "1/5", "3/10"), "1", "K_A"),
             generate("saturated_vapor", ("3/10", "1/5", "1/5", "3/10"), "0", "either"),
             generate("two_phase_feed", ("3/10", "1/5", "1/5", "3/10"), "1/2", "either"),
             generate("bottom_bottleneck", ("1/20", "1/20", "2/5", "1/2"), "1", "K_C"),
             generate("rational_middle_root", ("3/10", "1/5", "1/5", "3/10"), "3/20", "either"),
             generate("superheated_parameter", ("3/10", "1/5", "1/5", "3/10"), "-1/2", "either"),
             generate("subcooled_parameter", ("3/10", "1/5", "1/5", "3/10"), "3/2", "either")]
    path = ROOT / "certificates" / "examples.json"
    path.write_text(json.dumps(dict(schema=1, cases=cases), indent=2) + "\n", encoding="utf-8")
    print("WROTE", path)
