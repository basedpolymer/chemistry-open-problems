#!/usr/bin/env python3
"""Exact audit behind the "already solved" record for the McClure-Shiu question.

Written from the statements of McClure-Shiu and of White's paper, without reading or running White's
scripts. Prints PASS or FAIL for each claim and exits non-zero on any failure.

Environment: Python 3.12, standard library only (fractions, itertools, random with a fixed seed).
No floating point anywhere. Runtime: about 5 s on one core.
Output: certificates/audit-record.log (run as `python scripts/audit_record.py | tee ...`).
"""
import itertools
import random
import sys
from fractions import Fraction as Q

FAILS = []


def check(name, ok, detail=""):
    print(("PASS " if ok else "FAIL ") + name + (f"  [{detail}]" if detail else ""))
    if not ok:
        FAILS.append(name)


# ---------------------------------------------------------------- exact polynomial tools
def trim(p):
    p = list(p)
    while p and p[-1] == 0:
        p.pop()
    return p


def pmul(p, q):
    out = [Q(0)] * (len(p) + len(q) - 1)
    for i, a in enumerate(p):
        for j, b in enumerate(q):
            out[i + j] += a * b
    return out


def pfrom_roots(lead, roots):
    p = [Q(lead)]
    for r in roots:
        p = pmul(p, [-Q(r), Q(1)])
    return p


def pval(p, x):
    v = Q(0)
    for a in reversed(p):
        v = v * x + a
    return v


def pderiv(p):
    return [i * p[i] for i in range(1, len(p))]


def prem(a, b):
    a = [Q(x) for x in a]
    while len(a) >= len(b) and a:
        f = a[-1] / b[-1]
        d = len(a) - len(b)
        for i, c in enumerate(b):
            a[i + d] -= f * c
        a = trim(a[:-1]) if a[-1] == 0 else a
    return trim(a)


def sign(x):
    return (x > 0) - (x < 0)


def changes(seq):
    s = [v for v in seq if v != 0]
    return sum(1 for u, v in zip(s, s[1:]) if u != v)


def distinct_pos_roots(p):
    """Number of distinct positive roots, by a Sturm chain. p must not be the zero polynomial."""
    p = trim(p)
    while p[0] == 0:
        p = p[1:]
    if len(p) == 1:
        return 0
    chain = [p, pderiv(p)]
    while True:
        r = prem(chain[-2], chain[-1])
        if not r:
            break
        chain.append([-c for c in r])
    at0 = [sign(c[0]) for c in chain]
    atinf = [sign(c[-1]) for c in chain]
    return changes(at0) - changes(atinf)


def coeffs_of_network(reactions, kappa, deg=None):
    """Steady-state polynomial sum_i kappa_i (b_i - a_i) x^{a_i}, reactions = [(a_i, b_i)]."""
    deg = deg if deg is not None else max(a for a, _ in reactions)
    c = [Q(0)] * (deg + 1)
    for (a, b), k in zip(reactions, kappa):
        c[a] += Q(k) * (b - a)
    return c


# ---------------------------------------------------------------- A. the source
def part_a():
    print("--- A. statements of the source, exact")
    star = [(1, 0), (1, 2), (2, 1), (2, 3), (3, 2), (3, 4)]  # kappa_1L,1R,2L,2R,3L,3R
    ku = [Q(v) for v in (1, 3, 4, 1, 1, 2)]
    kv = [Q(v) for v in (3, 1, 1, 4, 2, 1)]
    cu, cv = coeffs_of_network(star, ku), coeffs_of_network(star, kv)
    target = pmul([Q(0), Q(1)], pfrom_roots(1, [1, 2]))
    check("A1 Prop 5.1 witnesses: f_U = x(x-1)(x-2), f_V = -f_U",
          trim(cu) == target and trim(cv) == [-t for t in target])
    ok = True
    for t in (Q(i, 16) for i in range(17)):
        kt = [(1 - t) * a + t * b for a, b in zip(ku, kv)]
        c = coeffs_of_network(star, kt)
        ok &= all(v > 0 for v in kt)
        ok &= c == [(1 - 2 * t) * v for v in cu]
        if t != Q(1, 2):
            ok &= distinct_pos_roots(c) == 2
    check("A2 segment kappa(U)--kappa(V): f_t = (1-2t) f_U, two roots for t != 1/2", ok)
    mid = [(a + b) / 2 for a, b in zip(ku, kv)]
    check("A3 midpoint (2,2,5/2,5/2,3/2,3/2): f = 0 identically, so every x > 0 is a steady state",
          mid == [2, 2, Q(5, 2), Q(5, 2), Q(3, 2), Q(3, 2)] and not trim(coeffs_of_network(star, mid)))
    check("A4 Lemma 3.5(1) needs f != 0: {A->0, A->2A} with equal constants has f = 0",
          not trim(coeffs_of_network([(1, 0), (1, 2)], [1, 1])))
    quoted = [(0, 1), (1, 0), (2, 3), (3, 2)]
    c = coeffs_of_network(quoted, (2, 5, 4, 1))
    check("A5 quoted network, kappa=(2,5,4,1): f = -(x-1)^2(x-2), two distinct roots",
          c == pfrom_roots(-1, [1, 1, 2]) and distinct_pos_roots(c) == 2)
    up = distinct_pos_roots(coeffs_of_network(quoted, (2, 5, 4 + Q(1, 100), 1)))
    check("A6 that point is not interior: raising kappa_3 by 1/100 leaves one positive root", up == 1)
    p = pmul(pfrom_roots(1, [1, 1]), pfrom_roots(1, [2, 2]))
    ok = p == [4, -12, 13, -6, 1] and distinct_pos_roots(p) == 2
    check("A7 brief's 'even' start: (x-1)^2(x-2)^2 has 4 sign changes, 2 roots, f >= 0 (it is a square)",
          ok)


# ---------------------------------------------------------------- B. combinatorics
def achievable(t):
    return list(itertools.product(*[(1, -1) if v == 0 else (v,) for v in t]))


def sc(sig):
    return sum(1 for u, v in zip(sig, sig[1:]) if u != v)


def orthant_graph_components(t):
    """White's Definition 2.6, reimplemented from the text."""
    verts = [p for p in achievable(t) if sc(p) >= 2]
    parent = {v: v for v in verts}

    def find(v):
        while parent[v] != v:
            parent[v] = parent[parent[v]]
            v = parent[v]
        return v

    for u, v in itertools.combinations(verts, 2):
        rest = tuple(a for a, b in zip(u, v) if a == b)
        if sc(rest) >= 2:
            parent[find(u)] = find(v)
    return len({find(v) for v in verts})


def has_alternating(t, length, first):
    """Is there an achievable alternating subsequence of the given length starting with sign `first`?"""
    need, got = first, 0
    for v in t:
        if v == 0 or v == need:
            got += 1
            need = -need
            if got == length:
                return True
    return False


def part_b():
    print("--- B. combinatorics, exhaustive")
    bad = []
    n = 0
    for s in range(1, 9):
        for t in itertools.product((1, -1, 0), repeat=s):
            n += 1
            comps = orthant_graph_components(t)
            exceptional = t == (0, 0, 0)
            if comps != (2 if exceptional else min(comps, 1)):
                bad.append(t)
    check("B1 Prop 5.2: orthant graph connected for every signature with s <= 8, except (0,0,0)",
          not bad, f"{n} signatures")

    bad = []
    n = 0
    for s in range(1, 13):
        for t in itertools.product((1, -1, 0), repeat=s):
            n += 1
            both = has_alternating(t, 3, 1) and has_alternating(t, 3, -1)
            four = has_alternating(t, 4, 1) or has_alternating(t, 4, -1)
            if both and not four and not (s == 3 and t == (0, 0, 0)):
                bad.append(t)
    check("B2 both (+,-,+) and (-,+,-) achievable => achievable 4-alternation, except (0,0,0)",
          not bad, f"{n} signatures, s <= 12")

    # Every signature that r <= 5 reactions can realise: s pure/mixed sources with s + b <= 5.
    rows = {}
    for s in range(1, 6):
        for t in itertools.product((1, -1, 0), repeat=s):
            b = t.count(0)
            if s + b <= 5:
                rows.setdefault((s, b), 0)
                rows[(s, b)] += 1
    excep = [(s, b) for (s, b) in rows if s == 3 and b == 3]
    check("B3 no signature realisable with r <= 5 is exceptional (three mixed sources need r >= 6)",
          not excep, f"(s,b) counts {sorted(rows.items())}")

    pairs = [(a, b) for a in range(5) for b in range(5) if a != b]

    def signature(net):
        src = sorted({a for a, _ in net})
        out = []
        for a0 in src:
            d = {(b > a) for a, b in net if a == a0}
            out.append(0 if len(d) == 2 else (1 if True in d else -1))
        return tuple(out)

    tally = {}
    for r in (4, 5, 6):
        tally[r] = [0, 0]
        for net in itertools.combinations(pairs, r):
            tally[r][0] += 1
            if signature(net) == (0, 0, 0):
                tally[r][1] += 1
    check("B4 all networks on complexes 0A..4A: none exceptional for r = 4, 5; some for r = 6",
          tally[4][1] == 0 and tally[5][1] == 0 and tally[6][1] > 0, str(tally))
    star = ((1, 0), (1, 2), (2, 1), (2, 3), (3, 2), (3, 4))
    check("B5 the source's six-reaction network is exceptional", signature(star) == (0, 0, 0))


# ---------------------------------------------------------------- C. White's Lemma 3.1
def section(a, sig, tau):
    """Formulas (14)-(18) of White's paper, with exact rationals."""
    s = len(a)
    w = [tau ** e for e in a]
    ch = [j for j in range(s - 1) if sig[j] != sig[j + 1]]
    p, q, r = ch[0], ch[0] + 1, ch[1] + 1
    K = [Q(0)] * s
    K[p], K[q], K[r] = (sig[p] * (w[r] - w[q]), sig[p] * (w[p] - w[r]), sig[p] * (w[q] - w[p]))
    others = [k for k in range(s) if k not in (p, q, r)]
    lam = 1 + (sum(abs(w[k] - w[r]) for k in others) / (w[r] - w[q])
               + sum(abs(w[p] - w[k]) for k in others) / (w[q] - w[p])) / (w[r] - w[p])
    out = [lam * v for v in K]
    for k in others:
        out[k] += sig[k]
        out[p] += sig[k] * (w[k] - w[r]) / (w[r] - w[p])
        out[r] += sig[k] * (w[p] - w[k]) / (w[r] - w[p])
    return out, w


def part_c():
    print("--- C. White's Lemma 3.1, section s(tau), exact rationals")
    rng = random.Random(20261031)
    n, ok = 0, True
    for _ in range(4000):
        s = rng.randint(3, 8)
        a = sorted(rng.sample(range(0, 14), s))
        sig = tuple(rng.choice((1, -1)) for _ in range(s))
        if sc(sig) < 2:
            continue
        tau = Q(rng.randint(5, 40), 4)
        if tau <= 1:
            continue
        d, w = section(a, sig, tau)
        ok &= sum(d) == 0 and sum(x * y for x, y in zip(d, w)) == 0
        ok &= all(u * v > 0 for u, v in zip(d, sig))
        n += 1
    check("C1 s(tau) lies in F(tau): both roots and the prescribed sign pattern", ok, f"{n} samples")


# ---------------------------------------------------------------- D. the short route
def box_contains(c, signature):
    return all(t == 0 or u * t > 0 for u, t in zip(c, signature))


def part_d():
    print("--- D. fibres in the whole box, not in one orthant")
    # Five reactions, signature (0,+,-,+): {A->0, A->2A, 2A->3A, 3A->2A, 4A->5A}.
    net = [(1, 0), (1, 2), (2, 3), (3, 2), (4, 5)]
    sig = (0, 1, -1, 1)

    def lift(cs, big):
        return [big, cs[0] + big, cs[1], -cs[2], cs[3]]

    def poly(a):
        return pmul([Q(0), Q(1)], pmul([Q(a), Q(1)], pfrom_roots(1, [1, 2])))

    ok, signs = True, set()
    for i in range(29):
        a = Q(1, 2) - Q(i, 8)                  # a runs over [-3, 1/2]
        c = poly(a)
        cs = c[1:]
        ok &= box_contains(cs, sig)
        ok &= pval(c, Q(1)) == 0 and pval(c, Q(2)) == 0
        ok &= distinct_pos_roots(c) >= 2
        k = lift(cs, Q(10))
        ok &= all(v > 0 for v in k) and coeffs_of_network(net, k) == c
        signs.add(sign(cs[0]))
    check("D1 r = 5: segment in rate space from pattern (+,+,-,+) to (-,+,-,+) with roots 1, 2 throughout",
          ok and signs == {1, 0, -1})

    quoted = [(0, 1), (1, 0), (2, 3), (3, 2)]
    ok = True
    for i in range(33):
        a = 1 + Q(i, 16)                       # a runs over [1, 3]
        c = pfrom_roots(-1, [1, 2, a])
        k = (c[0], -c[1], c[2], -c[3])
        ok &= all(v > 0 for v in k) and coeffs_of_network(quoted, k) == c
        ok &= distinct_pos_roots(c) == (2 if a in (1, 2) else 3)
    check("D2 r = 4: segment from a three-root point to the degenerate point (2,5,4,1), roots 1, 2 kept", ok)

    # Exceptional case: for exponents a1 < a2 < a3 the fibre F(tau) is a line through 0.
    rng = random.Random(7)
    ok = True
    for _ in range(300):
        a = sorted(rng.sample(range(0, 12), 3))
        tau = Q(rng.randint(5, 30), 4)
        w = [tau ** e for e in a]
        d = [w[2] - w[1], w[0] - w[2], w[1] - w[0]]      # the kernel of [1,1,1; w1,w2,w3]
        ok &= sum(d) == 0 and sum(x * y for x, y in zip(d, w)) == 0
        ok &= d[0] > 0 > d[1] and d[2] > 0
    check("D3 three sources: F(tau) is a line, its two rays have patterns (+,-,+) and (-,+,-)", ok)

    p = pmul(pfrom_roots(1, [1, 1, 1]), pfrom_roots(1, [2]))
    check("D4 nondegenerate reading: (x-1)^3(x-2) has two sign changes but one simple root",
          p == [2, -7, 9, -5, 1] and distinct_pos_roots(p) == 2 and pval(pderiv(pderiv(p)), Q(1)) == 0)


def main():
    part_a()
    part_b()
    part_c()
    part_d()
    print("ALL PASS" if not FAILS else f"FAILED: {FAILS}")
    return 1 if FAILS else 0


if __name__ == "__main__":
    sys.exit(main())
