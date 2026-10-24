#!/usr/bin/env python3
"""Exhaustive check (in a finite box) of the classification of case-9 and case-10 networks used in result.md, Section 2.

Environment: Python 3.12/3.13 + numpy (no other dependency). Single thread. Runtime: about 1 minute (WSL/container) for the default
box T = 20 (all targets with both coefficients <= 20, i.e. 440^3 ~ 85 million networks per case).

Networks (BBH arXiv:2304.02303v1 eq. (13); notation c_i, d_i = X-, Y-components of the i-th reaction vector):
  case 9 : sources 2X, X+Y, Y      ; case 10 : sources 2X, X+Y, 0
  target i = source i + (c_i, d_i), any target with nonnegative coefficients <= T, target != source.

Every network is put in exactly one class, in this order:
  rank<2      : rank of Gamma = [c; d] is < 2                         (no periodic orbit: trajectories lie on lines)
  trivial     : rank 2 and u = c x d has a zero entry or entries of both signs (dynamically trivial)
  saddle      : u < 0 (every positive equilibrium is a saddle)
  c1<=0       : u > 0, c1 <= 0  (Dulac 1/(xy), BBH Lemma 4 arXiv = Lemma 3 published)
  LVA         : u > 0, c1 > 0, d1 = 0 (case 9 only; Dulac 1/(xy) gives div > 0)
  HOPF        : BBH Theorem 9 (arXiv) = Theorem 3 (published) Andronov-Hopf condition, implemented verbatim
  9a, 9b, 10a : the three families of result.md (Liénard argument)
Claims checked (PASS/FAIL), for each case:
  C1  no network is left unclassified (every network with u>0, c1>0, d1>0 that is not BBH-Hopf is in 9a/9b/10a);
  C2  BBH-Hopf  <=>  u>0, c1>0, d1>0, A := c1 u1 - c3 u3 > 0   (the trace criterion of result.md, Lemma 2.2);
  C3  every network in 9a/9b/10a has k := c3 u3/(c1 u1) >= 1 and d2 in {-1, 0} (hypotheses of Theorem 3.1);
  C4  LVA networks are exactly 2X -> (2+c1)X, X+Y -> (1+d2)Y, Y -> 0 with c1, d2 >= 1 (case 9), none in case 10;
  C5  the BBH example 2X -> 3X+2Y, X+Y -> 0, Y -> X+Y is in family 9a with k = 1.
Exits with status 1 on any FAIL.
"""
import sys
import time
import numpy as np

T = int(sys.argv[1]) if len(sys.argv) > 1 else 20
fails = 0


def check(cond, msg):
    global fails
    print(("PASS " if cond else "FAIL ") + msg, flush=True)
    fails += (not cond)


def vectors(src, T):
    a, b = src
    out = []
    for p in range(T + 1):
        for q in range(T + 1):
            if (p, q) != (a, b):
                out.append((p - a, q - b))
    return np.array(out, dtype=np.int64)


def classify(case, T):
    srcs = [(2, 0), (1, 1), (0, 1)] if case == 9 else [(2, 0), (1, 1), (0, 0)]
    V1, V2, V3 = (vectors(s, T) for s in srcs)
    counts = {}
    examples = {}
    c2v, d2v = V2[:, 0][:, None], V2[:, 1][:, None]
    c3v, d3v = V3[:, 0][None, :], V3[:, 1][None, :]
    bad = {"C1": 0, "C2": 0, "C3": 0, "C4": 0}
    for (c1, d1) in V1:
        c1 = int(c1); d1 = int(d1)
        c2, d2, c3, d3 = np.broadcast_arrays(c2v, d2v, c3v, d3v)
        u1 = c2 * d3 - c3 * d2
        u2 = c3 * d1 - c1 * d3
        u3 = c1 * d2 - c2 * d1
        rank2 = (u1 != 0) | (u2 != 0) | (u3 != 0)  # rank 2 iff c x d != 0
        pos = (u1 > 0) & (u2 > 0) & (u3 > 0)
        neg = (u1 < 0) & (u2 < 0) & (u3 < 0)
        A = c1 * u1 - c3 * u3
        # BBH Hopf condition, verbatim (exact integer arithmetic; c1>0, c2=-1, c3>0 makes all denominators nonzero)
        if case == 9:
            # (1/2)(d3/c3 + d1/c1) < d2/c2 < d1/c1 with c2 = -1, multiply by 2 c1 c3 > 0, resp. c1 > 0
            br1 = (c1 > 0) & (c2 == -1) & (c3 > 0) & (d1 > 0) & (d2 >= -1) & (d3 >= -1) \
                & (c1 * d3 + c3 * d1 < -2 * c1 * c3 * d2) & (-d2 * c1 < d1)
            br2 = (c1 > 0) & (c2 == -1) & (c3 == 0) & (d1 > 0) & (d2 >= -1) & (d3 == -1) & (-d2 * c1 < d1)
            hopf = br1 | br2
        else:
            hopf = (c1 > 0) & (c2 == -1) & (c3 > 0) & (d1 > 0) & (d2 == -1) & (d3 >= 0) \
                & (c1 * d3 + c3 * d1 < -2 * c1 * c3 * d2) & (-d2 * c1 < d1)
        cls = np.full(c2.shape, "", dtype="<U8")
        cls[~rank2] = "rank<2"
        m = rank2 & ~pos & ~neg & (cls == ""); cls[m] = "trivial"
        m = neg & (cls == ""); cls[m] = "saddle"
        m = pos & (c1 <= 0) & (cls == ""); cls[m] = "c1<=0"
        m = pos & (c1 > 0) & (d1 == 0) & (cls == ""); cls[m] = "LVA"
        # C4: LVA networks shape
        mL = cls == "LVA"
        if case == 10 and mL.any():
            bad["C4"] += int(mL.sum())
        if case == 9 and mL.any():
            okL = (c2 == -1) & (c3 == 0) & (d3 == -1) & (d2 >= 1)
            bad["C4"] += int((mL & ~okL).sum())
        # Hopf networks must not have been classified above (they have u > 0, c1 > 0, d1 > 0)
        bad["C2"] += int((hopf & (cls != "")).sum())
        m = hopf & (cls == ""); cls[m] = "HOPF"
        # C2: Hopf <=> u>0, c1>0, d1>0, A>0
        crit = pos & (c1 > 0) & (d1 > 0) & (A > 0)
        bad["C2"] += int((crit != hopf).sum())
        if case == 9:
            f9a = (c2 == -1) & (d2 == -1) & (c1 >= 1) & (c3 >= 1) & (d1 > c1) & (d3 >= -1) & (d3 < c3) \
                & (2 * c1 * c3 - c1 * d3 - c3 * d1 <= 0)
            f9b = (c2 == -1) & (d2 == 0) & (d3 == -1) & (c1 >= 1) & (c3 >= 1) & (d1 >= 1) & (c3 * d1 >= c1)
            m = f9a & (cls == ""); cls[m] = "9a"
            m = f9b & (cls == ""); cls[m] = "9b"
        else:
            f10a = (c2 == -1) & (d2 == -1) & (c1 >= 1) & (c3 >= 1) & (d1 > c1) & (d3 >= 0) & (d3 < c3) \
                & (2 * c1 * c3 - c1 * d3 - c3 * d1 <= 0)
            m = f10a & (cls == ""); cls[m] = "10a"
        bad["C1"] += int((cls == "").sum())
        # C3: in the families, k = c3 u3/(c1 u1) >= 1 (i.e. c3 u3 >= c1 u1 > 0) and d2 in {-1,0}
        fam = (cls == "9a") | (cls == "9b") | (cls == "10a")
        okk = (c1 * u1 > 0) & (c3 * u3 >= c1 * u1) & ((d2 == -1) | (d2 == 0)) & pos
        bad["C3"] += int((fam & ~okk).sum())
        for key in np.unique(cls):
            n = int((cls == key).sum())
            counts[key] = counts.get(key, 0) + n
            if key not in examples:
                i, j = np.argwhere(cls == key)[0]
                examples[key] = ((c1, d1), (int(c2[i, j]), int(d2[i, j])), (int(c3[i, j]), int(d3[i, j])))
    return counts, examples, bad


def fmt_net(case, ex):
    srcs = [(2, 0), (1, 1), (0, 1)] if case == 9 else [(2, 0), (1, 1), (0, 0)]
    def cplx(p, q):
        parts = []
        if p: parts.append(("%dX" % p) if p > 1 else "X")
        if q: parts.append(("%dY" % q) if q > 1 else "Y")
        return " + ".join(parts) if parts else "0"
    return ",  ".join("%s -> %s" % (cplx(*s), cplx(s[0] + v[0], s[1] + v[1])) for s, v in zip(srcs, ex))


t0 = time.time()
print(f"box: all target coefficients in 0..{T}")
for case in (9, 10):
    counts, examples, bad = classify(case, T)
    total = sum(counts.values())
    print(f"\n=== case {case}: {total} networks")
    for key in ["rank<2", "trivial", "saddle", "c1<=0", "LVA", "HOPF", "9a", "9b", "10a", ""]:
        if key in counts:
            print(f"   {key or 'UNCLASSIFIED':8s} {counts[key]:9d}   e.g. {fmt_net(case, examples[key])}")
    check(bad["C1"] == 0, f"C1 case {case}: every network is classified (unclassified: {bad['C1']})")
    check(bad["C2"] == 0, f"C2 case {case}: BBH-Hopf <=> (u>0, c1>0, d1>0, c1u1-c3u3>0) (mismatches: {bad['C2']})")
    check(bad["C3"] == 0, f"C3 case {case}: families 9a/9b/10a have k>=1 and d2 in {{-1,0}} (violations: {bad['C3']})")
    check(bad["C4"] == 0, f"C4 case {case}: LVA class has the stated shape (violations: {bad['C4']})")
# C5 the example: c = (1,-1,1), d = (2,-1,0)
c1, c2, c3, d1, d2, d3 = 1, -1, 1, 2, -1, 0
u1, u2, u3 = c2 * d3 - c3 * d2, c3 * d1 - c1 * d3, c1 * d2 - c2 * d1
in9a = (d1 > c1) and (-1 <= d3 < c3) and (2 * c1 * c3 - c1 * d3 - c3 * d1 <= 0)
check(in9a and c3 * u3 == c1 * u1, f"C5 BBH example is in family 9a with k = c3u3/(c1u1) = {c3*u3}/{c1*u1} = 1")
print(f"\nruntime {time.time()-t0:.1f} s")
print("ALL PASS" if not fails else f"{fails} FAILURES")
sys.exit(1 if fails else 0)
