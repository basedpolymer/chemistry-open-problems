"""Tests of agent Q3's candidates (round2/ec/candidates.md) on exhaustive EC families.

Environment: Python 3.12 (WSL Fedora 44), standard library only, exact rational arithmetic.  Networks: agent P's
generator (petri/scripts/enumerate_ec.gen_networks, imported unchanged): every EC network up to species
permutation, connected, (H1), (H2), with the given species count, molecularity, reaction count, ordinary flag.
Usage (from round2/computation/): python scripts/q3_tests.py N MAXMOL RMAX ORDINARY [NPROC]

Claims (WF-EC = EC + connected + (H1) + (H2) + s = m - 1):
  Q3-C    every WF-EC network has no critical siphon (minimal siphons = supports of extreme conservation laws);
  Q3-Uc   sigma(S) cap sigma(D^perp) = {0} (D = span of source differences), decided twice: by exhaustive sign
          vectors with exact LPs (ms_wfec.sign_vectors_common) and by the maximal-minor criterion (all products
          det S_I det D_I of the same sign, not all zero); the two must agree and both must say 'injective';
  Q3-Ub   at two sampled steady states (x rational, v > 0 in ker Gamma, k_t = v_t / x^{y_t}) the characteristic
          polynomial of the Jacobian is t^(n-s) q(t) with q Hurwitz (exact Hurwitz determinants);
  Q3-T*'  N* is absolutely consistent (exact log-cone verdict of k12_abscons.verdict; OPEN counted separately);
  Q3-T*cex (report) EC, connected, (H1), (H2), s = m: how many have an absolutely consistent dual, split by
          balanced / unbalanced joins (Q3 says a balanced one cannot).
"""
import sys
import time
import random
import platform
from fractions import Fraction
from itertools import permutations, combinations
from multiprocessing import Pool

import e2lib
from e2lib import CNet, crn, ecore
import ms_wfec
import k12_abscons

sys.path.insert(0, e2lib.PSCRIPTS)
import enumerate_ec as EE  # noqa: E402


def det(M):
    M = [[Fraction(x) for x in row] for row in M]
    m = len(M)
    d = Fraction(1)
    for c in range(m):
        p = next((i for i in range(c, m) if M[i][c] != 0), None)
        if p is None:
            return Fraction(0)
        if p != c:
            M[c], M[p] = M[p], M[c]
            d = -d
        d *= M[c][c]
        for i in range(c + 1, m):
            if M[i][c] != 0:
                f = M[i][c] / M[c][c]
                M[i] = [x - f * y for x, y in zip(M[i], M[c])]
    return d


def charpoly(J):
    """Faddeev-LeVerrier: coefficients [1, c1, ..., cn] of det(tI - J)."""
    n = len(J)
    I = [[Fraction(int(i == j)) for j in range(n)] for i in range(n)]
    M = [[Fraction(0)] * n for _ in range(n)]
    coeffs = [Fraction(1)]
    for k in range(1, n + 1):
        # M_k = J M_{k-1} + c_{k-1} I
        AM = [[sum(J[i][l] * M[l][j] for l in range(n)) for j in range(n)] for i in range(n)]
        M = [[AM[i][j] + coeffs[-1] * I[i][j] for j in range(n)] for i in range(n)]
        JM = [[sum(J[i][l] * M[l][j] for l in range(n)) for j in range(n)] for i in range(n)]
        c = -sum(JM[i][i] for i in range(n)) / k
        coeffs.append(c)
    return coeffs


def hurwitz(q):
    """q = [1, a1, ..., as] monic; True iff all Hurwitz determinants are > 0."""
    s = len(q) - 1
    if s == 0:
        return True
    a = lambda k: q[k] if 0 <= k <= s else Fraction(0)  # noqa: E731
    H = [[a(2 * (j + 1) - (i + 1)) for j in range(s)] for i in range(s)]
    return all(det([row[:k] for row in H[:k]]) > 0 for k in range(1, s + 1))


def minor_criterion(Sb, Db, n):
    s = len(Sb)
    sg = set()
    for I in combinations(range(n), s):
        p = det([[Sb[l][i] for i in I] for l in range(s)]) * det([[Db[l][i] for i in I] for l in range(s)])
        if p > 0:
            sg.add(1)
        elif p < 0:
            sg.add(-1)
    return len(sg) == 1  # injective (no common nonzero sign vector)


def analyse(n, rx, seed):
    rng = random.Random(seed)
    N = CNet(n, rx)
    r = len(rx)
    s = N.rank()
    m = e2lib.n_sources(rx)
    out = dict(n=n, rx=rx, s=s, m=m)
    bal = all(len({x for x in a if x > 0}) <= 1 for a, b in rx)
    out["bal"] = bal
    vD = k12_abscons.verdict(r, e2lib.dual_rx(n, rx))
    out["vD"] = vD[0]
    if s != m - 1:
        return out
    out["wf"] = True
    # Q3-C
    sy = ecore.siphon_masks(N)
    mins = [ecore.to_set(S) for S in ecore.minimal_masks(sy)]
    sf = set(ecore.semiflow_supports(N))
    out["C"] = set(mins) == sf
    # Q3-Uc
    G = N.Gamma
    srcs = sorted({a for a, b in rx})
    Db = [[Fraction(c[i] - srcs[0][i]) for i in range(n)] for c in srcs[1:]]
    Sb = ms_wfec.basis([[Fraction(G[i][j]) for i in range(n)] for j in range(r)]) if hasattr(ms_wfec, "basis") \
        else None
    if Sb is None:
        Sb = []
        for j in range(r):
            col = [Fraction(G[i][j]) for i in range(n)]
            if crn.rank(Sb + [col]) > len(Sb):
                Sb.append(col)
    common = ms_wfec.sign_vectors_common(Sb, Db, n)
    inj_minor = minor_criterion(Sb, Db, n)
    out["Uc_lp"] = common is None
    out["Uc_minor"] = inj_minor
    # Q3-Ub
    rays = N.flux_rays()
    ok_b = True
    for _ in range(2):
        x = [Fraction(rng.randint(1, 9), rng.randint(1, 9)) for _ in range(n)]
        v = [sum(Fraction(rng.randint(1, 5)) * Fraction(z[j]) for z in rays) for j in range(r)]
        J = [[sum(Fraction(G[i][t]) * v[t] * rx[t][0][j] / x[j] for t in range(r)) for j in range(n)]
             for i in range(n)]
        cp = charpoly(J)
        if any(c != 0 for c in cp[s + 1:]):
            ok_b = False
            continue
        if not hurwitz(cp[:s + 1]):
            ok_b = False
    out["Ub"] = ok_b
    return out


def worker(args):
    n, m, rmax, ordinary, pidx = args
    perms = list(permutations(range(n)))
    seen = set()
    out = []
    for rx in EE.gen_networks(n, m, rmax, ordinary, pidx):
        rx = tuple(sorted(rx))
        can = EE.canonical(rx, perms)
        if can != rx or can in seen:
            continue
        seen.add(can)
        if not e2lib.petri_connected(n, rx):
            continue
        N = CNet(n, rx)
        if not N.is_consistent() or not N.is_conservative():
            continue
        s, mm = N.rank(), e2lib.n_sources(rx)
        if s not in (mm - 1, mm):
            continue
        out.append(analyse(n, rx, hash(rx) & 0xffffffff))
    return out


def main():
    n, m, rmax, ordinary = (int(a) for a in sys.argv[1:5])
    nproc = int(sys.argv[5]) if len(sys.argv) > 5 else 4
    t0 = time.time()
    print(f"# q3_tests.py n={n} molecularity<={m} r<={rmax} ordinary={ordinary}; Python {platform.python_version()} "
          f"on {platform.system()}; exact arithmetic")
    nparts = sum(1 for _ in EE.set_partitions(list(range(n))))
    res = []
    with Pool(nproc) as pool:
        for chunk in pool.imap_unordered(worker, [(n, m, rmax, ordinary, p) for p in range(nparts)]):
            res.extend(chunk)
    wf = [q for q in res if q.get("wf")]
    sm = [q for q in res if not q.get("wf")]
    ok = True

    def claim(name, rows, pred, desc):
        nonlocal ok
        bad = [q for q in rows if not pred(q)]
        ok &= not bad
        print(f"{'PASS' if not bad else 'FAIL'} {name}: {desc}  ({len(rows) - len(bad)}/{len(rows)})")
        bad.sort(key=lambda q: (len(q["rx"]), q["rx"]))
        for q in bad[:2]:
            print(f"     counterexample: {e2lib.describe(q['n'], q['rx'])}")
    print(f"# WF-EC networks: {len(wf)} (unbalanced joins: {sum(not q['bal'] for q in wf)}); EC with s = m: {len(sm)}")
    claim("Q3-C", wf, lambda q: q["C"], "no critical siphon (minimal siphons = extreme semiflow supports)")
    claim("Q3-Uc", wf, lambda q: q["Uc_lp"] and q["Uc_minor"], "no common sign vector (LP enumeration and maximal "
          "minors agree, injective)")
    claim("Q3-Uc'", wf, lambda q: q["Uc_lp"] == q["Uc_minor"], "the two injectivity tests agree")
    claim("Q3-Ub", wf, lambda q: q["Ub"], "Jacobian at 2 sampled steady states: t^(n-s) q(t), q Hurwitz")
    dec = [q for q in wf if q["vD"] != "OPEN"]
    claim("Q3-T*'", dec, lambda q: q["vD"].startswith("AC"), "N* absolutely consistent (exact verdict; "
          f"{len(wf) - len(dec)} OPEN not counted)")
    for lab, rows in (("balanced", [q for q in sm if q["bal"]]), ("unbalanced", [q for q in sm if not q["bal"]])):
        cnt = {}
        for q in rows:
            cnt[q["vD"]] = cnt.get(q["vD"], 0) + 1
        print(f"# Q3-T*cex report, EC with s = m, {lab} joins: {len(rows)} networks; verdicts for N*: "
              f"{dict(sorted(cnt.items()))}")
        acs = sorted([q for q in rows if q["vD"].startswith("AC")], key=lambda q: (len(q["rx"]), q["rx"]))
        for q in acs[:2]:
            print(f"     N* absolutely consistent: N = {e2lib.describe(q['n'], q['rx'])}")
    print(f"{'PASS' if ok else 'FAIL'}: all claims")
    print(f"# runtime {time.time() - t0:.1f} s")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
