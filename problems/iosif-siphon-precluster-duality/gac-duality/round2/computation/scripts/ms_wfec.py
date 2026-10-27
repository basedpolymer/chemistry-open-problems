"""Multistationarity of well-formed equal-conflict (WF-EC) networks, decided exactly (task 3 of BRIEF-E2; Q3).

Environment: Python 3.12 (WSL Fedora 44), standard library only, exact rational arithmetic.  Networks are generated
by agent P's enumerator (../../../petri/scripts/enumerate_ec.py, gen_networks, imported unchanged): every EC network
up to species permutation with the given number of species, molecularity and number of reactions, connected,
conservative (H1) and consistent (H2).  Kept: WF-EC = s = m - 1 (m = number of source complexes).
Usage (from round2/computation/): python scripts/ms_wfec.py N MAXMOL RMAX ORDINARY [NPROC]

Criterion (proved in notes.md, Prop. E2.3; it is the classical sign criterion for toric steady states, Pérez Millán-
Dickenstein-Shiu-Conradi 2012, Müller-Feliu-Regensburger-Conradi-Shiu-Dickenstein 2016).  For a WF-EC network, by
agent P's Theorem T the positive steady states for rates k are E_k = {x > 0 : x^{y_c} = lambda u_c(k)}, i.e.
log E_k = log x* + L^perp with L = span{y_c - y_c'} (sources), and every u > 0 occurs as u(k) for some k.  Hence:
some k and some compatibility class carry two distinct positive steady states  <=>  there are d in S = im Gamma and
w in L^perp, nonzero, with sign(d) = sign(w).  If L = S (Birch case) this never happens.
Witness (exact, rational): w integer in L^perp with the common sign vector, rho_i = 2^{w_i}, x2_i = d_i/(rho_i - 1)
(x2_i = 1 where w_i = 0), x1 = rho * x2, k_t = v_t / x1^{y_c(t)} with v > 0 in ker Gamma.  Claims (PASS/FAIL):
  M1 every witness satisfies f_k(x1) = f_k(x2) = 0 exactly, x1 != x2, x1 - x2 in S;
  M2 Birch case L = S => no common sign vector (consistency check of the code);
  M3 (agent P, Theorem T) for every WF-EC network ker A_k is one-dimensional and positive at k = v / x1^{y} ... checked
     through M1 on the witnesses.
Reported: counts (WF-EC, L = S, L != S, multistationary), smallest multistationary network with its witness.
"""
import os
import sys
import time
import platform
from fractions import Fraction
from itertools import permutations, product
from multiprocessing import Pool

import e2lib
from e2lib import CNet, crn

sys.path.insert(0, e2lib.PSCRIPTS)
import enumerate_ec as EE  # noqa: E402


def sign_vectors_common(S_basis, Lperp_eqs, n):
    """Find sigma != 0 with d in span(S_basis), w in {w: e.w = 0 for e in Lperp_eqs}, sign d = sign w = sigma.
    Returns (sigma, d, w) or None.  Exhaustive over sigma up to global sign."""
    k = len(S_basis)
    seen = set()
    for sigma in product((-1, 0, 1), repeat=n):
        if not any(sigma):
            continue
        neg = tuple(-x for x in sigma)
        if neg in seen:
            continue
        seen.add(sigma)
        # w: variables w_0..w_{n-1}
        eqs = [(list(e), 0) for e in Lperp_eqs]
        ineqs = []
        for i, sg in enumerate(sigma):
            unit = [int(t == i) for t in range(n)]
            if sg == 0:
                eqs.append((unit, 0))
            elif sg > 0:
                ineqs.append(([-x for x in unit], -1))
            else:
                ineqs.append((unit, -1))
        w = e2lib.lp_point(eqs, ineqs, n)
        if w is None:
            continue
        # d = sum_l lam_l S_basis[l]
        eqs2, ineqs2 = [], []
        for i, sg in enumerate(sigma):
            row = [S_basis[l][i] for l in range(k)]
            if sg == 0:
                eqs2.append((row, 0))
            elif sg > 0:
                ineqs2.append(([-x for x in row], -1))
            else:
                ineqs2.append((row, -1))
        lam = e2lib.lp_point(eqs2, ineqs2, k)
        if lam is None:
            continue
        d = [sum(lam[l] * S_basis[l][i] for l in range(k)) for i in range(n)]
        return sigma, d, w
    return None


def integer_scale(v):
    from math import lcm
    den = 1
    for x in v:
        den = lcm(den, Fraction(x).denominator)
    return [int(Fraction(x) * den) for x in v]


def analyse(n, rx):
    N = CNet(n, rx)
    G = N.Gamma
    r = len(rx)
    srcs = sorted({a for a, b in rx})
    s = N.rank()
    Lgen = [[Fraction(c[i] - srcs[0][i]) for i in range(n)] for c in srcs[1:]]
    Scols = [[Fraction(G[i][j]) for i in range(n)] for j in range(r)]
    birch = crn.rank(Lgen + Scols) == s and crn.rank(Lgen) == s
    # basis of S
    Sb = []
    for col in Scols:
        if crn.rank(Sb + [col]) > len(Sb):
            Sb.append(col)
    res = sign_vectors_common(Sb, Lgen, n)
    out = dict(birch=birch, s=s, m=len(srcs), multi=res is not None)
    if res is None:
        return out
    sigma, d, w = res
    wi = integer_scale(w)
    x2, x1 = [], []
    for i in range(n):
        if wi[i] == 0:
            x2.append(Fraction(1))
            x1.append(Fraction(1))
        else:
            rho = Fraction(2) ** wi[i]
            x2.append(Fraction(d[i]) / (rho - 1))
            x1.append(rho * x2[-1])
    # rates: v > 0 in ker Gamma (sum of flux rays), k_t = v_t / x1^{y_c(t)}
    rays = N.flux_rays()
    v = [sum(Fraction(z[j]) for z in rays) for j in range(r)]
    kk = []
    for j, (a, b) in enumerate(rx):
        mon = Fraction(1)
        for i in range(n):
            mon *= x1[i] ** a[i]
        kk.append(v[j] / mon)

    def f(x):
        out = [Fraction(0)] * n
        for j, (a, b) in enumerate(rx):
            rate = kk[j]
            for i in range(n):
                rate *= x[i] ** a[i]
            for i in range(n):
                out[i] += rate * (b[i] - a[i])
        return out
    diff = [x1[i] - x2[i] for i in range(n)]
    ok = (all(t == 0 for t in f(x1)) and all(t == 0 for t in f(x2)) and any(diff)
          and crn.rank(Sb + [diff]) == len(Sb) and all(t > 0 for t in x1 + x2))
    out.update(sigma=sigma, x1=x1, x2=x2, k=kk, witness_ok=ok)
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
        if N.rank() != e2lib.n_sources(rx) - 1:
            continue
        res = analyse(n, rx)
        res["rx"] = rx
        res["H3"] = N.rows_distinct()
        res["WR"] = N.weakly_reversible()
        out.append(res)
    return out


def main():
    n, m, rmax, ordinary = (int(a) for a in sys.argv[1:5])
    nproc = int(sys.argv[5]) if len(sys.argv) > 5 else 10
    t0 = time.time()
    print(f"# ms_wfec.py n={n} molecularity<={m} r<={rmax} ordinary={ordinary}; Python {platform.python_version()} "
          f"on {platform.system()}; exact arithmetic; generator: agent P's enumerate_ec.gen_networks")
    nparts = sum(1 for _ in EE.set_partitions(list(range(n))))
    allres = []
    with Pool(nproc) as pool:
        for chunk in pool.imap_unordered(worker, [(n, m, rmax, ordinary, p) for p in range(nparts)]):
            allres.extend(chunk)
    tot = len(allres)
    birch = sum(r["birch"] for r in allres)
    multi = [r for r in allres if r["multi"]]
    m2fail = [r for r in allres if r["birch"] and r["multi"]]
    m1fail = [r for r in multi if not r["witness_ok"]]
    print(f"# WF-EC networks: {tot} (with (H3): {sum(r['H3'] for r in allres)}, WR: {sum(r['WR'] for r in allres)})")
    print(f"# Birch case L = S: {birch}; L != S: {tot - birch}")
    print(f"# multistationary for some rates (common sign vector): {len(multi)} "
          f"(with (H3): {sum(r['H3'] for r in multi)})")
    print(f"{'PASS' if not m1fail else 'FAIL'} M1: exact witnesses (two positive steady states in one class) verified "
          f"for {len(multi) - len(m1fail)}/{len(multi)}")
    print(f"{'PASS' if not m2fail else 'FAIL'} M2: Birch case never multistationary ({birch} networks)")
    multi.sort(key=lambda r: (len(r["rx"]), sum(map(sum, (a for a, b in r["rx"]))), r["rx"]))
    for r in multi[:3]:
        print(f"  multistationary: {e2lib.describe(n, r['rx'])}   s={r['s']} m={r['m']} H3={r['H3']}")
        print(f"     sign vector {r['sigma']}; x1 = {[str(x) for x in r['x1']]}; x2 = {[str(x) for x in r['x2']]}")
        print(f"     k = {[str(x) for x in r['k']]}")
    nonb = [r for r in allres if not r["birch"]]
    nonb.sort(key=lambda r: (len(r["rx"]), r["rx"]))
    for r in nonb[:2]:
        print(f"  smallest L != S (unique by the sign criterion: {not r['multi']}): {e2lib.describe(n, r['rx'])}")
    ok = not m1fail and not m2fail
    print(f"{'PASS' if ok else 'FAIL'}: all claims")
    print(f"# runtime {time.time() - t0:.1f} s")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
