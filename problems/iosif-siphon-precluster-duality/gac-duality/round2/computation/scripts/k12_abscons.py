"""K12 (agent P's question): absolute consistency of the dual, decided exactly by a log-cone criterion.

Environment: Python 3.12 (WSL Fedora 44), standard library only, exact rational arithmetic (e2lib -> ecore -> crn).
Usage (from round2/computation/): python scripts/k12_abscons.py [--nproc K]
Input: every first-round record file (~/gacE-data/*.pkl.gz), both sides (X = N and X = N*).

Criterion (proved in notes.md, Prop. E2.4).  For a network X (n species, r reactions, Ye, Gamma) with (H2), put
V = ker Gamma cap R^r_{>0} (open flux cone, rays v^1..v^p) and K = ker Ye (subset of R^r), d = dim K.  A positive
steady state for rates k exists iff some v in V has log v - log k in im Ye^T = K^perp.  Hence X is absolutely
consistent (a positive steady state for every k) iff h : V -> R^d, h(v) = (c_l . log v)_l (c_1..c_d a basis of K) is
onto.  Decisions:
  AC0     d = 0 (Ye of full column rank): absolutely consistent.
  ACB     K subset of ker Gamma: absolutely consistent (Birch: exp(log k + K^perp) meets v0 + K subset ker Gamma).
  d = 1, c spanning K.  Writing v = sum_k lambda_k v^k (lambda > 0), h = Phi(log lambda) + O(1) with the tropical
          function Phi(alpha) = sum_j c_j max_{k : v^k_j > 0} alpha_k.  h is onto iff Phi takes both signs, iff
          sum_j c_j != 0, or Phi(1_A) > 0 and Phi(1_B) < 0 for some sets A, B of rays
          (Phi(1_A) = sum of c_j over the reactions j covered by a ray of A).
          AC1: onto; NOT1: not onto (then h is bounded on one side and some k has no positive steady state).
  OPEN    d >= 2 (not decided here; necessary condition only).
Claims (PASS/FAIL): B every weakly reversible X with a decided verdict is AC (Boros 2019); T every EC, connected
X with (H1), (H2) and a decided verdict is AC iff s = m - 1 (agent P's Theorem T); S for EC with balanced joins
and (H1), (H2), N and N* get the same decided verdict (Corollary T*); U agent P's U is AC and U* is NOT1 (manuscript
network; P checked it with sympy).
Reported: K12(a) EC networks N with unbalanced joins, (H1), (H2), connected, s = m - 1 (absolutely consistent by
Theorem T): verdicts for N*; K12(b) weakly reversible N with (H1)-(H3): verdicts for N*; smallest NOT1 examples.
"""
import sys
import time
import platform
from fractions import Fraction
from itertools import combinations
from multiprocessing import Pool

import e2lib
from e2lib import CNet, crn


def verdict(n, rx, maxrays=18):
    X = CNet(n, rx)
    r = len(rx)
    Ye = X.Ye
    K = crn.nullspace([[Fraction(x) for x in row] for row in Ye], r)
    d = len(K)
    if d == 0:
        return "AC0", None
    G = X.Gamma
    if all(all(sum(G[i][j] * c[j] for j in range(r)) == 0 for i in range(n)) for c in K):
        return "ACB", None
    if d >= 2:
        return "OPEN", d
    c = K[0]
    if sum(c) != 0:
        return "AC1", "sum"
    rays = X.flux_rays()
    p = len(rays)
    if p > maxrays:
        return "OPEN", ("rays", p)
    cover = [frozenset(k for k in range(p) if rays[k][j] > 0) for j in range(r)]
    pos = neg = None
    for size in range(1, p):
        for A in combinations(range(p), size):
            A = set(A)
            val = sum(c[j] for j in range(r) if cover[j] & A)
            if val > 0 and pos is None:
                pos = tuple(sorted(A))
            if val < 0 and neg is None:
                neg = tuple(sorted(A))
            if pos is not None and neg is not None:
                return "AC1", (pos, neg)
    return "NOT1", dict(c=[str(x) for x in c], bounded="above" if pos is None else "below", rays=p)


def info(rec):
    n, rx = rec["n"], rec["rx"]
    out = {}
    for side in ("N", "D"):
        if side == "N":
            m_, rx_ = n, rx
        else:
            m_, rx_ = len(rx), e2lib.dual_rx(n, rx)
        R = rec[side]
        if not R["consist"]:
            out[side] = None
            continue
        out[side] = verdict(m_, rx_)
    return out


def work(rec):
    n, rx = rec["n"], rec["rx"]
    v = info(rec)
    Nn = rec["N"]
    D = rec["D"]
    ec = e2lib.is_EC(rx)
    drx = e2lib.dual_rx(n, rx)
    ecd = e2lib.is_EC(drx)
    bal = all(len({x for x in a if x > 0}) <= 1 for a, b in rx)
    bald = all(len({x for x in a if x > 0}) <= 1 for a, b in drx)
    conn = e2lib.petri_connected(n, rx)
    return dict(key=(n, rx), v=v, H=e2lib.H(rec), WR=Nn["WR"], WRd=D["WR"], ec=ec, ecd=ecd, bal=bal, bald=bald,
                conn=conn, s=rec["s"], m=e2lib.n_sources(rx), md=len({a for a, b in drx}),
                h12=Nn["cons"] and Nn["consist"], h12d=D["cons"] and D["consist"])


def main():
    t0 = time.time()
    nproc = int(sys.argv[sys.argv.index("--nproc") + 1]) if "--nproc" in sys.argv else 4
    print(f"# k12_abscons.py  Python {platform.python_version()} on {platform.system()}; exact arithmetic; nproc={nproc}")
    recs = e2lib.load(e2lib.E_FILES)
    print(f"# {len(recs)} distinct networks (all first-round record files)")
    with Pool(nproc) as pool:
        res = pool.map(work, recs, chunksize=200)
    fails = {"B": [], "T": [], "S": [], "U": []}
    napp = {"B": 0, "T": 0, "S": 0, "U": 1}
    dec = lambda x: x is not None and x[0] != "OPEN"  # noqa: E731
    isac = lambda x: x[0].startswith("AC")  # noqa: E731
    tally = {}
    a_rows, b_rows = [], []
    for q in res:
        n, rx = q["key"]
        for side in ("N", "D"):
            x = q["v"][side]
            if x is None:
                continue
            tally.setdefault(side, {}).setdefault(x[0], 0)
            tally[side][x[0]] += 1
            wr = q["WR"] if side == "N" else q["WRd"]
            if wr and dec(x):
                napp["B"] += 1
                if not isac(x):
                    fails["B"].append((side, rx))
            ecx = q["ec"] if side == "N" else q["ecd"]
            h12 = q["h12"] if side == "N" else q["h12d"]
            sx = q["s"]
            mx = q["m"] if side == "N" else q["md"]
            if ecx and h12 and q["conn"] and dec(x):
                napp["T"] += 1
                if isac(x) != (sx == mx - 1):
                    fails["T"].append((side, rx))
        if q["ec"] and q["bal"] and q["h12"] and q["conn"] and dec(q["v"]["N"]) and q["v"]["D"] is not None \
                and dec(q["v"]["D"]):
            napp["S"] += 1
            if isac(q["v"]["N"]) != isac(q["v"]["D"]):
                fails["S"].append(rx)
        if q["ec"] and not q["bal"] and q["h12"] and q["conn"] and q["s"] == q["m"] - 1:
            a_rows.append(q)
        if q["H"] and q["WR"]:
            b_rows.append(q)
    # U and U*
    U = (((0, 0, 0, 1), (0, 1, 1, 0)), ((0, 0, 1, 1), (1, 0, 0, 0)), ((1, 1, 0, 0), (0, 0, 0, 2)))
    vU = verdict(4, U)
    vUs = verdict(3, e2lib.dual_rx(4, U))
    if not (vU[0].startswith("AC") and vUs[0] == "NOT1"):
        fails["U"].append((vU, vUs))
    print()
    for side in ("N", "D"):
        print(f"# verdicts, side {side} (networks with (H2) on that side): {dict(sorted(tally[side].items()))}")
    print()
    names = {"B": "weakly reversible => absolutely consistent (Boros 2019), on decided cases",
             "T": "EC, connected, (H1), (H2): AC iff s = m - 1 (Theorem T), on decided cases",
             "S": "balanced EC, (H1), (H2): same verdict for N and N* (Corollary T*)",
             "U": f"U is AC ({vU[0]}) and U* is NOT1 ({vUs[0]}, {vUs[1]})"}
    ok = True
    for k in "BTSU":
        st = "PASS" if not fails[k] else "FAIL"
        ok &= not fails[k]
        print(f"{st} {k}: {names[k]}  ({napp[k]} applicable, {len(fails[k])} failures)")
        for f in fails[k][:3]:
            print(f"     failure: {f}")
    print()

    def summary(rows, label):
        cnt = {}
        for q in rows:
            x = q["v"]["D"]
            key = "no (H2) for N*" if x is None else x[0]
            cnt[key] = cnt.get(key, 0) + 1
        print(f"{label}: {len(rows)} networks; verdicts for N*: {dict(sorted(cnt.items()))}")
        bad = [q for q in rows if q["v"]["D"] is not None and q["v"]["D"][0] == "NOT1"]
        bad.sort(key=lambda q: (q["key"][0] + len(q["key"][1]), len(q["key"][1]), q["key"][1]))
        for q in bad[:3]:
            n, rx = q["key"]
            print(f"   N* not absolutely consistent: N = {e2lib.describe(n, rx)}")
            print(f"      N* = {e2lib.describe(len(rx), e2lib.dual_rx(n, rx), [f'k{j+1}' for j in range(len(rx))])}"
                  f"   certificate {q['v']['D'][1]}")
        return cnt

    summary(a_rows, "K12(a) EC, unbalanced joins, connected, (H1), (H2), s = m - 1 (N absolutely consistent)")
    summary(b_rows, "K12(b) weakly reversible with (H1)-(H3) (N absolutely consistent, Boros 2019)")
    print()
    print(f"{'PASS' if ok else 'FAIL'}: all claims")
    print(f"# runtime {time.time() - t0:.1f} s")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
