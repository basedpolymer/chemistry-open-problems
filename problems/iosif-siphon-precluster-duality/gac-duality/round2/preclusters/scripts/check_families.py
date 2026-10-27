"""Corrected preclusters on the first-round families (agent Q2, Claude Opus 5.5, xhigh).

Environment: Python 3.12 (WSL Fedora 44), standard library only; exact arithmetic (crn.py).  Reads agent E's
records (pickled lists of dicts with keys 'rx', 'n') from a directory given on the command line (default
~/gacE-data).  Usage: python check_families.py FAMILY_FILE [FAMILY_FILE ...] [--data DIR] [--limit K]
Prints PASS/FAIL per claim and family, counts, smallest examples; exit status 1 on any FAIL.

Claims checked on every network N (records satisfy (H1), (H2)) and on its dual N* (which then satisfies them too):
  K1  chain  P <= Pcs <= pi <= pi+ <= CFup  (pi+ = sound lower bound, CFup = meet of kappa_V, upper bound);
  K2  pi = Pcs v M*  (M* = full-coupling classes);
  K3  M* <= pi, and  [pi = M*] <=> [CFup = M*] <=> [every doubling pair is fully coupled];
  K4  EC + connected  =>  CFup is one block  (Theorem 4.1: CF is one block);
  K5  weakly reversible with deficiency 0  =>  CFup = linkage classes and pi+ refines them  (Theorem 4.2);
  K6  isolation (every internal cycle inside a CF block, tested with pi+ when exact)  =>  dim ker A_V = #blocks of
      kappa_V at a generic V  (Prop. 5.2: PMDSC Condition 3.1);
  K7  Q3 for the corrected preclusters of N*: every block of CF(N*) lies in a minimal siphon of N  <=>  N has no twin
      pair and every MIPS of N lies in a minimal siphon (Prop. 6.3).
Reports (no claim): how often CF is determined (pi+ = CFup), how often P != Pcs, Pcs != pi, pi != CF, kappa_gen
!= CF; siphon relations Q1, Q2, Q4, Q5 and R23 for CF on both sides; toric (Condition 3.1) frequency.
"""
import os
import sys
import gzip
import time
import pickle
import random
import platform

import pclust as pc
from pclust import Network, refines, normalize_partition

sys.path.insert(0, pc.PETRI_DIR)
import petri  # noqa: E402

T0 = time.time()
FAILS = 0


def lc_partition(N):
    cs = N.complexes()
    idx = {c: i for i, c in enumerate(cs)}
    parent = list(range(len(cs)))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    for a, b in N.reactions:
        parent[find(idx[a])] = find(idx[b])
    blocks = {}
    for t, (a, b) in enumerate(N.reactions):
        blocks.setdefault(find(idx[a]), []).append(t)
    return normalize_partition(blocks.values())


def side(N, rng):
    o = {}
    o["P"] = pc.iosif_P(N)
    o["Pcs"] = pc.P_cs(N)
    o["pi"] = pc.flat_closure(N)
    o["M"] = pc.full_coupling(N)
    o["lo"] = pc.pi_plus(N)
    V = pc.random_flux(N, rng)
    o["kgen"] = pc.kappa(N, V)
    up, _ = pc.cf_upper(N, o["lo"], rng, ntrials=2)
    o["up"] = pc.meet(up, o["kgen"])
    X = pc.A_matrix(N, V)
    o["dimker"] = len(pc.nullspace(X, len(X[0]))) if X and X[0] else 0
    o["nk"] = len({frozenset(pc.source_classes(N)[i]) for i in range(len(pc.source_classes(N)))})
    # number of kappa_V blocks counted on source classes = number of reaction blocks
    o["nkappa"] = len(o["kgen"])
    return o


def q_relations(blocks, sets, union_ok=True):
    """Siphon-type relations between a partition `blocks` and a family `sets` (minimal siphons / T-traps)."""
    blocks = [frozenset(b) for b in blocks]
    sets = [frozenset(s) for s in sets]
    Q1 = set(blocks) == set(sets)
    Q2 = all(any(s <= b for b in blocks) for s in sets)
    Q3 = all(any(b <= s for s in sets) for b in blocks)
    return Q1, Q2, Q3


def main():
    global FAILS
    args = sys.argv[1:]
    data = os.path.expanduser("~/gacE-data")
    limit = None
    if "--data" in args:
        i = args.index("--data")
        data = args[i + 1]
        del args[i:i + 2]
    if "--limit" in args:
        i = args.index("--limit")
        limit = int(args[i + 1])
        del args[i:i + 2]
    print(f"# check_families.py; Python {platform.python_version()} on {platform.system()}; families {args}")
    rng = random.Random(271027)
    for fam in args:
        t1 = time.time()
        with gzip.open(os.path.join(data, fam), "rb") as fh:
            recs = pickle.load(fh)
        recs = [r for r in recs if r["N"]["cons"] and r["N"]["consist"]]
        if limit:
            recs = recs[:limit]
        cnt = {}
        ex = {}
        fails = {}

        def inc(k, rx=None, n=None):
            cnt[k] = cnt.get(k, 0) + 1
            if rx is not None and (k not in ex or (len(rx), n) < (len(ex[k][0]), ex[k][1])):
                ex[k] = (rx, n)

        def claim(k, cond, rx, n):
            if not cond:
                fails[k] = fails.get(k, 0) + 1
                if k not in ex:
                    ex["FAIL " + k] = (rx, n)
            inc("checked " + k)

        for rec in recs:
            rx, n = rec["rx"], rec["n"]
            N = Network(n, rx)
            D = N.dual()
            H3 = rec["N"]["rows_distinct"]
            inc("networks")
            if H3:
                inc("networks with (H3)")
            OO = {}
            for sd, M_ in (("N", N), ("D", D)):
                o = side(M_, rng)
                OO[sd] = o
                tag = sd
                claim("K1 " + tag, refines(o["P"], o["Pcs"]) and refines(o["Pcs"], o["pi"]) and
                      refines(o["pi"], o["lo"]) and refines(o["lo"], o["up"]), rx, n)
                claim("K2 " + tag, o["pi"] == pc.join(o["Pcs"], o["M"], M_.r), rx, n)
                allc = all(any({i, j} <= B for B in o["M"]) for i, j in M_.doubling_pairs())
                claim("K3 " + tag, refines(o["M"], o["pi"]) and (o["pi"] == o["M"]) == allc
                      and (o["up"] == o["M"]) == allc, rx, n)
                ec = petri.is_EC(M_) and petri.is_connected(M_)
                if ec:
                    claim("K4 " + tag, len(o["up"]) == 1, rx, n)
                    if len(o["lo"]) == 1:
                        inc(f"EC sides with pi+ one block ({tag})")
                if M_.weakly_reversible() and M_.deficiency() == 0:
                    lcp = lc_partition(M_)
                    claim("K5 " + tag, o["up"] == lcp and refines(o["lo"], lcp), rx, n)
                    if o["lo"] == lcp:
                        inc(f"WR0 sides with pi+ = linkage classes ({tag})")
                exact = o["lo"] == o["up"]
                if exact:
                    inc(f"CF determined ({tag})")
                else:
                    inc(f"CF undetermined ({tag})", rx, n)
                if o["P"] != o["Pcs"]:
                    inc(f"P != Pcs ({tag})", rx, n)
                if o["Pcs"] != o["pi"]:
                    inc(f"Pcs != pi ({tag})", rx, n)
                if exact and o["pi"] != o["lo"]:
                    inc(f"pi != CF ({tag})", rx, n)
                if exact and o["kgen"] != o["lo"]:
                    inc(f"kappa_gen != CF ({tag})", rx, n)
                ics = M_.internal_cycles()
                iso = exact and all(any(c <= B for B in o["lo"]) for c in ics)
                if iso:
                    inc(f"isolation (R23 for CF) ({tag})", rx, n)
                    claim("K6 " + tag, o["dimker"] == o["nkappa"], rx, n)
                if o["dimker"] == o["nkappa"]:
                    inc(f"Condition 3.1 at generic V ({tag})")
                else:
                    inc(f"Condition 3.1 fails at generic V ({tag})", rx, n)
            # siphon relations, species side: CF(N*) vs minimal siphons of N; reaction side: CF(N) vs MS(N*)
            oD, oN = OO["D"], OO["N"]
            msN = N.minimal_siphons()
            msD = D.minimal_siphons()
            twins = D.doubling_pairs()
            MN = N.mips_ray()
            q3_crit = (not twins) and all(any(B <= s for s in msN) for B in MN)
            claim("K7", q_relations(oD["lo"], msN)[2] == q3_crit and q_relations(oD["up"], msN)[2] == q3_crit, rx, n)
            if oD["lo"] == oD["up"]:
                q1, q2, q3 = q_relations(oD["lo"], msN)
                for nm_, val in (("Q1", q1), ("Q2", q2), ("Q3", q3)):
                    inc(f"species side {nm_} holds" if val else f"species side {nm_} fails",
                        None if val else rx, n)
                    if val and H3:
                        inc(f"species side {nm_} holds, (H3)")
            if oN["lo"] == oN["up"]:
                q1, q2, q3 = q_relations(oN["lo"], msD)
                for nm_, val in (("Q1", q1), ("Q2", q2), ("Q3", q3)):
                    inc(f"reaction side {nm_} holds" if val else f"reaction side {nm_} fails",
                        None if val else rx, n)
        print(f"## family {fam}: {len(recs)} networks with (H1),(H2); {time.time() - t1:.1f} s")
        for k in sorted(cnt):
            extra = ""
            if k in ex and not k.startswith("checked"):
                extra = f"   smallest: n={ex[k][1]} {ex[k][0]}"
            print(f"   {k}: {cnt[k]}{extra}")
        for k in sorted(c for c in cnt if c.startswith("checked ")):
            kk = k[len("checked "):]
            nf = fails.get(kk, 0)
            print(f"{'PASS' if nf == 0 else 'FAIL'} {kk} on {fam}: {cnt[k]} checked, {nf} failures"
                  + (f"; first {ex.get('FAIL ' + kk)}" if nf else ""))
            FAILS += nf > 0
    print(f"# runtime {time.time() - T0:.1f} s; {FAILS} FAIL")
    sys.exit(1 if FAILS else 0)


if __name__ == "__main__":
    main()
