# V3 copy of round3/dichotomy/scripts/certify_s2.py (agent D2), unchanged except: imports from the original
# folder and writes its list to verification/V3/logs/ instead of the author's logs folder.
"""Agent D2 (round 3): exact rate-independent certificates for persistence of the dual N* at s = 2, on the 18 861 WR
networks with (H1)-(H3) and s = 2 of the first-round records (E2's list, rebuilt by extract_undecided.py).

Environment: Python 3.12 (WSL Fedora 44), standard library only; exact rational arithmetic; LPs by exact
Fourier-Motzkin elimination (E2's e2lib, unchanged); Q1's s1lib (unchanged) for product duals.
Usage (from round3/dichotomy/scripts/): python certify_s2.py [--nproc P]
Input: ~/gacD2-data/s2_all.pkl.gz (written by extract_undecided.py).  Output: ../logs/certify_s2.log (stdout) and
../logs/undecided_s2.txt (the networks left undecided).

Claims (PASS/FAIL):
  W   no realizable siphon facet of any of the 18 861 duals is DOM (Theorem W, proved for every WR network);
  K   no network satisfies both a NONE certificate and the ALL conditions (internal consistency);
  E2a every network certified ALL by E2 (no critical siphon) gets ALL here;
  E2b no network certified NONE by E2 gets ALL here;
  F   (Conjecture F, checked on the family) no network with a realizable OUT siphon facet is certified ALL;
  D   (dichotomy, checked) every decided network is ALL or NONE; the counts on E2's 14 825 undecided networks.
"""
import os
import sys
import gzip
import time
import pickle
import platform
from collections import Counter
from multiprocessing import Pool

HERE = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..", "round3", "dichotomy", "scripts"))
MYLOGS = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "logs"))
sys.path.insert(0, HERE)
import d2cert  # noqa: E402
import d2lib  # noqa: E402


def work(q):
    n, rx = q["n"], q["rx"]
    r = d2cert.analyse(n, rx)
    N = d2lib.Net(n, rx)
    F = r["F"] if r["F"] else d2cert.facets(N)
    sip_types = sorted({f["type"] for f in F if f["siphon"]})
    kinds = sorted({c[0] for c in r["none"]})
    return {"n": n, "rx": rx, "e2": q["cat"], "verdict": r["verdict"], "kinds": kinds, "sip": sip_types,
            "uv": r["uv"], "uf": r["uf"], "cycle": r["cycle"], "none": r["none"]}


def size(o):
    return (o["n"] + len(o["rx"]), len(o["rx"]), sum(sum(a) for a, b in o["rx"]), o["rx"])


def main():
    t0 = time.time()
    argv = sys.argv[1:]
    nproc = int(argv[argv.index("--nproc") + 1]) if "--nproc" in argv else 10
    print(f"# certify_s2.py  Python {platform.python_version()} on {platform.system()}; exact rational arithmetic")
    recs = pickle.load(gzip.open(os.path.expanduser("~/gacD2-data/s2_all.pkl.gz")))
    print(f"# {len(recs)} WR networks with (H1)-(H3) and s = 2 (E2's list); E2 undecided: "
          f"{sum(q['cat'] == 'UNDECIDED' for q in recs)}")
    with Pool(nproc) as P:
        out = P.map(work, recs, chunksize=10)
    fails = 0

    def claim(name, ok, msg):
        nonlocal fails
        print(f"{'PASS' if ok else 'FAIL'}: {name}: {msg}")
        fails += 0 if ok else 1
    ndom = sum("DOM" in o["sip"] for o in out)
    claim("W", ndom == 0, f"DOM siphon facets in {ndom} of {len(out)} duals")
    nconf = sum(o["verdict"] == "CONFLICT" for o in out)
    claim("K", nconf == 0, f"{nconf} networks with both a NONE certificate and the ALL conditions")
    e2all = [o for o in out if o["e2"] == "ALL"]
    claim("E2a", all(o["verdict"] == "ALL" for o in e2all), f"{sum(o['verdict'] == 'ALL' for o in e2all)} of "
          f"{len(e2all)} E2-ALL networks certified ALL")
    e2none = [o for o in out if o["e2"].startswith("NONE")]
    claim("E2b", not any(o["verdict"] == "ALL" for o in e2none), f"E2-NONE networks: "
          f"{dict(Counter(o['verdict'] for o in e2none))}")
    without = [o for o in out if "OUT" in o["sip"]]
    claim("F", not any(o["verdict"] == "ALL" for o in without),
          f"networks with a realizable OUT siphon facet: {len(without)}, verdicts "
          f"{dict(Counter(o['verdict'] for o in without))}")
    und = [o for o in out if o["e2"] == "UNDECIDED"]
    cnt = Counter(o["verdict"] for o in und)
    claim("D", set(cnt) <= {"ALL", "NONE", "UNDECIDED"},
          f"on E2's {len(und)} undecided networks: ALL {cnt['ALL']}, NONE {cnt['NONE']}, SOME 0 (none certified), "
          f"undecided {cnt['UNDECIDED']}")
    kc = Counter()
    for o in und:
        if o["verdict"] == "NONE":
            kc["+".join(o["kinds"])] += 1
    print(f"# NONE certificates used (on E2's undecided networks; several may apply, all listed): {dict(sorted(kc.items()))}")
    first = Counter()
    order = ["B", "E", "H", "CM", "L", "S"]
    for o in und:
        if o["verdict"] == "NONE":
            first[next(k for k in order if k in o["kinds"])] += 1
    print(f"# NONE by first applicable certificate (order B, E, H, CM, L, S): {dict(first)}")
    allc = Counter(o["verdict"] for o in out)
    print(f"# all {len(out)} networks: {dict(allc)}")
    print("# smallest examples (E2's undecided networks):")
    for v in ("ALL", "NONE", "UNDECIDED"):
        for kind in ([None] if v != "NONE" else order):
            cand = [o for o in und if o["verdict"] == v and (kind is None or (o["kinds"] and
                    next(k for k in order if k in o["kinds"]) == kind))]
            if not cand:
                continue
            o = min(cand, key=size)
            N = d2lib.Net(o["n"], o["rx"])
            print(f"  {v}{'' if kind is None else ' [' + kind + ']'}: N = {d2lib.describe(o['n'], o['rx'])}")
            print(f"      N* = {d2lib.describe(N.r, N.dual_rx(), [f'k{j+1}' for j in range(N.r)])}")
            if v == "NONE":
                print(f"      witness: {o['none'][0]}")
            if v == "UNDECIDED":
                print(f"      open features: facets {o['uf']}, vertices {o['uv']}, cycle possible {o['cycle']}")
    rest = sorted([o for o in und if o["verdict"] == "UNDECIDED"], key=size)
    cat = Counter((tuple(sorted({f'{t}/{fl}' for Z, t, fl, mu in o['uf']})), bool(o["uv"]), bool(o["cycle"]))
                  for o in rest)
    print(f"# undecided by open feature (facet types/flows, vertex open, cycle possible): {dict(cat)}")
    with open(os.path.join(MYLOGS, "rerun_D2_undecided_s2.txt"), "w") as fh:
        fh.write(f"# {len(rest)} networks left undecided by certify_s2.py (smallest first): N | open features\n")
        for o in rest:
            fh.write(f"{d2lib.describe(o['n'], o['rx'])} | facets {o['uf']} vertices {o['uv']} cycle {o['cycle']}\n")
    print(f"# list written to verification/V3/logs/rerun_D2_undecided_s2.txt")
    print(f"# runtime {time.time() - t0:.1f} s")
    print("ALL CLAIMS PASS" if fails == 0 else f"{fails} FAIL")
    sys.exit(1 if fails else 0)


if __name__ == "__main__":
    main()
