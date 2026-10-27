"""Beyond equal conflict (round2/ec/notes.md §7.1; agent Q3, Claude Opus 5.5, effort xhigh).

On the general first-round families of agent E (WSL ~/gacE-data/all_*.pkl.gz: every network with (H1), (H2) in the
families (2,2,5), (2,2,6), (2,3,4), (3,2,4), (3,3,3), (4,1,12), (4,2,3)), select the networks N that are connected,
satisfy (H1)-(H3), and the 'rank equation' s = m - 1 (m = number of distinct source complexes), and report, for EC,
TEFC-not-EC (sources sharing a species have the same support) and not-TEFC networks (kinds 'EC', 'TEFC', 'other');
the same with the rank equation counted on Petri clusters, s = #clusters - 1 (kinds '.../cl'; 'BH/cl' =
block-homogeneous: TEFC with every source a multiple of a 0/1 vector; Corollary BH predicts no failure there):
  critN  : N has a critical siphon (record field siphons_eq_SF false; E proved the equivalence);
  critD  : N* has a critical siphon;
  multi  : (numerical, floating point) two distinct positive steady states found in one class of N for one of the
           sampled rate vectors (multistart Newton; NOT a proof - candidates for exact re-check);
smallest example of each.  Environment: Python 3.12 + numpy/scipy (WSL).  Usage: python beyond_ec.py CAP FILES...  (CAP = maximal number
of non-EC networks per kind on which the multistart search is run, smallest first).
"""
import sys
import time
import gzip
import pickle
import platform
import warnings

import numpy as np
from scipy.optimize import root

import eclib
from eclib import Network

warnings.filterwarnings("ignore")
np.seterr(all="ignore")


def load(files):
    seen = {}
    for fn in files:
        with gzip.open(fn, "rb") as fh:
            for r in pickle.load(fh):
                seen.setdefault((r["n"], tuple(r["rx"])), r)
    return list(seen.values())


def kind(N):
    srcs = sorted({a for a, b in N.reactions})
    supp = [frozenset(i for i in range(N.n) if y[i]) for y in srcs]
    ec = all(not (supp[i] & supp[j]) for i in range(len(srcs)) for j in range(i + 1, len(srcs)))
    if ec:
        return "EC", len(srcs)
    tefc = all(not (supp[i] & supp[j]) or supp[i] == supp[j] for i in range(len(srcs)) for j in range(i + 1, len(srcs)))
    return ("TEFC" if tefc else "other"), len(srcs)


def multistart(N, rng, nk=3, ncl=3, nst=20):
    G = np.array(N.Gamma, dtype=float)
    Ye = np.array(N.Ye, dtype=float)
    U, sv, Vt = np.linalg.svd(G.T)          # G.T is r x n; left kernel of G = conservation laws
    s = int((sv > 1e-9).sum())
    L = Vt[s:]                               # rows: basis of conservation laws (n - s)
    Sb = np.linalg.svd(G)[0][:, :s]          # basis of S
    for _ in range(nk):
        k = np.exp(rng.uniform(np.log(0.1), np.log(10), N.r))
        for _ in range(ncl):
            x0 = np.exp(rng.uniform(np.log(0.1), np.log(10), N.n))

            def F(q):
                x = np.exp(q)
                rates = k * np.exp(Ye.T @ q)
                return np.concatenate([L @ (x - x0), Sb.T @ (G @ rates)])
            sols = []
            for _ in range(nst):
                res = root(F, np.log(x0) + rng.normal(0, 1.5, N.n), method="hybr", options={"xtol": 1e-13})
                if not res.success or not np.all(np.isfinite(res.x)):
                    continue
                if np.max(np.abs(F(res.x))) > 1e-9:
                    continue
                x = np.exp(res.x)
                if np.min(x) < 1e-6 * np.max(x) or np.max(x) > 1e6:
                    continue                      # (near-)boundary or diverging point: not a positive steady state
                if all(np.max(np.abs(x - y) / (np.abs(y) + 1e-12)) > 1e-5 for y in sols):
                    sols.append(x)
            if len(sols) >= 2:
                return True, (list(k), list(x0), [list(x) for x in sols])
    return False, None


def main():
    t0 = time.time()
    cap = int(sys.argv[1])
    files = sys.argv[2:]
    print(f"# beyond_ec.py; Python {platform.python_version()}, numpy {np.__version__} on {platform.system()}; "
          f"multistart on at most {cap} networks per kind (smallest first)", flush=True)
    recs = load(files)
    rng = np.random.default_rng(2026)
    sel = []
    for r in recs:
        if not (r["N"]["cons"] and r["N"]["consist"] and r["N"]["rows_distinct"]):
            continue
        N = Network(r["n"], r["rx"])
        if not eclib.petri.is_connected(N):
            continue
        kd, m = kind(N)
        ncl = len([c for c in eclib.petri.clusters(N) if c[1]])
        rk = N.rank()
        if rk == ncl - 1:
            srcs = {a for a, b in N.reactions}
            bh = all(len({c for c in y if c}) == 1 for y in srcs)   # every source a multiple of a 0/1 vector
            if kd == "TEFC" and bh:
                kd = "BH"
            sel.append((N.n + N.r, kd + "/cl", r, N))
        if rk == m - 1:
            sel.append((N.n + N.r, kd if kd != "BH" else "TEFC", r, N))
    sel.sort(key=lambda t: (t[0], t[2]["rx"]))
    stats, smallest, tried = {}, {}, {}
    for size, kd, r, N in sel:
        props = {"all": True, "critN": not r["N"]["siphons_eq_SF"], "critD": not r["D"]["siphons_eq_SF"]}
        if kd != "EC" and tried.get(kd, 0) < cap:
            tried[kd] = tried.get(kd, 0) + 1
            ok, data = multistart(N, rng, nk=2, ncl=2, nst=12)
            props["multi"] = ok
            props["multi_tried"] = True
            if ok:
                props["_multi"] = data
                print(f"MULTI [{kd}] {N.describe()}  data={data}", flush=True)
        for p_, v in props.items():
            if p_.startswith("_") or not v:
                continue
            stats[(kd, p_)] = stats.get((kd, p_), 0) + 1
            if (kd, p_) not in smallest:
                smallest[(kd, p_)] = (size, N.describe())
    for key in sorted(stats):
        sz, desc = smallest[key]
        print(f"REPORT  [{key[0]}] {key[1]}: {stats[key]}  smallest: {desc}")
    print(f"# runtime {time.time() - t0:.1f} s")


if __name__ == "__main__":
    main()
