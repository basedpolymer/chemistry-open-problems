"""Q3-T*m: EC networks N (connected, (H1), (H2)) with s = m and unbalanced joins: is N* always absolutely
consistent?  Exact log-cone verdict (k12_abscons.verdict) of N* on an exhaustive EC family.

Environment: Python 3.12 (WSL Fedora 44), standard library only, exact.  Generator: agent P's
enumerate_ec.gen_networks (imported unchanged).  Usage (from round2/computation/): python scripts/tstar_m.py N MAXMOL
RMAX [NPROC].  A search (no claim): prints the verdict counts and the smallest network whose dual is NOT1 (not
absolutely consistent, certified) or OPEN (d >= 2, undecided here).
"""
import sys
import time
import platform
from itertools import permutations
from multiprocessing import Pool

import e2lib
from e2lib import CNet
import k12_abscons

sys.path.insert(0, e2lib.PSCRIPTS)
import enumerate_ec as EE  # noqa: E402


def worker(args):
    n, m, rmax, pidx = args
    perms = list(permutations(range(n)))
    seen = set()
    out = []
    for rx in EE.gen_networks(n, m, rmax, 0, pidx):
        rx = tuple(sorted(rx))
        can = EE.canonical(rx, perms)
        if can != rx or can in seen:
            continue
        seen.add(can)
        if all(len({x for x in a if x > 0}) <= 1 for a, b in rx):
            continue
        if not e2lib.petri_connected(n, rx):
            continue
        N = CNet(n, rx)
        if N.rank() != e2lib.n_sources(rx) or not N.is_consistent() or not N.is_conservative():
            continue
        v = k12_abscons.verdict(len(rx), e2lib.dual_rx(n, rx))
        out.append((n, rx, v[0], v[1]))
    return out


def main():
    n, m, rmax = (int(a) for a in sys.argv[1:4])
    nproc = int(sys.argv[4]) if len(sys.argv) > 4 else 4
    t0 = time.time()
    print(f"# tstar_m.py n={n} molecularity<={m} r<={rmax}; Python {platform.python_version()} on {platform.system()}")
    nparts = sum(1 for _ in EE.set_partitions(list(range(n))))
    with Pool(nproc) as pool:
        res = [q for ch in pool.imap_unordered(worker, [(n, m, rmax, p) for p in range(nparts)]) for q in ch]
    cnt = {}
    for q in res:
        cnt[q[2]] = cnt.get(q[2], 0) + 1
    print(f"# EC, unbalanced, connected, (H1), (H2), s = m: {len(res)} networks; verdicts for N*: {cnt}")
    for lab in ("NOT1", "OPEN"):
        bad = sorted([q for q in res if q[2] == lab], key=lambda q: (len(q[1]), q[1]))
        for q in bad[:3]:
            print(f"  {lab}: N = {e2lib.describe(q[0], q[1])}; N* = "
                  f"{e2lib.describe(len(q[1]), e2lib.dual_rx(q[0], q[1]), [f'k{j+1}' for j in range(len(q[1]))])}; {q[3]}")
    print(f"# runtime {time.time() - t0:.1f} s")


if __name__ == "__main__":
    main()
