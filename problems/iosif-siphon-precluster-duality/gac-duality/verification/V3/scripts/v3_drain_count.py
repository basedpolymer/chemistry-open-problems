"""Verifier V3 (Claude Opus 5.5, effort xhigh), 27 October 2026: independent count of duals with a drain species
(E's Prop. E.4 hypotheses, own implementation) among the 18 861 WR networks with (H1)-(H3), s = 2, of
~/gacD2-data/s2_all.pkl.gz (E2's list, rebuilt by D2).  E2 reports 199 (round2/computation/notes.md, section 5.1).
Usage: python v3_drain_count.py.  Environment: Python 3.12, standard library."""
import os, sys, gzip, pickle, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import v3lib as L

def drain_species(n, R):
    out = []
    for i in range(n):
        for s, p in R:
            if s[i] >= 1 and all(s[l] == 0 for l in range(n) if l != i) and p[i] < s[i]:
                if all(s2[i] > s[i] for s2, p2 in R if p2[i] > s2[i]):
                    out.append(i)
    return sorted(set(out))

t0 = time.time()
recs = pickle.load(gzip.open(os.path.expanduser("~/gacD2-data/s2_all.pkl.gz")))
cnt = 0
for q in recs:
    nd, Rd = L.dual(q["n"], q["rx"])
    if drain_species(nd, Rd):
        cnt += 1
ok = cnt == 199
print("networks %d; duals with a drain species %d" % (len(recs), cnt))
print(("PASS" if ok else "FAIL") + ": equals E2's count 199")
print("# runtime %.1f s" % (time.time() - t0))
sys.exit(0 if ok else 1)
