"""Finite exact checks, not proofs, on a precisely specified WR family.

Environment: Python 3.12+, standard library and repository crn.py.
Expected runtime: <30 s. Families: three species, six degree-two complexes;
all unions of 1, 2 or 3 distinct reversible pairs (15+105+455 networks),
and every directed simple 3- and 4-cycle, with a fixed smallest first vertex.
PASS/FAIL per aggregated universal claim; exact witnesses written as JSON.
"""
from itertools import combinations, permutations
from pathlib import Path
import json
import platform
import sys
import time

sys.dont_write_bytecode=True
from check_gac_duality import Network, criterion, incidence, rank

ROOT=Path(__file__).resolve().parents[1]

def family():
    cs=[(2,0,0),(0,2,0),(0,0,2),(1,1,0),(1,0,1),(0,1,1)]
    pairs=list(combinations(cs,2))
    for size in (1,2,3):
        for chosen in combinations(pairs,size):
            yield Network(3,[edge for a,b in chosen for edge in [(a,b),(b,a)]])
    for size in (3,4):
        for chosen in combinations(range(len(cs)),size):
            for tail in permutations(chosen[1:]):
                cycle=(chosen[0],)+tail
                yield Network(3,[(cs[cycle[i]],cs[cycle[(i+1)%size]]) for i in range(size)])

def main():
    print(f"Environment: Python {platform.python_version()}, {platform.platform()}; exact integer/Fraction arithmetic; expected runtime <30 s")
    start=time.perf_counter()
    total=eligible=0
    row_ok=cycle_ok=adls_ok=True
    counts={"both_criteria":0,"primal_only":0,"dual_only":0,"neither":0,"dual_wr":0}
    witnesses={}
    for N in family():
        total+=1
        if not N.rows_distinct():
            continue
        eligible+=1
        D=N.dual()
        row_ok &= D.weakly_reversible()==(sorted(map(tuple,N.Ye))==sorted(map(tuple,N.Yp)))
        cycle_ok &= N.n-N.rank()==D.r-rank(incidence(D))+D.deficiency()
        c,d=criterion(N),criterion(D)
        key="both_criteria" if c and d else "primal_only" if c else "dual_only" if d else "neither"
        counts[key]+=1
        counts["dual_wr"]+=int(D.weakly_reversible())
        if c!=d and key not in witnesses:
            witnesses[key]={"primal":N.reactions,"dual":D.reactions,"primal_siphons":[sorted(s) for s in N.minimal_siphons()],"dual_siphons":[sorted(s) for s in D.minimal_siphons()],"primal_conservation_rays":N.conservation_rays(),"dual_conservation_rays":D.conservation_rays()}
    results=[("family size 705", total==705), ("row-permutation dual-WR equivalence",row_ok), ("nullity decomposition",cycle_ok)]
    for label,ok in results:
        print(("PASS " if ok else "FAIL ")+label)
    print(f"Family total {total}; H3-eligible {eligible}; counts {json.dumps(counts,sort_keys=True)}")
    data={"family_total":total,"eligible":eligible,"counts":counts,"witnesses":witnesses}
    (ROOT/"certificates"/"finite_family.json").write_text(json.dumps(data,indent=2)+"\n",encoding="utf-8")
    print(f"Witnesses: {json.dumps(witnesses,sort_keys=True)}")
    print(f"Runtime: {time.perf_counter()-start:.3f} s")
    return 0 if all(ok for label,ok in results) else 1

if __name__=="__main__":
    sys.exit(main())
