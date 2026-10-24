"""Build an exact rational solution file from a list of matches.

usage: python3 build_solution.py <instance.dat> '<json list of [i,j]>' <out.json> [source-note]
Prunes unnecessary matches (exact max-flow), computes an exact flow, verifies it with check_solution, writes JSON.
Environment: Python >= 3.10, standard library only.  Runtime: < 10 s per instance.
"""
import sys, json, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from hen_instance import parse_dat
from hen_exact import flow_for_matches, check_solution, save_solution, prune_matches

inst = parse_dat(sys.argv[1])
M = [tuple(p) for p in json.loads(sys.argv[2])]
M = prune_matches(inst, M)
ok, q = flow_for_matches(inst, M)
assert ok
used = sorted(set((i, j) for (i, j, t) in q))
errs = check_solution(inst, used, q)
print(inst.name, "matches:", len(used), "errors:", errs[:3])
assert not errs
save_solution(sys.argv[3], inst, used, q, {"source": sys.argv[4] if len(sys.argv) > 4 else ""})
