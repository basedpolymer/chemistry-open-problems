"""Count N = 5 regular-column configurations left open in the permissive model after: sign rule (signcheck),
merge-free theorem (Agent G, Thm G3).  A (configuration, j) pair is open iff it is sign-hard AND its mixed
sub-network has a merge node (a mixed stream produced by two mixed columns).  Uses theory/data/signhard_N5.json
(written by signcheck.dump(5, ...)).  Environment: WSL Fedora, Python 3.  Runtime: seconds.  Counting only, no PASS/FAIL."""
import sys, json
sys.path.insert(0, '../../computation/scripts')
from hsmodel import all_configurations
N = 5
A = all_configurations(N)
d = json.load(open('../data/signhard_N5.json'))
open_pairs = []
for j in range(1, N):
    for cid in d['hard'][str(j)]:
        c = A[cid]
        mixed = {s for s in c['streams'] if s[0] <= j < s[1]}
        merge = any(len([p for (t, p) in c['producers'].get(s, []) if p in mixed]) >= 2 for s in mixed)
        if merge:
            open_pairs.append((cid, j))
cfgs = sorted(set(c for c, j in open_pairs))
basics = sorted(set(A[c]['basic_id'] for c in cfgs))
bj = sorted(set((A[c]['basic_id'], j) for c, j in open_pairs))
print('open (variant, j) pairs:', len(open_pairs), '; configurations with an open j:', len(cfgs), 'of', len(A),
      '; basic configurations:', len(basics), '; (basic, j) pairs:', len(bj))
