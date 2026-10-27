"""A2 targeted exact search for rate-dependent stable linear vertices.

Environment: Python 3.14 standard library; exact Fraction arithmetic.
Expected runtime: under 10 minutes for 5000 seeded networks, measured at exit.
This is a finite search, not a proof of the dichotomy.
"""
from pathlib import Path
import sys
import time
import random
from fractions import Fraction as F
from itertools import product
from collections import Counter

HERE = Path(__file__).resolve().parent
BASE = HERE.parents[2]
sys.dont_write_bytecode = True
sys.path.insert(0, str(BASE / 'round3' / 'dichotomy' / 'scripts'))
sys.path.insert(0, str(BASE.parent / 'scripts'))
import crn
import d2cert as D
from d2lib import Net


def random_networks(count):
    rng = random.Random(721087)
    for _ in range(count):
        m = rng.choice([2,3,4])
        pts = [(a,b) for a in range(m+1) for b in range(m-a+1)]
        coords = rng.sample(pts, rng.choice([3,4,5,6]))
        n = rng.choice([4,5,6])
        extra = rng.sample([(a,b) for a in range(-2,3) for b in range(-2,3)
                            if (a,b) not in [(0,0),(1,0),(0,1),(-1,-1)]], n-3)
        offsets = [-min(a*x+b*y for x,y in coords)+rng.choice([0,0,1]) for a,b in extra]
        cxs = [tuple([x,y,m-x-y]+[c+a*x+b*y for c,(a,b) in zip(offsets,extra)]) for x,y in coords]
        rx = [(cxs[i], cxs[(i+1)%len(cxs)]) for i in range(len(cxs))]
        if crn.rank(crn.Network(n,rx).Gamma) != 2:
            continue
        yield Net(n,rx)


def matrix(N, v, rates):
    J = [[F(0),F(0)],[F(0),F(0)]]
    for i, pp in v['loc'].items():
        for j,power in enumerate([(1,0),(0,1)]):
            c = pp.get(power,0)
            J[0][j] += rates[i]*c*N.G[i][v['e1']]
            J[1][j] += rates[i]*c*N.G[i][v['e2']]
    return J


def main():
    start = time.perf_counter()
    print('Environment: Python standard library; exact Fraction arithmetic.', flush=True)
    print('Runtime: measured at exit; expected under 10 minutes.', flush=True)
    count = int(sys.argv[1]) if len(sys.argv)>1 else 5000
    stats = Counter()
    for N in random_networks(count):
        stats['networks'] += 1
        for v in D.vertices(N):
            if not v['eq']:
                continue
            active = [i for i,P in v['loc'].items() if (1,0) in P or (0,1) in P]
            if len(active) < 3:
                continue
            stats['linear_vertices'] += 1
            if D.vertex_hurwitz_cert(N,v):
                continue
            stats['no_tropical_certificate'] += 1
            for smallrates in product((1,16), repeat=len(active)):
                rates = [1]*N.n
                for i,c in zip(active,smallrates):
                    rates[i] = c
                J = matrix(N,v,rates)
                tr = J[0][0]+J[1][1]
                det = J[0][0]*J[1][1]-J[0][1]*J[1][0]
                if tr < 0 and det > 0:
                    stats['candidate'] += 1
                    print('CANDIDATE', N.rx, 'zero', v['Zv'], 'rates', rates, 'J', J, flush=True)
                    break
        if stats['networks'] % 500 == 0:
            print('COUNTS', dict(stats), flush=True)
    ok=stats['candidate']==0
    print(('PASS' if ok else 'FAIL')+': LV no sampled Hurwitz vertex lacks the earlier tropical certificate.',dict(stats))
    print('Runtime seconds:', round(time.perf_counter()-start,3))
    return 0 if ok else 1


if __name__ == '__main__':
    raise SystemExit(main())
