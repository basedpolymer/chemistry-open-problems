"""A2 general rank-two exact search. Python standard library, Fraction only.

Runtime: measured at exit; about 8 minutes per 1000 input networks, with a
15-minute time cap. A time-capped run reports its finite prefix explicitly.
Imports all earlier geometry unchanged. Reports candidates, not a theorem.
"""
from pathlib import Path
import sys
import time
import random
import json
from fractions import Fraction as F
from collections import Counter, defaultdict

HERE = Path(__file__).resolve().parent
sys.dont_write_bytecode = True
sys.path.insert(0,str(HERE))
from search_linear_vertices import random_networks
from search_rates import fixed_repeller
import d2cert as D


def structural_check(N,v,cert):
    kind,lam,mu = cert
    terms = defaultdict(F)
    for i,pp in v['loc'].items():
        g1,g2 = N.G[i][v['e1']],N.G[i][v['e2']]
        for (p,q),c in pp.items():
            if kind == 'linear':
                terms[i,p,q] += c*(lam*g1+mu*g2)
            else:
                terms[i,p-1,q] += c*lam*g1
                terms[i,p,q-1] += c*mu*g2
    return D.dominated([((p,q),1 if c>0 else -1) for (i,p,q),c in terms.items() if c])


def main():
    start = time.perf_counter()
    print('Environment: Python standard library; exact rational arithmetic.', flush=True)
    print('Runtime: measured at exit; about 8 minutes per 1000 inputs; cap 15 minutes.', flush=True)
    count = int(sys.argv[1]) if len(sys.argv)>1 else 1000
    stats = Counter()
    rng = random.Random(9203)
    hits = []
    capped=False
    for N in random_networks(count):
        if time.perf_counter()-start>900:
            capped=True
            break
        stats['networks'] += 1
        facets = D.facets(N)
        if any(f['siphon'] and f['type'] != 'IN' and f['fl'] == 'zero' for f in facets):
            continue
        verts = [v for v in D.vertices(N) if v['eq']]
        bad = [v for v in verts if not D.vertex_repelling(N,v)]
        if not bad:
            continue
        stats['eligible'] += 1
        rates_list = [(1,)*N.n] + [tuple(rng.choice([1,16]) for _ in range(N.n)) for _ in range(20)]
        for rates in rates_list:
            certs = [fixed_repeller(N,v,rates) for v in bad]
            if not all(certs):
                continue
            if all(structural_check(N,v,c) for v,c in zip(bad,certs)):
                stats['structural_success'] += 1
                break
            stats['rate_candidate'] += 1
            rec = {'rx':N.rx,'rates':rates,'corners':[(v['Zv'],c) for v,c in zip(bad,certs)],
                   'facets':[(f['Z'],f['type'],f['fl'],f['mu']) for f in facets if f['siphon']]}
            hits.append(rec)
            print('RATE-CANDIDATE',rec,flush=True)
            break
        if stats['networks'] % 1000 == 0:
            print('COUNTS',dict(stats),flush=True)
    (HERE.parent / 'logs').mkdir(exist_ok=True)
    (HERE.parent / 'logs' / 'general_rate_candidates.json').write_text(json.dumps(hits,default=str,indent=2)+'\n')
    print('PASS: GR finite exact prefix recorded; time cap reached='+str(capped)+'.',dict(stats))
    print('Runtime seconds:',round(time.perf_counter()-start,3))


if __name__ == '__main__':
    main()
