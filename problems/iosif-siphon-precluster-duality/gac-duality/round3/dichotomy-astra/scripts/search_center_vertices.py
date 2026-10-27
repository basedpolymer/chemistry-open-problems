"""A2 exact center-manifold coefficient search; not a universal certificate.

Environment: Python 3.14 / standard library, exact Fraction arithmetic.
Expected runtime: under 5 minutes for 5000 seeded networks; measured at exit.
Earlier geometry is imported unchanged. Reported negative coefficients need
the hypotheses and persistence counterpart discussed in notes.md.
"""
from pathlib import Path
import sys
import time
from itertools import product
from collections import Counter
from fractions import Fraction as F
HERE=Path(__file__).resolve().parent
sys.dont_write_bytecode=True
sys.path.insert(0,str(HERE))
from search_linear_vertices import random_networks
import d2cert as D


def main():
    start=time.perf_counter()
    print('Environment: standard library; exact rationals. Runtime measured at exit.',flush=True)
    stats=Counter()
    count=int(sys.argv[1]) if len(sys.argv)>1 else 5000
    for N in random_networks(count):
        stats['networks']+=1
        for v in D.vertices(N):
            if not v['eq']:
                continue
            for side in (0,1):
                ea,eb=(v['e1'],v['e2']) if not side else (v['e2'],v['e1'])
                loc={i:{((p,q) if not side else (q,p)):c for (p,q),c in P.items()} for i,P in v['loc'].items()}
                if any((1,0) in P for P in loc.values()):
                    continue
                lin=[i for i,P in loc.items() if (0,1) in P]
                if not lin or not any(N.G[i][eb]<0 for i in lin):
                    continue
                axes=sorted({p for i,P in loc.items() for p,q in P if q==0 and (N.G[i][ea] or N.G[i][eb])})
                if not axes:
                    continue
                m=axes[0]
                if m<2:
                    continue
                active=sorted(set(lin)|{i for i,P in loc.items() if (m,0) in P})
                stats['center_configurations']+=1
                got=False
                for kk in product((1,16),repeat=len(active)):
                    k=dict(zip(active,kk))
                    delta=sum(k[i]*loc[i][0,1]*N.G[i][eb] for i in lin)
                    if delta>=0:
                        continue
                    beta=sum(k[i]*loc[i][0,1]*N.G[i][ea] for i in lin)
                    fm=sum(k[i]*loc[i].get((m,0),0)*N.G[i][ea] for i in active)
                    gm=sum(k[i]*loc[i].get((m,0),0)*N.G[i][eb] for i in active)
                    cc=fm-beta*gm/delta
                    if cc<0:
                        got=True
                        if D.vertex_cm_cert(N,v) is None and D.vertex_hurwitz_cert(N,v) is None:
                            stats['candidate']+=1
                            print('CANDIDATE',N.rx,'Z',v['Zv'],'side',side,'m',m,'rates',k,'delta',delta,'coefficient',cc,flush=True)
                        break
                if got:
                    stats['sampled_negative']+=1
        if stats['networks']%1000==0:
            print('COUNTS',dict(stats),flush=True)
    ok=stats['candidate']==0
    print(('PASS' if ok else 'FAIL')+': CM no sampled negative center coefficient lacks an earlier CM/H certificate.',dict(stats))
    print('Runtime seconds:',round(time.perf_counter()-start,3))
    return 0 if ok else 1


if __name__=='__main__':
    raise SystemExit(main())
