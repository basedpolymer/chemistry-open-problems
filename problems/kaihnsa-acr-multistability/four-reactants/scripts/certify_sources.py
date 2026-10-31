#!/usr/bin/env python3
"""Enumerate all four distinct at-most-bimolecular sources, exactly.
Python 3.12+, NumPy; expected runtime < 1 minute. No floating-point arithmetic.
Species absent from every source are removed before enumeration.
"""
import itertools as it
import json
import time
from pathlib import Path
import numpy as np
from certify_subspaces import complexes, det3, check, SOURCES

BASE=Path(__file__).resolve().parents[1]
EXPECTED={
 'S1':[(0,0,0,0,0),(0,0,0,0,1),(0,0,1,1,0),(1,1,0,0,0)],
 'S2':[(0,0,0,0,0),(0,0,0,0,2),(0,0,1,1,0),(1,1,0,0,0)],
 'S4':[(0,0,0,0,1),(0,0,0,1,0),(0,0,1,0,1),(1,1,0,0,0)],
 'S5':[(0,0,0,0,1),(0,0,0,1,1),(0,0,1,0,1),(1,1,0,0,0)],
 'S6':[(0,0,0,0,1),(0,0,0,1,1),(0,0,1,1,0),(1,1,0,0,0)],
 'S8':[(0,0,0,0,1),(0,0,0,1,1),(0,1,1,0,0),(1,0,1,0,0)],
 **SOURCES}

def canonical(y):
 # Source permutations (24) and species permutations (sort columns).
 return min(tuple(sorted(tuple(y[r][j] for r in p) for j in range(len(y[0]))))
            for p in it.permutations(range(4)))

start=time.monotonic();result={}
for n in range(4,9):
 seen={};checked=0
 for y in it.combinations(complexes(n),4):
  if any(all(row[j]==0 for row in y) for j in range(n)): continue
  a=np.array([[y[r][j]-y[3][j] for j in range(n)] for r in range(3)],dtype=np.int64)
  nonzero=[p for p in it.combinations(range(n),3) if det3(a[:,p])]
  if not nonzero: continue
  acr=set.intersection(*(set(p) for p in nonzero))
  if not acr: continue
  checked+=1
  if n==4:
   v=[(-1)**j*int(det3(a[:,[i for i in range(n) if i!=j]])) for j in range(n)]
   if len(set(v))>3: check(False,'n=4 has at most three distinct toric exponents, including zero')
  seen[canonical(y)]=y
 expected={canonical(y) for y in EXPECTED.values() if len(y[0])==n}
 if n==4: check(len(seen)==26,'n=4: all 26 source types have at most three distinct exponents')
 else: check(set(seen)==expected,f'n={n}: exactly {len(expected)} specified source types')
 result[str(n)]={'labelled_sets':checked,'types':len(seen),'representatives':list(seen.values())}
out=BASE/'certificates';out.mkdir(exist_ok=True)
(out/'sources.json').write_text(json.dumps(result,indent=2)+'\n')
check(True,f'exhaustive source classification completed in {time.monotonic()-start:.1f}s')
