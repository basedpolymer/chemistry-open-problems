#!/usr/bin/env python3
"""Exact, exhaustive certificates for the three remaining four-source structures.

Python 3.12+, NumPy, SymPy. Integer arrays are used only for bounded integer
operations (entries initially in [-2,2]); no floating-point rank or LP is used.
Run from any directory. Expected runtime: a few minutes on one CPU.
Every assertion prints PASS/FAIL; a failed assertion exits nonzero.
"""
import itertools as it
import json
import math
import time
from collections import Counter
from pathlib import Path
import numpy as np
import sympy as sp

BASE=Path(__file__).resolve().parents[1]
SOURCES={
 'S1': [(0,0,0,0,0),(0,0,0,0,1),(0,0,1,1,0),(1,1,0,0,0)],
 'S2': [(0,0,0,0,0),(0,0,0,0,2),(0,0,1,1,0),(1,1,0,0,0)],
 'S3': [(0,0,0,0,1),(0,0,0,0,2),(0,0,1,1,0),(1,1,0,0,0)],
 'S7': [(0,0,0,0,1),(0,0,0,1,1),(0,0,2,0,0),(1,1,0,0,0)],
 'S9': [(0,0,0,0,0,1),(0,0,0,0,1,1),(0,0,1,1,0,0),(1,1,0,0,0,0)]}
NORMAL_FORMS={
 'S3': [
 [[0,1,0,0],[0,0,1,1]],
 [[1,0,-2,2],[0,1,1,0]],
 [[1,0,-1,1],[0,1,'1/2','1/2']],
 [[1,0,-1,1],[0,1,1,0]],
 [[1,0,-1,1],[0,1,1,1]],
 [[1,0,-1,1],[0,1,2,0]],
 [[1,0,'-1/2',1],[0,1,'1/2',0]],
 [[1,0,0,0],[0,1,1,1]],
 [[1,0,0,1],[0,1,1,0]],
 [[1,0,'1/2','1/2'],[0,1,'1/2','1/2']],
 [[1,0,1,1],[0,1,1,1]]],
 'S7': [
 [[0,1,0,1],[0,0,1,-1]],
 [[1,0,-1,2],[0,1,1,0]],
 [[1,0,'-1/2','1/2'],[0,1,'1/2','1/2']],
 [[1,0,0,-1],[0,1,0,0]],
 [[1,1,0,1],[0,0,1,0]]],
 'S9': [
 [[0,1,0,0,1],[0,0,1,0,-1],[0,0,0,1,0]],
 [[1,0,'-1/2',0,'1/2'],[0,1,'1/2',0,'1/2'],[0,0,0,1,0]],
 [[1,0,0,0,0],[0,1,0,0,0],[0,0,1,1,1]]]
}

def check(ok,msg):
 print(('PASS ' if ok else 'FAIL ')+msg,flush=True)
 if not ok: raise AssertionError(msg)

def complexes(n):
 z=(0,)*n
 out=[z]
 for i in range(n):
  e=list(z);e[i]=1;out.append(tuple(e))
 for i in range(n):
  for j in range(i,n):
   e=list(z);e[i]+=1;e[j]+=1;out.append(tuple(e))
 return out

def det3(a):
 return (a[...,0,0]*(a[...,1,1]*a[...,2,2]-a[...,1,2]*a[...,2,1])
 -a[...,0,1]*(a[...,1,0]*a[...,2,2]-a[...,1,2]*a[...,2,0])
 +a[...,0,2]*(a[...,1,0]*a[...,2,1]-a[...,1,1]*a[...,2,0]))

def adj3(m):
 return np.stack([np.cross(m[1],m[2]),np.cross(m[2],m[0]),np.cross(m[0],m[1])],axis=1)

def spaces(g):
 """All rank-three spans, keyed by their primitive exterior coordinates."""
 unique=np.unique(g,axis=0)
 inds=np.array(list(it.combinations(range(len(unique)),3)),dtype=np.int64)
 triples=unique[inds]
 cols=list(it.combinations(range(g.shape[1]),3))
 wedges=np.stack([det3(triples[:,:,c]) for c in cols],axis=1)
 divisor=np.gcd.reduce(np.abs(wedges),axis=1)
 good=divisor!=0
 triples=triples[good];wedges=wedges[good]//divisor[good,None]
 first=np.argmax(wedges!=0,axis=1)
 wedges*=np.sign(wedges[np.arange(len(wedges)),first])[:,None]
 _,where=np.unique(wedges,axis=0,return_index=True)
 return triples[where]

def consistent(c):
 """The columns/rows c positively span R^3 iff its dual cone is zero.

 Any nonzero pointed dual cone has an extreme ray perpendicular to two
 independent rows of c. Exact cross products thus suffice to reject it.
 """
 c=np.unique(c,axis=0)
 ij=np.array(list(it.combinations(range(len(c)),2)),dtype=np.int64)
 rays=np.cross(c[ij[:,0]],c[ij[:,1]])
 rays=rays[np.any(rays!=0,axis=1)]
 dots=rays@c.T
 return not np.any(np.all(dots>=0,axis=1)|np.all(dots<=0,axis=1))

def kinetic_det(g,sourceids):
 G=sp.Matrix(g.T.tolist())
 K=sp.Matrix.hstack(*G.nullspace())
 rows=G.T.rref()[1]
 mats=[sp.zeros(3,K.cols) for _ in range(3)]
 for j,a in enumerate(sourceids):
  if a<3: mats[a]+=G[rows,j]*K[j,:]
 # A nonzero evaluation suffices to avoid the exclusion; it is never used
 # to assert that a zero polynomial is nonzero, or the converse.
 for seed in range(4):
  t=sp.Matrix([((j+1)*(seed+3))%11-5 for j in range(K.cols)])
  d=sp.Matrix.hstack(*(m*t for m in mats)).det()
  if d: return False
 t=sp.Matrix(sp.symbols('z:'+str(K.cols)))
 d=sp.Poly(sp.expand(sp.Matrix.hstack(*(m*t for m in mats)).det()),*t)
 return d.is_zero

def canon(name,R):
 if name=='S7': perms=[(0,1,2,3),(1,0,2,3)]
 else:
  perms=[(a,b,c,d) for a,b in [(0,1),(1,0)] for c,d in [(2,3),(3,2)]]
  perms+= [(p[2],p[3],p[0],p[1]) for p in perms.copy()]
  if name=='S9': perms=[p+(4,) for p in perms]
 return min((R[:,p].rref()[0] for p in perms),key=lambda m:tuple(m))

def run(name,sources):
 n=len(sources[0]);k=n-3
 reactions=[(a,p) for a,r in enumerate(sources) for p in complexes(n) if p!=r]
 g=np.array([[p[j]-sources[a][j] for j in range(n)] for a,p in reactions],dtype=np.int64)
 ids=np.array([a for a,p in reactions],dtype=np.int64)
 bases=spaces(g)
 check(len(bases)>0,f'{name}: enumerated {len(bases)} distinct rank-three reaction spans')
 Y=sp.Matrix(sources);A=Y[:3,:]-sp.ones(3,1)*Y[3,:]
 V=sp.Matrix.hstack(*A.nullspace())
 acr=next(i for i in range(n) if not any(V[i,:]))
 formname='S3' if name in ('S1','S2') else name
 normal=[sp.Matrix([[sp.Rational(v) for v in row] for row in R]) for R in NORMAL_FORMS[formname]]
 hist=Counter();records=[];observed=set()
 for b in bases:
  piv=next(c for c in it.combinations(range(n),3) if det3(b[:,c]))
  m=b[:,piv];det=det3(m);adj=adj3(m)
  mask=np.all((g[:,piv]@adj)@b==det*g,axis=1)
  if set(ids[mask])!={0,1,2,3}: hist['missing-source']+=1;continue
  c=g[mask][:,piv]
  if not consistent(c): hist['no-positive-flux']+=1;continue
  W=sp.Matrix.hstack(*sp.Matrix(b.tolist()).nullspace()).T.rref()[0]
  terms=[W[:,I].det()*V[I,:].det() for I in it.combinations(range(n),k)]
  signs={sp.sign(t) for t in terms if t}
  label=None;form=None
  if not signs: label='class-degenerate'
  elif len(signs)==1: label='injective'
  elif kinetic_det(g[mask],ids[mask]): label='kinetic-degenerate'
  else:
   R=canon(formname,W[:,[i for i in range(n) if i!=acr]])
   if R not in normal:
    print(name,'UNCOVERED',W.tolist(),R.tolist(),flush=True)
    check(False,'every surviving conservation space has a proved normal form')
   form=normal.index(R);observed.add(form);label='normal-form'
  hist[label]+=1
  records.append({'basis':b.tolist(),'W':[[str(v) for v in row] for row in W.tolist()],
                  'classification':label,'normal_form':form})
 check(sum(hist.values())==len(bases),f'{name}: every enumerated span classified exactly once')
 if name in ('S1','S2'):
  check(not observed,f'{name}: no noninjective nondegenerate cases remain')
 else:
  check(observed==set(range(len(normal))),f'{name}: exactly {len(normal)} remaining normal forms')
 check(all(abs(int(v))<1000000 for row in g for v in row),f'{name}: bounded exact integer arithmetic')
 print(name,dict(hist),flush=True)
 return {'sources':sources,'spans':len(bases),'counts':dict(hist),'records':records}

if __name__=='__main__':
 start=time.monotonic()
 results={name:run(name,sources) for name,sources in SOURCES.items()}
 out=BASE/'certificates';out.mkdir(exist_ok=True)
 (out/'subspaces.json').write_text(json.dumps(results,indent=2)+'\n')
 check(True,f'completed in {time.monotonic()-start:.1f}s; certificates/subspaces.json written')
