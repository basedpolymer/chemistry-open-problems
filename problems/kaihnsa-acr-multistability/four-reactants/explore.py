"""Exploration only. No theorem is certified by this script."""
import json
import itertools
import random
from collections import Counter
from pathlib import Path
import sympy as s

BASE = Path(__file__).resolve().parent
OLD = BASE.parent / 'examples/scripts/explore'
STRUCT = {
 'S3': [(0,0,0,0,1),(0,0,0,0,2),(0,0,1,1,0),(1,1,0,0,0)],
 'S7': [(0,0,0,0,1),(0,0,0,1,1),(0,0,2,0,0),(1,1,0,0,0)],
 'S9': [(0,0,0,0,0,1),(0,0,0,0,1,1),(0,0,1,1,0,0),(1,1,0,0,0,0)]
}

for name, sources in STRUCT.items():
 n=len(sources[0]); k=n-3
 Y=s.Matrix(sources)
 A=Y[:3,:]-s.ones(3,1)*Y[3,:]
 V=s.Matrix.hstack(*A.nullspace())
 hist=Counter(); hard=[]
 for entry in json.loads((OLD/f'valid_{name}.json').read_text()):
  W=s.Matrix([[s.Rational(v) for v in row] for row in entry['W']])
  terms=[W[:,inds].det()*V[inds,:].det() for inds in itertools.combinations(range(n),k)]
  signs={s.sign(t) for t in terms if t}
  if not signs: hist['class-degenerate']+=1; continue
  if len(signs)==1: hist['injective']+=1; continue
  rx=entry['rx']
  G=s.Matrix.hstack(*(s.Matrix(p)-s.Matrix(r) for r,p in rx))
  K=s.Matrix.hstack(*G.nullspace())
  inds=G.T.rref()[1]
  sourceids=[sources.index(tuple(r)) for r,p in rx]
  M=[]
  for a in range(3):
   M.append(s.zeros(3,K.cols))
   for j,rid in enumerate(sourceids):
    if rid==a: M[-1]+=G[inds,j]*K[j,:]
  rng=random.Random(0)
  nonzero=False
  for _ in range(4):
   z=s.Matrix([rng.randrange(-4,5) for _ in range(K.cols)])
   if s.Matrix.hstack(*(m*z for m in M)).det(): nonzero=True;break
  if not nonzero:
   z=s.Matrix(s.symbols('t:'+str(K.cols)))
   det=s.expand(s.Matrix.hstack(*(m*z for m in M)).det())
   if det==0: hist['kinetic-degenerate']+=1;continue
  # Remove the ACR column and compute the affine plane in remaining species.
  acr=next(i for i in range(n) if not any(V[i,:]))
  ids=[i for i in range(n) if i!=acr]
  R=W[:,ids].rref()[0]
  key=tuple(tuple(row) for row in R.tolist())
  hist['mixed']+=1
  hard.append({'W':entry['W'],'reduced':[[str(v) for v in row] for row in R.tolist()]})
 (BASE/f'hard_{name}.json').write_text(json.dumps(hard,indent=2)+'\n')
 print(name,dict(hist),'distinct reduced',len({str(e['reduced']) for e in hard}),flush=True)
 print('Reduced mixed:',Counter(str(e['reduced']) for e in hard),flush=True)
