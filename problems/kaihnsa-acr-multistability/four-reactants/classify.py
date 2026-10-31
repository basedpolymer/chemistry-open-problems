"""Exploratory classification of the remaining conservation spaces."""
import json,itertools
from pathlib import Path
import sympy as s
BASE=Path(__file__).resolve().parent
perms={
 'S3': [(a,b,c,d) for a,b in [(0,1),(1,0)] for c,d in [(2,3),(3,2)]]+[(c,d,a,b) for a,b in [(0,1),(1,0)] for c,d in [(2,3),(3,2)]],
 'S7':[(0,1,2,3),(1,0,2,3)],
 'S9':[(a,b,c,d,4) for a,b in [(0,1),(1,0)] for c,d in [(2,3),(3,2)]]+[(c,d,a,b,4) for a,b in [(0,1),(1,0)] for c,d in [(2,3),(3,2)]]}
for name in perms:
 classes={}
 for entry in json.loads((BASE/f'hard_{name}.json').read_text()):
  W=s.Matrix([[s.Rational(v) for v in row] for row in entry['reduced']])
  vals=[W[:,p].rref()[0] for p in perms[name]]
  R=min(vals,key=lambda m:tuple(m))
  classes[tuple(R)]=R
 print(name,len(classes))
 for i,R in enumerate(sorted(classes.values(),key=lambda m:tuple(m))): print(i,R.tolist())
