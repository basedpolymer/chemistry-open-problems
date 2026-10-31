#!/usr/bin/env python3
"""Five-source upper bound, derived exactly from the reaction list.

Python 3.12+, SymPy 1.14; runtime below 2 seconds. All arithmetic is exact.
No numerical search output is read. Mass action uses k*x^y.
"""
import json
from pathlib import Path
import sympy as s
from certify_subspaces import check

a,b,c,r=s.symbols('a b c r', positive=True)
x=s.Matrix([a,b,c,r])
# Source, product, rate; coordinates A,B,C,R.
reactions=[
 ((0,0,0,0),(0,0,1,0),64),
 ((1,1,0,0),(0,2,0,0),3),
 ((1,0,1,0),(0,1,1,0),5),
 ((0,1,1,0),(1,0,0,0),64),
 ((0,0,0,0),(0,0,0,1),1),
 ((0,0,0,1),(0,0,0,0),1),
]
G=s.Matrix.hstack(*(s.Matrix(v)-s.Matrix(u) for u,v,k in reactions))
rates=s.Matrix([k*s.prod(xi**ui for xi,ui in zip(x,u)) for u,v,k in reactions])
f=s.expand(G*rates)
check(len(set(u for u,v,k in reactions))==5,'exactly five distinct reactant complexes')
check(len(set(u for u,v,k in reactions)|set(v for u,v,k in reactions))==8,'exactly eight complexes')
check(all(sum(u)<=2 and sum(v)<=2 and k>0 for u,v,k in reactions),'all sources and products at most bimolecular; all rates positive')
check(G.rank()==3,'stoichiometric dimension is three')
check(G.T.nullspace()==[s.Matrix([1,1,0,0])],'the only independent conservation law is a+b')
check(f==s.Matrix([-3*a*b-5*a*c+64*b*c,3*a*b+5*a*c-64*b*c,64-64*b*c,1-r]),'mass-action equations derived from the reactions')
check(f[3]==1-r,'global ACR: every positive equilibrium has r=1')

reduced=s.Matrix([f[0],f[2],f[3]]).subs(b,9-a)
J=reduced.jacobian([a,c,r])
scalar=s.factor(reduced[0].subs({c:1/(9-a),r:1}))
check(s.cancel(scalar+3*(a-4)*(a-6)*(a-8)/(9-a))==0,'all equilibria in a+b=9: cubic factors exactly')
states=[s.Matrix([4,5,s.Rational(1,5),1]),s.Matrix([6,3,s.Rational(1,3),1]),s.Matrix([8,1,1,1])]
expected=[(s.Rational(-1684,5),s.Integer(1536)),(-206,-768),(-112,1536)]
records=[]
lam=s.symbols('lambda')
for i,(state,(tr,de)) in enumerate(zip(states,expected)):
 check(f.subs(dict(zip(x,state)))==s.zeros(4,1),f'equilibrium {i+1}: exact residual zero')
 check(state[0]+state[1]==9 and all(z>0 for z in state),f'equilibrium {i+1}: positive and in the same class')
 jac=J.subs(dict(zip(x,state)))
 core=jac[:2,:2]
 check(core.trace()==tr and core.det()==de,f'equilibrium {i+1}: exact planar trace {tr} and determinant {de}')
 poly=s.factor(jac.charpoly(lam).as_expr())
 check(s.expand(poly-(lam+1)*(lam**2-tr*lam+de))==0,f'equilibrium {i+1}: third eigenvalue is -1')
 stable=tr<0 and de>0
 check(bool(stable)==(i!=1),f'equilibrium {i+1}: '+('exponentially stable' if stable else 'saddle'))
 records.append({'state':[str(v) for v in state],'restricted_jacobian':[[str(v) for v in row] for row in jac.tolist()],
                 'characteristic_polynomial':str(poly),'stable':bool(stable)})

# An ACR species can also be coupled to the module without adding a source:
# replace 0->R and R->0 by 0->C+R and R->C, and lower the original C input to 62.
# The equilibrium C input is still 64 and the restricted Jacobian is triangular
# in blocks, so the same certificate applies. This avoids relying on a completely
# dynamically independent regulator, although independence is allowed by ACR.
coupled=list(reactions)
coupled[0]=((0,0,0,0),(0,0,1,0),62)
coupled[4]=((0,0,0,0),(0,0,1,1),1)
coupled[5]=((0,0,0,1),(0,0,1,0),1)
Gc=s.Matrix.hstack(*(s.Matrix(v)-s.Matrix(u) for u,v,k in coupled))
fc=s.expand(Gc*s.Matrix([k*s.prod(xi**ui for xi,ui in zip(x,u)) for u,v,k in coupled]))
check(fc==f+s.Matrix([0,0,r-1,0]),'optional coupled variant: c-dot gains r-1 only')
check(len(set(u for u,v,k in coupled))==5 and Gc.rank()==3,'coupled variant retains five sources and dimension three')
Jc=s.Matrix([fc[0],fc[2],fc[3]]).subs(b,9-a).jacobian([a,c,r])
check(Jc.charpoly(lam).as_expr()==J.charpoly(lam).as_expr(),'coupling preserves every restricted characteristic polynomial')

out=Path(__file__).resolve().parents[1]/'certificates'
out.mkdir(exist_ok=True)
(out/'five-example.json').write_text(json.dumps({'reactions':reactions,'equilibria':records,'coupled_variant':coupled},indent=2)+'\n')
check(True,'five-source upper bound certified in exact arithmetic')
