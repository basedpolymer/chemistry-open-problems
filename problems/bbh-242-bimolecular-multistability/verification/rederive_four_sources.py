#!/usr/bin/env python3
"""Re-derive the four-source classification without importing the solver's modules.

Exact circuit determinants give the positive kernel cone; cofactor minors give
the affine source relation. SymPy isolates polynomial roots over the rationals.
Python 3.12, SymPy 1.14; about a minute on one core. Exit nonzero on failure.
"""
from itertools import combinations, product
from math import gcd
from functools import reduce
import sympy as S

C = tuple((a,b) for a in range(3) for b in range(3-a))
def det(v,w): return v[0]*w[1]-v[1]*w[0]
def primitive(v):
    g=reduce(gcd,(abs(x) for x in v))
    return tuple(x//g for x in v)
def canonical(net):
    def key(n):
        return tuple(sorted((s,primitive((p[0]-s[0],p[1]-s[1]))) for s,p in n))
    return min(key(net),key([((s[1],s[0]),(p[1],p[0])) for s,p in net]))
def circuits(v):
    rays=set()
    for i,j in combinations(range(4),2):
        if det(v[i],v[j]) == 0 and sum(a*b for a,b in zip(v[i],v[j])) < 0:
            z=[0]*4
            z[i]=sum(abs(x) for x in v[j]); z[j]=sum(abs(x) for x in v[i])
            rays.add(primitive(z))
    for i,j,k in combinations(range(4),3):
        w=(det(v[j],v[k]),-det(v[i],v[k]),det(v[i],v[j]))
        if all(x>0 for x in w) or all(x<0 for x in w):
            z=[0]*4
            for ix,x in zip((i,j,k),w): z[ix]=abs(x)
            rays.add(primitive(z))
    return sorted(rays)
all_reactions=[(s,p) for s in C for p in C if s!=p]
total=0
for net in combinations(all_reactions,4):
    v=[(p[0]-s[0],p[1]-s[1]) for s,p in net]
    total+=any(det(x,y) for x,y in combinations(v,2))
assert total==27195,total
print('PASS: 27195 rank-two networks from all 30 choose 4 reaction sets',flush=True)
classes={}
for sources in combinations(C,4):
    for targets in product(*(tuple(p for p in C if p!=s) for s in sources)):
        net=tuple(zip(sources,targets)); v=[(p[0]-s[0],p[1]-s[1]) for s,p in net]
        if not any(det(x,y) for x,y in combinations(v,2)): continue
        rays=circuits(v)
        if not all(any(r[i]>0 for r in rays) for i in range(4)): continue
        assert len(rays)==2
        classes.setdefault(canonical(net),(net,rays))
assert len(classes)==838,len(classes)
print('PASS: 838 nontrivial four-source classes from positive circuits',flush=True)
t=S.symbols('t'); constant=fold=0
def roots(p):
    p=S.Poly(p,t)
    return int(p.count_roots(0,1))-int(p.eval(0)==0)-int(p.eval(1)==0)
for net,(u,v) in classes.values():
    source=S.Matrix([[*s,1] for s,p in net])
    c=[(-1)**i*S.Matrix([list(source.row(j)) for j in range(4) if j!=i]).det() for i in range(4)]
    assert any(c)
    h=[u[i]+t*(v[i]-u[i]) for i in range(4)]
    N=S.Poly(sum(c[i]*(v[i]-u[i])*S.prod(h[j] for j in range(4) if j!=i) for i in range(4)),t)
    J=S.zeros(2)
    for i,(s,p) in enumerate(net):
        J+=S.Matrix([p[0]-s[0],p[1]-s[1]])*S.Matrix([[s[0],s[1]]])*h[i]
    Q=S.Poly(J.det(),t)
    if N.is_zero:
        assert Q.is_zero
        constant+=1; continue
    assert not Q.is_zero
    assert all(roots(f)==0 for f,m in S.sqf_list(Q*N)[1] if m%2)
    nr=roots(Q)
    assert nr in (0,1)
    if nr:
        assert all(roots(f)==0 for f,m in S.sqf_list(Q)[1] if m>1)
        fold+=1
assert (constant,fold)==(9,30),(constant,fold)
print('PASS: 9 constant families, 30 simple-fold classes; Q*N never changes sign on (0,1)')
print('ALL PASS')
