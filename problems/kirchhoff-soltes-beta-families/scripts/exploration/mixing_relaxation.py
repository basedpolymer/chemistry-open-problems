"""mixing_relaxation.py -- FLOATING-POINT EXPLORATION: relaxation of Section 5 (classes grouped by rho only, the
condition on c = delta - a dropped): best single-class efficiency and best balanced pair (a1 > 0 > a2, equal rho).
Usage: python3 mixing_relaxation.py classes.txt  (output of class_pairs.py dump)."""
# relaxation: classes grouped by rho only (ignore c); single eff and balanced pairs (a1>0>a2, same rho)
import sys
from collections import defaultdict
byrho=defaultdict(list)
for line in open(sys.argv[1]):
    s,x,g,m,r,d,a=line.split(); g=int(g); m=int(m); r=float(r); a=float(a)
    R=a+m*r
    byrho[round(r,7)].append((g,m,a,R,s,int(x)))
best=[]
for rho,L in byrho.items():
    # best single
    for (g,m,a,R,s,x) in L:
        e = g/m if a<=1e-12 else g*rho/R
        best.append((e,'single',rho,(g,m,a,s,x)))
    pos=[t for t in L if t[2]>1e-9]; neg=[t for t in L if t[2]<-1e-9]
    # for pairs, only need Pareto-best: maximize over pos (g,m,a) and neg
    for p in pos:
        for q in neg:
            a1,a2=p[2],-q[2]
            fr=(a2*p[0]+a1*q[0])/(a2*p[1]+a1*q[1])
            best.append((fr,'pair',rho,(p[0],p[1],p[2],p[4],p[5]),(q[0],q[1],q[2],q[4],q[5])))
best.sort(key=lambda t:-t[0])
for t in best[:12]: print('%.4f'%t[0], t[1:])
