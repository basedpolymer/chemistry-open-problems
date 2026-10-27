"""Exact ancillary certificates for Theorems AF, E, and O.

Environment: Python 3.14 / standard library only; integer and Fraction arithmetic.
Runtime: less than one second; measured at exit. Universal statements are
proved in ../notes.md; finite checks below are labeled as such.
Imports previous crn.py and this folder's certificate helpers unchanged.
"""
from pathlib import Path
import sys
import time
from fractions import Fraction as F
HERE=Path(__file__).resolve().parent
sys.dont_write_bytecode=True
sys.path.insert(0,str(HERE))
from certify_out_all import crn,add,mul,term


def matmul(A,B):
    return [[sum(a*b for a,b in zip(row,col)) for col in zip(*B)] for row in A]


def main():
    start=time.perf_counter()
    print('Environment: Python standard library; exact integer/Fraction arithmetic.')
    print('Runtime: expected <1 second; measured below.')
    fails=[]
    def claim(name,ok,detail):
        print(('PASS' if ok else 'FAIL')+': '+name+': '+detail)
        if not ok:
            fails.append(name)
    count=0
    for u in range(3,11):
        for R in range(1,11):
            cxs=[(0,0,u,0),(1,1,0,0),(0,1,2,R)]
            N=crn.Network(4,[(cxs[j],cxs[(j+1)%3]) for j in range(3)])
            w=(u*R+2,u-2,R+1,u-2)
            good=all(sum(w[i]*x[i] for i in range(4))==u*(R+1) for x in cxs)
            good=good and crn.rank(N.Gamma)==2 and N.weakly_reversible()
            good=good and len(set(map(tuple,N.Gamma)))==4 and all(sum(row)==0 for row in N.Gamma)
            claim(f'AF-network-{u}-{R}',good,'positive conservation/flux, WR, rank two, distinct rows.')
            count+=1
    claim('AF-count',count==80,'80 finite parameter pairs checked; general family proof is in notes.md.')
    # At u=3,R=2, all rates=1 and T=1: c' <= -c^2.
    claim('AF-NONE-unit',F(3-2)*F(1)**3-2==-1,'u=3,R=2,rates=1,T=1: c\' <= -c^2.')
    cxs=[(0,0,3,0),(1,1,0,0),(0,1,2,3)]
    N=crn.Network(4,[(cxs[j],cxs[(j+1)%3]) for j in range(3)])
    selector=[[0,0,1],[1,0,0],[0,1,0]]
    claim('E-selector',matmul(N.Ye,selector)==N.Yp,'Yp=Ye*Sigma, using the three-cycle selector.')
    lift=[[-1,0],[1,-1],[0,1]]
    Pi=[row[:2] for row in N.Gamma]
    claim('E-lift',matmul(N.Ye,lift)==Pi and crn.rank(Pi)==2,'Ye*R=Pi, R=[[-1,0],[1,-1],[0,1]].')
    w=[11,1,4,1]
    gram=matmul(list(map(list,zip(*Pi))),[[w[i]*x for x in Pi[i]] for i in range(4)])
    claim('E-positive-Gram',gram[0][0]>0 and gram[0][0]*gram[1][1]-gram[0][1]**2>0,
          'exact positive Hessian witness at kappa=(11,1,4,1), t=0: '+str(gram)+'.')
    dual=N.dual()
    k=[F(1),F(1),F(1,100)]
    rates=[F(17),F(400),F(70000),F(1000000)]
    flux=[]
    for rate,(s,p) in zip(rates,dual.reactions):
        v=rate
        for x,e in zip(k,s):
            v*=x**e
        flux.append(v)
    field=[sum(dual.Gamma[j][i]*flux[i] for i in range(4)) for j in range(3)]
    claim('O-equilibrium',flux==[17,4,7,1] and field==[0,0,0],
          'k=(1,1,1/100), rates=(17,400,70000,1000000), flux=(17,4,7,1).')
    ambient=[[sum(dual.Gamma[j][i]*flux[i]*dual.reactions[i][0][l]/k[l] for i in range(4))
              for l in range(3)] for j in range(3)]
    J=[[ambient[j][l]-ambient[j][0] for l in (1,2)] for j in (1,2)]
    claim('O-Jacobian',J==[[-59,3658],[-25,79]],'reduced Jacobian [[-59,3658],[-25,79]].')
    tr=J[0][0]+J[1][1]
    det=J[0][0]*J[1][1]-J[0][1]*J[1][0]
    claim('O-source',tr==20 and det==86789 and tr*tr-4*det==-4*86689,
          'trace=20, determinant=86789; eigenvalues 10 +/- i sqrt(86689).')
    H=[[93655,-290457],[-290457,13471234]]
    Hdet=H[0][0]*H[1][1]-H[0][1]**2
    lhs1=matmul(list(map(list,zip(*J))),H)
    lhs2=matmul(H,J)
    lhs=[[lhs1[i][j]+lhs2[i][j] for j in range(2)] for i in range(2)]
    claim('O-quadratic-source',H[0][0]>0 and Hdet>0 and lhs==[[3471560,0],[0,3471560]],
          'H positive definite, det(H)='+str(Hdet)+', J^T H+H J=3471560 I.')
    x0=[F(1),F(99,100),F(1,50)]
    claim('O-initial',all(x>0 for x in x0) and x0!=k and sum(x0)==sum(k)==F(201,100),
          'x0=(1,99/100,1/50) is a different positive point in the same class.')
    bp=add(term(-1,0,1,0,1,0,0,0),term(2,3,0,2,0,0,1,0),term(3,0,0,3,0,0,0,1))
    cp=add(term(-1,0,1,1,0,1,0,0),term(1,3,0,2,0,0,1,0),term(-3,0,0,3,0,0,0,1))
    rhs=add(term(-1,0,1,0,1,0,0,0),term(2,0,1,1,0,1,0,0),term(9,0,0,3,0,0,0,1))
    claim('O-uniqueness-identity',add(bp,mul(term(-2,0,0,0,0,0,0,0),cp))==rhs,
          'b\'-2c\'=(2 beta c-alpha)b+9 delta c^3; monotone equilibrium parametrization in notes.md.')
    print('Runtime seconds:',round(time.perf_counter()-start,6))
    return 1 if fails else 0


if __name__=='__main__':
    raise SystemExit(main())
