"""A2 candidate from three linked-direction three-cycles.

Python standard library, exact arithmetic. Runtime: printed at exit.
This script explores a candidate; it does not certify persistence by itself.
"""
from pathlib import Path
import sys
import time
HERE=Path(__file__).resolve().parent
BASE=HERE.parents[2]
sys.dont_write_bytecode=True
sys.path.insert(0,str(BASE.parent/'scripts'))
sys.path.insert(0,str(BASE/'round3'/'dichotomy'/'scripts'))
import crn
import d2cert as D
from d2lib import Net

G=[(-2,-1,0,0,2,1),(2,1,-2,-1,0,0),(0,0,2,1,-2,-1)]
S=[(2,1,0,0,1,2),(1,2,2,1,0,0),(0,0,1,2,2,1)]


def candidate():
    rx=[]
    for j,order in enumerate([(0,2,1),(1,0,2),(2,1,0)]):
        z=S[j]
        for d in order:
            zz=tuple(a+b for a,b in zip(z,G[d]))
            rx.append((z,zz))
            z=zz
        assert z==S[j]
    return crn.Network(6,rx)


def main():
    start=time.perf_counter()
    print('Environment: Python standard library, exact rationals. Runtime measured at exit.',flush=True)
    N=candidate()
    ok=N.weakly_reversible() and crn.rank(N.Gamma)==2 and len(set(map(tuple,N.Gamma)))==6
    ok=ok and all(sum(z)==6 for a,b in N.reactions for z in (a,b))
    print(('PASS' if ok else 'FAIL')+': C0 WR, H1 ones, H2 cycle fluxes, H3, rank two.',flush=True)
    print('RX',N.reactions,flush=True)
    print('DUAL',N.dual().reactions,flush=True)
    result=D.analyse(N.n,N.reactions)
    print('CLASSIFICATION',result['verdict'],'first NONE certificate',result['none'][0] if result['none'] else None,flush=True)
    geometry=Net(N.n,N.reactions)
    p=[1]*9
    p[0]=p[2]=0
    z=[0,-1,0,-2,0,0]
    direction=[sum(N.Gamma[i][j]*z[i] for i in range(6)) for j in range(9)]
    from fractions import Fraction as F
    positive=[F(x)+F(d,4) for x,d in zip(p,direction)]
    claim=all(x>0 for x in positive) and direction[0]==direction[2]==1
    print(('PASS' if claim else 'FAIL')+': C1 vertex p=(0,1,0,1,1,1,1,1,1) lies in a positive class.',flush=True)
    ok=ok and claim
    v=next(v for v in D.vertices(geometry) if v['Zv']==(0,2))
    # Coordinates are explicitly (k_1,k_3), irrespective of the library ordering.
    J=[[[F(0)]*6 for _ in range(2)] for _ in range(2)]
    for i in range(6):
        for column in range(2):
            if (N.Ye[i][0],N.Ye[i][2])==((1,0) if column==0 else (0,1)):
                for row,j in enumerate((0,2)):
                    J[row][column][i]=F(N.Gamma[i][j])
    expected=[[[0,-1,0,0,0,0],[0,0,0,0,0,0]],
              [[0,1,0,0,0,0],[0,0,0,-1,0,0]]]
    claim=J==expected and v['eq']
    print(('PASS' if claim else 'FAIL')+': C2 J=[[-kappa_2,0],[kappa_2,-kappa_4]]; Hurwitz for every positive rate vector.',flush=True)
    ok=ok and claim
    print('Candidate SOME status: REFUTED; this network is NONE, by the hyperbolic boundary sink.',flush=True)
    print('Runtime seconds:',round(time.perf_counter()-start,3))
    return 0 if ok else 1


if __name__=='__main__':
    raise SystemExit(main())
