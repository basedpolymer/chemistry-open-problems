"""NON-RIGOROUS shooting exploration. scipy/numpy, one thread; about 10 s.

t=s^2/4, lambda=4*t_zero. Integrate w and its first two alpha derivatives.
No floating point output from this script is a proof certificate.
"""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
import json
from pathlib import Path
import time
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import root, brentq

def initial(a,e,t):
    d=1+e*a; f=np.exp(a/d); fp=f/d**2; fpp=f*(1-2*e*d)/d**4
    return [a-f*t,-f,1-fp*t,-fp,-fpp*t,-fpp]

def rhs(t,y,e):
    w,p,v,q,z,j=y
    d=1+e*w; f=np.exp(w/d); fp=f/d**2; fpp=f*(1-2*e*d)/d**4
    return [p,(-p-f)/t,q,(-q-fp*v)/t,j,(-j-fp*z-fpp*v*v)/t]

def shoot(a,e,tend=None,tol=2e-12):
    t0=1e-9/np.exp(a/(1+e*a))
    def event(t,y): return y[0]
    event.terminal=True; event.direction=-1
    sol=solve_ivp(lambda t,y:rhs(t,y,e),(t0,tend or a+1),initial(a,e,t0),method='DOP853',rtol=tol,atol=tol*.01,events=None if tend else event)
    if not sol.success: raise RuntimeError(sol.message)
    return sol.t[-1],sol.y[:,-1]

def main():
    start=time.monotonic()
    def F(x):
        a,e=x
        _,y=shoot(a,e)
        return y[[2,4]]
    result=root(F,[6,0.24],tol=1e-10)
    assert result.success and np.linalg.norm(F(result.x))<1e-8
    a,e=result.x
    t,_=shoot(a,e)
    output={'warning':'NON-RIGOROUS: floating point only','cusp':{'alpha':a,'epsilon':e,'lambda':4*t,'residual':list(F(result.x))},'folds':[]}
    print('PASS numerical cusp residual (not a proof):',output['cusp'])
    for eps in [0.05,0.1,0.15,0.2,0.22,0.23,0.235,0.238,0.239,0.24,0.241,0.242,0.245,0.249]:
        aa=np.geomspace(.3,max(100,3/eps**2),140)
        vv=[shoot(x,eps)[1][2] for x in aa]
        folds=[]
        for lo,hi,vl,vh in zip(aa[:-1],aa[1:],vv[:-1],vv[1:]):
            if vl*vh<0:
                af=brentq(lambda x:shoot(x,eps)[1][2],lo,hi,xtol=1e-10)
                tf,yf=shoot(af,eps)
                folds.append({'alpha':af,'lambda':4*tf,'v':yf[2],'z':yf[4]})
        output['folds'].append({'epsilon':eps,'folds':folds})
        print('PASS sampled fold search (not exhaustive):',eps,folds)
    Path('certificates/exploration.json').write_text(json.dumps(output,indent=2)+'\n')
    print('Runtime seconds:',round(time.monotonic()-start,2))

if __name__=='__main__':main()
