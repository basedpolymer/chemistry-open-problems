"""Non-rigorous exploration of the parameter-free transformed shooting.

Python/scipy, one thread, seconds. No output here is a proof certificate.
"""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
import json
import time
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp

def rhs(x,y):
    q,v,z,c,p,j,k,l=y
    if q<=0: f=fp=fpp=fppp=0.
    else:
        f=np.exp(-1/q); fp=f/q**2
        fpp=f*(1-2*q)/q**4
        fppp=f*(1-6*q+6*q*q)/q**6
    return [p,j,k,l,-np.exp(x)*f,-np.exp(x)*fp*v,
            -np.exp(x)*(fp*z+fpp*v*v),
            -np.exp(x)*(fp*c+3*fpp*v*z+fppp*v**3)]

def shoot(a,which='v'):
    def event(x,y):return y[1] if which=='v' else y[0]-.02
    event.terminal=True;event.direction=-1
    def zero(x,y):return y[0]
    zero.terminal=True;zero.direction=-1
    sol=solve_ivp(rhs,(-35,12),[a,1,0,0,0,0,0,0],rtol=2e-12,atol=1e-13,
                  method='DOP853',events=[event,zero])
    if not sol.success:raise RuntimeError('FAIL numerical integration')
    return sol.t[-1],sol.y[:,-1]

if __name__=='__main__':
    start=time.monotonic();out=[]
    for a in [.5,.55,.58,.59,.5904,.6,.65,.7,.8,1,1.2,1.4,1.5,1.52,1.53,1.54,1.6,2]:
        x,y=shoot(a)
        row={'a':a,'x':float(x),'epsilon':float(y[0]),'v':float(y[1]),'z':float(y[2]),'c':float(y[3])}
        out.append(row);print('PASS numerical sample (not a proof):',row,flush=True)
    runtime=round(time.monotonic()-start,3)
    Path('certificates/reduced-exploration.json').write_text(json.dumps({'warning':'NON-RIGOROUS','samples':out,'runtime_seconds':runtime},indent=2)+'\n')
    print('Runtime seconds:',runtime,flush=True)
