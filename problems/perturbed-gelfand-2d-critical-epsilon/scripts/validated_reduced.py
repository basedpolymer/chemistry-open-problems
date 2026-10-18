"""Validated parameter-free shooting for phase 2.

Python 3, python-flint 0.9, scipy for proposals only. One thread per worker.
Arb 160 bits; Taylor order 14 with an independent interval remainder.
Full cover: about 145 seconds with six workers; individual boxes 5-10 s.
Commands and resource limits are documented in README.md.
"""
import argparse
import json
import math
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
os.environ['MKL_NUM_THREADS']='1'
from decimal import Decimal
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import time
from flint import arb, arb_series, ctx

ctx.prec=160
ORDER=14
ctx.cap=ORDER+3
ZERO=arb(0)

def ball(c,r):return arb(c,arb(r).abs_upper())
def hull(a,b):return arb(a).union(arb(b))
def bounds(x):return [x.lower().str(30),x.upper().str(30)]
def check(name,ok):
    print(('PASS ' if ok else 'FAIL ')+name,flush=True)
    if not ok:raise AssertionError(name)

def rhs(x,y):
    q,v,z,c,p,j,k,l=y
    f=(-1/q).exp();fp=f/q**2
    fpp=f*(1-2*q)/q**4
    fppp=f*(1-6*q+6*q*q)/q**6
    return [p,j,k,l,-x.exp()*f,-x.exp()*fp*v,
            -x.exp()*(fp*z+fpp*v*v),
            -x.exp()*(fp*c+3*fpp*v*z+fppp*v**3)]

def coeff(x,y,n):
    cs=[[a] for a in y]
    for k in range(n):
        ctx.cap=k+2
        rr=rhs(arb_series([x,arb(1)],prec=k+2),
               [arb_series(a,prec=k+2) for a in cs])
        for j in range(8):cs[j].append(rr[j][k]/(k+1))
    ctx.cap=ORDER+3
    return cs

def tube(x,y,h):
    X=x+hull(0,h);f0=rhs(x,y)
    floor=h*h*(1+sum((v.abs_upper() for v in f0),ZERO))/10000+sum((v.rad() for v in y),ZERO)/1000000+arb('1e-40')
    rad=[(h*v).abs_upper()*2+floor for v in f0]
    for _ in range(20):
        B=[v+ball(0,r) for v,r in zip(y,rad)]
        if not B[0]>0:raise ArithmeticError('positive-q tube unavailable')
        f=rhs(X,B)
        if all(b.contains_interior(v+hull(0,h)*g) for b,v,g in zip(B,y,f)):return X,B
        rad=[r.max((h*v).abs_upper()*2+floor) for r,v in zip(rad,f)]
    raise ArithmeticError('Picard inclusion failed')

def flat_forces(X,B):
    # C-infinity extension f(q)=0 for q<=0; q<=1/40 throughout the tube.
    if not B[0]<arb(1)/40:raise ArithmeticError('flat bound domain')
    v,z,c=B[1:4];t=X.exp();f=arb(-40).exp()
    fp=hull(0,f*40**2)
    fpp=ball(0,f*(40**4+2*40**3))
    fppp=ball(0,f*(40**6+6*40**5+6*40**4))
    return [-t*hull(0,f),-t*fp*v,-t*(fp*z+fpp*v*v),
            -t*(fp*c+3*fpp*v*z+fppp*v**3)]

def flat_step(x,y,h,window=False):
    H=hull(0,h);X=x+H
    # Fixed broad tube, followed by strict inclusion. No Taylor analyticity
    # across q=0 is assumed or used.
    B=[v+H*p+ball(0,'0.00001') for v,p in zip(y[:4],y[4:])]+[v+ball(0,'0.00001') for v in y[4:]]
    force=flat_forces(X,B)
    Z=[v+H*p for v,p in zip(y[:4],B[4:])]+[v+H*f for v,f in zip(y[4:],force)]
    if not all(b.contains_interior(z) for b,z in zip(B,Z)):raise ArithmeticError('flat Picard inclusion')
    H=H if window else h
    return [v+H*p+H*H*f/2 for v,p,f in zip(y[:4],y[4:],force)]+[v+H*f for v,f in zip(y[4:],force)]

def step(x,y,h,window=False):
    if y[0]<arb('.02') and y[4]<0:return flat_step(x,y,h,window)
    X,B=tube(x,y,h)
    if window:return B
    C=coeff(x,y,ORDER);R=coeff(X,B,ORDER+1)
    H=hull(0,h) if window else h
    out=[]
    for c,r in zip(C,R):
        v=c[-1]
        for a in c[-2::-1]:v=v*H+a
        remainder=r[-1]*H**(ORDER+1)
        if not remainder.abs_upper()<arb('1e-14'):raise ArithmeticError('remainder too wide')
        out.append(v+remainder)
    return out

def initial(a):
    if not (a>=arb('.49') and a<=arb('1.61')):raise ArithmeticError('initialization domain')
    t=arb(-50).exp()
    err=[t*math.factorial(k)*8**k for k in range(4)]
    return [c+ball(0,r) for c,r in zip([a,arb(1),ZERO,ZERO],err)]+[ball(0,r) for r in err]

def recenter(y,yc,da):
    # Integral mean-value form in the one parameter a. Sensitivities are part
    # of the validated state, so this intersection introduces no assumption.
    for j in [2,1,0,6,5,4]:
        y[j]=y[j].intersection(yc[j]+da*y[j+1])
    return y

def integrate(a,target,rx=ZERO):
    x=arb(-50);y=initial(a);yc=initial(a.mid());da=a-a.mid();n=0;flat=0
    target=target-rx
    while x<target:
        h=arb(1)/(64 if x>-4 else 2)
        if not x+h<target:h=target-x
        for _ in range(18):
            try:
                ycn=step(x,yc,h)
                yn=recenter(step(x,y,h),ycn,da)
                if all(v.is_finite() for v in yn):break
            except ArithmeticError:pass
            h=h/2
        else:raise ArithmeticError('step failed at '+str(x)+' state '+str([str(v) for v in y]))
        if y[0]<arb('.02'):flat+=1
        x=x+h;y=yn;yc=ycn;n+=1
    left=y
    if rx>0:
        end=x+2*rx;whole=list(y)
        while x<end:
            h=arb(1)/128
            if not x+h<end:h=end-x
            for _ in range(18):
                try:
                    ycn=step(x,yc,h)
                    yn=recenter(step(x,y,h),ycn,da)
                    W=recenter(step(x,y,h,True),step(x,yc,h,True),da)
                    break
                except ArithmeticError:h=h/2
            else:raise ArithmeticError('event-window step failed')
            whole=[u.union(v) for u,v in zip(whole,W)]
            x=x+h;y=yn;yc=ycn
        right=y
    else:whole=right=y
    return left,whole,right,n,flat

def proposal(a):
    # Floating output only chooses a rational window; all conclusions are Arb.
    from explore_reduced import rhs as nrhs
    from scipy.integrate import solve_ivp
    def event(x,y):return y[1]
    event.terminal=True;event.direction=-1
    sol=solve_ivp(nrhs,(-35,4),[float(a),1,0,0,0,0,0,0],rtol=1e-10,atol=1e-12,
                  method='DOP853',events=event)
    if not sol.success or not len(sol.t_events[0]):raise ArithmeticError('proposal failed')
    return arb(format(sol.t[-1],'.12f'))

def certify(lo,hi,kind):
    a=hull(arb(lo),arb(hi));xc=proposal(a.mid())
    ra=(arb(hi)-arb(lo))/2
    rx=arb('.00001')+2*ra
    L,Y,R,n,nflat=integrate(a,xc,rx)
    if not (L[1]>0 and R[1]<0 and Y[5]<0):raise ArithmeticError('event bracket '+str([str(L[1]),str(R[1]),str(Y[5])]))
    if kind=='negative' and not Y[2]<arb('-.01'):raise ArithmeticError('z not below -0.01')
    if kind=='positive' and not Y[2]>arb('.01'):raise ArithmeticError('z not above 0.01')
    derivative=Y[3]-Y[6]*Y[2]/Y[5]
    if kind=='cusp' and not derivative>1:raise ArithmeticError('event z derivative not above 1')
    if kind in ('negative','cusp') and not Y[0]>0:raise ArithmeticError('positive fold level')
    if not Y[4]<0:raise ArithmeticError('q_x not negative')
    return {'a':[str(lo),str(hi)],'kind':kind,'x_center':xc.str(25),'x_radius':rx.str(25),
            'v_left':bounds(L[1]),'v_right':bounds(R[1]),'q_event':bounds(Y[0]),
            'z_event':bounds(Y[2]),'v_x_event':bounds(Y[5]),'q_x_event':bounds(Y[4]),
            'dz_da_event':bounds(derivative),'steps':n,'flat_steps':nflat}

def cover_piece(job):
    lo,hi,kind,depth=job
    try:return [certify(lo,hi,kind)],0
    except (ArithmeticError,ValueError) as exc:
        if depth>=12:raise RuntimeError('FAIL cover '+str(job)+': '+str(exc)) from None
        mid=str((Decimal(lo)+Decimal(hi))/2)
        L,nl=cover_piece((lo,mid,kind,depth+1))
        R,nr=cover_piece((mid,hi,kind,depth+1))
        return L+R,1+nl+nr

def full_cover(workers):
    start=time.monotonic();path=Path('certificates/reduced-cover.json')
    path.unlink(missing_ok=True)
    jobs=[]
    for lo,hi,kind in [('0.5','0.58','negative'),('0.58','0.6','cusp'),('0.6','1.6','positive')]:
        x=Decimal(lo)
        while x<Decimal(hi):
            y=min(x+Decimal('.01'),Decimal(hi));jobs.append((str(x),str(y),kind,0));x=y
    records=[];splits=0
    with ProcessPoolExecutor(max_workers=workers) as pool:
        futures={pool.submit(cover_piece,j):j for j in jobs}
        for f in futures:
            rec,n=f.result();records.extend(rec);splits+=n
            print('PASS cover block',futures[f][:3],'leaves',len(rec),flush=True)
    records.sort(key=lambda r:Decimal(r['a'][0]))
    check('cover starts at 1/2',Decimal(records[0]['a'][0])==Decimal('.5'))
    check('cover ends at 8/5',Decimal(records[-1]['a'][1])==Decimal('1.6'))
    check('no gaps between cover intervals',all(Decimal(L['a'][1])==Decimal(R['a'][0]) for L,R in zip(records[:-1],records[1:])))
    for r in records:
        lo,hi=map(Decimal,r['a']);k=r['kind']
        check_domain=(k=='negative' and Decimal('.5')<=lo<hi<=Decimal('.58')) or (k=='cusp' and Decimal('.58')<=lo<hi<=Decimal('.6')) or (k=='positive' and Decimal('.6')<=lo<hi<=Decimal('1.6'))
        if not check_domain:raise AssertionError('cover kind/domain mismatch')
    check('all cover kinds have their stated domains',True)
    points=[]
    for a,kind in [('0.5','negative'),('1.5','positive'),('1.6','positive')]:
        rec=certify(a,a,kind);points.append(rec)
        qlo=arb(rec['q_event'][0]);qhi=arb(rec['q_event'][1])
        if a=='0.5':check('E(1/2) > 0.238',qlo>arb('.238'))
        if a=='1.5':check('E(3/2) > 0.0078',qlo>arb('.0078'))
        if a=='1.6':check('E(8/5) < 0',qhi<0)
    data=json.loads(Path('certificates/cusp-input.json').read_text())
    aa=ball(data['centers'][0],data['radii'][0]);ee=ball(data['centers'][1],data['radii'][1])
    ac=ee+ee*ee*aa
    check('certified local cusp overlaps reduced cusp strip',ac>arb('.58') and ac<arb('.6'))
    # This is an overlap calculation, not a substitute for rerunning the local certificate.
    out={'status':'PASS reduced global cover (combine with local cusp certificate and analytic lemmas)',
         'precision_bits':ctx.prec,'order':ORDER,'workers':workers,'leaf_count':len(records),
         'subdivisions':splits,'records':records,'points':points,'transformed_cusp_amplitude':bounds(ac),
         'runtime_seconds':round(time.monotonic()-start,3)}
    path.write_text(json.dumps(out,indent=2)+'\n')
    print('PASS full reduced cover; leaves:',len(records),'runtime seconds:',out['runtime_seconds'],flush=True)

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--test',nargs=3,metavar=('LO','HI','KIND'))
    ap.add_argument('--full',action='store_true');ap.add_argument('--workers',type=int,default=6)
    args=ap.parse_args();start=time.monotonic()
    if args.test:
        result=certify(*args.test);check('test interval '+str(args.test),True);print(json.dumps(result,indent=2))
    if args.full:
        check('resource limit: 1 to 8 single-thread workers',1<=args.workers<=8)
        full_cover(args.workers)
    print('Runtime seconds:',round(time.monotonic()-start,3),flush=True)

if __name__=='__main__':
    try:main()
    except Exception as exc:
        print('FAIL',type(exc).__name__,str(exc),flush=True)
        raise SystemExit(1) from None
