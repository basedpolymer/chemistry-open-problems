"""Rigorous local cusp certificate, Python 3 + python-flint 0.9 (Arb).

One thread, below 1 GB; expected runtime 1-5 minutes. All analytic assumptions
and the Taylor remainder construction are proved in paper.tex. No scipy is used.
Run from the problem folder. Every assertion prints PASS/FAIL and aborts on FAIL.
"""
import json
import math
from pathlib import Path
import time
from flint import arb, arb_series, arb_mat, ctx

ctx.prec=160
ORDER=16
ctx.cap=ORDER+3
ZERO=arb(0)

def check(name, condition):
    print(('PASS ' if condition else 'FAIL ')+name,flush=True)
    if not condition: raise SystemExit(1)

def ball(c,r):return arb(c,arb(r).abs_upper())
def hull(a,b):return arb(a).union(arb(b))

def rhs(x,y,e):
    w,v,z,c,b,d,h=y[:7]
    D=1+e*w
    F=(w/D).exp()
    lw=1/D**2; lww=-2*e/D**3; lwww=6*e**2/D**4
    le=-w*w/D**2; lwe=-2*w/D**3; lwwe=(-2+4*e*w)/D**4
    fw=F*lw; fww=F*(lw*lw+lww)
    fwww=F*(lw**3+3*lw*lww+lwww)
    fe=F*le; fwe=F*(lw*le+lwe)
    fwwe=F*((lw*lw+lww)*le+2*lw*lwe+lwwe)
    force=[F,fw*v,fw*z+fww*v*v,fw*c+3*fww*v*z+fwww*v**3,
           fw*b+fe,fw*d+fww*b*v+fwe*v,
           fw*h+fww*b*z+fwe*z+2*fww*v*d+fwww*b*v*v+fwwe*v*v]
    return y[7:]+[-x.exp()*a for a in force]

def coeff(x,y,e,n):
    cs=[[a] for a in y]
    for k in range(n):
        ctx.cap=k+2
        ys=[arb_series(a,prec=k+2) for a in cs]
        xs=arb_series([x,arb(1)],prec=k+2)
        rr=rhs(xs,ys,e)
        for j in range(14):cs[j].append(rr[j][k]/(k+1))
    ctx.cap=ORDER+3
    return cs

def tube(x,y,e,h):
    # Strict Picard inclusion is a validated existence/enclosure test.
    X=x+hull(0,h)
    f0=rhs(x,y,e)
    floor=h*h*(1+sum([a.abs_upper() for a in f0],ZERO))+sum([a.rad() for a in y],ZERO)/1000000+arb('1e-40')
    rad=[(a*h).abs_upper()*2+floor for a in f0]
    for _ in range(25):
        B=[a+ball(0,r) for a,r in zip(y,rad)]
        if not 1+e*B[0]>0:
            raise ArithmeticError('Nonlinearity denominator not positive on tube')
        f=rhs(X,B,e)
        Z=[a+hull(0,h)*b for a,b in zip(y,f)]
        if all(b.contains_interior(z) for b,z in zip(B,Z)):
            return X,B
        rad=[r.max((h*a).abs_upper()*2+floor) for r,a in zip(rad,f)]
    raise ArithmeticError('Picard inclusion failed')

def step(x,y,e,h,end_window=False):
    X,B=tube(x,y,e,h)
    C=coeff(x,y,e,ORDER)
    R=coeff(X,B,e,ORDER+1)
    H=hull(0,h) if end_window else h
    out=[]
    for c,r in zip(C,R):
        a=c[-1]
        for cj in c[-2::-1]:a=a*H+cj
        out.append(a+r[-1]*H**(ORDER+1))
    return out

def initial(a,e):
    # Complex Volterra contraction and Cauchy estimates; see manuscript.
    t=arb(-60).exp()
    ks=[(0,0),(1,0),(2,0),(3,0),(0,1),(1,1),(2,1)]
    errors=[10000*t*math.factorial(k)*100**j for k,j in ks]
    return [c+ball(0,r) for c,r in zip([a,arb(1),ZERO,ZERO,ZERO,ZERO,ZERO],errors)]+[ball(0,r) for r in errors]

def integrate(a,e,xc,rx=ZERO,verbose=False):
    x=arb(-60); y=initial(a,e)
    target=xc-rx
    # Exact rational endpoints (dyadic except the final step).
    n=0
    while x<target:
        h=arb(1)/(128 if x> -4 else 8)
        if not x+h<target:h=target-x
        for attempt in range(12):
            try:
                yn=step(x,y,e,h)
                if all(v.is_finite() for v in yn):break
            except ArithmeticError:pass
            h=h/2
        else:
            print('FAIL state at x',x,'values:',[strb(v) for v in y],flush=True)
            raise ArithmeticError('Step validation failed')
        y=yn
        x=x+h;n+=1
        if verbose and n%100==0:print('Validated steps:',n,'x:',x,'max radius:',max(float(v.rad()) for v in y),flush=True)
    if rx>0:y=step(x,y,e,2*rx,True)
    return y,n

def strb(x):return x.str(35)
def endpoints(x):return [x.lower().str(35),x.upper().str(35)]

def main():
    start=time.monotonic()
    Path('certificates/cusp-certificate.json').unlink(missing_ok=True)
    # Rational decimal centres; radii chosen independently of numerical accuracy.
    data=json.loads(Path('certificates/cusp-input.json').read_text())
    centers=data['centers']; radii=data['radii']
    c=list(map(arb,centers)); r=list(map(arb,radii))
    check('all box radii are positive',all(a>0 for a in r))
    check('parameter box within analytic initialization domain',c[0]-r[0]>arb('5.94') and c[0]+r[0]<arb('5.95') and c[1]-r[1]>arb('.242') and c[1]+r[1]<arb('.243'))
    Y0,n=integrate(c[0],c[1],c[2],verbose=True)
    F=arb_mat([[Y0[i]] for i in range(3)])
    print('Centre residual:',[strb(Y0[i]) for i in range(3)],flush=True)
    Y,n=integrate(ball(c[0],r[0]),ball(c[1],r[1]),c[2],r[2],verbose=True)
    # Columns: alpha, epsilon, logarithmic endpoint x.
    J=arb_mat([[Y[1],Y[4],Y[7]],[Y[2],Y[5],Y[8]],[Y[3],Y[6],Y[9]]])
    J0=arb_mat([[Y0[1].mid(),Y0[4].mid(),Y0[7].mid()],
                [Y0[2].mid(),Y0[5].mid(),Y0[8].mid()],
                [Y0[3].mid(),Y0[6].mid(),Y0[9].mid()]])
    # Midpoints of an approximate inverse are exact dyadic rationals.
    inv=J0.inv(); C=arb_mat([[inv[i,j].mid() for j in range(3)] for i in range(3)])
    I=arb_mat([[int(i==j) for j in range(3)] for i in range(3)])
    M=I-C*J
    correction=-C*F
    K=[]; norms=[]
    for i in range(3):
        radius=sum([M[i,j].abs_upper()*r[j] for j in range(3)],ZERO)
        K.append(correction[i,0]+ball(0,radius))
        norms.append(radius/r[i])
        check('Krawczyk strict inclusion coordinate '+str(i),K[-1].abs_upper()<r[i])
        check('Krawczyk offset below 1 percent of radius '+str(i),K[-1].abs_upper()<r[i]/100)
    check('weighted contraction norm < 1',all(q<1 for q in norms))
    check('weighted contraction norm < 1e-6',all(q<arb('1e-6') for q in norms))
    check('first zero transversality w_x < 0',Y[7]<0)
    # Eliminate x through w=0, then evaluate derivatives of the boundary v.
    Be=Y[5]-Y[8]*Y[4]/Y[7]
    Baa=Y[3]  # at the zero of (w,v,z), x_a=x_aa=0.
    check('positive unfolding derivative B_epsilon',Be>0)
    check('positive quadratic derivative B_alpha_alpha',Baa>0)
    lambda_box=4*ball(c[2],r[2]).exp()
    check('reported lambda enclosure',lambda_box>arb('3.006301478866') and lambda_box<arb('3.006301478873'))
    # An explicit local rectangle: all folds in this amplitude window are counted.
    ra=arb('0.002'); re=arb('0.000000002'); rx=arb('0.0001')
    check('full rectangle within initialization domain',c[0]-ra>arb('5.94') and c[0]+ra<arb('5.95') and c[1]-re>arb('.242') and c[1]+re<arb('.243'))
    V,_=integrate(ball(c[0],ra),ball(c[1],re),c[2],rx,verbose=True)
    check('rectangle endpoint transversality',V[7]<0)
    # Parameter mean-value theorem and time monotonicity bracket every first zero.
    margin=(-V[7]).lower()*rx-Y0[0].abs_upper()-V[1].abs_upper()*ra-V[4].abs_upper()*re
    check('every rectangle shooting zero is bracketed',margin>0)
    xa=-V[1]/V[7]
    X=ball(c[2],rx); ebox=ball(c[1],re)
    ff=(V[0]/(1+ebox*V[0])).exp()
    px=-X.exp()*ff
    qx=-X.exp()*ff/(1+ebox*V[0])**2*V[1]
    xaa=-(V[2]+2*V[8]*xa+px*xa**2)/V[7]
    Baa_rect=V[3]+2*V[9]*xa+qx*xa**2+V[8]*xaa
    Be_rect=V[5]-V[8]*V[4]/V[7]
    check('strict convexity B_alpha_alpha on full rectangle',Baa_rect>0)
    check('positive B_epsilon on full rectangle',Be_rect>0)
    check('rectangle bounds 0.005 < B_alpha_alpha < 0.011',Baa_rect>arb('.005') and Baa_rect<arb('.011'))
    check('rectangle bounds 2 < B_epsilon < 4',Be_rect>2 and Be_rect<4)
    edge_margin=Baa_rect.lower()*(ra-r[0])**2/2-Be_rect.abs_upper()*(re+r[1])
    check('B positive at both amplitude edges for every epsilon',edge_margin>0)
    output={'status':'PASS local cusp only','precision_bits':ctx.prec,'order':ORDER,'centers':centers,'radii':radii,
            'center_residual':[strb(Y0[i]) for i in range(3)],
            'jacobian':[[strb(J[i,j]) for j in range(3)] for i in range(3)],
            'krawczyk_offsets':[strb(a) for a in K],'weighted_norm_bounds':[strb(a) for a in norms],
            'B_epsilon':strb(Be),'B_alpha_alpha':strb(Baa),'w_x':strb(Y[7]),
            'lambda_box':strb(lambda_box),'steps_per_integration':n,
            'local_rectangle':{'alpha_radius':strb(ra),'epsilon_radius':strb(re),'x_radius':strb(rx),
                               'event_bracket_margin':strb(margin),'B_alpha_alpha':endpoints(Baa_rect),
                               'B_epsilon':endpoints(Be_rect),'edge_margin':strb(edge_margin)},
            'runtime_seconds':round(time.monotonic()-start,3)}
    Path('certificates/cusp-certificate.json').write_text(json.dumps(output,indent=2)+'\n')
    print('PASS local ordinary cusp; global uniqueness is NOT certified',flush=True)
    print('Runtime seconds:',output['runtime_seconds'],flush=True)

if __name__=='__main__':main()
