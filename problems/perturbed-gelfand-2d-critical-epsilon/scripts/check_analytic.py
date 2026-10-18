"""Exact symbolic and Arb checks of analytic identities. Runtime < 5 s.

Python 3, sympy, python-flint; one thread, below 1 GB. This checks the
displayed identities and inequalities, not the logical proof in the manuscript.
"""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
import json
from pathlib import Path
import time
import sympy as s
from flint import arb,ctx

def check(name,yes):
    print(('PASS ' if yes else 'FAIL ')+name)
    if not yes:raise SystemExit(1)

def main():
    start=time.monotonic();ctx.prec=160
    Path('certificates/analytic-checks.json').unlink(missing_ok=True)
    w,e,t,a=s.symbols('w e t a',positive=True)
    D=1+e*w; l=w/D; f=s.exp(l)
    lw=1/D**2; lww=-2*e/D**3; lwww=6*e**2/D**4
    le=-w*w/D**2; lwe=-2*w/D**3; lwwe=(-2+4*e*w)/D**4
    formulas={(1,0):f*lw,(2,0):f*(lw**2+lww),(3,0):f*(lw**3+3*lw*lww+lwww),
              (0,1):f*le,(1,1):f*(lw*le+lwe),(2,1):f*((lw*lw+lww)*le+2*lw*lwe+lwwe)}
    for (i,j),formula in formulas.items():
        check(f'exact derivative f_({i},{j})',s.simplify(s.diff(f,w,i,e,j)-formula)==0)
    explicit=a-2*s.log(1+s.exp(a)*t/2)
    check('exact unperturbed Gelfand solution',s.simplify(t*s.diff(explicit,t,2)+s.diff(explicit,t)+s.exp(explicit))==0)
    check('exact sublinearity identity',s.simplify(w/D**2-(1/(4*e)-(1-e*w)**2/(4*e*D**2)))==0)
    check('complex initial denominator exceeds one',1+arb('.242')*arb('5.94')-arb('.243')*2-arb('.01')*arb('7.95')>1)
    check('complex forcing bound exp(8) < 10000',arb(8).exp()<10000)
    check('initial contraction 10000 exp(-60) < 1',10000*arb(-60).exp()<1)
    q=s.symbols('q',positive=True);g=s.exp(-1/q)
    reduced=[g/q**2,g*(1-2*q)/q**4,g*(1-6*q+6*q*q)/q**6]
    for k,formula in enumerate(reduced,1):
        check('exact reduced derivative order '+str(k),s.simplify(s.diff(g,q,k)-formula)==0)
    check('exact transformed exponent',s.simplify(w/(1+e*w)-(1/e-1/(e+e**2*w)))==0)
    check('reduced complex half-plane bound (exact rationals)',s.Rational(49,100)-s.Rational(1,4)==s.Rational(6,25))
    check('reduced complex Lipschitz bound',1/arb('.24')**2<18)
    check('reduced origin contraction and ball',18*arb(-50).exp()<1 and arb(-50).exp()<arb(1)/8)
    check('flat derivative envelope monotonicity through order three',arb(1)/40<arb(1)/6)
    check('alpha >= 24 and epsilon >= .24 imply transformed a > 8/5',arb('.24')+24*arb('.24')**2>arb('1.6'))
    data={'status':'PASS','reduced_initial_lipschitz_upper':'18',
          'uniform_tail_conversion':{'epsilon_lower':'0.24','alpha_lower':'24','transformed_a_lower':'1.6224'},
          'tail_exclusion_also_requires_reduced_cover_and_limiting_problem_theorem':True,
          'runtime_seconds':round(time.monotonic()-start,3)}
    Path('certificates/analytic-checks.json').write_text(json.dumps(data,indent=2)+'\n')
    print('Runtime seconds:',data['runtime_seconds'])

if __name__=='__main__':main()
