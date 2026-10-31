#!/usr/bin/env python3
"""Exact determinant and concavity identities used in the proof.
Python 3.12+, SymPy 1.14; runtime about 3 seconds. No floating-point arithmetic.
"""
import json
from pathlib import Path
import sympy as s
from certify_subspaces import NORMAL_FORMS, check

t,p,q,a,V=s.symbols('t p q a V',positive=True)
T,U=s.symbols('T U',real=True)
x=s.symbols('A B C D E',positive=True)
records=[]

def verify(name,i,curve,F,expected,relation=None):
 R=s.Matrix([[s.Rational(v) for v in row] for row in NORMAL_FORMS[name][i]])
 xx=x[:R.cols]
 if name=='S3': fs=[xx[0]*xx[1]-p,xx[2]*xx[3]-q]
 elif name=='S7': fs=[xx[0]*xx[1]-p*xx[3],xx[2]**2-q*xx[3]]
 else: fs=[xx[0]*xx[1]-p*xx[4],xx[2]*xx[3]-q*xx[4]]
 substitution=dict(zip(xx,curve))
 conserved=R*s.Matrix(curve)
 check(all(s.simplify(s.diff(v,t))==0 for v in conserved),f'{name}/{i}: scalar curve remains in its compatibility class')
 residuals=[s.factor(v.subs(substitution)) for v in fs]
 check(sum(v==0 for v in residuals)==1,f'{name}/{i}: one binomial vanishes identically along the scalar curve')
 residual=next(v for v in residuals if v!=0)
 factor=s.factor(s.cancel(residual/F))
 check(factor in (s.Integer(1),t**-2),f'{name}/{i}: remaining equation has positive scalar factor {factor}')
 D=R.col_join(s.Matrix(fs).jacobian(xx)).det()
 D=s.factor(D.subs(dict(zip(xx,curve))))
 target=s.diff(F,t) if relation is None else relation
 ratio=s.factor(s.cancel(D/target))
 check(s.simplify(ratio-expected)==0 or s.simplify(ratio+expected)==0,
       f'{name}/{i}: determinant factor {ratio} (up to fixed orientation)')
 records.append({'structure':name,'form':i,'determinant_factor':str(ratio),'scalar':str(F)})

# S3: all scalar curves satisfy the conservation laws and the other binomial.
verify('S3',0,[p/T,T,t,U-t],t*(U-t)-q,-T)
for i in range(1,7):
 R=s.Matrix(NORMAL_FORMS['S3'][i])
 alpha=-s.Rational(R[0,2]);beta=s.Rational(R[0,3])
 gamma=s.Rational(R[1,2]);delta=s.Rational(R[1,3])
 A=T+alpha*t-beta*q/t;B=U-gamma*t-delta*q/t
 F=A*B-p
 verify('S3',i,[A,B,t,q/t],F,-t)
 check(alpha>0 and beta>0 and gamma>0 and delta>=0,f'S3/{i}: concavity coefficients')
 check(s.diff(A,t,2)==-2*beta*q/t**3 and s.diff(B,t,2)==-2*delta*q/t**3,
       f'S3/{i}: exact second derivatives')
verify('S3',7,[T,p/T,t,U-p/T-t],t*(U-p/T-t)-q,T)
F=(U-p/t)*(T-t)-q
verify('S3',8,[t,p/t,U-p/t,T-t],F,t)
check(s.Poly(s.expand(t*F),t).degree()==2,'S3/8: positive denominator leaves a quadratic')
for i,h in [(9,s.Rational(1,2)),(10,s.Integer(1))]:
 verify('S3',i,[a,p/a,t,V-t],t*(V-t)-q,h*(a+p/a))

# S7: coordinates are A,B,C,E; C^2=qE, AB=pE.
verify('S7',0,[p*(t-U)/(T-t+U),T-t+U,t,t-U],t**2-q*(t-U),T-t+U)
for i in [1,2]:
 R=s.Matrix(NORMAL_FORMS['S7'][i])
 alpha=-s.Rational(R[0,2]);beta=s.Rational(R[0,3])
 gamma=s.Rational(R[1,2]);delta=s.Rational(R[1,3])
 A=T+alpha/t-beta/(q*t**2);B=U-gamma/t-delta/(q*t**2)
 H=s.factor((t*A)*(t*B)-p/q)
 # The determinant identity holds off the scalar zero set too:
 # D=q*(H'-2H/t), so on H=0 its multiplier is the positive q.
 verify('S7',i,[A,B,1/t,1/(q*t**2)],H,q,s.diff(H,t)-2*H/t)
 check(beta>0 and delta>=0,f'S7/{i}: concavity coefficients')
 check(s.simplify(s.diff(t*A,t,2)+2*beta/(q*t**3))==0 and
       s.simplify(s.diff(t*B,t,2)+2*delta/(q*t**3))==0,
       f'S7/{i}: exact second derivatives')
verify('S7',3,[T+t**2/q,U,t,t**2/q],T*U+(U-p)*t**2/q,q)
verify('S7',4,[t,T-U**2/q-t,U,U**2/q],t*(T-U**2/q-t)-p*U**2/q,-q)

# S9: coordinates are A,B,C,D,F; AB=pF and CD=qF.
verify('S9',0,[p*t/(T-t),T-t,t+U,V,t],V*(t+U)-q*t,-(T-t))
A=T+(q/V-1)*t/2;B=U-(q/V+1)*t/2
verify('S9',1,[A,B,q*t/V,V,t],A*B-p*t,-V)
verify('S9',2,[T,U,t,V-T*U/p-t,T*U/p],t*(V-T*U/p-t)-q*T*U/p,-p)

out=Path(__file__).resolve().parents[1]/'certificates'
out.mkdir(exist_ok=True)
(out/'normal-forms.json').write_text(json.dumps(records,indent=2)+'\n')
check(len(records)==19,'all 19 normal forms have exact determinant identities')
