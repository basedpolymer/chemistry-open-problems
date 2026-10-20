#!/usr/bin/env python3
"""Exact algebra checks for the stationary proof.

Environment: Python 3 with sympy 1.14 (tested). Runtime: about 1-10 seconds.
One computational thread, well below 1 GB RAM. Run from this problem folder.
These checks verify identities and SIX coefficient signs, not the infinite
analytic continuation proof. A failed assertion prints FAIL and exits nonzero.
"""
import os
for key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[key] = "1"
import json
from pathlib import Path
import sympy as s

ROOT = Path(__file__).resolve().parents[1]
claims = []


def check(name, value):
    good = bool(value)
    print(("PASS " if good else "FAIL ") + name, flush=True)
    claims.append({"claim": name, "pass": good})
    if not good:
        raise AssertionError(name)


def zero(name, expression):
    check(name, s.cancel(expression) == 0)


def main():
    z, m, p, v, r, rp, R = s.symbols("z m p v r rp R")
    a = 2*m+1
    N = m*(1-m)
    h = N-(1-z)*(N+a*p-p*p)/(1+z)
    zero("H recovery satisfies implicit generating equation",
         (1-z)*(a*p-p*p)-2*z*N+(1+z)*h)
    B = z*(1-z*z)
    D = 2*z*(m-p)*(p+1-m)-(1-z)**2*p
    # Differentiate recovered H and impose H'=P/z.
    diff = s.diff(h, z)+s.diff(h, p)*v-p/z
    zero("differentiation gives the P ordinary differential equation",
         diff.subs(v, D/(B*(a-2*p))))
    target = (m-r+2*z*(z+2)*r*r)/(z*(1+z)*(1+2*(1-z)*r))
    expr = (B*(a-2*p)*v-D).subs({p:m-(1-z)*r, v:r-(1-z)*rp})
    zero("r=(m-P)/(1-z) regularization", expr.subs(rp, target))
    zero("coefficient at the removable endpoint",
         target.subs(z, 1)-(m-r+6*r*r)/2)
    zero("large-r endpoint coefficient is 3",
         s.limit((target/r**2).subs(z, 1), r, s.oo)-3)
    zero("first-contact numerator at P=m", D.subs(p, m)+m*(1-z)**2)
    pstar = m+s.Rational(1,2)
    zero("turning-point numerator is strictly negative for z>1",
         D.subs(p, pstar)+s.Rational(3,2)*z+(1-z)**2*pstar)
    inverse = B*(a-2*p)/D
    zero("inverse has a quadratic turning point",
         s.diff(inverse,p).subs(p,pstar)+2*B/D.subs(p,pstar))
    bsq = ((R-1)**2*pstar+s.Rational(3,2)*R)/(R*R-1)
    zero("square-root amplitude formula",
         R*D.subs({z:R,p:pstar})/B.subs(z,R)-bsq)
    zero("coefficient transfer amplitude",
         s.Rational(2,3)/s.gamma(-s.Rational(3,2))-1/(2*s.sqrt(s.pi)))
    rho = [s.Integer(0)]
    tail = N
    coefficient_data = []
    for n in range(1,7):
        conv = sum(i*(n-i)*rho[i]*rho[n-i] for i in range(1,n))
        c = s.cancel((2*tail+conv)/(a*n+1))
        zero(f"recursion coefficient {n}", (a*n+1)*c-2*tail-conv)
        rho.append(c)
        tail = s.cancel(tail-c)
        num, den = s.fraction(s.cancel(c/(m*(1-m))))
        np = s.Poly(num,m)
        dp = s.Poly(den,m)
        # Positive monomial coefficients suffice for m>0; m(1-m)>0.
        positive = (all(x>=0 for x in np.all_coeffs()) and
                    any(x>0 for x in np.all_coeffs()) and
                    all(x>=0 for x in dp.all_coeffs()) and
                    any(x>0 for x in dp.all_coeffs()))
        check(f"rho_{n}>0 for every 0<m<1", positive)
        coefficient_data.append({"n":n, "numerator":str(s.factor(num)),
                                 "denominator":str(s.factor(den))})
    certificate = {"arithmetic":"exact symbolic rational arithmetic",
                   "scope":"Identities and first six coefficient signs only",
                   "claims":claims, "rho_n_divided_by_m_times_one_minus_m":coefficient_data}
    out=ROOT/"certificates"/"algebra.json"
    out.write_text(json.dumps(certificate,indent=2)+"\n",encoding="utf-8")
    print("PASS all exact algebra checks; certificate: certificates/algebra.json")


if __name__=="__main__":
    main()
