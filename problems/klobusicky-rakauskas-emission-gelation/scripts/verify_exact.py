#!/usr/bin/env python3
"""Exact checks for emission coagulation. Python 3.12+, SymPy 1.14.
Run: python scripts/verify_exact.py. Typical runtime: under one minute.
These checks supplement the all-parameter analytic proof; finite checks alone
do not certify an infinite system. No floating-point arithmetic is used.
"""
from fractions import Fraction as Q
from math import comb
import json
from pathlib import Path
import sympy as s

checks = []


def check(name, value):
    value = bool(value)
    checks.append({"claim": name, "pass": value})
    print(("PASS " if value else "FAIL ") + name, flush=True)
    if not value:
        raise AssertionError(name)


def choose(x, r):
    out = Q(1)
    for i in range(r):
        out *= (x-i)/Q(i+1)
    return out


def coefficient(ell, k, j):
    n = ell+(k-ell)*j
    p = Q(2*k-ell, ell)
    return Q(ell**(j-1)*k**j, j*n)*choose(p*j, j-1)


# Reconstruct all coefficients by the triangular ODE, independently of inversion.
for ell, k in [(e, k) for e in range(1, 8) for k in range(e+1, e+9)]:
    values = [Q(0), Q(1)]
    for j in range(2, 41):
        values.append(sum(Q((ell+(k-ell)*i)*(ell+(k-ell)*(j-i)))
                          *values[i]*values[j-i] for i in range(1,j))/Q(j-1))
    check(f"recurrence and positivity ell={ell}, k={k}, j=1..40",
          all(values[j] == coefficient(ell,k,j) > 0 for j in range(1,41)))

# General characteristic equations: logarithmic derivatives, exact symbols.
h, a, tau, z, ell, k, t = s.symbols('h a tau z ell k t', positive=True)
p = a+1
W = h/(ell*tau*(1-h))
hz = h*(1-h)/(z*(1-p*h))
htau = h*(1-h)/(tau*(1-p*h))
Wt = s.diff(W,tau)+s.diff(W,h)*htau
Wz = s.diff(W,h)*hz
check("characteristic PDE W_tau = ell W^2 + ell a z W W_z",
      s.factor(Wt-ell*W**2-ell*a*z*W*Wz) == 0)
check("log derivative of h(1-h)^a",
      s.factor(1/h-a/(1-h)-(1-p*h)/(h*(1-h))) == 0)

d = k-ell
M = 1-ell*t/k
h0 = ell*t/k
p0 = (2*k-ell)/ell
mass = (M**2*W).subs({tau:t/k**2,h:h0})
check("self-consistency: sum n u_n = k-ell t", s.factor(mass-(k-ell*t))==0)
m2 = (k-ell*t)*(ell+d/(1-p0*h0))
m2_expected = (k-ell*t)*(k**2-ell*(2*k-ell)*t)/(k-(2*k-ell)*t)
check("second moment exact formula", s.factor(m2-m2_expected)==0)
check("second moment differential identity",
      s.factor(s.diff(m2,t)-(2*m2**2/(k-ell*t)**2-4*ell*m2/(k-ell*t)+ell**2))==0)
check("critical time equation", s.factor((1-p0*h0).subs(t,k/(2*k-ell)))==0)
check("critical tail constant squared",
      s.factor(((k*(2*d/ell)**2/p0)/((2*d/ell)*d))**2
               *p0/(2*(2*d/ell))-k**2/(d*(2*k-ell)))==0)

# Direct substitution in the original (ungauged) dimer equations.
tt = s.symbols('t')
u = {n: s.Rational(comb(3*n-3,n-2), n*(n-1)*8**(n-1))
        *tt**(n-2)*(2-tt)**(2*n) for n in range(2,19)}
for n in u:
    gain = sum(i*(n+1-i)*u[i]*u[n+1-i] for i in range(2,n))
    residual = s.diff(u[n],tt)-gain/(2-tt)**2+2*n*u[n]/(2-tt)
    check(f"original dimer ODE n={n}",s.cancel(residual)==0)
check("initial dimer data",all(u[n].subs(tt,0)==(1 if n==2 else 0) for n in u))

# Separate exact mass enclosure: coefficients positive; H's branch is known.
# Bound mass tail by evaluating its positive series at a rational z>1.
for time in [Q(1,10),Q(1,3),Q(1,2),Q(3,5)]:
    # Solve h(1-h)^2 = q*z with an exact rational supersolution H<1/3.
    hbase = time/2
    H = (hbase+Q(1,3))/2
    q = hbase*(1-hbase)**2
    radius_witness = H*(1-H)**2/q
    N = 100
    partial_mass = sum(Q(n*comb(3*n-3,n-2), n*(n-1)*8**(n-1))
                       *time**(n-2)*(2-time)**(2*n) for n in range(2,N+1))
    mass_at_witness = (2-time)**2/time*H/(1-H)
    tail_bound = mass_at_witness/radius_witness**N
    check(f"exact positive mass tail enclosure at t={time}",
          radius_witness>1 and 0 < (2-time)-partial_mass < tail_bound)

out = Path(__file__).resolve().parents[1]/'certificates'/'exact-checks.json'
out.parent.mkdir(exist_ok=True)
out.write_text(json.dumps({"arithmetic":"exact rational and symbolic", "checks":checks,
                           "all_pass":all(c['pass'] for c in checks)}, indent=2)+'\n')
print(f"PASS ALL ({len(checks)} checks); report: {out}")
