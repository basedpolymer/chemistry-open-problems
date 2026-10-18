"""Verifier's NON-RIGOROUS comparison of dimensions 2 and 3 (about 1 minute, one thread).

Same Du-Lou reduction in dimension n:  q'' + (n-1) q'/s + exp(-1/q) = 0, q(0) = a.  E_n(a) is the value
of q at the first (downward) zero of v = q_a; the critical secondary parameter is max_a E_n(a).
For n = 2 this reproduces 0.2421061656 (certified elsewhere); for n = 3 it gives about 0.2387971, the
three-dimensional value of Minamoto-Nakao (0.238797186), and is close to the "0.2385" printed in the
published version of Korman-Li-Ouyang.  Floating point (scipy DOP853), evidence only.
"""
import math, numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import minimize_scalar
def g(q,k):
    if q<=0: return 0.0
    e=math.exp(-1/q)
    return [e, e/q**2, e*(1-2*q)/q**4][k]
def E(a,n):
    def rhs(s,y):
        q,qp,v,vp,z,zp=y
        return [qp,-(n-1)/s*qp-g(q,0),vp,-(n-1)/s*vp-g(q,1)*v,zp,-(n-1)/s*zp-g(q,1)*z-g(q,2)*v*v]
    s0=1e-6
    y0=[a-g(a,0)*s0**2/(2*n),-g(a,0)*s0/n,1-g(a,1)*s0**2/(2*n),-g(a,1)*s0/n,-g(a,2)*s0**2/(2*n),-g(a,2)*s0/n]
    ev=lambda s,y:y[2]; ev.terminal=True; ev.direction=-1
    sol=solve_ivp(rhs,(s0,1e4),y0,method='DOP853',rtol=1e-12,atol=1e-14,events=ev)
    return sol.y_events[0][0][0] if sol.t_events[0].size else float('nan')
for n in [2,3]:
    grid=np.linspace(0.3,1.2,46)
    vals=[E(a,n) for a in grid]
    i=int(np.nanargmax(vals))
    r=minimize_scalar(lambda a:-E(a,n),bracket=(grid[i-1],grid[i],grid[i+1]),tol=1e-10)
    print(f"n={n}: a* ~ {r.x:.8f}, eps0 ~ {-r.fun:.10f}")
