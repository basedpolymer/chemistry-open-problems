#!/usr/bin/env python3
"""
EXPLORATION ONLY (floating point).  Neumann Galerkin on the half-sector T of the regular n-gon (apothem 1)
with Dubiner polynomials; torsion function from harmonic-polynomial boundary fit.
Returns eigenpairs lam, u_i, T_ij and helpers for G, H, K', K.
"""
import numpy as np, sys, math
from scipy.special import eval_jacobi, roots_legendre, roots_jacobi
from scipy.linalg import eigh

def jac(n, a, b, x):
    return eval_jacobi(n, a, b, x)

def djac(n, a, b, x):
    if n == 0:
        return np.zeros_like(x)
    return 0.5*(n+a+b+1)*eval_jacobi(n-1, a+1, b+1, x)

class Tri:
    def __init__(self, n, p, Nq=None):
        self.n, self.p = n, p
        self.t = math.tan(math.pi/n)
        Nq = Nq or (2*p+20)
        vs, wv = roots_legendre(Nq)
        vs = 0.5*(vs+1); wv = wv*0.5
        ss, ws = roots_jacobi(Nq, 0, 1)    # weight (1+x) on [-1,1]  -> s ds
        ss = 0.5*(ss+1); ws = ws*0.25
        S, V = np.meshgrid(ss, vs, indexing='ij')
        WW = np.outer(ws, wv)*self.t      # dA = t s ds dv  (x = s, y = t s v)
        self.s, self.v, self.w = S.ravel(), V.ravel(), WW.ravel()
        self.x = self.s; self.y = self.t*self.s*self.v
        idx = []
        for i in range(p+1):
            for j in range(i+1):
                idx.append((i, j))
        self.idx = idx
        m = len(idx)
        Q = len(self.s)
        phi = np.zeros((m, Q)); px = np.zeros((m, Q)); py = np.zeros((m, Q))
        s, v = self.s, self.v
        for k, (i, j) in enumerate(idx):
            a = jac(j, 0, 0, 2*v-1)
            da = 2*djac(j, 0, 0, 2*v-1)
            b0 = jac(i-j, 0, 2*j+1, 2*s-1)
            db0 = 2*djac(i-j, 0, 2*j+1, 2*s-1)
            b = s**j*b0
            db = j*s**(j-1)*b0 + s**j*db0 if j > 0 else db0
            bos = s**(j-1)*b0 if j > 0 else 0*s    # b/s (only needed when a' != 0, j>0)
            phi[k] = a*b
            px[k] = a*db - da*v*bos
            py[k] = da*bos/self.t
        # orthonormalise wrt mass (should be diagonal)
        M = (phi*self.w)@phi.T
        d = np.sqrt(np.diag(M))
        self.phi = (phi.T/d).T; self.px = (px.T/d).T; self.py = (py.T/d).T
        self.area = self.w.sum()
        self.M = (self.phi*self.w)@self.phi.T
        self.Kmat = (self.px*self.w)@self.px.T + (self.py*self.w)@self.py.T

def torsion(n, K=None, M=4000):
    """ut = (1-r^2)/2 + sum a_k Re((w/R)^{nk}), LSQ fit on side x=1, 0<=tau<=t.  Return function."""
    t = math.tan(math.pi/n); R = math.sqrt(1+t*t)
    K = K or max(8, min(40, int(240/n)))
    tau = t*0.5*(1-np.cos(np.linspace(0, np.pi, M)))
    w = (1+1j*tau)/R
    A = np.array([np.real(w**(n*k)) for k in range(K+1)]).T
    b = -(1-(1+tau**2))/2
    co, *_ = np.linalg.lstsq(A, b, rcond=1e-14)
    res = A@co-b
    def f(x, y):
        wz = (x+1j*y)/R
        val = (1-x*x-y*y)/2
        for k, c in enumerate(co):
            val = val + c*np.real(wz**(n*k))
        return val
    return f, np.abs(res).max()

class Polygon:
    def __init__(self, n, p):
        self.n = n
        tr = Tri(n, p)
        self.tr = tr
        f, res = torsion(n)
        self.res = res
        ut = f(tr.x, tr.y)
        avg = (ut*tr.w).sum()/tr.area
        self.c = avg
        self.u = ut-avg
        lam, vec = eigh(tr.Kmat/tr.area, tr.M/tr.area)  # normalised so that <phi^2>=1 w.r.t. average
        self.lam = lam
        self.vec = vec
        # eigenfunctions at quadrature points
        self.E = vec.T@tr.phi    # shape (m, Q); E[k] is k-th eigenfunction
        self.Kmat = tr.Kmat/tr.area
        w = tr.w/tr.area
        self.w = w
        self.ui = self.E@(self.u*w)
        self.Tij = (self.E*(self.u*w))@self.E.T
        self.u2 = (self.u**2*w).sum()
        self.u3 = (self.u**3*w).sum()
        self.lamK = None

    def proj_err(self):
        # L2 error of the projection of u
        return math.sqrt(max(self.u2-(self.ui**2).sum(), 0))

def kernels(lam, tau):
    """I_ij(tau) = int_0^tau (1-e^{-la a})(1-e^{-lb (tau-a)}) da   for la, lb > 0  (vectorised)"""
    la = lam[:, None]; lb = lam[None, :]
    # = tau - (1-e^{-la tau})/la - (1-e^{-lb tau})/lb + E, E = int_0^tau e^{-la a - lb (tau-a)} da
    d = la-lb
    with np.errstate(divide='ignore', invalid='ignore'):
        E = (np.exp(-lb*tau)-np.exp(-la*tau))/d
    i = np.abs(d) < 1e-9*(la+lb)
    E = np.where(i, tau*np.exp(-la*tau), E)
    return tau-(1-np.exp(-la*tau))/la-(1-np.exp(-lb*tau))/lb+E

def K_of_tau(lam, ui, Tij, taus, skip=1):
    """K(tau)=2 int_0^tau <u A(a) A(tau-a)> da, A(s)=sum (u_i/lam_i)(1-e^{-lam_i s}) phi_i"""
    l = lam[skip:]; g = ui[skip:]/l; T = Tij[skip:, skip:]
    out = []
    for tau in taus:
        I = kernels(l, tau)
        out.append(2*(g@(T*I)@g))
    return np.array(out)

def Kp_of_tau(lam, ui, Tij, taus, skip=1):
    """K'(s)=2 iint_{a+b<=s} G = 2 sum_ij u_i u_j T_ij J_ij(s), J=iint_{a+b<=s} e^{-la a - lb b}"""
    l = lam[skip:]; g = ui[skip:]/l; T = Tij[skip:, skip:]
    # J_ij*l_i*l_j: d/ds of I_ij ...  I'(s) = e^{...}: compute I' analytically
    la = l[:, None]; lb = l[None, :]
    out = []
    for s in taus:
        d = la-lb
        with np.errstate(divide='ignore', invalid='ignore'):
            E = (np.exp(-lb*s)-np.exp(-la*s))/d
            Ep = (-lb*np.exp(-lb*s)+la*np.exp(-la*s))/d
        i = np.abs(d) < 1e-9*(la+lb)
        E = np.where(i, s*np.exp(-la*s), E)
        Ep = np.where(i, (1-la*s)*np.exp(-la*s), Ep)
        Ip = 1-np.exp(-la*s)-np.exp(-lb*s)+Ep
        out.append(2*(g@(T*Ip)@g))
    return np.array(out)

def H_of_tau(lam, ui, Tij, taus, skip=1):
    l = lam[skip:]; a = ui[skip:]; T = Tij[skip:, skip:]
    la = l[:, None]; lb = l[None, :]
    out = []
    for s in taus:
        d = la-lb
        with np.errstate(divide='ignore', invalid='ignore'):
            E = (np.exp(-lb*s)-np.exp(-la*s))/d
        i = np.abs(d) < 1e-9*(la+lb)
        E = np.where(i, s*np.exp(-la*s), E)
        out.append(a@(T*E)@a)
    return np.array(out)

def G_ab(lam, ui, Tij, a, b, skip=1):
    l = lam[skip:]; c = ui[skip:]; T = Tij[skip:, skip:]
    return (c*np.exp(-l*a))@T@(c*np.exp(-l*b))

if __name__ == '__main__':
    p = int(sys.argv[1]) if len(sys.argv) > 1 else 24
    for n in [3, 4, 5, 6, 8, 12, 24]:
        P = Polygon(n, p)
        print(n, 'res', P.res, 'lam1..4', P.lam[1:5], 'c', P.c, 'u2', P.u2, 'u3', P.u3, 'E', P.c**2-P.u2, 'pe', P.proj_err(), flush=True)
        print('   T11', P.Tij[1, 1], 'u1', P.ui[1], 'T22', P.Tij[2,2], 'T12', P.Tij[1,2])
        taus = np.array([1e-4,3e-4,1e-3,3e-3,1e-2,3e-2,0.1,0.2,0.3,0.5,1,2,5,10])
        print('   K/tau^4', K_of_tau(P.lam, P.ui, P.Tij, taus)/taus**4)
        print('   Kp/tau^3', Kp_of_tau(P.lam, P.ui, P.Tij, taus)/taus**3)
        print('   Kp', Kp_of_tau(P.lam, P.ui, P.Tij, taus))
        print('   H', H_of_tau(P.lam, P.ui, P.Tij, taus))
