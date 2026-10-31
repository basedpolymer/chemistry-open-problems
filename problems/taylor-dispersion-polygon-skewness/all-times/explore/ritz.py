#!/usr/bin/env python3
"""EXPLORATION (floats): Ritz in D_n-invariant polynomials f(p,q); residual and tail budget."""
import numpy as np, math, sys
from numpy.polynomial import legendre as Lg
from scipy.special import roots_legendre, roots_jacobi
from scipy.linalg import eigh

class Inv:
    def __init__(self, n, D, Nq=None):
        self.n=n; self.t=math.tan(math.pi/n); self.R=math.sqrt(1+self.t**2)
        self.D=D
        idx=[]
        for k in range(0, D//n+1):
            for a in range(0,(D-n*k)//2+1):
                idx.append((a,k))
        self.idx=idx
        Nq=Nq or max(60, int(1.2*(D+n)))
        v,wv=roots_legendre(Nq); v=.5*(v+1); wv=.5*wv
        s,ws=roots_jacobi(Nq,0,1); s=.5*(s+1); ws=.25*ws
        S,V=np.meshgrid(s,v,indexing='ij')
        self.w=(np.outer(ws,wv)*self.t).ravel()
        self.area=self.w.sum()
        self.w/=self.area    # average on T
        x=S.ravel(); y=(S*V).ravel()*self.t
        self.x,self.y=x,y
        tau,wt=roots_legendre(4*Nq); tau=.5*(tau+1)*self.t; wt=.5*wt  # average over side (wt sums to 1)
        self.tau=tau; self.wt=wt
    def bas(self,x,y,deriv=False):
        n=self.n; R2=self.R**2; Rn=self.R**n
        z=x+1j*y
        p=(x*x+y*y); P=2*p/R2-1
        wn=z**n
        q=wn.real/Rn
        m=len(self.idx)
        amax=max(a for a,k in self.idx); kmax=max(k for a,k in self.idx)
        # Legendre tables
        Lv=np.zeros((amax+1,)+p.shape); Lp=Lv.copy(); Lpp=Lv.copy()
        for a in range(amax+1):
            c=np.zeros(a+1); c[a]=1
            Lv[a]=Lg.legval(P,c)
            Lp[a]=Lg.legval(P,Lg.legder(c)) *2/R2 if a>0 else 0
            Lpp[a]=Lg.legval(P,Lg.legder(c,2))*(2/R2)**2 if a>1 else 0
        qk=[q**k for k in range(kmax+1)]
        f=np.zeros((m,)+p.shape); fp=f.copy(); fq=f.copy(); fpp=f.copy(); fqq=f.copy(); fpq=f.copy()
        for i,(a,k) in enumerate(self.idx):
            f[i]=Lv[a]*qk[k]
            fp[i]=Lp[a]*qk[k]
            fpp[i]=Lpp[a]*qk[k]
            if k>0:
                fq[i]=Lv[a]*k*qk[k-1]; fpq[i]=Lp[a]*k*qk[k-1]
            if k>1:
                fqq[i]=Lv[a]*k*(k-1)*qk[k-2]
        # derivatives wrt raw variables: q_raw = Rn q  =>  d/dq_raw = (1/Rn) d/dq
        grad=dict(fp=fp,fq=fq/Rn)
        lap=4*fp+4*p*fpp+4*n*wn.real*fpq/Rn+n**2*p**(n-1)*fqq/Rn**2
        # gradient components
        wn1=z**(n-1)
        qx=n*wn1.real; qy=-n*wn1.imag
        gx=fp*2*x+grad['fq']*qx
        gy=fp*2*y+grad['fq']*qy
        return f,gx,gy,lap

def torsion(n,K,M=6000):
    t=math.tan(math.pi/n); R=math.sqrt(1+t*t)
    tau=t*.5*(1-np.cos(np.linspace(0,np.pi,M)))
    w=(1+1j*tau)/R
    wts=np.gradient(tau); W=np.sqrt(wts)
    A=np.array([np.real(w**(n*k)) for k in range(K+1)]).T
    b=-(1-(1+tau**2))/2
    co,*_=np.linalg.lstsq(A*W[:,None],b*W,rcond=1e-13)
    res=A@co-b
    def ut(x,y):
        z=(x+1j*y)/R
        val=(1-x*x-y*y)/2
        for k,c in enumerate(co):
            val=val+c*np.real(z**(n*k))
        return val
    return ut,co,np.abs(res).max(),math.sqrt((res**2*wts).sum()/t)

def run(n,D,Mcut,KL=None,verbose=True):
    I=Inv(n,D)
    KL=KL or max(4,int(120/n))
    ut,co,dinf,dl2=torsion(n,KL)
    f,gx,gy,lap=I.bas(I.x,I.y)
    w=I.w
    Mm=(f*w)@f.T
    Kk=(gx*w)@gx.T+(gy*w)@gy.T
    # conditioning: scale
    d=np.sqrt(np.diag(Mm)); f=(f.T/d).T;gx=(gx.T/d).T;gy=(gy.T/d).T;lap=(lap.T/d).T
    Mm=(f*w)@f.T; Kk=(gx*w)@gx.T+(gy*w)@gy.T
    ev,U=np.linalg.eigh(Mm)
    keep=ev>1e-13
    T=U[:,keep]/np.sqrt(ev[keep])
    Kt=T.T@Kk@T
    mu,V=np.linalg.eigh(Kt)
    C=T@V                     # columns: modes
    psi=C.T@f; 
    # mode 0: constant
    i0=np.argmin(mu); 
    uv=ut(I.x,I.y); c0=(uv*w).sum(); u=uv-c0
    u2=(u*u*w).sum()
    # boundary
    fb,gxb,gyb,lapb=I.bas(np.ones_like(I.tau),I.tau)
    fb=(fb.T/d).T;gxb=(gxb.T/d).T;lapb=(lapb.T/d).T
    psib=C.T@fb; bk=C.T@gxb
    Rint=-(C.T@lap)-mu[:,None]*psi
    cc=psi@(u*w)
    # keep modes
    order=[i for i in np.argsort(mu) if mu[i]>1e-6]
    res=[]
    for i in order[:Mcut]:
        Rk=Rint[i]; Rk=Rk-(Rk*w).sum()
        nR=math.sqrt((Rk**2*w).sum()); nb=math.sqrt(2*(bk[i]**2*I.wt).sum())
        res.append((mu[i],cc[i],nR,nb))
    res=np.array(res)
    CP=I.R/math.pi; tr=math.sqrt(2*CP**2+2*I.R*CP)
    ell=CP*res[:,2]+tr*res[:,3]
    dk=res[:,1]/res[:,0]
    L1=(np.abs(dk)*ell/np.sqrt(2*res[:,0])).sum()
    sel=np.array(order[:Mcut])
    gh=(dk[:,None]*psi[sel]).sum(axis=0)
    # exact-ish g: Galerkin over all modes
    allm=np.array(order)
    gfull=((cc[allm]/mu[allm])[:,None]*psi[allm]).sum(axis=0)
    rho=math.sqrt(((gfull-gh)**2*w).sum())
    uh=(res[:,1][:,None]*psi[sel]).sum(axis=0)
    tailu=math.sqrt(((u-uh)**2*w).sum())
    tailufull=math.sqrt(max(u2-(cc[allm]**2).sum(),0))
    if verbose:
        print(f"n={n} D={D} dim={len(I.idx)} rank={keep.sum()} KL={KL} res_inf={dinf:.2e} res_l2={dl2:.2e} c={c0:.6f} u2={u2:.6f}  mu1..={res[:6,0].round(3)}")
        print(f"   mu_M={res[-1,0]:.1f}  L1={L1:.3e} rho(tail g)={rho:.3e} tail u={tailu:.3e} (full space {tailufull:.2e})")
        print(f"   ell_k first: {ell[:5]}  last: {ell[-3:]}")
        print(f"   nR  : {res[[0,1,2,-2,-1],2]}")
        print(f"   nb  : {res[[0,1,2,-2,-1],3]}")
    return I,res,ell,dk,psi[sel],u
if __name__=='__main__':
    for n,D,Mc in [(4,60,30),(5,60,25),(6,60,20),(8,60,15),(12,60,12),(16,64,10),(32,80,10),(64,128,10)]:
        run(n,D,Mc)
