import numpy as np, math, sys
import ritz
from ritz import Inv, torsion
def runb(n,A0,A1,Kq,Lam,KL=8):
    I=Inv.__new__(Inv)
    I.n=n; I.t=math.tan(math.pi/n); I.R=math.sqrt(1+I.t**2)
    D=max(2*A0, n*Kq+2*A1); I.D=D
    I.idx=[(a,0) for a in range(A0+1)]+[(a,k) for k in range(1,Kq+1) for a in range(A1+1)]
    Nq=max(80,int(0.8*D)+40)
    from scipy.special import roots_legendre, roots_jacobi
    v,wv=roots_legendre(Nq); v=.5*(v+1); wv=.5*wv
    s,ws=roots_jacobi(Nq,0,1); s=.5*(s+1); ws=.25*ws
    S,V=np.meshgrid(s,v,indexing='ij')
    I.w=(np.outer(ws,wv)*I.t).ravel(); I.area=I.w.sum(); I.w/=I.area
    I.x=S.ravel(); I.y=(S*V).ravel()*I.t
    tau,wt=roots_legendre(3*Nq); tau=.5*(tau+1)*I.t; wt=.5*wt
    I.tau=tau; I.wt=wt
    ut,co,dinf,dl2=torsion(n,KL)
    f,gx,gy,lap=I.bas(I.x,I.y); w=I.w
    Mm=(f*w)@f.T; d=np.sqrt(np.diag(Mm)); f=(f.T/d).T;gx=(gx.T/d).T;gy=(gy.T/d).T;lap=(lap.T/d).T
    Mm=(f*w)@f.T; Kk=(gx*w)@gx.T+(gy*w)@gy.T
    ev,U=np.linalg.eigh(Mm); keep=ev>1e-14
    T=U[:,keep]/np.sqrt(ev[keep]); mu,Vv=np.linalg.eigh(T.T@Kk@T); C=T@Vv
    psi=C.T@f
    uv=ut(I.x,I.y); c0=(uv*w).sum(); u=uv-c0
    fb,gxb,gyb,lapb=I.bas(np.ones_like(I.tau),I.tau); gxb=(gxb.T/d).T
    bk=C.T@gxb; Rint=-(C.T@lap)-mu[:,None]*psi
    cc=psi@(u*w)
    order=[i for i in np.argsort(mu) if mu[i]>1e-6]
    sel=[i for i in order if mu[i]<=Lam]
    CP=I.R/math.pi; tr=math.sqrt(2*CP**2+2*I.R*CP)
    rows=[]
    for i in sel:
        Rk=Rint[i]-(Rint[i]*w).sum()
        nR=math.sqrt((Rk**2*w).sum()); nb=math.sqrt(2*(bk[i]**2*I.wt).sum())
        rows.append((mu[i],cc[i],nR,nb))
    r=np.array(rows); ell=CP*r[:,2]+tr*r[:,3]; dk=r[:,1]/r[:,0]
    contrib=np.abs(dk)*ell/np.sqrt(2*r[:,0]); L1=contrib.sum()
    gh=(dk[:,None]*psi[sel]).sum(axis=0)
    allm=np.array(order); gfull=((cc[allm]/mu[allm])[:,None]*psi[allm]).sum(axis=0)
    rho=math.sqrt(((gfull-gh)**2*w).sum())
    print(f"n={n} dim={len(I.idx)} rank={keep.sum()} M={len(sel)} Nq={Nq} L1={L1:.3e} rho={rho:.3e} ell1={ell[0]:.2e} contrib top {np.sort(contrib)[-3:]}",flush=True)
    return
cases=[(4,36,30,8,700),(5,36,24,8,700),(6,36,20,8,700),(8,34,14,8,700),(12,34,8,8,700),(16,34,6,8,700),(24,34,5,8,700),(32,34,4,8,700),(48,34,3,8,700),(64,34,3,8,700),(96,34,2,8,700)]
for c in cases:
    runb(*c)
