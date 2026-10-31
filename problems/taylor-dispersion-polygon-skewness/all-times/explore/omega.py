import numpy as np, math, sys
def fit(n,K,J=None,M=None):
    J=J or K
    t=math.tan(math.pi/n); R=math.sqrt(1+t*t)
    M=M or max(4000,40*n*K)
    # Chebyshev-like nodes on [0,t] with weights
    th=np.linspace(0,np.pi,M); tau=t*.5*(1-np.cos(th)); wt=np.gradient(tau); W=np.sqrt(wt)
    w=(1+1j*tau)
    lw=np.log(w/R)
    def rew(m): return np.real(np.exp(m*lw))
    def dxw(m): return np.real(m*np.exp((m-1)*lw))/R   # d/dx Re (w/R)^m on x=1
    # torsion
    A=np.array([rew(n*k) for k in range(K+1)]).T
    b=-(1-(1+tau**2))/2
    a,*_=np.linalg.lstsq(A*W[:,None],b*W,rcond=1e-14)
    res=A@a-b
    L2ut=math.sqrt(2*(res**2*wt).sum()/t)          # ||ut_N||_dP (avg normalisation (1/|P|) int_dP)
    Einf=np.abs(res).max()
    # particular solution: -Lap g = u_N = (1-r^2)/2 - m + sum a_k Re(w/R)^{nk}
    # mean of ut_N: from (1-r^2)/2 and harmonic terms
    r2=lambda: None
    # average of Re(w/R)^{m}: homogeneous degree m: (1/(m+2)) * int_0^t Re((1+i tau)^m)/R^m dtau / (t/2)
    def avg_h(m):
        return (2/t)*(1/(m+2))*(rew(m)*wt).sum()
    avr2=(1+t*t/3)/2
    mu=.5-avr2/2+sum(a[k]*avg_h(n*k) for k in range(K+1))
    # g_p = r^4/32 - (1/2-mu)... -Lap(r^2)=-4 => -Lap(c r^2)=-4c ; (const part C=1/2-mu + a_0) ;  -Lap(r^4)=-16 r^2
    C=.5-mu+a[0]
    # solve -Lap g_p = C - r^2/2 + sum_{k>=1} a_k Re(w/R)^{nk}
    # -Lap(-C r^2/4) = C ; -Lap(r^4/32) = -r^2/2 ; -Lap(r^2 Re w^m c) = -c 4 (1+m) Re w^m
    # flux on x=1: d/dx
    r2v=1+tau**2
    dxg=-C/2*1 + (4*r2v/32)*2*1/2*1  # placeholder overwritten below
    # d/dx (r^2)=2x=2; d/dx r^4 = 4 r^2 x = 4 r2v
    dxg=(-C/4)*2+(1/32)*4*r2v
    for k in range(1,K+1):
        m=n*k; c=-a[k]/(4*(1+m))
        dxg=dxg+c*(2*rew(m)+r2v*dxw(m))
    B=np.array([dxw(n*j) for j in range(1,J+1)]).T
    bb,*_=np.linalg.lstsq(B*W[:,None],-dxg*W,rcond=1e-14)
    e=dxg+B@bb
    L2e=math.sqrt(2*(e**2*wt).sum()/t)
    CP=R/math.pi; tr=math.sqrt(2*CP**2+2*R*CP)
    Enorm=math.sqrt((1+R)/2)*L2ut
    eta=CP*Enorm+tr*L2e
    om=CP*eta
    return Einf,L2ut,L2e,eta,om
for n in [4,5,6,8,12,16,24,32,48,64,100]:
    for K in [4,8,16,32,64]:
        if n*K>6000: continue
        print(n,K,'Einf %.2e L2ut %.2e L2e %.2e eta %.2e omega %.2e'%fit(n,K))
