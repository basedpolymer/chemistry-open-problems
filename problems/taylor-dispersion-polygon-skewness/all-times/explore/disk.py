import numpy as np, math
from scipy.special import jn_zeros
M=4000
j=jn_zeros(1,M); lam=j**2
# modes phi_k(r)=J0(j r)/J0(j), <phi^2>=1, phi(1)=1; u=-2 sum phi_k/lam_k
c=-2/lam; d=c/lam
def T(M):
    l=lam[:M]
    Tm=-2*(l[:,None]+l[None,:])/np.where(np.eye(M)>0,1,(l[:,None]-l[None,:])**2)
    Tm[np.arange(M),np.arange(M)]=1/12.
    return Tm
Mm=1500
Tm=T(Mm); l=lam[:Mm]; cc=c[:Mm]; dd=d[:Mm]
print('check <u^2>=1/48', (cc**2).sum(), 1/48, ' <u^3>', cc@Tm@cc, ' <u g1^2>=1/184320', dd@Tm@dd, 1/184320., ' <u g>', (cc*dd).sum(),1/768)
def G(a,b):
    return (cc*np.exp(-l*a))@Tm@(cc*np.exp(-l*b))
for s in [1e-4,3e-4,1e-3,3e-3,5e-3,1e-2,2e-2,3e-2,5e-2,.1,.2,.5]:
    print(s,'G(s,0)',G(s,0),'G(s/2,s/2)',G(s/2,s/2),'lb',s/12-.752*s**1.5, ' (s,0)/s', G(s,0)/s)
def kern(tau,M):
    l=lam[:M]
    la=l[:,None]; lb=l[None,:]; dl=la-lb
    with np.errstate(divide='ignore',invalid='ignore'):
        E=(np.exp(-lb*tau)-np.exp(-la*tau))/dl
    i=np.arange(M); E[i,i]=tau*np.exp(-l*tau)
    return tau-(1-np.exp(-la*tau))/la-(1-np.exp(-lb*tau))/lb+E
for M2 in [5,8,12,20,50,200,1000]:
    dd2=d[:M2]; Tm2=Tm[:M2,:M2]
    out=[]
    for tau in [1e-3,3e-3,5e-3,1e-2,2e-2,5e-2,.1,.5,1,5]:
        out.append(2*dd2@(Tm2*kern(tau,M2))@dd2)
    print(M2,np.array(out))
print('K/tau^4 (M=1000)')
for tau in [1e-4,1e-3,3e-3,5e-3,1e-2,1.5e-2,2e-2,3e-2,5e-2,.1,.2,.5,1,2]:
    K=2*d[:1000]@(Tm[:1000,:1000]*kern(tau,1000))@d[:1000]
    L=1.504;E=1/24;c0=.25
    KLB=E*tau**4/3+c0*tau**5/15-8/35*c0*L*tau**4.5-128/3465*L*tau**5.5
    print(tau,K,K/tau**4,KLB/tau**4, 'K/(tau*1.085e-5)',K/(tau/92160))
