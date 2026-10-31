import numpy as np, sys, math
from galerkin import *
p = int(sys.argv[1]) if len(sys.argv)>1 else 26
ns = [int(x) for x in sys.argv[2:]] or [4,5,6,7,8,10,12,16]
grid = np.concatenate([[0.0], np.geomspace(1e-5, 3.0, 61)])
for n in ns:
    P = Polygon(n, p)
    l = P.lam[1:]; c = P.ui[1:]; T = P.Tij[1:,1:]
    E = np.exp(-np.outer(grid, l))*c   # (ng, m)
    Gm = E@T@E.T
    # normalise by leading decays to look at sign robustly
    print('n',n,'T11',T[0,0],'u1^2 T11',c[0]**2*T[0,0], 'u1<u^2 phi1>', c[0]*(P.u**2*P.E[1]*P.w).sum(), 'G(0,0)',Gm[0,0],' u3', P.u3)
    i,j = np.unravel_index(np.argmin(Gm), Gm.shape)
    print('   min G', Gm[i,j], 'at', grid[i], grid[j], ' #neg', (Gm<0).sum(), 'of', Gm.size)
    # H on fine grid
    s = np.geomspace(1e-5,5,200)
    Hv = H_of_tau(P.lam,P.ui,P.Tij,s)
    print('   min H/s^2 ', (Hv/s**2).min(), 'H<0 count', (Hv<0).sum(), ' H/(s e^{-l1 s}) at end', Hv[-1]/(s[-1]*np.exp(-l[0]*s[-1])))
    # sign pattern for negative G
    neg = Gm<0
    if neg.any():
        ii = np.argwhere(neg)
        print('   neg region a range', grid[ii[:,0]].min(), grid[ii[:,0]].max(), 'b', grid[ii[:,1]].min(), grid[ii[:,1]].max())
        ratio = (Gm/np.exp(-l[0]*(grid[:,None]+grid[None,:])))
        print('   ', np.min(ratio))
