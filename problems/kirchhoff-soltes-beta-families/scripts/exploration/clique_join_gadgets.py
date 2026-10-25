"""FLOATING-POINT EXPLORATION: gadgets K_p + H rooted at a vertex of H (all graphs H on <= 7 vertices, p <= 39); the clique class has a = delta + q/(p(p-1)) + n^2 Gamma_xx - (n-1)^2 Gamma1_xx; lists classes with a <= 0. Runtime about 2 min.
Environment: Fedora 44 (WSL), Python 3.12, numpy, numba, networkx, nauty geng."""
# gadget K_p + H rooted at x in H, class = clique; a = R - m rho; scan all graphs H on q<=7 vertices (any, incl. disconnected)
import sys, numpy as np, networkx as nx, subprocess
best=[]
for q in range(1,8):
    out=subprocess.run(['geng','-q',str(q)],capture_output=True,text=True).stdout.split()
    for s in out:
        H=nx.from_graph6_bytes(s.encode()); LH=nx.laplacian_matrix(H,nodelist=range(q)).toarray().astype(float)
        for p in range(2,40):
            n=p+q
            G=np.linalg.inv(LH+p*np.eye(q)); G1=np.linalg.inv(LH+(p-1)*np.eye(q))
            trG, trG1 = np.trace(G), np.trace(G1)
            delta = 1 + q/(p*(p-1)) + n*trG - (n-1)*trG1
            for x in range(q):
                a = delta + q/(p*(p-1)) + n*n*G[x,x] - (n-1)**2*G1[x,x]
                frac = p/(p+q-1)
                if a <= 1e-9:
                    best.append((frac, s, q, p, x, a, delta))
best.sort(key=lambda t:-t[0])
for t in best[:15]: print('%.4f H=%s q=%d p=%d x=%d a=%.6f delta=%.6f'%t)
