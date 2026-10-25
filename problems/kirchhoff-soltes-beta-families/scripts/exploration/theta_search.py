"""FLOATING-POINT EXPLORATION: theta graphs with e threads of length 1, p of length 2, t of length 3, q of length 4; prints those with a good class. Runtime about 1 min.
Environment: Fedora 44 (WSL), Python 3.12, numpy, numba, networkx, nauty geng."""
import numpy as np, networkx as nx, itertools
def kf_np(G):
    L=nx.laplacian_matrix(G).toarray().astype(float); ev=np.linalg.eigvalsh(L)
    return G.number_of_nodes()*np.sum(1/ev[1:])
def theta(lens):
    G=nx.Graph(); G.add_nodes_from(['u','w']); reps={}
    for i,l in enumerate(lens):
        path=['u']+[(i,k) for k in range(1,l)]+['w']
        nx.add_path(G,path)
        if l>=2: reps.setdefault(l,[]).append([(i,k) for k in range(1,l)])
    return G,reps
res=[]
for e in [0,1]:
  for p in range(0,25):
    for t in range(1,30):
      for q in range(0,6):
        lens=[1]*e+[2]*p+[3]*t+[4]*q
        if len(lens)<2: continue
        G,reps=theta(lens); K=kf_np(G); n=G.number_of_nodes()
        good=0; names=[]
        for l,threads in reps.items():
            th=threads[0]
            for k in range(1,(l+1)//2+1):  # by u<->w symmetry positions k and l-k equivalent
                v=th[k-1]; H=G.copy(); H.remove_node(v)
                if abs(kf_np(H)-K)<1e-9*(1+K):
                    cnt=len(threads)*(2 if 2*k!=l else 1); good+=cnt; names.append((l,k))
        if good: res.append((good/n,e,p,t,q,names,good,n)); print('%.4f'%(good/n),e,p,t,q,names,good,n,flush=True)
