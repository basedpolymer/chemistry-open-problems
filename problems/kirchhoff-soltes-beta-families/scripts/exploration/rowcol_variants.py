"""FLOATING-POINT EXPLORATION: integer-grid search (a, b <= A, B) for exactly good orbits in row x column graphs with S_a x S_b symmetry (row/column cliques, row/column hubs, hub cliques, super hubs). Usage: python3 rowcol_variants.py 9 9. Runtime about 2 min.
Environment: Fedora 44 (WSL), Python 3.12, numpy, numba, networkx, nauty geng."""
# integer-grid search for exact good classes in row x column structures with S_a x S_b symmetry (float screen)
import sys, itertools, numpy as np, networkx as nx
def kf_np(G):
    L=nx.laplacian_matrix(G).toarray().astype(float); ev=np.linalg.eigvalsh(L)
    if ev[1] < 1e-9: return None
    return G.number_of_nodes()*np.sum(1/ev[1:])
def build(a,b,rowE,colE,rh,ch,rhE,chE,zr,zc):
    G=nx.Graph()
    V=[(i,j) for i in range(a) for j in range(b)]
    G.add_nodes_from(V)
    if rowE:
        for i in range(a):
            for j in range(b):
                for jj in range(j+1,b): G.add_edge((i,j),(i,jj))
    if colE:
        for j in range(b):
            for i in range(a):
                for ii in range(i+1,a): G.add_edge((i,j),(ii,j))
    for i in range(a):
        for t in range(rh):
            for j in range(b): G.add_edge(('r',i,t),(i,j))
    for j in range(b):
        for t in range(ch):
            for i in range(a): G.add_edge(('c',j,t),(i,j))
    if rhE and rh:
        hs=[('r',i,t) for i in range(a) for t in range(rh)]
        for x,y in itertools.combinations(hs,2): G.add_edge(x,y)
    if chE and ch:
        hs=[('c',j,t) for j in range(b) for t in range(ch)]
        for x,y in itertools.combinations(hs,2): G.add_edge(x,y)
    if zr and rh:
        for i in range(a):
            for t in range(rh): G.add_edge('zr',('r',i,t))
    if zc and ch:
        for j in range(b):
            for t in range(ch): G.add_edge('zc',('c',j,t))
    return G
def classes(G,a,b,rh,ch,zr,zc):
    reps={'rook':(0,0)}
    if rh: reps['rowhub']=('r',0,0)
    if ch: reps['colhub']=('c',0,0)
    if zr and rh: reps['zr']='zr'
    if zc and ch: reps['zc']='zc'
    sizes={'rook':a*b,'rowhub':a*rh,'colhub':b*ch,'zr':1,'zc':1}
    return reps,sizes
A=int(sys.argv[1]); B=int(sys.argv[2])
opts=[]
for rowE,colE in [(1,1),(1,0),(0,1),(0,0)]:
    for rh in [0,1,2]:
        for ch in [0,1,2]:
            for rhE in ([0,1] if rh else [0]):
                for chE in ([0,1] if ch else [0]):
                    for zr in ([0,1] if rh else [0]):
                        for zc in ([0,1] if ch else [0]):
                            opts.append((rowE,colE,rh,ch,rhE,chE,zr,zc))
best=[]
for o in opts:
    rowE,colE,rh,ch,rhE,chE,zr,zc=o
    for a in range(2,A+1):
        for b in range(2,B+1):
            if not rowE and not colE and rh==0 and ch==0: continue
            G=build(a,b,*o)
            if not nx.is_connected(G): continue
            K=kf_np(G)
            if K is None: continue
            reps,sizes=classes(G,a,b,rh,ch,zr,zc)
            good=0; gl=[]
            for name,v in reps.items():
                H=G.copy(); H.remove_node(v)
                if not nx.is_connected(H): continue
                K2=kf_np(H)
                if K2 is not None and abs(K2-K)<1e-8*(1+K):
                    good+=sizes[name]; gl.append(name)
            n=G.number_of_nodes()
            if good>0:
                best.append((good/n,o,a,b,gl,good,n))
                print('%.4f'%(good/n),o,a,b,gl,good,n,flush=True)
