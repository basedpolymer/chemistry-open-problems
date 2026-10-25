"""FLOATING-POINT EXPLORATION: sign tables of Delta for rook graphs with row hubs and for cycles with hubs. Runtime about 1 min.
Environment: Fedora 44 (WSL), Python 3.12, numpy, numba, networkx, nauty geng."""
import numpy as np, networkx as nx
def kf_np(G):
    L=nx.laplacian_matrix(G).toarray().astype(float); ev=np.linalg.eigvalsh(L)
    return G.number_of_nodes()*np.sum(1/ev[1:])
def dl(G,v):
    H=G.copy(); H.remove_node(v); 
    if not nx.is_connected(H): return float('nan')
    return kf_np(G)-kf_np(H)
def rook_rowhubs(a,b,hubclique=False):
    G=nx.cartesian_product(nx.complete_graph(a),nx.complete_graph(b))
    for i in range(a):
        for j in range(b): G.add_edge(('h',i),(i,j))
    if hubclique:
        for i in range(a):
            for k in range(i+1,a): G.add_edge(('h',i),('h',k))
    return G
def cycle_hubs(m,t,hubclique=False):
    G=nx.cycle_graph(m)
    for i in range(m): G.add_edge(('h',i%t),i)
    if hubclique:
        for i in range(t):
            for k in range(i+1,t): G.add_edge(('h',i),('h',k))
    return G
print('rook K_a x K_b + row hubs: Delta(rook vertex) [hub vertex]')
for a in range(2,8):
    print(a, ' '.join('%7.2f'%dl(rook_rowhubs(a,b),(0,0)) for b in range(2,10)))
print('same with hubs forming a clique')
for a in range(2,8):
    print(a, ' '.join('%7.2f'%dl(rook_rowhubs(a,b,True),(0,0)) for b in range(2,10)))
print('cycle C_m with t hubs (i mod t): Delta(rim), m = t*s')
for t in range(1,7):
    print(t, ' '.join('%8.2f'%dl(cycle_hubs(t*s,t),0) for s in range(2,12)))
print('cycle hubs with hub clique')
for t in range(2,7):
    print(t, ' '.join('%8.2f'%dl(cycle_hubs(t*s,t,True),0) for s in range(2,12)))
