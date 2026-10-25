"""FLOATING-POINT EXPLORATION: for vertex-transitive H, real roots r* of Delta_H(r) for the joins K_r + H and Kbar_r + H (continuous r), and the proportion h/(h+r*). Runtime about 1 min.
Environment: Fedora 44 (WSL), Python 3.12, numpy, numba, networkx, nauty geng."""
import numpy as np, networkx as nx
from scipy.optimize import brentq
def S(L, r):  # sum over nonprincipal Laplacian eigenvalues of 1/(r+mu) (all eigenvalues minus one zero)
    ev=np.linalg.eigvalsh(L); ev=np.sort(ev)[1:]; return np.sum(1/(r+ev))
def deltaH(H, r, clique=True):
    # Delta for v in H, G = K_r + H (clique) or Kbar_r + H
    h=H.number_of_nodes(); L=nx.laplacian_matrix(H).toarray().astype(float)
    v=list(H)[0]; Hv=H.copy(); Hv.remove_node(v); Lv=nx.laplacian_matrix(Hv).toarray().astype(float)
    n=r+h
    if clique:
        kf = r + n*S(L, r); kfv = r + (n-1)*S(Lv, r)
    else:
        # Kbar_r + H: eigenvalues 0, n, h (r-1 times), r + mu_j
        kf = 1 + n*((r-1)/h + S(L, r)); kfv = 1 + (n-1)*((r-1)/(h-1) + S(Lv, r))
    return kf - kfv
fams = {'C%d'%m: nx.cycle_graph(m) for m in range(5,16)}
fams.update({'K%dx%d'%(a,b): nx.cartesian_product(nx.complete_graph(a),nx.complete_graph(b)) for a in range(2,6) for b in range(a,7)})
fams.update({'Pet': nx.petersen_graph(), 'Q3': nx.hypercube_graph(3), 'Q4': nx.hypercube_graph(4), 'K33': nx.complete_bipartite_graph(3,3),
             'K44': nx.complete_bipartite_graph(4,4), 'Prism6': nx.circular_ladder_graph(6), 'Prism8': nx.circular_ladder_graph(8),
             'Dodec': nx.dodecahedral_graph(), 'Heawood': nx.heawood_graph(), '3K3': nx.disjoint_union_all([nx.complete_graph(3)]*3),
             '4K2': nx.disjoint_union_all([nx.complete_graph(2)]*4), 'C5xC5': nx.cartesian_product(nx.cycle_graph(5),nx.cycle_graph(5))})
for name,H in fams.items():
    h=H.number_of_nodes()
    for cl in (True, False):
        rs=np.linspace(0.05, 6*h, 600); vals=[deltaH(H,r,cl) for r in rs]
        roots=[]
        for i in range(len(rs)-1):
            if np.sign(vals[i])!=np.sign(vals[i+1]):
                roots.append(brentq(lambda r: deltaH(H,r,cl), rs[i], rs[i+1]))
        print('%-8s h=%3d %s roots r*=%s  fraction h/(h+r*)=%s' % (name, h, 'K_r+H' if cl else 'Kbar_r+H', ['%.3f'%x for x in roots], ['%.3f'%(h/(h+x)) for x in roots]))
