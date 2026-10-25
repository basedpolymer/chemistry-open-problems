"""FLOATING-POINT EXPLORATION: sign tables of Delta(rook vertex) for rows K_b joined through column hubs (and variants). Runtime about 2 min.
Environment: Fedora 44 (WSL), Python 3.12, numpy, numba, networkx, nauty geng."""
import numpy as np, networkx as nx
from hub_tables import dl
def rows_colhubs(a,b,rowhub=False,colclique=False):
    G=nx.Graph()
    for i in range(a):
        for j in range(b):
            for jj in range(j+1,b): G.add_edge((i,j),(i,jj))
            G.add_edge((i,j),('c',j))
            if rowhub: G.add_edge((i,j),('h',i))
    if colclique:
        for j in range(b):
            for jj in range(j+1,b): G.add_edge(('c',j),('c',jj))
    return G
for name,kw in [('rows K_b + column hubs',{}),('+ row hubs',{'rowhub':True}),('column hubs form clique',{'colclique':True})]:
    print(name)
    print('a\b '+' '.join('%7d'%b for b in range(2,10)))
    for a in range(2,13):
        print('%2d  '%a + ' '.join('%7.3f'%dl(rows_colhubs(a,b,**kw),(0,0)) for b in range(2,10)))
