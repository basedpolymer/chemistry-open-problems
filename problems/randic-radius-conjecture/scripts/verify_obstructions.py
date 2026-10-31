"""Certify obstructions to attempted lemmas, NOT to the Randić-radius conjecture.

Environment: WSL scientific Python, NetworkX; standard-library rational intervals.
Runtime: normally below 5 seconds, one core. No floating-point comparisons.
Run from repository root: python problems/randic-radius-conjecture/scripts/verify_obstructions.py
"""
from collections import Counter
from fractions import Fraction as F
from math import isqrt
import networkx as nx

SCALE = 10**30

def invsqrt(q):
    """Rigorous rational enclosure, obtained using an integer square root."""
    k=isqrt(q*SCALE*SCALE)
    if k*k==q*SCALE*SCALE:
        return F(SCALE,k), F(SCALE,k)
    return F(SCALE,k+1), F(SCALE,k)

def randic(g):
    counts=Counter(g.degree(u)*g.degree(v) for u,v in g.edges())
    lo=hi=F(0)
    for q,m in counts.items():
        a,b=invsqrt(q)
        lo+=m*a;hi+=m*b
    return lo,hi

def incident(g,v):
    vals=[invsqrt(g.degree(v)*g.degree(u)) for u in g[v]]
    return sum((x[0] for x in vals),F(0)),sum((x[1] for x in vals),F(0))

def quotient(g,v):
    """Simple quotient obtained by identifying the closed neighborhood of v."""
    s={v,*g[v]}
    h=nx.Graph()
    h.add_nodes_from(set(g)-s|{v})
    for a,b in g.edges():
        a=v if a in s else a
        b=v if b in s else b
        if a!=b:h.add_edge(a,b)
    return h

def check(label,condition):
    print(('PASS' if condition else 'FAIL')+': '+label,flush=True)
    if not condition:raise SystemExit(1)

def unique_maximum_example():
    g=nx.Graph()
    g.add_edges_from((0,i) for i in range(1,13))
    for group in range(2):
        for j in range(10):
            b=13+10*group+j
            g.add_edges_from((a,b) for a in range(1+6*group,7+6*group))
            g.add_edges_from((b,33+5*(b-13)+t) for t in range(5))
    return g

def main():
    # All vertices of K_(4,1224) are forced into any connected dominating set:
    # each has a private pendant neighbor. The core itself is connected and dominates.
    g=nx.complete_bipartite_graph(4,1224)
    core=set(g)
    nc=len(core)
    g.add_edges_from((v,nc+v) for v in core)
    check('corona core is a connected dominating set of size 1228',
          nc==1228 and nx.is_connected(g.subgraph(core)) and
          all(list(g[nc+v])==[v] for v in core))
    lo,hi=randic(g)
    check('corona has 2R < 1228 (so the proposed internal-vertex bound fails)',2*hi<1228)
    # Four automorphism orbits, justified by permutations inside the bipartition classes.
    check('corona orbit eccentricities are 3,3,4,4',
          [nx.eccentricity(g,v) for v in [0,4,nc,nc+4]]==[3,3,4,4])
    check('corona satisfies the target weak inequality',lo>2)

    g=unique_maximum_example()
    h=quotient(g,0)
    check('133-vertex construction has a unique maximum-degree vertex 0',
          len(g)==133 and g.number_of_edges()==232 and g.degree(0)==12 and
          max(g.degree(v) for v in g if v!=0)==11)
    check('contracting its unique maximum-degree neighborhood INCREASES R',randic(g)[1]<randic(h)[0])
    check('construction and quotient have radii 3 and 2',nx.radius(g)==3 and nx.radius(h)==2)
    check('133-vertex construction satisfies the target weak inequality',randic(g)[0]>2)

    g=nx.from_graph6_bytes(b'EJy?')
    weights={v:incident(g,v) for v in g}
    # Identify the unique maximizer by its neighbor degrees, then compare intervals.
    maxima=[v for v in g if sorted(g.degree(u) for u in g[v])==[1,3]]
    check('the unique incident-weight maximizer in EJy? is vertex 0',
          maxima==[0] and all(g.degree(v)==2 for v in maxima) and
          min(weights[v][0] for v in maxima)>max(weights[v][1] for v in g if v not in maxima))
    check('every maximum-incident-weight choice loses less than one unit of R',
          all(randic(g)[1]-randic(quotient(g,v))[0]<1 for v in maxima))
    print('EJy? edges:',sorted(g.edges()))

    g=nx.from_graph6_bytes(b'E^mG')
    check('E^mG is 2-connected, with radius 2 and diameter 3',
          nx.is_biconnected(g) and nx.radius(g)==2 and nx.diameter(g)==3)
    check('E^mG disproves the auxiliary inequality R >= diameter',randic(g)[1]<3)
    check('E^mG satisfies the target weak inequality',randic(g)[0]>1)
    print('E^mG edges:',sorted(g.edges()))

    # Check path exceptions directly; this finite check corroborates the formula in PROGRESS.
    check('P2 has R = radius = 1',randic(nx.path_graph(2))==(F(1),F(1)))
    check('even paths P4 through P20 violate R >= radius but satisfy R >= radius-1',
          all((lambda z,r:z[1]<r and z[0]>r-1)(randic(nx.path_graph(n)),n//2) for n in range(4,21,2)))
    print('PASS: all listed obstruction certificates; no target-conjecture counterexample claimed.')

if __name__=='__main__':main()
