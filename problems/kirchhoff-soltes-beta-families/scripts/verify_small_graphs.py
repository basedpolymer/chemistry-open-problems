"""verify_small_graphs.py -- Proposition 4.1 (small graphs) from the logs of scripts/run_goodcount.sh.

Input: certificates/search/goodcount_n<N>.log (histogram of the rigorous modular upper bound pg(G) >= #good(G) over
all connected graphs on N vertices, nauty geng -c) and goodcount_n<N>.out (graphs with pg >= 2N/5), N = 3..NMAX.
Claims:
  (G0) the counts of connected graphs agree with OEIS A001349 (2, 6, 21, 112, 853, 11117, 261080, 11716571, ...);
  (G1) for every listed graph, the exact number of good vertices (exact rationals) equals pg; hence the maximum
       number of good vertices among connected graphs on N vertices equals max pg, printed per N;
  (G2) C5 is the only connected graph with 3 <= N <= NMAX having at least 2N/3 good vertices;
  (G3) identification of the extremal graphs (isomorphism with G_s = K_s + K_{s+1,1} and with splices of K4 - e).
Usage: python3 verify_small_graphs.py [NMAX]   (default 10; 11 once the n = 11 slices are complete).
Runtime: a few seconds (Fedora 44 / WSL, python-flint 0.6, networkx 3).
"""
import sys, os, re
import networkx as nx
from fractions import Fraction
from kfexact import good_vertices, splice, Checker

NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 10
HERE = os.path.dirname(os.path.abspath(__file__))
SEARCH = os.path.join(HERE, '..', 'certificates', 'search')
A001349 = {3: 2, 4: 6, 5: 21, 6: 112, 7: 853, 8: 11117, 9: 261080, 10: 11716571, 11: 1006700565}
C = Checker()


def read_hist(nn):
    files = [os.path.join(SEARCH, 'goodcount_n%d.log' % nn)]
    if not os.path.exists(files[0]):
        files = [os.path.join(SEARCH, 'goodcount_n%d_part%d.log' % (nn, r)) for r in range(4)]
    hist = {}
    for f in files:
        for line in open(f):
            m = re.match(r'n=(\d+) graphs=(\d+) pg-histogram:(.*)', line.strip())
            if m and int(m.group(1)) == nn:
                for tok in m.group(3).split():
                    g, c = tok.split(':')
                    hist[int(g)] = hist.get(int(g), 0) + int(c)
    return hist


def read_out(nn):
    files = [os.path.join(SEARCH, 'goodcount_n%d.out' % nn)]
    if not os.path.exists(files[0]):
        files = [os.path.join(SEARCH, 'goodcount_n%d_part%d.out' % (nn, r)) for r in range(4)]
    rows = []
    for f in files:
        for line in open(f):
            s, n_, pg = line.split()
            rows.append((s, int(n_), int(pg)))
    return rows


D = nx.Graph([('x', 'a'), ('x', 'b'), ('a', 'b'), ('a', 'y'), ('b', 'y')])


def G_s(s):
    """K_s + K_{s+1,1} = clique K_{s+1} joined to an independent set of size s+1 (complete split graph)."""
    G = nx.complete_graph(range(s + 1))
    G.add_edges_from((i, j) for i in range(s + 1) for j in range(s + 1, 2 * s + 2))
    return G


known = {'C5': nx.cycle_graph(5)}
for s in range(1, 6):
    known['G_%d = K_%d + K_{%d,1}' % (s, s, s + 1)] = G_s(s)
for k in range(2, 4):
    known['%d copies of K4-e spliced at a degree-2 vertex' % k] = splice([(D, 'x')] * k)

maxima = {}
ok0 = ok1 = True
two_thirds = []
for nn in range(3, NMAX + 1):
    hist = read_hist(nn)
    tot = sum(hist.values())
    if tot != A001349[nn]:
        ok0 = False
    maxima[nn] = max(g for g, c in hist.items() if c)
    for g, c in hist.items():
        if c and 3 * g >= 2 * nn:
            two_thirds.append((nn, g, c))
    for s, n_, pg in read_out(nn):
        G = nx.from_graph6_bytes(s.encode())
        ng = len(good_vertices(G))
        name = [k for k, H in known.items() if H.number_of_nodes() == n_ and nx.is_isomorphic(G, H)]
        print('   n=%d %-12s pg=%d exact good=%d  %s' % (n_, s, pg, ng, name[0] if name else ''))
        if ng != pg:
            ok1 = False
C.check(ok0, '(G0) numbers of connected graphs agree with OEIS A001349 for 3 <= n <= %d' % NMAX)
C.check(ok1, '(G1) every listed graph has exactly pg good vertices (exact rationals)')
for nn in range(3, NMAX + 1):
    print('   n=%d: maximum number of good vertices = %d' % (nn, maxima[nn]))
C.check(two_thirds == [(5, 5, 1)], '(G2) pg >= 2n/3 occurs only for n = 5 (one graph, C5): %s' % two_thirds)
C.finish()
