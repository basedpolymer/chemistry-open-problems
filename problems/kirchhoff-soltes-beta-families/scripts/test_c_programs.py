"""test_c_programs.py -- validation of the two C programs against the exact Python library kfexact.py.

  (T1) goodcount: for every connected graph with 3 <= n <= 8 (nauty geng -c), the modular upper bound pg(G) equals
       the exact number of good vertices (computed with exact rationals); in particular pg is never smaller.
  (T2) gadgetbound -d: for every connected graph with 3 <= n <= 7 and every root, the multiset of exact classes
       (equal (rho, delta)) with their efficiencies eff = (g/m) min(1, rho m / R) (0 if rho <= 0) agrees with the
       values recomputed with exact rationals.
Requires the binaries scripts/goodcount and scripts/gadgetbound (gcc -O3 -march=native -o ... ...c) and geng.
Runtime: about 3 min (Fedora 44 / WSL).
"""
import os, sys, subprocess
from collections import Counter
import networkx as nx
from flint import fmpq
from kfexact import kf, rt, delta, minus, good_vertices, Checker

HERE = os.path.dirname(os.path.abspath(__file__))
GB = sys.argv[1] if len(sys.argv) > 1 else 'gadgetbound'     # binary name (default scripts/gadgetbound)
C = Checker()


def geng(n):
    return subprocess.run(['geng', '-c', '-q', str(n)], capture_output=True, text=True).stdout


ok = True
cnt = 0
for n in range(3, 9):
    out = subprocess.run([os.path.join(HERE, 'goodcount'), '0', '1'], input=geng(n), capture_output=True,
                         text=True).stdout
    pgs = {}
    for line in out.split('\n'):
        if line.strip():
            s, n_, pg = line.split()
            pgs[s] = int(pg)                      # graphs that are not printed have pg = 0
    for s in geng(n).split():
        G = nx.from_graph6_bytes(s.encode())
        cnt += 1
        if len(good_vertices(G)) != pgs.get(s, 0):
            ok = False
            print('   mismatch', s, pgs.get(s, 0))
C.check(ok, '(T1) goodcount pg = exact number of good vertices for all %d connected graphs with 3 <= n <= 8' % cnt)


def exact_classes(G, x):
    K, R = kf(G), rt(G, x)
    m = G.number_of_nodes() - 1
    cl = {}
    for v in G.nodes():
        if v == x:
            continue
        d = delta(G, v, K)
        if d is None:
            continue
        rho = rt(minus(G, v), x) - R
        cl.setdefault((rho, d), []).append(v)
    res = []
    for (rho, d), vs in cl.items():
        g = len(vs)
        if rho <= 0:
            e = fmpq(0)
        else:
            e = fmpq(g, m) * min(fmpq(1), rho * m / R)
        res.append((g, e))
    return Counter(res)


ok = True
cnt = 0
for n in range(3, 8):
    out = subprocess.run([os.path.join(HERE, GB), '-d'], input=geng(n), capture_output=True,
                         text=True).stdout
    got = {}
    for line in out.split('\n'):
        if line.startswith('D '):
            s, x, g, en, ed = line.split()[1:6]
            got.setdefault((s, int(x)), Counter())[(int(g), fmpq(int(en), int(ed)))] += 1
    for s in geng(n).split():
        G = nx.from_graph6_bytes(s.encode())
        for x in G.nodes():
            cnt += 1
            if exact_classes(G, x) != got.get((s, x), Counter()):
                ok = False
                print('   mismatch', s, x)
C.check(ok, '(T2) gadgetbound classes and efficiencies = exact recomputation for %d rooted graphs, 3 <= n <= 7' % cnt)
C.finish()
