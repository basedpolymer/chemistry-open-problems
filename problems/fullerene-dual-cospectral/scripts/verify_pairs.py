#!/usr/bin/env python3
"""verify_pairs.py -- independent re-check, from the certificate files alone, of every claimed cospectral pair.

Uses only the graph6 strings stored in certificates/*.json (no buckygen needed for V1-V7; V8 optionally
re-runs buckygen for the minimality claim n <= 42).  Every check is exact (integers / python-flint).

Claims (PASS/FAIL; exit code 1 on any failure):
  V1  every graph6 string is the dual T_n of a fullerene C_n: simple planar graph with m = n/2+2 vertices and
      3m-6 edges (hence a triangulation of the sphere), all degrees 5 or 6, exactly twelve of degree 5; the
      rotation system of its (unique up to mirror) planar embedding has only triangular faces, Euler char. 2;
  V2  inside each class the integer characteristic polynomials of the stated matrix (A(T), A+D, 2A+D or A(T6))
      are identical, and equal to the polynomial recorded in the certificate;
  V3  inside each class the duals are pairwise NON-isomorphic: their canonical face spirals differ (the canonical
      spiral = lexicographically least pentagon-position vector over all starts and both orientations is an
      isomorphism invariant), and networkx's exact VF2 test also reports non-isomorphic;
  V4  the canonical spirals recorded in the certificate are recomputed exactly;
  V5  (T6 classes) the T6 graphs recorded are the induced subgraphs on the degree-6 vertices, and the recorded
      isomorphic / non-isomorphic status of the T6 graphs is confirmed with VF2;
  V6  the counts table certificates/counts_per_n.csv agrees with the class files, the isomer counts agree with
      OEIS A007894, and every n in 60..NMAX has at least one A(T)-cospectral class (VOR Conjecture 1(b) on the
      scanned range), while no n >= 54 has a T6-cospectral class (VOR Conjecture 1(c) on the scanned range);
  V7  literal VOR Conjecture 1(a): T6 of C34:2 is C6 + K1, cospectral with the subdivided claw S(2,2,2);
  V8  (optional, needs buckygen: --exhaustive-small) exhaustive EXACT re-check for 20 <= n <= 42: no two
      non-isomorphic duals are A(T)-cospectral for n < 42, exactly one class (2 isomers) for n = 42.

Usage: python3 -I verify_pairs.py [--exhaustive-small]
Environment: python3.13, python-flint 0.9, networkx 3.x.  Runtime: ~2-10 min (most of it in V3/V4 spirals).
"""
import csv
import json
import os
import subprocess
import sys
import time

import flint
import networkx as nx

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
from planarcode import from_graph6, read_planar_code, check_triangulation_dual_of_fullerene  # noqa: E402
from spiral import canonical_spiral  # noqa: E402

CERT = os.path.join(ROOT, 'certificates')
FAILED = []


def report(name, ok, detail=''):
    print('%s  %s%s' % ('PASS' if ok else 'FAIL', name, ('  ' + detail) if detail else ''), flush=True)
    if not ok:
        FAILED.append(name)


def adj_from_graph6_embedded(g6):
    """adjacency lists in clockwise rotation order of the planar embedding (networkx), from graph6."""
    a = from_graph6(g6)
    G = nx.Graph()
    G.add_nodes_from(range(len(a)))
    G.add_edges_from((v, w) for v in range(len(a)) for w in a[v])
    ok, emb = nx.check_planarity(G)
    if not ok:
        return None
    return [list(emb.neighbors_cw_order(v)) for v in range(len(a))]


def is_fullerene_dual(g6, n):
    adj = adj_from_graph6_embedded(g6)
    if adj is None:
        return False, None
    m = len(adj)
    E = sum(len(x) for x in adj) // 2
    ok = (m == n // 2 + 2) and (E == 3 * m - 6) and not check_triangulation_dual_of_fullerene(adj)
    return ok, adj


def mat(adj, kind):
    m = len(adj)
    deg = [len(x) for x in adj]
    if kind == 'A(T6)':
        verts = [v for v in range(m) if deg[v] == 6]
    else:
        verts = list(range(m))
    loc = {v: i for i, v in enumerate(verts)}
    k = len(verts)
    M = [[0] * k for _ in range(k)]
    for v in verts:
        for w in adj[v]:
            if w in loc:
                M[loc[v]][loc[w]] = 2 if kind == '2A+D(T)' else 1
        if kind in ('A+D(T)', '2A+D(T)'):
            M[loc[v]][loc[v]] = deg[v]
    return M, verts


def charpoly(M):
    return [str(int(c)) for c in flint.fmpz_mat(M).charpoly().coeffs()] if M else ['1']


def nxgraph(adj, verts=None):
    G = nx.Graph()
    vs = range(len(adj)) if verts is None else verts
    S = set(vs)
    G.add_nodes_from(vs)
    G.add_edges_from((v, w) for v in vs for w in adj[v] if w in S and v < w)
    return G


def check_class_file(fname, expect_kind):
    path = os.path.join(CERT, fname)
    data = json.load(open(path))
    ok1 = ok2 = ok3 = ok4 = ok5 = True
    bad = []
    for c in data:
        n = c['n']
        kind = c['matrix']
        if kind != expect_kind:
            ok2 = False
        adjs = []
        for g6 in c['graph6_dual']:
            ok, adj = is_fullerene_dual(g6, n)
            ok1 &= ok
            adjs.append(adj)
        if not all(adjs):
            bad.append(n)
            continue
        cps = [charpoly(mat(a, kind)[0]) for a in adjs]
        if not (all(cp == cps[0] for cp in cps) and cps[0] == c['charpoly_coeffs_low_to_high']):
            ok2 = False
            bad.append(('cp', n))
        sps = [list(canonical_spiral(a)) for a in adjs]
        if len({tuple(s) for s in sps}) != len(sps):
            ok3 = False
            bad.append(('spiral-equal', n))
        Gs = [nxgraph(a) for a in adjs]
        for i in range(len(Gs)):
            for j in range(i + 1, len(Gs)):
                if nx.is_isomorphic(Gs[i], Gs[j]):
                    ok3 = False
                    bad.append(('iso', n))
        if sps != c['canonical_spirals']:
            ok4 = False
            bad.append(('spiral-record', n))
        if kind == 'A(T6)':
            subs = [nxgraph(a, mat(a, kind)[1]) for a in adjs]
            rec = [nxgraph(from_graph6(s)) for s in c['T6_graph6']]
            same = all(nx.is_isomorphic(x, y) for x, y in zip(subs, rec))
            alliso = all(nx.is_isomorphic(subs[0], s) for s in subs[1:])
            ok5 &= same and (alliso == c['T6_graphs_all_isomorphic'])
    tag = '%s (%d classes)' % (fname, len(data))
    report('V1 %s: all graphs are fullerene duals' % tag, ok1)
    report('V2 %s: equal exact characteristic polynomials of %s' % (tag, expect_kind), ok2)
    report('V3 %s: members pairwise non-isomorphic (spirals + VF2)' % tag, ok3)
    report('V4 %s: recorded canonical spirals recomputed' % tag, ok4)
    if expect_kind == 'A(T6)':
        report('V5 %s: T6 subgraphs and their isomorphism status confirmed' % tag, ok5)
    if bad:
        print('   problems:', bad[:20])
    return data


def main():
    t0 = time.time()
    T = check_class_file('cospectral_pairs_T.json', 'A(T)')
    Q = check_class_file('cospectral_pairs_signless_laplacian.json', 'A+D(T)')
    S = check_class_file('cospectral_pairs_2A_plus_D.json', '2A+D(T)') if \
        json.load(open(os.path.join(CERT, 'cospectral_pairs_2A_plus_D.json'))) else []
    H = check_class_file('T6_cospectral_classes.json', 'A(T6)')
    # smallest A(T) counterexample in the files
    nmin = min(c['n'] for c in T)
    c42 = [c for c in T if c['n'] == nmin]
    report('V2b smallest n with an A(T)-cospectral class in the certificate is 42 (one class of size 2)',
           nmin == 42 and len(c42) == 1 and len(c42[0]['graph6_dual']) == 2,
           'atlas numbers %s' % (c42[0].get('atlas_numbers'),))
    q44 = [c for c in Q if c['n'] == 44]
    report('V2c signless-Laplacian cospectral pair at n = 44 (Atlas 72, 80)',
           len(q44) == 1 and sorted(q44[0].get('atlas_numbers') or []) == [72, 80])
    # V6 counts table
    rows = list(csv.DictReader(open(os.path.join(CERT, 'counts_per_n.csv'))))
    ok = True
    for r in rows:
        n = int(r['n'])
        ok &= r['counts_match'] == 'True' and r['isomers_buckygen'] == r['isomers_OEIS_A007894']
        ok &= int(r['T_cospectral_classes']) == sum(1 for c in T if c['n'] == n)
        ok &= int(r['T6_cospectral_classes']) == sum(1 for c in H if c['n'] == n)
        ok &= int(r['signless_A_plus_D_classes'] or 0) == sum(1 for c in Q if c['n'] == n)
    ns = sorted(int(r['n']) for r in rows)
    nmax = max(ns)
    expected = [20] + list(range(24, nmax + 1, 2))
    report('V6a counts table consistent with class files and with OEIS A007894 for every scanned n (20..%d)' % nmax,
           ok and ns == expected)
    byn = {int(r['n']): r for r in rows}
    report('V6b every n in 60..%d has >= 1 A(T)-cospectral class (Conjecture 1(b) of the VOR, finite range)' % nmax,
           all(int(byn[n]['T_cospectral_classes']) >= 1 for n in range(60, nmax + 1, 2)))
    report('V6c no n in 54..%d has a T6-cospectral pair of isomers (Conjecture 1(c) of the VOR, finite range)' % nmax,
           all(int(byn[n]['T6_cospectral_classes']) == 0 for n in range(54, nmax + 1, 2)))
    report('V6d every T6-cospectral class found (n <= %d) consists of isomers with ISOMORPHIC T6 graphs' % nmax,
           all(c['T6_graphs_all_isomorphic'] for c in H))
    # V7 literal Conjecture 1(a): C6 + K1 versus subdivided claw
    c6k1 = nx.cycle_graph(6)
    c6k1.add_node(6)
    claw = nx.Graph([(0, 1), (1, 2), (0, 3), (3, 4), (0, 5), (5, 6)])
    A1 = nx.to_numpy_array(c6k1, nodelist=range(7)).astype(int).tolist()
    A2 = nx.to_numpy_array(claw, nodelist=range(7)).astype(int).tolist()
    report('V7 C6+K1 and the subdivided claw S(2,2,2) are cospectral and non-isomorphic',
           charpoly(A1) == charpoly(A2) and sorted(d for _, d in c6k1.degree()) != sorted(d for _, d in claw.degree()))
    t6file = os.path.join(CERT, 't6_ds_all_graphs_n24-38.json')
    if os.path.exists(t6file):
        rep = json.load(open(t6file))
        e = [x for x in rep if x['n'] == 34 and x['atlas'] == 2][0]
        adj = from_graph6(e['T6_graph6'])
        G = nxgraph(adj)
        first = min(x['n'] for x in rep if x['cospectral_mates_graph6'])
        report('V7b T6(C34:2) is C6+K1 (VF2) and n = 34 is the first n with a T6 having a cospectral mate among all '
               'graphs (geng-exhaustive record for n <= 38)', nx.is_isomorphic(G, c6k1) and first == 34)
    if '--exhaustive-small' in sys.argv:
        from scan_driver import BUCKYGEN
        okx = True
        cnt42 = None
        for n in [20] + list(range(24, 44, 2)):
            data = subprocess.run([BUCKYGEN, '-q', '%dd' % n], capture_output=True, check=True).stdout
            seen = {}
            for g in read_planar_code(data):
                cp = tuple(charpoly(mat(g, 'A(T)')[0]))
                seen.setdefault(cp, []).append(g)
            multi = [v for v in seen.values() if len(v) > 1]
            if n < 42:
                okx &= not multi
            else:
                cnt42 = [len(v) for v in multi]
        report('V8 exhaustive exact check: no A(T)-cospectral isomers for n < 42; one pair at n = 42', okx and cnt42 == [2])
    print('runtime %.1f s' % (time.time() - t0))
    if FAILED:
        print('FAILED:', FAILED)
        sys.exit(1)
    print('ALL PASS')


if __name__ == '__main__':
    main()
