#!/usr/bin/env python3
"""tube_family.py -- the infinite families of cospectral fullerene duals obtained by rotating the C3v cap I of
the thinnest zigzag tube by one step (60 degrees).

Construction of G(Z, k, delta) (all adjacency, no embedding needed; validity is re-checked):
  * Z: a cap of the 6-row tube, given as patch graph6 with the cap vertices first and then its boundary row
    c_0..c_5 in cyclic order (certificates/caps_6row_tube_n52.json, produced by find_caps.py); this row is row 1;
  * rows 2..k: row r+1 vertex j is adjacent to row r vertices j-1, j, and to row r+1 vertices j +- 1;
  * the cap I(delta) on top of row k (= S): a row O_0..O_5 with O_j ~ S_{j-1}, S_j, O_j ~ O_{j+-1}; O_j is a
    pentagon iff j = delta (mod 2); a triangle t_0 t_1 t_2 with t_i ~ O_{p-1}, O_p, O_{p+1}, p = 2i + delta.
Then I(0) and I(1) differ by a rotation of the cap by one step.

Claims checked here (PASS/FAIL, exit 1 on failure), for every cap Z in the file and every 1 <= k <= KMAX:
  F1  G(Z,k,0), G(Z,k,1) are fullerene duals (planar, 3m-6 edges, degrees 5/6, twelve 5s);
  F2  STRONG CAP CONDITION for I (independent of Z, k): det(xI - A_I) is the same for both attachments and
      W_delta(x) = B_delta^T adj(xI - A_I) B_delta (6 x 6 over Z[x]) is the same for delta = 0, 1, exactly
      (python-flint fmpz_poly arithmetic); hence by the Schur complement lemma
      det(xI - A(G(Z,k,0))) = det(xI - A(G(Z,k,1))) for EVERY rest -- re-checked directly here (exact charpoly);
  F3  the isomorphism status of G(Z,k,0) vs G(Z,k,1) (nauty canonical forms) is recorded; the theorem in
      result.md predicts: non-isomorphic iff the cap Z has no automorphism acting on the row as an odd rotation
      and k has the parity singled out by Z's reflections (see result.md, Prop. 3); F3 checks this prediction.
Usage: python3 tube_family.py [KMAX]   (default 8)
Environment: python3.13, python-flint 0.9, networkx, nauty labelg. Runtime ~1-2 min.
"""
import json
import os
import subprocess
import sys

import flint
import networkx as nx

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
from planarcode import from_graph6, to_graph6  # noqa: E402
from scan_driver import LABELG  # noqa: E402
from verify_pairs import adj_from_graph6_embedded  # noqa: E402
from planarcode import check_triangulation_dual_of_fullerene  # noqa: E402

CERT = os.path.join(ROOT, 'certificates')
FAILED = []


def report(name, ok, detail=''):
    print('%s  %s%s' % ('PASS' if ok else 'FAIL', name, ('  ' + detail) if detail else ''), flush=True)
    if not ok:
        FAILED.append(name)


def cap_I(delta):
    """adjacency of I(delta) glued to an abstract row S (labels 0..5); I's vertices 6..14: O_j = 6+j, t_i = 12+i."""
    E = set()
    for j in range(6):
        E.add((6 + j, (j - 1) % 6))
        E.add((6 + j, j))
        E.add((6 + j, 6 + (j + 1) % 6))
    for i in range(3):
        p = 2 * i + delta
        for q in (p - 1, p, p + 1):
            E.add((12 + i, 6 + q % 6))
        E.add((12 + i, 12 + (i + 1) % 3))
    return E


def build(capg6, k, delta):
    h = from_graph6(capg6)
    nz = len(h) - 6
    adj = [set(a) for a in h]
    row = list(range(nz, nz + 6))
    # sanity: row is a 6-cycle in the given cyclic order
    assert all(row[(j + 1) % 6] in adj[row[j]] for j in range(6))
    for r in range(2, k + 1):
        new = [len(adj) + j for j in range(6)]
        adj.extend(set() for _ in range(6))
        for j in range(6):
            for w in (row[(j - 1) % 6], row[j], new[(j + 1) % 6], new[(j - 1) % 6]):
                adj[new[j]].add(w)
                adj[w].add(new[j])
        row = new
    base = len(adj)
    adj.extend(set() for _ in range(9))
    for (a, b) in cap_I(delta):
        va = row[a] if a < 6 else base + a - 6
        vb = row[b] if b < 6 else base + b - 6
        adj[va].add(vb)
        adj[vb].add(va)
    return [sorted(a) for a in adj]


def is_dual(adj):
    g6 = to_graph6(adj) if len(adj) <= 62 else None
    G = nx.Graph()
    G.add_nodes_from(range(len(adj)))
    G.add_edges_from((v, w) for v in range(len(adj)) for w in adj[v] if v < w)
    ok, emb = nx.check_planarity(G)
    if not ok:
        return False
    rot = [list(emb.neighbors_cw_order(v)) for v in range(len(adj))]
    m = len(adj)
    return G.number_of_edges() == 3 * m - 6 and not check_triangulation_dual_of_fullerene(rot)


def charpoly(adj):
    m = len(adj)
    return flint.fmpz_mat(m, m, [1 if w in adj[v] else 0 for v in range(m) for w in range(m)]).charpoly()


def strong_condition_I():
    """exact polynomial identity: p_I(x) = det(xI - A_I) is the same for delta = 0, 1 and
    W_delta(x) = B_delta^T adj(xI - A_I) B_delta agree.  Every entry of W is a polynomial of degree <= 8 (adjugate
    of a 9 x 9 matrix), p_I has degree 9; so agreement at 12 distinct integer points is an exact identity.
    Values are computed exactly over Q (python-flint fmpq_mat), adj = det * inverse at points that are not
    eigenvalues (|x0| >= 10 > spectral radius 5)."""
    data = []
    for delta in (0, 1):
        E = cap_I(delta)
        A = [[0] * 9 for _ in range(9)]
        B = [[0] * 6 for _ in range(9)]
        for (a, b) in E:
            for (u, v) in ((a, b), (b, a)):
                if u >= 6 and v >= 6:
                    A[u - 6][v - 6] = 1
                elif u >= 6 > v:
                    B[u - 6][v] = 1
        data.append((A, B))
    pts = list(range(10, 22))
    ok = True
    for x0 in pts:
        vals = []
        for A, B in data:
            M = flint.fmpq_mat(9, 9, [(x0 if i == j else 0) - A[i][j] for i in range(9) for j in range(9)])
            d = M.det()
            Bm = flint.fmpq_mat(9, 6, [B[i][j] for i in range(9) for j in range(6)])
            vals.append((d, Bm.transpose() * M.inv() * Bm * d))
        ok &= vals[0][0] == vals[1][0] and vals[0][1] == vals[1][1]
    A0 = data[0][0]
    pI = flint.fmpz_mat(9, 9, [A0[i][j] for i in range(9) for j in range(9)]).charpoly()
    return ok, pI, None


def row_action(capg6):
    """automorphisms of the patch H = (cap + its row) and their action on the row (positions 0..5, cyclic):
    returns lists of rotation amounts t (pi(j) = j + t) and reflection constants c (pi(j) = c - j)."""
    h = from_graph6(capg6)
    nz = len(h) - 6
    H = nx.Graph()
    H.add_nodes_from(range(len(h)))
    H.add_edges_from((i, j) for i in range(len(h)) for j in h[i] if i < j)
    rots, refl = [], []
    for iso in nx.algorithms.isomorphism.GraphMatcher(H, H).isomorphisms_iter():
        pi = [iso[nz + j] - nz for j in range(6)]
        assert all(0 <= x < 6 for x in pi), 'automorphism must preserve the row (degree-4 vertices)'
        if len({(pi[j] - j) % 6 for j in range(6)}) == 1:
            rots.append((pi[0]) % 6)
        elif len({(pi[j] + j) % 6 for j in range(6)}) == 1:
            refl.append(pi[0] % 6)
        else:
            raise AssertionError('not a dihedral action')
    return rots, refl


def I_patch_g6():
    """the cap I(0) together with its row S, in the find_caps format (cap vertices first, then the row)."""
    E = cap_I(0)
    # relabel: cap vertices 6..14 -> 0..8, row 0..5 -> 9..14
    lab = lambda v: v - 6 if v >= 6 else v + 9
    adj = [set() for _ in range(15)]
    for (a, b) in E:
        adj[lab(a)].add(lab(b))
        adj[lab(b)].add(lab(a))
    for j in range(6):
        adj[9 + j].add(9 + (j + 1) % 6)
        adj[9 + (j + 1) % 6].add(9 + j)
    return to_graph6([sorted(a) for a in adj])


def predicted_isomorphic(rots, refl, is_I, k):
    """Prop. B of result.md: G(Z,k,0) ~ G(Z,k,1) iff an automorphism of the rest acts on row k as a
    parity-swapping dihedral element (odd rotation, or reflection with constant c + k - 1 even), or (Z = I) an
    end-swapping isomorphism exists (iff k is even)."""
    return any(t % 2 for t in rots) or any((c + k - 1) % 2 == 0 for c in refl) or (is_I and k % 2 == 0)


def main():
    kmax = int(sys.argv[1]) if len(sys.argv) > 1 else 8
    ok, pI, W0 = strong_condition_I()
    report('F2a strong cap condition for I: det(xI-A_I) and W(x) identical for the two attachments (exact, Z[x])', ok,
           'det(xI-A_I) = %s' % pI)
    caps = json.load(open(os.path.join(CERT, 'caps_6row_tube_n52.json')))
    base = [c for c in caps if c['size'] in (7, 8, 9, 10, 12)]
    table = []
    okF1 = okF2 = okF3 = True
    Ican = subprocess.run([LABELG, '-q'], input=(I_patch_g6() + '\n').encode(), capture_output=True,
                          check=True).stdout.decode().split()[0]
    for c in base:
        pattern = []
        rots, refl = row_action(c['patch_graph6'])
        mycan = subprocess.run([LABELG, '-q'], input=(c['patch_graph6'] + '\n').encode(), capture_output=True,
                               check=True).stdout.decode().split()[0]
        is_I = mycan == Ican
        print('cap |Z|=%d: automorphisms act on its row as rotations %s and reflections c = %s; Z = I: %s'
              % (c['size'], sorted(rots), sorted(refl), is_I))
        for k in range(1, kmax + 1):
            g0, g1 = build(c['patch_graph6'], k, 0), build(c['patch_graph6'], k, 1)
            okF1 &= is_dual(g0) and is_dual(g1)
            okF2 &= charpoly(g0) == charpoly(g1)
            r = subprocess.run([LABELG, '-q'], input=(to_graph6(g0) + '\n' + to_graph6(g1) + '\n').encode(),
                               capture_output=True, check=True).stdout.decode().split()
            noniso = r[0] != r[1]
            okF3 &= noniso == (not predicted_isomorphic(rots, refl, is_I, k))
            n = 2 * (len(g0) - 2)
            pattern.append((k, n, noniso))
            table.append({'cap_size': c['size'], 'cap_aut': c['aut'], 'k': k, 'n': n, 'nonisomorphic': noniso,
                          'graph6_delta0': to_graph6(g0), 'graph6_delta1': to_graph6(g1)})
        print('cap |Z|=%d |Aut|=%d: ' % (c['size'], c['aut']) +
              ' '.join('k=%d(n=%d):%s' % (k, n, 'PAIR' if ni else 'iso') for k, n, ni in pattern))
    report('F1 all constructed graphs are fullerene duals', okF1)
    report('F3 nauty (non)isomorphism of G(Z,k,0), G(Z,k,1) agrees with the symmetry prediction of Prop. B '
           '(all caps, k <= %d)' % kmax, okF3)
    report('F2b G(Z,k,0) and G(Z,k,1) have identical integer characteristic polynomials (all Z, all k <= %d)' % kmax,
           okF2)
    with open(os.path.join(CERT, 'tube_family_members.json'), 'w') as f:
        json.dump(table, f, indent=1)
    if FAILED:
        sys.exit(1)


if __name__ == '__main__':
    main()
