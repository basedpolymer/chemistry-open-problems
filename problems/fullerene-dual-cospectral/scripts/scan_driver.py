#!/usr/bin/env python3
"""scan_driver.py -- exhaustive scan, per n, of all fullerene isomers C_n for spectral collisions of their duals.

Pipeline for each even n (C_n has a dual triangulation T_n with m = n/2 + 2 vertices):
  1. buckygen -q <n>d | dualspec32 MASK K   (hash of tr(M^k) mod 2^32, k <= K, for every isomer and every
     matrix type M; complete: cospectral => same key; prints all key collisions);
  2. buckygen -q <n>d | pcselect idx        (planar codes of every graph involved in a collision);
  3. EXACT settlement with python-flint: integer characteristic polynomials of the relevant matrices; a
     collision is genuine iff the polynomials are equal; non-isomorphism of the two duals is certified by
     nauty's canonical labelling (labelg), and for T6 collisions we also test whether the two T6 graphs
     are isomorphic (labelg on the induced subgraphs);
  4. writes certificates/scan/n<NNN>.json (counts, classes with graph6 codes, timings).

Matrix types: 0 A(T), 1 A(T6), 2 A(T5) (hash counts only), 3 2A+D on T, 4 A+D on T.
Usage: python3 scan_driver.py NMIN NMAX [MASK] [K]
Environment: python3.13, python-flint 0.9, buckygen 1.1, nauty 2.8.8 (labelg); paths below.
Runtime: n <= 100: minutes; n = 150: ~1 h on one (loaded) core.
"""
import json
import os
import subprocess
import sys
import time

import flint

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
from planarcode import read_planar_code, to_graph6, check_triangulation_dual_of_fullerene  # noqa: E402

SRC = os.path.join(ROOT, 'certificates', 'sources')
BUCKYGEN = os.path.join(SRC, 'tools', 'bg', 'buckygen-1.1', 'buckygen')
DUALSPEC = os.path.join(SRC, 'build', 'dualspec32')
PCSELECT = os.path.join(SRC, 'build', 'pcselect')
LABELG = os.path.join(SRC, 'tools', 'nauty2_8_8', 'labelg')
OUTDIR = os.path.join(ROOT, 'certificates', 'scan')
RAWDIR = os.path.join(SRC, 'scanraw')   # large raw files stay git-ignored
OEIS = os.path.join(SRC, 'oeis', 'b007894.txt')
TYPENAMES = {0: 'A(T)', 1: 'A(T6)', 2: 'A(T5)', 3: '2A+D(T)', 4: 'A+D(T)', 5: 'A+tP(T)'}


def oeis_counts():
    d = {}
    with open(OEIS) as f:
        for line in f:
            p = line.split()
            if len(p) == 2 and not line.startswith('#'):
                d[2 * int(p[0])] = int(p[1])
    return d


def matrix(adj, ty):
    m = len(adj)
    deg = [len(a) for a in adj]
    if ty == 1:
        verts = [v for v in range(m) if deg[v] == 6]
    elif ty == 2:
        verts = [v for v in range(m) if deg[v] == 5]
    else:
        verts = list(range(m))
    loc = {v: i for i, v in enumerate(verts)}
    k = len(verts)
    M = [[0] * k for _ in range(k)]
    for v in verts:
        i = loc[v]
        for w in adj[v]:
            if w in loc:
                M[i][loc[w]] = 2 if ty == 3 else 1
        if ty in (3, 4):
            M[i][i] = deg[v]
    return M, verts


def charpoly(adj, ty):
    M, _ = matrix(adj, ty)
    if not M:
        return (1,)
    return tuple(int(c) for c in flint.fmpz_mat(M).charpoly().coeffs())


def induced_adj(adj, verts):
    loc = {v: i for i, v in enumerate(verts)}
    return [[loc[w] for w in adj[v] if w in loc] for v in verts]


def labelg_canon(g6list):
    """canonical graph6 strings from nauty labelg (same order)."""
    if not g6list:
        return []
    r = subprocess.run([LABELG, '-q'], input=('\n'.join(g6list) + '\n').encode(), capture_output=True, check=True)
    out = r.stdout.decode().split()
    assert len(out) == len(g6list)
    return out


def run_n(n, mask, K, oe):
    os.makedirs(OUTDIR, exist_ok=True)
    os.makedirs(RAWDIR, exist_ok=True)
    m = n // 2 + 2
    t0 = time.time()
    rawpath = os.path.join(RAWDIR, 'n%03d.txt' % n)
    with open(rawpath, 'w') as fout:
        p1 = subprocess.Popen(['nice', '-n', '10', BUCKYGEN, '-q', '%dd' % n], stdout=subprocess.PIPE,
                              stderr=subprocess.DEVNULL)
        p2 = subprocess.Popen(['nice', '-n', '10', DUALSPEC, str(mask), str(K)], stdin=p1.stdout, stdout=fout)
        p1.stdout.close()
        p2.wait()
        p1.wait()
        assert p2.returncode == 0 and p1.returncode == 0, (p1.returncode, p2.returncode)
    t_hash = time.time() - t0
    coll = []
    summ = {}
    with open(rawpath) as f:
        for line in f:
            if line.startswith('C '):
                _, ty, i, j = line.split()
                coll.append((int(ty), int(i), int(j)))
            elif line.startswith('SUMMARY'):
                parts = dict(kv.split('=', 1) for kv in line.replace('SUMMARY', ' ').split() if '=' in kv)
                summ[int(parts['type'])] = parts
    graphs = int(next(iter(summ.values()))['graphs'])
    # exact settlement
    idx = sorted({i for c in coll for i in c[1:]})
    sel = {}
    t1 = time.time()
    if idx:
        idxfile = os.path.join(RAWDIR, 'n%03d_idx.txt' % n)
        with open(idxfile, 'w') as f:
            f.write('\n'.join(map(str, idx)) + '\n')
        p1 = subprocess.Popen([BUCKYGEN, '-q', '%dd' % n], stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
        p2 = subprocess.run([PCSELECT, idxfile], stdin=p1.stdout, capture_output=True)
        p1.stdout.close()
        p1.wait()
        assert p2.returncode == 0
        gl = list(read_planar_code(p2.stdout))
        assert len(gl) == len(idx)
        sel = dict(zip(idx, gl))
    res = {'n': n, 'm': m, 'isomers': graphs, 'oeis_A007894': oe.get(n), 'count_matches_oeis': oe.get(n) == graphs,
           'mask': mask, 'K': K, 'hash_seconds': round(t_hash, 1), 'types': {}}
    for g in sel.values():
        errs = check_triangulation_dual_of_fullerene(g)
        assert not errs and len(g) == m, errs
    canon_T = {}
    if sel:
        il = sorted(sel)
        cs = labelg_canon([to_graph6(sel[i]) for i in il])
        canon_T = dict(zip(il, cs))
    for ty in sorted(summ):
        s = summ[ty]
        entry = {'name': TYPENAMES[ty], 'distinct_keys': int(s['distinct_keys']),
                 'keys_with_multiple_graphs': int(s['keys_with_multiple_graphs']),
                 'graphs_in_multiple_keys': int(s.get('graphs_in_multiple_keys', -1))}
        if ty == 2:
            entry['note'] = 'hash counts only (traces mod 2^32, k <= min(K,12)); not settled exactly'
            res['types'][str(ty)] = entry
            continue
        cl = [(i, j) for (t, i, j) in coll if t == ty]
        # key groups: dualspec reports (first graph with this key, new graph); partition each group EXACTLY
        groups = {}
        for (i, j) in cl:
            groups.setdefault(i, {i}).add(j)
        cp = {}
        classes = []
        nonsplit = 0
        for first, members in sorted(groups.items()):
            parts = {}
            for x in sorted(members):
                if x not in cp:
                    cp[x] = charpoly(sel[x], ty)
                parts.setdefault(cp[x], []).append(x)
            real = [sorted(v) for v in parts.values() if len(v) > 1]
            classes.extend(real)
            if len(parts) > 1:
                nonsplit += 1
        cls_out = []
        for mem in sorted(classes):
            cans = [canon_T[x] for x in mem]
            item = {'buckygen_indices': mem, 'graph6_dual': [to_graph6(sel[x]) for x in mem],
                    'duals_pairwise_nonisomorphic_labelg': len(set(cans)) == len(cans),
                    'charpoly_degree': len(cp[mem[0]]) - 1,
                    'charpoly_coeffs_low_to_high': [str(c) for c in cp[mem[0]]]}
            if ty == 1:
                sub = [to_graph6(induced_adj(sel[x], matrix(sel[x], 1)[1])) for x in mem]
                sc = labelg_canon(sub)
                item['T6_graph6'] = sub
                item['T6_graphs_all_isomorphic_labelg'] = len(set(sc)) == 1
                item['T6_graphs_pairwise_nonisomorphic_labelg'] = len(set(sc)) == len(sc)
            cls_out.append(item)
        entry.update({'hash_collisions': len(cl), 'key_groups': len(groups),
                      'key_groups_not_a_single_spectrum': nonsplit,
                      'cospectral_classes_exact': len(cls_out),
                      'isomers_in_cospectral_classes': sum(len(c['buckygen_indices']) for c in cls_out),
                      'classes': cls_out})
        res['types'][str(ty)] = entry
    res['settle_seconds'] = round(time.time() - t1, 1)
    with open(os.path.join(OUTDIR, 'n%03d.json' % n), 'w') as f:
        json.dump(res, f, indent=1)
    line = 'n=%d isomers=%d oeis_ok=%s' % (n, graphs, res['count_matches_oeis'])
    for ty, e in res['types'].items():
        if ty == '2':
            line += ' | T5 multi=%d' % e['keys_with_multiple_graphs']
        else:
            line += ' | %s: groups=%d mixed=%d classes=%d' % (e['name'], e['key_groups'],
                                                             e['key_groups_not_a_single_spectrum'],
                                                             e['cospectral_classes_exact'])
    line += ' | %.0fs+%.0fs' % (t_hash, res['settle_seconds'])
    print(line, flush=True)
    return res


def main():
    nmin, nmax = int(sys.argv[1]), int(sys.argv[2])
    mask = int(sys.argv[3]) if len(sys.argv) > 3 else 31
    K = int(sys.argv[4]) if len(sys.argv) > 4 else 20
    oe = oeis_counts()
    for n in range(nmin, nmax + 1, 2):
        if n == 22 or n < 20:
            continue
        if os.path.exists(os.path.join(OUTDIR, 'n%03d.json' % n)) and '--force' not in sys.argv:
            continue
        run_n(n, mask, K, oe)


if __name__ == '__main__':
    main()
