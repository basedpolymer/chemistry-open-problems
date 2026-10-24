#!/usr/bin/env python3
"""make_certificates.py -- collect the per-n scan results (certificates/scan/nNNN.json) into the certificate files:

  certificates/counts_per_n.csv            one row per n: isomer count (buckygen) vs OEIS A007894, and the
                                           exact numbers of cospectral classes for A(T), A(T6), A+D, 2A+D,
                                           and the (hash) number of non-unique T5 spectra
  certificates/cospectral_pairs_T.json     every class of non-isomorphic C_n isomers whose duals T_n have the same
                                           adjacency characteristic polynomial (graph6 of the duals, canonical
                                           face spirals, Atlas numbers for n <= ATLAS_MAX, characteristic polynomial)
  certificates/cospectral_pairs_signless_laplacian.json   same for the signless Laplacian A + D of T_n
  certificates/cospectral_pairs_2A_plus_D.json            same for 2A + D (empty if none)
  certificates/T6_cospectral_classes.json                 classes of isomers with cospectral T_n^6

Usage: python3 make_certificates.py [ATLAS_MAX]   (default 70; Atlas numbers need canonical spirals of all isomers)
Environment: python3.13, python-flint 0.9, networkx; buckygen only for the Atlas numbering. Runtime ~5 min.
"""
import csv
import glob
import json
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
from planarcode import read_planar_code, to_graph6  # noqa: E402
from spiral import canonical_spiral  # noqa: E402
from scan_driver import BUCKYGEN  # noqa: E402
from verify_pairs import adj_from_graph6_embedded  # noqa: E402

CERT = os.path.join(ROOT, 'certificates')


def atlas_map(n):
    """canonical spiral -> Atlas number (rank in lexicographic order) for all isomers of C_n."""
    data = subprocess.run([BUCKYGEN, '-q', '%dd' % n], capture_output=True, check=True).stdout
    sps = sorted(canonical_spiral(g) for g in read_planar_code(data))
    return {sp: r + 1 for r, sp in enumerate(sps)}


def main():
    atlas_max = int(sys.argv[1]) if len(sys.argv) > 1 else 70
    files = sorted(glob.glob(os.path.join(CERT, 'scan', 'n*.json')))
    rows = []
    outT, outQ, out2, out6 = [], [], [], []
    for fn in files:
        d = json.load(open(fn))
        n = d['n']
        t = d['types']
        row = {'n': n, 'dual_vertices': d['m'], 'isomers_buckygen': d['isomers'], 'isomers_OEIS_A007894': d['oeis_A007894'],
               'counts_match': d['count_matches_oeis'],
               'T_cospectral_classes': t['0']['cospectral_classes_exact'],
               'T_isomers_in_classes': t['0']['isomers_in_cospectral_classes'],
               'T6_cospectral_classes': t['1']['cospectral_classes_exact'],
               'T6_isomers_in_classes': t['1']['isomers_in_cospectral_classes'],
               'T6_classes_with_nonisomorphic_T6': sum(1 for c in t['1']['classes']
                                                       if not c['T6_graphs_all_isomorphic_labelg']),
               'T5_nonunique_spectra_hash': t['2']['keys_with_multiple_graphs'],
               'signless_A_plus_D_classes': t['4']['cospectral_classes_exact'] if '4' in t else '',
               'twoA_plus_D_classes': t['3']['cospectral_classes_exact'] if '3' in t else '',
               'K': d['K'], 'hash_seconds': d['hash_seconds']}
        rows.append(row)
        amap = None
        need = any(t[k]['classes'] for k in ('0', '1', '3', '4') if k in t)
        if need and n <= atlas_max:
            amap = atlas_map(n)
        for key, out in (('0', outT), ('4', outQ), ('3', out2), ('1', out6)):
            if key not in t:
                continue
            for c in t[key]['classes']:
                sps = []
                for g6 in c['graph6_dual']:
                    adj = adj_from_graph6_embedded(g6)
                    sps.append(list(canonical_spiral(adj)))
                item = {'n': n, 'matrix': t[key]['name'], 'graph6_dual': c['graph6_dual'],
                        'canonical_spirals': sps,
                        'atlas_numbers': [amap[tuple(s)] for s in sps] if amap else None,
                        'charpoly_coeffs_low_to_high': c['charpoly_coeffs_low_to_high']}
                if key == '1':
                    item['T6_graph6'] = c['T6_graph6']
                    item['T6_graphs_all_isomorphic'] = c['T6_graphs_all_isomorphic_labelg']
                out.append(item)
        print('n=%d done' % n, flush=True)
    with open(os.path.join(CERT, 'counts_per_n.csv'), 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    for name, out in (('cospectral_pairs_T.json', outT), ('cospectral_pairs_signless_laplacian.json', outQ),
                      ('cospectral_pairs_2A_plus_D.json', out2), ('T6_cospectral_classes.json', out6)):
        with open(os.path.join(CERT, name), 'w') as f:
            json.dump(out, f, indent=1)
        print(name, len(out), 'classes')


if __name__ == '__main__':
    main()
