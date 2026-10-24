#!/usr/bin/env python3
"""Exact embedded-complex certificates for all 84+37 catalogue entries.

Default: verify saved integer data (Python, numpy, python-flint).
--build: regenerate from the downloaded CGD catalogue (also scipy/gemmi).
Runtime: about a minute for verification, a few minutes to build.
This proves geometric/topological feasibility, not any energy lower bound.
The source papers, not this program, establish enumeration completeness.
"""
from collections import Counter, defaultdict
from pathlib import Path
import argparse
import hashlib
import json
import time

import numpy as np
from audit_certificates import canon, cross, dot, require, sub, surface

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT/'certificates'/'catalogue_topology.npz'
MANIFEST = ROOT/'certificates'/'catalogue_topology_manifest.json'
K = 44
S = 1 << K


def build():
    from dfm import cgd, topology, registry
    data, entries = {}, []
    for part, count, folder in ((1, 84, 'ListDATApaper'), (2, 37, 'ListDATApaperStr')):
        paths = sorted((Path(registry.SRC)/folder).glob('*.cgd'))
        paths = [p for p in paths if 'total' not in p.name]
        require(len(paths) == count, f'missing catalogue inputs: expected {count}, found {len(paths)}')
        for path in paths:
            name = ('I_' if part == 1 else 'II_') + path.stem
            block = cgd.parse_cgd(str(path))
            L, P, coordination, edges, missing, sg = cgd.expand_cgd(block)
            require(missing == 0, 'missing CGD edge endpoint')
            degree = cgd.contact_degree(len(P), edges)
            require(np.all(degree == coordination) or np.all(2*degree == coordination), 'coordination mismatch')
            tets = topology.tets_from_graph(len(P), edges)
            data[name+'_points'] = np.rint(P*S).astype(np.int64)
            data[name+'_tets'] = np.array([[(i, *s) for i, s in tet] for tet in tets], dtype=np.int32)
            entries.append({'id': name, 'group': block['group'], 'source': folder+'/'+path.name,
                            'source_sha256': hashlib.sha256(path.read_bytes()).hexdigest(),
                            'ncell': len(P)})
            print('BUILD', name, len(P), len(tets), flush=True)
    np.savez_compressed(DATA, **data)
    MANIFEST.write_text(json.dumps({'grid_exponent': K, 'archive_md5': '582371cee665c62e01a2f47a73dd2c9d',
                                  'npz_sha256': hashlib.sha256(DATA.read_bytes()).hexdigest(), 'entries': entries}, indent=2)+'\n')


def oriented_key(vertices):
    vertices = tuple(vertices)
    key = canon(vertices)
    # Translate to the frame used by the canonical (sorted) key.
    for o in vertices:
        translated = tuple((v[0], *(v[k]-o[k] for k in range(1, 4))) for v in vertices)
        if tuple(sorted(translated)) == key:
            indices = [key.index(v) for v in translated]
            inversions = sum(indices[i] > indices[j] for i in range(len(indices)) for j in range(i+1, len(indices)))
            return key, (-1 if inversions % 2 else 1)
    raise ValueError('canonical orientation failure')


def verify(points, tetrahedra):
    require(points.dtype.kind in 'iu' and tetrahedra.dtype.kind in 'iu', 'integer data required')
    P = [tuple(map(int, row)) for row in points]
    links = [[] for _ in P]
    oriented_faces, face_counts, edges = Counter(), Counter(), set()
    volume6, min_det, seen = 0, None, set()
    for tet_arr in tetrahedra:
        tet = [tuple(map(int, row)) for row in tet_arr]
        require(all(0 <= v[0] < len(P) for v in tet), 'invalid cell index')
        require(len(set(tet)) == 4, 'repeated tetrahedron vertex')
        key = canon(tet)
        require(key not in seen, 'duplicate periodic tetrahedron')
        seen.add(key)
        positions = [tuple(P[v[0]][i]+S*v[i+1] for i in range(3)) for v in tet]
        det = dot(sub(positions[1], positions[0]), cross(sub(positions[2], positions[0]), sub(positions[3], positions[0])))
        require(det != 0, 'degenerate tetrahedron')
        if det < 0:
            tet[0], tet[1] = tet[1], tet[0]
            det = -det
        volume6 += det
        min_det = det if min_det is None else min(min_det, det)
        for i in range(4):
            face = [tet[j] for j in range(4) if j != i]
            if i % 2:
                face[0], face[1] = face[1], face[0]
            fkey, sign = oriented_key(face)
            oriented_faces[fkey] += sign
            face_counts[fkey] += 1
            center = tet[i]
            link_face = tuple((v[0], *(v[k]-center[k] for k in range(1, 4))) for v in face)
            links[center[0]].append(link_face)
            for j in range(i+1, 4):
                edges.add(canon((tet[i], tet[j])))
    require(all(v == 0 for v in oriented_faces.values()), 'oriented tetrahedron faces do not cancel')
    require(all(v == 2 for v in face_counts.values()), 'triangle incidence is not two')
    require(volume6 == 6*S**3, 'degree is not exactly one')
    require(len(P)-len(edges)+len(face_counts)-len(tetrahedra) == 0, 'nonzero torus Euler characteristic')
    histogram = Counter()
    for ts in links:
        surface(ts)
        adj = defaultdict(set)
        for a, b, c in ts:
            adj[a].update((b, c))
            adj[b].update((a, c))
            adj[c].update((a, b))
        z = len(adj)
        require(z in (12, 14, 15, 16), 'non-FK coordination')
        require(all(len(nb) in (5, 6) for nb in adj.values()), 'non-pentagon/hexagon in dual cell')
        require(sum(len(nb) == 5 for nb in adj.values()) == 12, 'wrong pentagon count')
        require(all(not (len(adj[a]) == len(adj[b]) == 6) for a in adj for b in adj[a]), 'non-isolated hexagons')
        histogram[z] += 1
    return {'cells': len(P), 'tetrahedra': len(tetrahedra), 'degree': 1,
            'smallest_positive_determinant': str(min_det), 'coordination_histogram': dict(sorted(histogram.items())),
            'oriented_faces_cancel': True, 'vertex_links_are_spheres': True, 'isolated_hexagons': True}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--build', action='store_true')
    a = ap.parse_args()
    if a.build:
        build()
    manifest = json.loads(MANIFEST.read_text())
    require(manifest['grid_exponent'] == K, 'unexpected grid')
    require(hashlib.sha256(DATA.read_bytes()).hexdigest() == manifest['npz_sha256'], 'archive hash mismatch')
    entries = manifest['entries']
    require(len(entries) == 121 and len({e['id'] for e in entries}) == 121, 'expected 121 distinct entries')
    require(sum(e['id'].startswith('I_') for e in entries) == 84, 'expected 84 part-I entries')
    rows = []
    with np.load(DATA, allow_pickle=False) as data:
        require(set(data.files) == {e['id']+'_'+s for e in entries for s in ('points', 'tets')}, 'unexpected certificate keys')
        for e in entries:
            started = time.time()
            try:
                r = verify(data[e['id']+'_points'], data[e['id']+'_tets'])
                require(r['cells'] == e['ncell'], 'cell count mismatch')
                rows.append({'id': e['id'], 'status': 'PASS', **r})
                print('PASS', e['id'], 'exact embedded FK complex', r['cells'], 'cells', f'{time.time()-started:.2f}s', flush=True)
            except Exception as ex:
                rows.append({'id': e['id'], 'status': 'FAIL', 'error': str(ex)})
                print('FAIL', e['id'], str(ex), flush=True)
    (ROOT/'certificates'/'catalogue_exact_report.json').write_text(json.dumps(rows, indent=2)+'\n')
    fails = sum(r['status'] != 'PASS' for r in rows)
    print(len(rows)-fails, 'PASS;', fails, 'FAIL')
    raise SystemExit(bool(fails))


if __name__ == '__main__':
    main()
