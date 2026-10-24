#!/usr/bin/env python3
"""Recheck topology, exact energies and intended boxes of the six saved meshes.

Python 3.12+, numpy, python-flint; scipy/gemmi for the optional identity check.
No optimizer or legacy energy/certificate code is used to evaluate energies.
Runtime: a few minutes, one process. Output: PASS/FAIL and a JSON audit report.
The lower endpoint concerns THIS mesh, never the infimum of its topology.
"""
from collections import Counter, defaultdict
from fractions import Fraction as Q
from pathlib import Path
import argparse
import hashlib
import json
import math
import platform
import time

import numpy as np
from flint import arb, arb_mat, ctx

ROOT = Path(__file__).resolve().parents[1]
ctx.prec = 192


def require(condition, message):
    if not condition:
        raise ValueError(message)


def sub(a, b):
    return tuple(x - y for x, y in zip(a, b))


def cross(a, b):
    return (a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0])


def dot(a, b):
    return sum(x*y for x, y in zip(a, b))


def matvec(m, v):
    return tuple(dot(row, v) for row in m)


def connected(adj):
    if not adj:
        return False
    seen, todo = set(), [next(iter(adj))]
    while todo:
        v = todo.pop()
        if v not in seen:
            seen.add(v)
            todo.extend(adj[v] - seen)
    return len(seen) == len(adj)


def surface(triangles, disk=False):
    """Exact combinatorial 2-manifold test, including every vertex link."""
    directed, links, graph = Counter(), defaultdict(dict), defaultdict(set)
    count = 0
    for a, b, c in triangles:
        require(len({a, b, c}) == 3, 'repeated vertex in triangle')
        count += 1
        for v, u, w in ((a, b, c), (b, c, a), (c, a, b)):
            directed[v, u] += 1
            graph[v].update((u, w))
            lk = links[v]
            lk.setdefault(u, set()).add(w)
            lk.setdefault(w, set()).add(u)
    require(connected(graph), 'disconnected surface')
    boundary = defaultdict(set)
    for (v, u), multiplicity in directed.items():
        require(multiplicity == 1, 'duplicate oriented edge')
        if (u, v) not in directed:
            require(disk, 'unpaired boundary edge')
            boundary[v].add(u)
            boundary[u].add(v)
    for v, lk in links.items():
        require(connected(lk), 'disconnected vertex link')
        degrees = Counter(map(len, lk.values()))
        if v in boundary:
            require(degrees[1] == 2 and set(degrees) <= {1, 2}, 'non-path boundary link')
        else:
            require(set(degrees) == {2}, 'non-circular interior link')
    edges = (len(directed) + sum(map(len, boundary.values()))//2)//2
    require(len(graph) - edges + count == (1 if disk else 2), 'wrong Euler characteristic')
    if disk:
        require(connected(boundary) and all(len(x) == 2 for x in boundary.values()), 'boundary is not one circle')
    return boundary


def canon(vertices):
    """Periodic simplex key with labels (cell, sx, sy, sz)."""
    vertices = tuple(vertices)
    return min(tuple(sorted((v[0], *(v[k]-o[k] for k in range(1, 4))) for v in vertices))
               for o in vertices)


def topology_checks(c, name, compare_identity=True):
    tri, shifts, tf = c['tri'], c['tsh'], c['tface']
    fi, fj, fsj = c['fi'], c['fj'], c['fsj']
    n = int(c['ncell'])
    require(len({canon(((int(i), 0, 0, 0), (int(j), *map(int, s))))
                 for i, j, s in zip(fi, fj, fsj)}) == len(fi), 'duplicate periodic cell interface')
    incident = [set() for _ in c['X']]
    cell_tris = [[] for _ in range(n)]
    face_tris = [[] for _ in fi]
    for t, f in enumerate(tf):
        f = int(f)
        vs = [tuple([int(tri[t, k])] + [int(x) for x in shifts[t, k]]) for k in range(3)]
        face_tris[f].append(tuple(vs))
        cell_tris[int(fi[f])].append(tuple(vs))
        other = tuple((v[0], *(v[k+1]-int(fsj[f, k]) for k in range(3))) for v in vs[::-1])
        cell_tris[int(fj[f])].append(other)
        for v in vs:
            incident[v[0]].add((int(fi[f]), *(-s for s in v[1:])))
            incident[v[0]].add((int(fj[f]), *(int(fsj[f, k])-v[k+1] for k in range(3))))
    for ts in cell_tris:
        surface(ts)
    require(all(len(s) in (2, 3, 4) for s in incident), 'invalid cell incidence at vertex')
    corners, fsize = [], []
    foam_edge_faces, foam_edge_cells = defaultdict(set), defaultdict(set)
    corner_faces = defaultdict(set)
    for f, ts in enumerate(face_tris):
        boundary = surface(ts, disk=True)
        for v, neighbours in boundary.items():
            for u in neighbours:
                if v > u:
                    continue
                ek = canon((v, u))
                for o in (v, u):
                    if tuple(sorted((w[0], *(w[k]-o[k] for k in range(1, 4))) for w in (v, u))) == ek:
                        break
                foam_edge_faces[ek].add(f)
                foam_edge_cells[ek].add((int(fi[f]), *(-o[k] for k in range(1, 4))))
                foam_edge_cells[ek].add((int(fj[f]), *(int(fsj[f, k])-o[k+1] for k in range(3))))
        cv = [v for v in boundary if len(incident[v[0]]) == 4]
        require(all(len(incident[v[0]]) in (3, 4) for v in boundary), 'invalid foam edge')
        fsize.append(len(cv))
        corners.append(cv)
        for v in cv:
            corner_faces[v[0]].add(f)
    require(all(len(fs) == 3 and len(foam_edge_cells[e]) == 3 for e, fs in foam_edge_faces.items()),
            'foam edge does not meet exactly three faces and cells')
    dual_tets = set()
    for v, cells in enumerate(incident):
        if len(cells) == 4:
            require(len(corner_faces[v]) == 6, 'foam vertex does not have six faces')
            dual_tets.add(canon(cells))
    require(len(dual_tets) == sum(len(x) == 4 for x in incident), 'duplicate periodic foam vertex')
    degree = Counter(map(int, fi)) + Counter(map(int, fj))
    hexcount = Counter()
    for f, size in enumerate(fsize):
        if name == 'BCC':
            require(size in (4, 6), 'BCC face is not square/hexagon')
        else:
            require(size in (5, 6), 'non-fullerene face')
        if size == 6:
            hexcount[int(fi[f])] += 1
            hexcount[int(fj[f])] += 1
    if name != 'BCC':
        require(all(degree[i] in (12, 14, 15, 16) and hexcount[i] == degree[i]-12 for i in range(n)), 'wrong FK composition')
        # At a corner of any cell at most one of its three faces may be a hexagon.
        # Count in each lifted cell frame, not just by representative cell label.
        at_corner = defaultdict(int)
        for f, cv in enumerate(corners):
            if fsize[f] != 6:
                continue
            for v in cv:
                for cell, off in ((int(fi[f]), (0, 0, 0)), (int(fj[f]), fsj[f])):
                    key = (cell, v[0], *(v[k+1]-int(off[k]) for k in range(3)))
                    at_corner[key] += 1
        require(max(at_corner.values(), default=0) <= 1, 'adjacent hexagons in an FK cell')
    identity = 'not requested'
    if compare_identity:
        from dfm import registry
        expected = {canon((int(i), *map(int, s)) for i, s in tet) for tet in registry.get(name)['tets']}
        require(dual_tets == expected, 'dual complex differs from registry topology')
        identity = 'exact equality of periodic labelled dual tetrahedra'
    return {'cells': n, 'faces': len(fi), 'foam_vertices': len(dual_tets),
            'cell_degree_histogram': dict(sorted(Counter(degree.values()).items())),
            'face_size_histogram': dict(sorted(Counter(fsize).items())),
            'cell_boundaries': 'connected oriented 2-spheres with circular vertex links',
            'faces_topology': 'disks with one boundary circle', 'identity': identity}


def decimal_root_bounds(lo, hi, digits=12):
    q = 10**digits
    a = math.floor(float(lo)**(1/3)*q)
    b = math.ceil(float(hi)**(1/3)*q)
    while Q(a, q)**3 > lo:
        a -= 1
    while Q(a+1, q)**3 <= lo:
        a += 1
    while Q(b, q)**3 < hi:
        b += 1
    while Q(b-1, q)**3 >= hi:
        b -= 1
    return Q(a, q), Q(b, q)


def exact_energy(c):
    """Integrate each positive cone using simplex barycentric moments."""
    L, X, centers = (c[key].astype(object).tolist() for key in ('L', 'X', 'C'))
    n = int(c['ncell'])
    volume6, first24, moment120 = [0]*n, [[0]*3 for _ in range(n)], [0]*n
    area_floor, area_ceil, smallest = 0, 0, None
    for ids, sh, f in zip(c['tri'], c['tsh'], c['tface']):
        f = int(f)
        points = [tuple(int(x)+s for x, s in zip(X[int(v)], matvec(L, list(map(int, off))))) for v, off in zip(ids, sh)]
        normal = cross(sub(points[1], points[0]), sub(points[2], points[0]))
        nn = dot(normal, normal)
        require(nn > 0, 'zero-area triangle')
        root = math.isqrt(nn)
        area_floor += root
        area_ceil += root + (root*root != nn)
        offset = matvec(L, list(map(int, c['fsj'][f])))
        for cell, sign, off in ((int(c['fi'][f]), 1, (0, 0, 0)), (int(c['fj'][f]), -1, offset)):
            origin = tuple(a+b for a, b in zip(centers[cell], off))
            a, b, d = [sub(p, origin) for p in points]
            det = sign * dot(a, cross(b, d))
            require(det > 0, 'non-positive cone determinant')
            smallest = det if smallest is None else min(smallest, det)
            total = tuple(a[k]+b[k]+d[k] for k in range(3))
            volume6[cell] += det
            for k in range(3):
                first24[cell][k] += det * total[k]
            moment120[cell] += det * (dot(total, total)+dot(a, a)+dot(b, b)+dot(d, d))
    box = abs(dot(L[0], cross(L[1], L[2])))
    require(box > 0 and sum(volume6) == 6*box, 'exact volume mismatch')
    moments = [Q(s, 120)-sum(Q(x, 24)**2 for x in m)/Q(v, 6)
               for v, m, s in zip(volume6, first24, moment120)]
    require(all(i > 0 for i in moments), 'nonpositive centered second moment')
    factor = Q(20, 27) * sum(moments) / box**3
    low, high = decimal_root_bounds(factor*Q(area_floor, 2)**2, factor*Q(area_ceil, 2)**2)
    return low, high, {'positive_cones': True, 'exact_volume_sum': True,
                       'smallest_cone_determinant': str(smallest),
                       'mesh_energy_lower': str(low), 'mesh_energy_upper': str(high)}


def intended_box(c, name, upper):
    """Transport to the exact caption box via an Arb-bounded affine map.

    If ||T-Id||_2 <= e < 1, F(T Pi) <= (1+e)^2/(1-e)^2 F(Pi).
    Frobenius norm bounds operator norm; all comparisons use rigorous balls.
    """
    if name == 'sigma':
        return upper, {'box': 'exact dyadic matrix in certificate', 'affine_factor_bound': '1'}
    targets = {'BCC': ('1', '1', '1'), 'P': ('1.91', '3.57', '1'),
               'M': ('1.89', '1', '3.3'), 'H': ('1', '3.88', '1')}
    if name == 'mu':
        B = arb_mat([[1, arb(-1)/2, 0], [0, arb(3).sqrt()/2, 0], [0, 0, arb(26)/5]])
    else:
        values = [Q(x) for x in targets[name]]
        B = arb_mat([[arb(values[i].numerator)/values[i].denominator if i == j else 0 for j in range(3)] for i in range(3)])
    L = arb_mat([[arb(int(x))/(1 << int(c['K'])) for x in row] for row in c['L']])
    T = B * L.inv()
    entries = [T[i, j]-int(i == j) for i in range(3) for j in range(3)]
    # Multiplication works for intervals crossing zero; older python-flint
    # versions may route **2 through real powers and return nan there.
    sq = sum((x*x for x in entries), arb(0))
    # A slightly wider rational radius avoids sqrt of a ball touching zero.
    epsilon = Q(1, 10**12)
    require(sq < arb(epsilon.numerator)**2/epsilon.denominator**2,
            'affine perturbation exceeds bound: squared norm '+str(sq))
    factor = (1+epsilon)**2/(1-epsilon)**2
    transformed = upper*factor
    return transformed, {'box': 'exact caption metric', 'frobenius_error_bound': str(epsilon),
                         'affine_factor_bound': str(factor)}


def audit(claim, compare_identity=True):
    started = time.time()
    path = ROOT/'certificates'/'meshes'/claim['file']
    with np.load(path, allow_pickle=False) as z:
        c = {key: z[key] for key in z.files}
    for key, arr in c.items():
        require(np.issubdtype(arr.dtype, np.integer), 'noninteger array '+key)
    topo = topology_checks(c, claim['name'], compare_identity)
    lo, hi, geo = exact_energy(c)
    transported, box = intended_box(c, claim['name'], hi)
    threshold = Q(claim['threshold'])
    require(transported < threshold, 'threshold not beaten in intended box')
    # The publication prints six decimals. Beat even the lower rounding endpoint.
    if claim['name'] != 'sigma':
        require(transported < threshold-Q(1, 2_000_000), 'margin smaller than publication rounding')
    q = 10**9
    upper9 = Q(-((-transported.numerator*q)//transported.denominator), q)
    return {'name': claim['name'], 'file': claim['file'], 'sha256': hashlib.sha256(path.read_bytes()).hexdigest(),
            'status': 'PASS', 'topology': topo, 'geometry': geo, 'box_transport': box,
            'intended_box_upper_9dp': f'{float(upper9):.9f}', 'threshold': str(threshold),
            'runtime_s': round(time.time()-started, 3)}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('names', nargs='*')
    ap.add_argument('--skip-identity', action='store_true')
    ap.add_argument('--output', type=Path, default=ROOT/'certificates'/'audit_report.json')
    a = ap.parse_args()
    claims = json.loads((ROOT/'certificates'/'claims.json').read_text())
    rows = []
    for c in claims:
        if a.names and c['name'] not in a.names:
            continue
        try:
            row = audit(c, not a.skip_identity)
            print('PASS', row['name'], 'F <=', row['intended_box_upper_9dp'], 'topology / exact energy / box', flush=True)
        except Exception as e:
            row = {'name': c['name'], 'status': 'FAIL', 'error': str(e)}
            print('FAIL', c['name'], str(e), flush=True)
        rows.append(row)
        a.output.write_text(json.dumps({'python': platform.python_version(), 'claims': rows}, indent=2)+'\n')
    require(bool(rows), 'no claims selected')
    raise SystemExit(0 if all(r['status'] == 'PASS' for r in rows) else 1)


if __name__ == '__main__':
    main()
