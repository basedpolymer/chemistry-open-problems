#!/usr/bin/env python3
"""Topology-wide DFM lower bounds and an exact ceiling on this proof method.

Python 3.12+, numpy, python-flint; --build also needs cvxpy/CLARABEL.
One process, about a minute to verify all 121 entries, several to build.
Default verification uses only saved integers, not an optimizer or source archive.
The ceiling limits THIS lower-bound construction, not physical DFM energies.
"""
from fractions import Fraction as Q
from itertools import combinations
from pathlib import Path
import argparse
import hashlib
import json
import math
import time

import numpy as np
from flint import arb, ctx

from audit_certificates import canon, cross, dot, sub, require, matvec
from certify_catalogue import verify as verify_complex, S

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT/'certificates'/'catalogue_topology.npz'
CERT = ROOT/'certificates'/'topology_calibrations.json'
REPORT = ROOT/'certificates'/'topology_calibration_report.json'
ctx.prec = 192


def det(matrix):
    return dot(matrix[0], cross(matrix[1], matrix[2]))


def edges(tets):
    return sorted({canon(pair) for tet in tets for pair in
                   combinations([tuple(map(int, v)) for v in tet], 2)})


def rounded_certificate(W, Z):
    """Round then shrink: every edge inequality is subsequently checked exactly."""
    return (np.rint(W*2**40).astype(np.int64).tolist(),
            np.rint(Z*2**40).astype(np.int64).tolist())


def edge_squared(W, Z, edge):
    a, b = edge
    shift = tuple(b[k]-a[k] for k in range(1, 4))
    ws = matvec(W, shift)
    v = tuple(Z[a[0]][k]-Z[b[0]][k]+ws[k] for k in range(3))
    return dot(v, v)


def lower_power9(W, denominator, n):
    d = Q(det(W), denominator**3)
    return arb(36*d.numerator**2)/(d.denominator**2*n*n)/arb.pi()**2


def build(data, ids, fresh=False):
    import cvxpy as cp
    rows = []
    digest = hashlib.sha256(DATA.read_bytes()).hexdigest()
    cached = {}
    if CERT.exists() and not fresh:
        previous = json.loads(CERT.read_text())
        require(previous['catalogue_sha256'] == digest, 'cached field/archive mismatch')
        cached = {c['id']: c for c in previous['certificates']}
    for ident in ids:
        if ident in cached:
            rows.append(cached[ident])
            continue
        n = len(data[ident+'_points'])
        E = edges(data[ident+'_tets'])
        Z = cp.Variable((n, 3))
        W = cp.Variable((3, 3), symmetric=True)
        constraints = [Z[0] == 0, W >> 0]
        for a, b in E:
            shift = np.array(b[1:])-np.array(a[1:])
            constraints.append(cp.norm(Z[a[0]]-Z[b[0]]+W@shift, 2) <= 1)
        problem = cp.Problem(cp.Maximize(cp.log_det(W)), constraints)
        try:
            problem.solve(solver='CLARABEL', tol_gap_abs=1e-9,
                          tol_gap_rel=1e-9, tol_feas=1e-9, max_iter=250)
        except cp.error.SolverError:
            # Only proposals: exact feasibility is enforced below and rechecked.
            try:
                problem.solve(solver='CLARABEL', tol_gap_abs=1e-6,
                              tol_gap_rel=1e-6, tol_feas=1e-6, max_iter=500)
            except cp.error.SolverError:
                pass
        if W.value is not None and Z.value is not None:
            wi, zi = rounded_certificate(W.value, Z.value)
            solver_status = problem.status
        else:
            # Any nonsingular periodic embedding gives a weaker feasible field
            # after shrinking. Solver success is never needed for certification.
            wi, zi = rounded_certificate(np.eye(3), -data[ident+'_points']/S)
            solver_status = 'fallback_catalogue_embedding'
        largest = max(edge_squared(wi, zi, e) for e in E)
        denominator = math.isqrt(largest)+1
        p9 = lower_power9(wi, denominator, n)
        target_int = math.floor(float(p9)**(1/9)*10**9)
        while not p9 > (arb(target_int)/10**9)**9:
            target_int -= 1
        rows.append(dict(id=ident, denominator=denominator, W=wi, Z=zi,
                         F_lower=f'{target_int//10**9}.{target_int%10**9:09d}',
                         proposal_solver_status=solver_status))
        print('BUILT', ident, rows[-1]['F_lower'], flush=True)
        CERT.write_text(json.dumps(dict(catalogue_sha256=digest,
                                       certificates=rows), indent=2)+'\n')


def verify_field(c, points, tets):
    W, Z, denominator = c['W'], c['Z'], c['denominator']
    n = len(points)
    require(type(denominator) is int and denominator > 0, 'invalid denominator')
    require(len(W) == 3 and all(len(r) == 3 for r in W), 'invalid W shape')
    require(len(Z) == n and all(len(r) == 3 for r in Z), 'invalid Z shape')
    require(all(type(v) is int for row in W+Z for v in row), 'integer fields required')
    require(det(W) > 0, 'nonpositive determinant')
    E = edges(tets)
    margin = denominator**2-max(edge_squared(W, Z, e) for e in E)
    require(margin >= 0, 'interface norm constraint failed')
    # Independent exact check of the signed-volume identity used in the ceiling.
    volume6 = 0
    for raw_tet in tets:
        tet = [tuple(map(int, v)) for v in raw_tet]
        old = [tuple(int(points[v[0], k])+S*v[k+1] for k in range(3)) for v in tet]
        olddet = dot(sub(old[1], old[0]), cross(sub(old[2], old[0]), sub(old[3], old[0])))
        new = [tuple(-Z[v[0]][k]+matvec(W, v[1:])[k] for k in range(3)) for v in tet]
        newdet = dot(sub(new[1], new[0]), cross(sub(new[2], new[0]), sub(new[3], new[0])))
        volume6 += (1 if olddet > 0 else -1)*newdet
    require(volume6 == 6*det(W), 'mapped signed volume identity failed')
    target = Q(c['F_lower'])
    require(target > 0, 'nonpositive claimed bound')
    require(lower_power9(W, denominator, n) >
            (arb(target.numerator)/target.denominator)**9, 'DFM lower bound failed')
    return dict(id=c['id'], status='PASS', cells=n, tetrahedra=len(tets),
                F_lower=c['F_lower'], min_squared_norm_margin=str(margin),
                det_W=str(Q(det(W), denominator**3)), signed_volume_identity=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--build', action='store_true')
    ap.add_argument('--fresh', action='store_true', help='recompute every proposal with --build')
    args = ap.parse_args()
    manifest = json.loads((ROOT/'certificates'/'catalogue_topology_manifest.json').read_text())
    digest = hashlib.sha256(DATA.read_bytes()).hexdigest()
    require(digest == manifest['npz_sha256'], 'catalogue archive hash mismatch')
    ids = sorted(e['id'] for e in manifest['entries'])
    require(len(ids) == len(set(ids)) == 121, 'expected all 121 topologies')
    with np.load(DATA, allow_pickle=False) as data:
        if args.build:
            build(data, ids, args.fresh)
        cert = json.loads(CERT.read_text())
        require(cert['catalogue_sha256'] == digest, 'field/archive mismatch')
        require(sorted(c['id'] for c in cert['certificates']) == ids, 'field set mismatch')
        rows = []
        for c in cert['certificates']:
            start = time.time()
            ident = c['id']
            p, t = data[ident+'_points'], data[ident+'_tets']
            verify_complex(p, t)
            row = verify_field(c, p, t)
            rows.append(row)
            print('PASS', ident, 'all shapes and boxes: F >', row['F_lower'],
                  f'{time.time()-start:.2f}s', flush=True)
    # A unit-edge tetrahedron has volume <= 1/(6 sqrt(2)).
    # Consequently det(W) <= t/(6 sqrt(2)); the generated lower number cannot
    # exceed [t^2/(2 pi^2 n^2)]^(1/9), even for the best possible feasible field.
    ratio = max(Q(r['tetrahedra'], r['cells']) for r in rows)
    require(ratio == Q(124, 21), 'unexpected maximum t/n')
    ceiling9 = arb(ratio.numerator)**2/(2*ratio.denominator**2*arb.pi()**2)
    ceiling = Q('1.065253040')
    require(ceiling9 < (arb(ceiling.numerator)/ceiling.denominator)**9,
            'method ceiling not below 1.065253040')
    summary = dict(catalogue_sha256=digest, claims=rows,
                   bounds_above_uniform_gap=sum(Q(r['F_lower']) > Q('1.00000025') for r in rows),
                   max_tetrahedra_per_cell=str(ratio),
                   method_ceiling_strictly_below='1.065253040',
                   ceiling_scope='linear divergence calibration plus separate ball moment bound; NOT an upper bound on DFM energies')
    REPORT.write_text(json.dumps(summary, indent=2)+'\n')
    print('PASS all 121 topology bounds;', summary['bounds_above_uniform_gap'], 'improve the uniform gap')
    print('PASS this proof method has a strict ceiling < 1.065253040 < 1.066 for every catalogue topology')


if __name__ == '__main__':
    try:
        main()
    except Exception as exc:
        print('FAIL topology calibration verification:', exc, flush=True)
        raise SystemExit(1)
