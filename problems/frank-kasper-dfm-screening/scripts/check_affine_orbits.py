#!/usr/bin/env python3
"""Global affine-orbit bounds for already audited integer partitions.

Python 3.12+, numpy/scipy/python-flint; several minutes in one process.
Integer integration computes the full centered moment tensor. Floating point
only proposes an SPD metric; Arb certifies the supporting-plane lower bound
for EVERY nonsingular affine deformation of the stored mesh. This is not a
lower bound for all partitions of its topology. See paper.tex for the proof.
"""
import argparse
from fractions import Fraction as Q
import hashlib
import json
import math
from pathlib import Path
import time

import numpy as np
from scipy.optimize import minimize
from scipy.linalg import eigvalsh
from flint import arb, arb_mat, ctx
from audit_certificates import sub, cross, dot, matvec, require

ctx.prec = 192
ROOT = Path(__file__).resolve().parents[1]
PAIRS = [(0, 0), (1, 1), (2, 2), (0, 1), (0, 2), (1, 2)]


def ball(x):
    x = Q(x)
    return arb(x.numerator)/x.denominator


def tensor_data(c):
    L, X, centers = (c[k].astype(object).tolist() for k in ('L', 'X', 'C'))
    n = int(c['ncell'])
    v6 = [0]*n
    m24 = [[0]*3 for _ in range(n)]
    s120 = [[0]*6 for _ in range(n)]
    normals = []
    for ids, sh, f in zip(c['tri'], c['tsh'], c['tface']):
        f = int(f)
        ps = [tuple(int(x)+s for x,s in zip(X[int(v)], matvec(L, list(map(int, off)))))
              for v,off in zip(ids, sh)]
        normal = cross(sub(ps[1], ps[0]), sub(ps[2], ps[0]))
        require(dot(normal, normal) > 0, 'zero normal')
        normals.append(normal)
        offset = matvec(L, list(map(int, c['fsj'][f])))
        for cell, sign, off in ((int(c['fi'][f]), 1, (0,0,0)),
                                (int(c['fj'][f]), -1, offset)):
            origin = tuple(a+b for a,b in zip(centers[cell], off))
            a,b,d = [sub(p, origin) for p in ps]
            det = sign*dot(a, cross(b,d))
            require(det > 0, 'nonpositive cone')
            total = tuple(a[k]+b[k]+d[k] for k in range(3))
            v6[cell] += det
            for i in range(3):
                m24[cell][i] += det*total[i]
            for k,(i,j) in enumerate(PAIRS):
                s120[cell][k] += det*(total[i]*total[j]+a[i]*a[j]+b[i]*b[j]+d[i]*d[j])
    volume = abs(dot(L[0], cross(L[1],L[2])))
    require(sum(v6) == 6*volume, 'volume mismatch')
    scale = 1 << int(c['K'])
    tensor = [[Q(0) for _ in range(3)] for _ in range(3)]
    for k,(i,j) in enumerate(PAIRS):
        tensor[i][j] = tensor[j][i] = sum(
            (Q(s120[cid][k],120)-Q(m24[cid][i]*m24[cid][j],96*v6[cid]))
            for cid in range(n))/scale**5
    return normals, tensor, Q(volume,scale**3), scale


def sym(z):
    return np.array([[z[0],z[3],z[4]], [z[3],z[1],z[5]], [z[4],z[5],z[2]]])


def six_gradient(m):
    return np.array([m[0,0],m[1,1],m[2,2],2*m[0,1],2*m[0,2],2*m[1,2]])


def propose(normals, tensor, scale):
    ns = np.array(normals,dtype=float)/scale**2
    qt = np.array(tensor,dtype=float)
    j0 = np.trace(qt)
    a0 = np.linalg.norm(ns,axis=1).sum()/2

    def values(z):
        g = sym(z)
        if np.linalg.eigvalsh(g)[0] <= 0:
            return 1e10, np.zeros(6)
        w = ns@np.linalg.inv(g)
        lengths = np.sqrt(np.einsum('ij,ij->i',ns,w))
        a = lengths.sum()/2
        kmat = (w.T/lengths)@w/2
        return a/a0, -six_gradient(kmat)/(2*a0)

    start = np.array([1.,1.,1.,0.,0.,0.])
    res = minimize(values,start,jac=True,method='SLSQP',
                   constraints=[{'type':'eq','fun':lambda z: np.sum(sym(z)*qt)/j0-1,
                                 'jac':lambda z: six_gradient(qt)/j0},
                                {'type':'ineq','fun':lambda z: np.linalg.eigvalsh(sym(z))-1e-6}],
                   options={'ftol':1e-14,'maxiter':200})
    g = sym(res.x)
    # Rational metric; optimizer convergence does not certify anything.
    metric = [[Q(round(float(x)*2**40),2**40) for x in row] for row in g]
    gf = np.array(metric,dtype=float)
    w = ns@np.linalg.inv(gf)
    lengths = np.sqrt(np.einsum('ij,ij->i',ns,w))
    kmat = (w.T/lengths)@w/2
    kappa = float(eigvalsh(kmat,qt)[-1])
    # Strict slack for the interval positive-definiteness check.
    kappa_q = Q(math.ceil(kappa*(1+1e-9)*2**40),2**40)
    return metric,kappa_q,str(res.message)


def positive_definite(m):
    return m[0,0] > 0 and (m[0,0]*m[1,1]-m[0,1]*m[1,0]) > 0 and m.det() > 0


def verify_data(normals,tensor,volume,scale,metric,kappa):
    gm = arb_mat([[ball(x) for x in row] for row in metric])
    qm = arb_mat([[ball(x) for x in row] for row in tensor])
    require(positive_definite(gm) and positive_definite(qm),'non-SPD metric or moment tensor')
    gi = gm.inv()
    area = arb(0)
    kval = [[arb(0) for _ in range(3)] for _ in range(3)]
    area0 = arb(0)
    for n in normals:
        nv = [arb(x)/scale**2 for x in n]
        w = [sum((gi[i,j]*nv[j] for j in range(3)),arb(0)) for i in range(3)]
        length = sum((nv[i]*w[i] for i in range(3)),arb(0)).sqrt()
        area += length/2
        area0 += sum((x*x for x in nv),arb(0)).sqrt()/2
        for i in range(3):
            for j in range(i,3):
                kval[i][j] += w[i]*w[j]/(2*length)
    for i in range(3):
        for j in range(i):
            kval[i][j] = kval[j][i]
    km = arb_mat(kval)
    diff = ball(kappa)*qm-km
    require(positive_definite(diff),'kappa Q - K not certified positive definite')
    moment = sum((gm[i,j]*qm[j,i] for i in range(3) for j in range(3)),arb(0))
    tangent = (3*area-ball(kappa)*moment)/2
    require(tangent > 0,'supporting plane gives no positive area bound')
    constant = arb(20)/(27*ball(volume)**3)
    lower_cube = constant*tangent*tangent*moment
    upper_cube = constant*area*area*moment
    lower = Q(math.floor(float(lower_cube.root(3))*10**9)-1,10**9)
    upper = Q(math.ceil(float(upper_cube.root(3))*10**9)+1,10**9)
    require(lower_cube > ball(lower)**3 and upper_cube < ball(upper)**3,'decimal bounds failed')
    baseline = (constant*area0*area0*sum((qm[i,i] for i in range(3)),arb(0))).root(3)
    return lower,upper,baseline


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('names',nargs='*')
    ap.add_argument('--rebuild',action='store_true',help='propose new metrics instead of reading saved ones')
    args = ap.parse_args()
    claims = json.loads((ROOT/'certificates'/'claims.json').read_text())
    audited = {r['name']:r for r in json.loads((ROOT/'certificates'/'audit_report.json').read_text())['claims']}
    outpath = ROOT/'certificates'/'affine_orbits_report.json'
    previous = json.loads(outpath.read_text())['claims'] if outpath.exists() else []
    saved = {r['name']:r for r in previous}
    rows = []
    for claim in claims:
        name = claim['name']
        if args.names and name not in args.names:
            continue
        started = time.time()
        path = ROOT/'certificates'/'meshes'/claim['file']
        sha = hashlib.sha256(path.read_bytes()).hexdigest()
        require(audited[name]['status']=='PASS' and audited[name]['sha256']==sha,'mesh not covered by prior topology audit')
        with np.load(path,allow_pickle=False) as z:
            c = {k:z[k] for k in z.files}
        print('Checking',name,'triangles',len(c['tri']),flush=True)
        ns,qt,v,scale = tensor_data(c)
        if name in saved and not args.rebuild:
            metric = [[Q(x) for x in row] for row in saved[name]['metric']]
            kappa = Q(saved[name]['kappa'])
            optimizer = 'not used; saved rational metric'
        else:
            metric,kappa,optimizer = propose(ns,qt,scale)
        lower,upper,baseline = verify_data(ns,qt,v,scale,metric,kappa)
        # Trace and normalization must recover the independently audited energy.
        geo = audited[name]['geometry']
        require(baseline > ball(Q(geo['mesh_energy_lower'])) and
                baseline < ball(Q(geo['mesh_energy_upper'])), 'baseline disagrees with exact scalar audit')
        row = {'name':name,'status':'PASS','mesh_sha256':sha,
               'scope':'all nonsingular affine images of this stored mesh, NOT its topology class',
               'metric':[[str(x) for x in r] for r in metric],'kappa':str(kappa),
               'affine_inf_lower_9dp':f'{float(lower):.9f}',
               'affine_inf_upper_9dp':f'{float(upper):.9f}',
               'original_mesh_energy':str(baseline),'precision_bits':ctx.prec,
               'proposal_status':optimizer,'runtime_s':round(time.time()-started,3)}
        rows.append(row)
        outpath.write_text(json.dumps({'status':'PASS','claims':rows},indent=2)+'\n')
        print('PASS',name,'affine infimum in',row['affine_inf_lower_9dp'],row['affine_inf_upper_9dp'],flush=True)
    require(bool(rows),'no selected certificates')


if __name__ == '__main__':
    main()
