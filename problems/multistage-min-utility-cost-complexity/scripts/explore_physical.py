#!/usr/bin/env python3
"""Floating-point exploration ONLY, not a proof/certificate.
WSL Python + pyscipopt 6.3.0 / SCIP; deterministic seed, 30 s per case.
Run: python scripts/explore_physical.py. Typical total runtime <5 minutes.
Prints PASS/FAIL for numerical expectations; nonzero exit on mismatch/unresolved.
The exact companion script, not this solver log, certifies mathematical examples.
"""
import json
from pathlib import Path
from pyscipopt import Model, quicksum

ROOT = Path(__file__).resolve().parents[1]


def solve(hot, cold, ell, version, limit=30):
    model = Model('utility_' + version)
    model.hideOutput()
    model.setRealParam('limits/time', limit)
    model.setIntParam('randomization/randomseedshift', 0)
    model.setIntParam('parallel/maxnthreads', 1)
    model.setRealParam('numerics/feastol', 1e-8)
    n, m = len(hot), len(cold)
    h = {(i,k): model.addVar(lb=ho, ub=hi, name=f'h_{i}_{k}') for i,(f,hi,ho) in enumerate(hot) for k in range(ell+1)}
    c = {(j,k): model.addVar(lb=ci, ub=co, name=f'c_{j}_{k}') for j,(f,ci,co) in enumerate(cold) for k in range(ell+1)}
    q = {(i,j,k): model.addVar(lb=0, ub=min(hot[i][0]*(hot[i][1]-hot[i][2]),cold[j][0]*(cold[j][2]-cold[j][1])), name=f'q_{i}_{j}_{k}') for i in range(n) for j in range(m) for k in range(ell)}
    for i,(fi,hi,ho) in enumerate(hot):
        model.addCons(h[i,0] == hi)
        for k in range(ell):
            model.addCons(fi*(h[i,k]-h[i,k+1]) == quicksum(q[i,j,k] for j in range(m)))
    for j,(fj,ci,co) in enumerate(cold):
        model.addCons(c[j,ell] == ci)
        for k in range(ell):
            model.addCons(fj*(c[j,k]-c[j,k+1]) == quicksum(q[i,j,k] for i in range(n)))
    if version == 'boundary':
        for i,j,k in q:
            z = model.addVar(vtype='B', name=f'z_{i}_{j}_{k}')
            cap = min(hot[i][0]*(hot[i][1]-hot[i][2]),cold[j][0]*(cold[j][2]-cold[j][1]))
            big = max(0, cold[j][2]-hot[i][2])
            model.addCons(q[i,j,k] <= cap*z)
            model.addCons(h[i,k]-c[j,k] >= -big*(1-z))
            model.addCons(h[i,k+1]-c[j,k+1] >= -big*(1-z))
    elif version == 'countercurrent':
        fh, fc = {}, {}
        max_t = max([x[1] for x in hot]+[x[2] for x in cold])
        for i,j,k in q:
            fi, hi, ho = hot[i]
            fj, ci, co = cold[j]
            fh[i,j,k] = model.addVar(lb=0, ub=fi, name=f'fh_{i}_{j}_{k}')
            fc[i,j,k] = model.addVar(lb=0, ub=fj, name=f'fc_{i}_{j}_{k}')
            a = model.addVar(lb=0, ub=max_t, name=f'a_{i}_{j}_{k}')
            b = model.addVar(lb=0, ub=max_t, name=f'b_{i}_{j}_{k}')
            model.addCons(q[i,j,k] == fh[i,j,k]*(h[i,k]-a))
            model.addCons(q[i,j,k] == fc[i,j,k]*(b-c[j,k+1]))
            model.addCons(q[i,j,k]*(h[i,k]-b) >= 0)
            model.addCons(q[i,j,k]*(a-c[j,k+1]) >= 0)
        for i in range(n):
            for k in range(ell):
                model.addCons(quicksum(fh[i,j,k] for j in range(m)) == hot[i][0])
        for j in range(m):
            for k in range(ell):
                model.addCons(quicksum(fc[i,j,k] for i in range(n)) == cold[j][0])
        # Stage heat balances plus branch heat equations and split sums imply mixing.
    else:
        raise ValueError(version)
    model.setObjective(quicksum(q.values()), 'maximize')
    model.optimize()
    result = {'hot': hot, 'cold': cold, 'stages': ell, 'version': version,
              'status': str(model.getStatus()), 'seconds': model.getSolvingTime(),
              'primal': model.getPrimalbound(), 'dual': model.getDualbound(),
              'nodes': model.getNNodes()}
    if model.getNSols():
        result['h'] = [[model.getVal(h[i,k]) for k in range(ell+1)] for i in range(n)]
        result['c'] = [[model.getVal(c[j,k]) for k in range(ell+1)] for j in range(m)]
        result['q'] = [[[model.getVal(q[i,j,k]) for j in range(m)] for i in range(n)] for k in range(ell)]
    return result


def main():
    cases = []
    # Small NO and YES instances for the boundary reduction.
    for items,B,label in [([4,4,4,4,4,6],13,'no'), ([4,4,5,4,4,5],13,'yes'), ([2,2,2],3,'partition_no')]:
        bins = sum(items)//B
        hot = [(a,bins,0) for a in items]
        cold = [(B,j,bins) for j in range(bins)]
        for ell,version in [(1,'boundary'),(1,'countercurrent'),(2,'boundary')]:
            result = solve(hot,cold,ell,version,10)
            result['label'] = label
            full = B*bins*(bins+1)/2
            expected_full = version == 'countercurrent' or ell >= bins or label == 'yes'
            ok = result['status']=='optimal' and ((abs(result['primal']-full)<1e-5) if expected_full else result['dual'] < full-1e-5)
            result['numerical_expectation_passed'] = ok
            print('PASS' if ok else 'FAIL', label, ell, version, result['status'], result['primal'], result['dual'], flush=True)
            cases.append(result)
    hot,cold = [(1,10,0)],[(1,8,10),(1,0,8)]
    for version in ['boundary','countercurrent']:
        for ell in [1,2]:
            result = solve(hot,cold,ell,version,10)
            result['label'] = 'stage_gap'
            expected = 10 if ell==2 else (8 if version=='boundary' else 8.4)
            ok = result['status']=='optimal' and abs(result['primal']-expected)<1e-5
            result['numerical_expectation_passed'] = ok
            print('PASS' if ok else 'FAIL', 'stage_gap',ell,version,result['status'],result['primal'],result['dual'],flush=True)
            cases.append(result)
    # Zero quadratic-loss construction: all cold pairs have complementary ranges.
    for items,B,label in [([2,2,2],3,'potential_no'),([1,2,1,2],3,'potential_yes'),
                           ([4,4,4,4,4,6],13,'potential_3partition_no')]:
        bins = sum(items)//B
        hi = bins+1
        hot = [(a,hi,0) for a in items]
        cold = [(B,j,hi) for j in range(1,bins+1)] + [(B,0,j) for j in range(1,bins+1)]
        for version in ['boundary','countercurrent']:
            result = solve(hot,cold,2,version,20)
            result['label'] = label
            full = sum(items)*hi
            ok = (result['status']=='optimal' and abs(result['primal']-full)<1e-5) if label=='potential_yes' else result['dual'] < full-1e-5
            result['numerical_expectation_passed'] = ok
            print('PASS' if ok else 'FAIL', label,2,version,result['status'],result['primal'],result['dual'],flush=True)
            cases.append(result)
    out = ROOT/'certificates/physical-exploration.json'
    out.write_text(json.dumps({'WARNING':'Floating point only; not a mathematical certificate','cases':cases},indent=2)+'\n')
    return int(any(not x['numerical_expectation_passed'] for x in cases))


if __name__ == '__main__':
    raise SystemExit(main())
