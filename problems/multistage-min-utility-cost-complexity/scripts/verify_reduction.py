#!/usr/bin/env python3
"""Exact corroboration of the proofs in paper.tex; no floating-point arithmetic.

Environment: Python 3.12+ and SymPy 1.14; tested in WSL scientific environment.
Runtime: typically <2 minutes. Run: python scripts/verify_reduction.py
Every claim prints PASS/FAIL; assertion failure gives a nonzero exit status.
Writes certificates/reduction-checks.json with exact rational witnesses.
The general reduction is proved on paper. Exhaustive finite checks alone do not
certify arbitrary nonlinear infeasibility or the complexity theorem.
"""
from fractions import Fraction as F
from itertools import combinations, combinations_with_replacement, product
import json
from pathlib import Path
import sympy as sp
from verify_model_audit import instance, blank, finish, violations, serialize

ROOT = Path(__file__).resolve().parents[1]
REPORT = []


def check(label, ok):
    print(('PASS ' if ok else 'FAIL ') + label, flush=True)
    REPORT.append({'claim': label, 'passed': bool(ok)})
    if not ok:
        raise AssertionError(label)


def gadget(items, B, ell=2, delta=F(0)):
    m = sum(items)//B
    assert sum(items) == m*B and ell >= 2
    H, p, A = m+1, ell-2, sum(items)
    hot = [(a, H+p+delta, delta) for a in items]
    cold = [(B,j,H) for j in range(1,m+1)] + [(B,0,j) for j in range(1,m+1)]
    cold += [(A,H+s-1,H+s) for s in range(1,p+1)]
    return instance(hot,cold,ell)


def partitions(items, B):
    """Independent subset-recursive exact bin partition, unordered bins."""
    def rec(indices):
        if not indices:
            return []
        first, rest = indices[0], indices[1:]
        for size in range(len(rest)+1):
            for sub in combinations(rest,size):
                group = (first,) + sub
                if sum(items[i] for i in group) == B:
                    other = tuple(i for i in rest if i not in sub)
                    tail = rec(other)
                    if tail is not None:
                        return [group] + tail
        return None
    return rec(tuple(range(len(items))))


def witness(items, B, groups, ell=2, delta=F(0)):
    inst = gadget(items,B,ell,delta)
    s = blank(inst)
    n, m, p = len(items), len(groups), ell-2
    H, nc = m+1, len(inst['cold'])
    labels = {i:j+1 for j,group in enumerate(groups) for i in group}
    for i in range(n):
        s['h'][i] = [F(H+p-k)+delta for k in range(p+1)] + [F(labels[i])+delta, delta]
    for j,(_,tin,tout) in enumerate(inst['cold']):
        active = p if j < m else (p+1 if j < 2*m else p-1-(j-2*m))
        s['c'][j] = [tout if k <= active else tin for k in range(ell+1)]
    for k in range(ell):
        for i,a in enumerate(items):
            j = 2*m+(p-1-k) if k < p else labels[i]-1+(m if k==p+1 else 0)
            s['fh'][k][i][j] = F(a)
            s['fc'][k][i][j] = F(a)
            s['q'][k][i][j] = a*(s['h'][i][k]-s['h'][i][k+1])
        for j,(cap,_,_) in enumerate(inst['cold']):
            if not any(s['q'][k][i][j] for i in range(n)):
                s['fc'][k][0][j] = cap  # idle cold flow, with zero temperature change
        for i in range(n):
            for j in range(nc):
                s['a'][k][i][j] = s['h'][i][k+1]
                s['b'][k][i][j] = s['c'][j][k]
    return inst, finish(inst,s)


def validate(inst,s,delta,version):
    err = violations(inst,s,'repaired')
    for k in range(inst['stages']):
        for i in range(len(inst['hot'])):
            for j in range(len(inst['cold'])):
                if s['q'][k][i][j] > 0:
                    if version == 'countercurrent':
                        ok = s['h'][i][k] >= s['b'][k][i][j]+delta and s['a'][k][i][j] >= s['c'][j][k+1]+delta
                    else:
                        ok = s['h'][i][k] >= s['c'][j][k]+delta and s['h'][i][k+1] >= s['c'][j][k+1]+delta
                    if not ok:
                        err.append(version+' approach')
    return err


def potential_identity(inst,s,delta):
    """Hot coordinates translated by -delta for these shifted constructions."""
    terms = []
    for k in range(inst['stages']):
        left = sum(fi*((s['h'][i][k]-delta)**2-(s['h'][i][k+1]-delta)**2)
                   for i,(fi,_,_) in enumerate(inst['hot']))
        left -= sum(fj*(s['c'][j][k]**2-s['c'][j][k+1]**2)
                    for j,(fj,_,_) in enumerate(inst['cold']))
        gap = vh = vc = F(0)
        for i in range(len(inst['hot'])):
            for j in range(len(inst['cold'])):
                q = s['q'][k][i][j]
                gap += q*(s['h'][i][k]-delta-s['b'][k][i][j]+s['a'][k][i][j]-delta-s['c'][j][k+1])
                vh += s['fh'][k][i][j]*(s['a'][k][i][j]-s['h'][i][k+1])**2
                vc += s['fc'][k][i][j]*(s['b'][k][i][j]-s['c'][j][k])**2
        assert left == gap+vh+vc
        terms.append([left,gap,vh,vc])
    return terms


def interval_enumeration(items,B):
    """Enumerate all equality-forced two-stage hot paths, without bin labels.

    Temperatures 0..H include idle and full-span candidates. Cold U/V paths
    are fixed by the equality lemma. Each hot interval must coincide with a
    positive cold interval; their capacity-weighted heat balances must agree.
    Completeness of this finite reduction follows from the paper, not SCIP.
    """
    m, H = sum(items)//B, sum(items)//B+1
    cold_intervals = [[(H,j) for j in range(1,m+1)],[(j,0) for j in range(1,m+1)]]
    count, first = 0, None
    for middle in product(range(H+1),repeat=len(items)):
        sums = [[0]*m for _ in range(2)]
        ok = True
        for a,t in zip(items,middle):
            for k,edge in enumerate(((H,t),(t,0))):
                if edge[0] == edge[1]:
                    continue
                if edge not in cold_intervals[k]:
                    ok = False
                    break
                j = cold_intervals[k].index(edge)
                sums[k][j] += a*(edge[0]-edge[1])
            if not ok:
                break
        demands = [[B*(u-v) for u,v in row] for row in cold_intervals]
        if ok and sums == demands:
            count += 1
            first = first or middle
    return count, first


def main():
    # Polynomial identities checked symbolically, with no division or sampling.
    fh,fc,h,a,b,c,q = sp.symbols('fh fc h a b c q')
    P = fh*(h*h-a*a)-fc*(b*b-c*c)-q*(h-b+a-c)
    R = (h+a)*(fh*(h-a)-q)-(b+c)*(fc*(b-c)-q)
    check('exchanger quadratic identity as a polynomial consequence of heat balances',sp.expand(P-R)==0)
    w1,w2,x1,x2,mu,cap = sp.symbols('w1 w2 x1 x2 mu cap')
    V = w1*(x1-mu)**2+w2*(x2-mu)**2-(w1*x1*x1+w2*x2*x2-cap*mu*mu)
    R = -2*mu*(w1*x1+w2*x2-cap*mu)+mu*mu*(w1+w2-cap)
    check('mixing variance identity as a polynomial consequence of mixing and split sums',sp.expand(V-R)==0)

    saved = {'symbolic_residuals':['0','0'],'witnesses':[],'enumerations':[]}
    examples = [([4,4,5,4,4,5],13),([4,4,4,4,4,6],13),([2,2,2],3),([1,2,1,2],3),([2,3,4,2,3,4,2,3,4],9)]
    for items,B in examples:
        groups = partitions(items,B)
        count,first = interval_enumeration(items,B)
        check(f'exact interval enumeration agrees with subset partition: {items}, B={B}',bool(count)==(groups is not None))
        saved['enumerations'].append({'items':items,'B':B,'partition':groups,'equality_solutions':count,'first_middle':first})
        for ell in [2,3,4,5]:
            inst = gadget(items,B,ell)
            energy = sum(f*(u-v) for f,u,v in inst['hot'])-sum(f*(v-u) for f,u,v in inst['cold'])
            square = sum(f*(u*u-v*v) for f,u,v in inst['hot'])-sum(f*(v*v-u*u) for f,u,v in inst['cold'])
            assert energy == square == 0
        if groups is not None:
            for ell in [2,3,4,5]:
                for delta in [F(0),F(1),F(3,2)]:
                    inst,s = witness(items,B,groups,ell,delta)
                    assert not validate(inst,s,delta,'countercurrent')
                    assert not validate(inst,s,delta,'boundary')
                    assert s['cost_unit_prices'] == 0
                    terms = potential_identity(inst,s,delta)
                    assert all(all(v==0 for v in row) for row in terms)
                    saved['witnesses'].append({'instance':inst,'delta':delta,'solution':s,'potential_terms':terms})
    check('all gadget energy and quadratic budgets balance for stages 2,3,4,5',True)
    check('36 complete rational zero-utility witnesses satisfy every equation and both approach versions',len(saved['witnesses'])==36)

    # Additional two-bin 3-PARTITION promise instances; exhaustive sorted multisets.
    checked = yes = no = 0
    for B in range(7,17):
        allowed = range(B//4+1,(B-1)//2+1)
        for items in combinations_with_replacement(allowed,6):
            if sum(items) != 2*B:
                continue
            groups = partitions(items,B)
            count,_ = interval_enumeration(items,B)
            assert bool(count) == (groups is not None)
            checked += 1
            yes += bool(count)
            no += not bool(count)
    check(f'exhaustive two-bin promise sweep: {checked} instances ({yes} YES, {no} NO)',checked>0 and yes>0 and no>0)
    saved['promise_sweep'] = {'B_range':[7,16],'instances':checked,'yes':yes,'no':no}

    # Exact suffix multiplier calculation used in the separate one-stage B proof.
    for m in range(2,13):
        lam = {k:F(m,(m-k)*(m-k+1)) for k in range(1,m)}
        assert all(v>0 for v in lam.values())
        assert all((m-r)*sum(lam[k] for k in range(1,r+1))==r for r in range(1,m))
        assert sum(lam[k]*sum(m-j for j in range(k,m)) for k in range(1,m)) == F(m*(m-1),2)
    check('one-stage boundary suffix identities for bin counts 2 through 12',True)
    saved['claims'] = REPORT
    (ROOT/'certificates/reduction-checks.json').write_text(json.dumps(serialize(saved),indent=2)+'\n',encoding='utf-8')


if __name__ == '__main__':
    main()
