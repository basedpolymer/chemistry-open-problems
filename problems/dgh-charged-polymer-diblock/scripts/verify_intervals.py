#!/usr/bin/env python3
"""Certified threshold and finite-volume DP (including all omitted-tail errors).
Environment: Python >=3.10, python-flint, mpmath; WSL sci environment.
Runtime: about 15-20 seconds measured (verifier rerun 16 s), one thread, <250 MB RAM, no network.
Default verifies stored certificates; --write creates/refreshes numerical data.
Every claim prints PASS/FAIL; any failure exits nonzero. See paper.tex, Section 9, for the bounds (notes.md is not shipped).
"""
from pathlib import Path
import argparse
import json
import time
import mpmath as mp
import flint
from flint import arb, ctx

ROOT = Path(__file__).resolve().parents[1]
PARAM = json.loads((ROOT/'certificates'/'parameters.json').read_text())
ctx.prec = PARAM['precision_bits']
ctx.threads = 1
mp.mp.dps = 105


def check(ok, name):
    print(('PASS ' if ok else 'FAIL ') + name, flush=True)
    if not ok:
        raise AssertionError(name)


def series_bounds(beta, z, cutoff):
    lower = sum(((-beta*k*k/2).exp()*z**k for k in range(1,cutoff+1)), arb(0))
    ratio = z*(-beta*(2*cutoff+3)/2).exp()
    assert ratio < 1
    tail = (-beta*(cutoff+1)**2/2).exp()*z**(cutoff+1)/(1-ratio)
    return lower, lower+tail, tail


def root_certificate():
    lo, hi = arb(PARAM['beta_lower']), arb(PARAM['beta_upper'])
    L = PARAM['root_series_cutoff']
    low_sum = series_bounds(lo,arb(1),L)[0]
    high_sum = series_bounds(hi,arb(1),L)[1]
    check(low_sum > 1, 'series at beta_lower is strictly above one')
    check(high_sum < 1, 'series plus rigorous tail at beta_upper is below one')
    pi_point = 2*arb.pi()/9
    check(series_bounds(pi_point,arb(1),L)[0] > 1, 'beta_c is strictly greater than 2*pi/9')
    print('INFO beta interval:', lo.union(hi))
    print('INFO beta_c - 2*pi/9:', lo.union(hi)-pi_point)
    return lo.union(hi)


def guess_r(beta):
    # This floating-point operation proposes endpoints only. All endpoints are
    # subsequently accepted or rejected using Arb and the infinite-series bound.
    b = mp.mpf(beta)
    weights = [mp.exp(-b*k*k/2) for k in range(1,81)]
    def f(z):
        return mp.fsum(w*z**k for k,w in enumerate(weights,1))-1
    return mp.findroot(f, (mp.mpf('.5'),mp.exp(b/2)))


def scaling_bounds(label, beta):
    if label == 'critical_interval':
        blo,bhi = PARAM['beta_lower'],PARAM['beta_upper']
    else:
        blo=bhi=label
    # Arb's union has an outward-rounded radius. Leave a wider margin for the
    # critical parameter interval so it also covers that radius rounding.
    digits=45 if label=='critical_interval' else 60
    unit=mp.mpf(10)**digits
    sl=mp.floor(guess_r(blo)*unit)-1
    su=mp.ceil(guess_r(bhi)*unit)+1
    s=arb(str(int(sl)))/arb(10)**digits
    upper=arb(str(int(su)))/arb(10)**digits
    L=PARAM['dp_cutoff']
    check(series_bounds(beta,s,L)[1] < 1, f'{label}: certified lower bound on renewal radius')
    check(series_bounds(beta,upper,L)[0] > 1, f'{label}: certified upper bound on renewal radius')
    return s,upper


def full_dp(beta, n):
    charges=[1]*n+[-1]*n
    prefix=[0]
    for w in charges:
        prefix.append(prefix[-1]+w)
    z=[arb(1)]
    for j in range(1,2*n+1):
        z.append(sum((z[k]*(-beta*(prefix[j]-prefix[k])**2/2).exp()
                      for k in range(j)),arb(0)))
    return z[-1]


def finite_volumes(label,beta):
    s, rhi=scaling_bounds(label,beta)
    L,D=PARAM['dp_cutoff'],PARAM['diagonal_cutoff']
    sizes=PARAM['half_lengths']
    N=max(sizes)
    p=[arb(0)]+[(-beta*k*k/2).exp()*s**k for k in range(1,L+1)]
    eps=series_bounds(beta,s,L)[2].upper()
    truncated=[arb(1)]
    for n in range(1,N+1):
        truncated.append(sum((p[k]*truncated[n-k] for k in range(1,min(L,n)+1)),arb(0)))
    # Resolvent identity: 0 <= s^n(A_n-A_n^{<=L}) <= (n+1)*eps.
    # Do not feed these tail enclosures back into the truncated recurrence.
    A=[]
    scale=arb(1)
    for n,a in enumerate(truncated):
        err=arb(0) if n<=L else (n+1)*eps
        A.append((a+arb(0,err))*scale)
        scale/=s
    kernel=[(-beta*k*k/2).exp() for k in range(D+2)]
    pair_sum, prefix_sum=arb(0),arb(0)
    rows=[]
    for n in range(1,N+1):
        j=n-1
        pair_sum += A[j]*A[j]+2*A[j]*sum((kernel[d]*A[j-d] for d in range(1,min(D,j)+1)),arb(0))
        prefix_sum += A[j]
        if n in sizes:
            tail=kernel[D+1]*prefix_sum**2
            total=A[n]**2+pair_sum+arb(0,tail)
            check(total > 0 and total.rel_accuracy_bits() > 105,
                  f'{label}, N={2*n}: positive partition enclosure, >105 accuracy bits')
            if n==sizes[0]:
                direct=full_dp(beta,n)
                check(total.overlaps(direct),f'{label}, N={2*n}: agrees with untruncated full-chain DP')
            f=beta/2+total.log()/(2*n)
            p_mid=A[n]**2/total
            f_alt=total.log()/(2*n)-arb(2).log()+arb(2).log()/(2*n)
            if label=='critical_interval':
                f_limit=beta/2
            elif rhi < 1:
                f_limit=beta/2-s.union(rhi).log()
            else:
                assert s>1
                f_limit=beta/2
            row={'beta':label,'half_length':n,'F_original':f.str(40),
                 'F_alternative_at_beta_over_2':f_alt.str(40),
                 'p_middle':p_mid.str(40),'F_limit':f_limit.str(40),
                 'log_reduced_partition':total.log().str(40)}
            rows.append(row)
    return rows


def main(write=False):
    start=time.perf_counter()
    critical=root_certificate()
    rows=[]
    for label in PARAM['beta_grid']:
        beta=critical if label=='critical_interval' else arb(label)
        rows.extend(finite_volumes(label,beta))
        last=rows[-1]
        print(f"INFO {label}, N={2*last['half_length']}: F={last['F_original']}; p_mid={last['p_middle']}",flush=True)
    record={'python_flint':flint.__version__,'precision_bits':ctx.prec,'rows':rows}
    path=ROOT/'certificates'/'finite_volumes.json'
    if write:
        path.write_text(json.dumps(record,indent=2)+'\n')
        check(True,'wrote finite-volume interval data')
    else:
        old=json.loads(path.read_text())
        check(len(old['rows'])==len(rows),'stored interval row count')
        for r0,r1 in zip(old['rows'],rows):
            assert (r0['beta'],r0['half_length'])==(r1['beta'],r1['half_length'])
            for field in ['F_original','F_alternative_at_beta_over_2','p_middle','F_limit','log_reduced_partition']:
                assert arb(r0[field]).contains(arb(r1[field]))
        check(True,'all stored finite-volume intervals agree with certified recomputation')
    # These are finite-volume checks, not substitutes for the asymptotic proofs.
    crit=next(r for r in rows if r['beta']=='critical_interval' and r['half_length']==max(PARAM['half_lengths']))
    ratio=3*crit['half_length']*arb(crit['p_middle'])
    check(ratio > arb('.999') and ratio < arb('1.001'),'finite-volume check of 3*n*p_middle near one at criticality')
    print(f'PASS all interval claims; elapsed={time.perf_counter()-start:.2f}s',flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--write',action='store_true')
    args=parser.parse_args()
    try:
        main(args.write)
    except Exception as exc:
        print(f'FAIL {type(exc).__name__}: {exc}',flush=True)
        raise
