"""Rigorous certification of the Hopf hypotheses (H1)-(H4) and of the sign of the direction
coefficient eta_2 at every grid point listed in certificates/grid_tintervals.csv.

For each (beta, nu) (exact rationals, enclosed in Arb balls) and the stored interval
T = [t_lo, t_hi] the routine certify_lib.certify_T checks, in ball arithmetic:
  C1  H(t_lo) > 0 > H(t_hi)                          (a zero t0 of H(t) = Re G(it) in T)
  C2  H' < 0 on T' = [t_lo - d, t_hi + d]            (t0 unique in T', H'(t0) < 0: transversality)
  C3  Im G(it) < 0 on T                              (vartheta0 > 0: a genuine crossing)
  C4  scan of (0, G(0)/vartheta0] \\ T'               (hypothesis (H4), which implies the
                                                     nonresonance hypothesis (H3))
  C5  Re ghat_1(1) > 0, L(2) != 0, eta_2 != 0         (enclosures written to the csv)
Writes certificates/grid_results.csv and prints PASS/FAIL per point; exit code 1 on any FAIL.

Runtime: about 15 minutes on 8 cores (NPROC environment variable, default 8).
Environment: WSL Fedora, Python 3, python-flint 0.9.  Deterministic.
"""
import sys, os, time
from multiprocessing import Pool
from fractions import Fraction
from flint import arb, ctx
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import bdhopf as B
import certify_lib as C

HERE = os.path.dirname(os.path.abspath(__file__))
CERT = os.path.join(HERE, '..', 'certificates')


def q2arb(q):
    return arb(q.numerator) / q.denominator


def work(row):
    bs, ns, tlo, thi = row
    ctx.prec = 128
    B.TOL_BITS = 70
    beta, nu = q2arb(Fraction(bs)), q2arb(Fraction(ns))
    o = C.certify_T(beta, nu, float(tlo), float(thi))
    q = o['q']
    f = lambda x: x.str(12, radius=True)
    rec = dict(beta=bs, nu=ns, t_lo=tlo, t_hi=thi,
               C1=o['C1'], C2=o['C2'], C3=o['C3'], C4=o['C4'] and o.get('C4_origin', False), C5=o['C5'],
               eta0=f(q['eta0']), kappa0=f(q['kappa0']), vartheta0=f(q['vartheta0']),
               eta2=f(q['eta2']), kappa2=f(q['kappa2']), Re_ghat1_1=f(q['ghat1_1'].real),
               abs_L2=f(abs(q['L2'])), Tmax='%.6g' % o['Tmax'], scan_evals=o['C4_evals'],
               eta2_sign=('-' if q['eta2'] < 0 else ('+' if q['eta2'] > 0 else '?')))
    rec['ok'] = all([rec['C1'], rec['C2'], rec['C3'], rec['C4'], rec['C5']])
    return rec


def main():
    t0 = time.time()
    rows = []
    with open(os.path.join(CERT, 'grid_tintervals.csv')) as fh:
        next(fh)
        for line in fh:
            if line.strip():
                rows.append(tuple(line.strip().split(',')))
    with Pool(int(os.environ.get('NPROC', '8'))) as p:
        res = p.map(work, rows, chunksize=2)
    keys = ['beta', 'nu', 't_lo', 't_hi', 'C1', 'C2', 'C3', 'C4', 'C5', 'ok', 'eta0', 'kappa0', 'vartheta0',
            'eta2', 'kappa2', 'eta2_sign', 'Re_ghat1_1', 'abs_L2', 'Tmax', 'scan_evals']
    with open(os.path.join(CERT, 'grid_results.csv'), 'w') as fh:
        fh.write(';'.join(keys) + '\n')
        for r in res:
            fh.write(';'.join(str(r[k]) for k in keys) + '\n')
    nfail = 0
    for r in res:
        if not r['ok']:
            nfail += 1
            print('FAIL beta=%s nu=%s' % (r['beta'], r['nu']), {k: r[k] for k in ['C1', 'C2', 'C3', 'C4', 'C5']})
    nneg = sum(1 for r in res if r['eta2_sign'] == '-')
    print('PASS' if nfail == 0 else 'FAIL', 'Hopf hypotheses (H1)-(H4) at %d/%d grid points' % (len(res) - nfail, len(res)))
    print('PASS' if nneg == len(res) else 'FAIL', 'eta_2 < 0 (supercritical direction) at %d/%d grid points' % (nneg, len(res)))
    print('runtime %.0f s' % (time.time() - t0))
    sys.exit(0 if (nfail == 0 and nneg == len(res)) else 1)


if __name__ == '__main__':
    main()
