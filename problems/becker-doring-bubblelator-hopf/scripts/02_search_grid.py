"""Exploration (NOT part of the proof): locate, in floating point, the first zero t0 of
H(t) = Re G_{beta,nu}(it) for every grid point, and write a tiny interval [t_lo, t_hi] around it
to certificates/grid_tintervals.csv.  The rigorous checks are done by 03_certify_grid.py.

Grid: nu in {1/9, 1/4, 3/7, 2/3, 1, 3/2, 7/3, 4, 9} (the values of Niethammer-Pego-Schlichting-
Velazquez, alpha = 0.1, ..., 0.9) plus nu in {1/2, 2} (alpha = 1/3, 2/3), beta = j/10, j = 0..100.
Runtime: about 10 minutes on 8 cores.  Environment: WSL Fedora, Python 3, python-flint 0.9.
"""
import sys, os, time
from multiprocessing import Pool
from fractions import Fraction
from flint import arb, ctx
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import bdhopf as B
import certify_lib as C

NUS = [Fraction(1, 9), Fraction(1, 4), Fraction(3, 7), Fraction(2, 3), Fraction(1), Fraction(3, 2),
       Fraction(7, 3), Fraction(4), Fraction(9), Fraction(1, 2), Fraction(2)]
BETAS = [Fraction(j, 10) for j in range(101)]


def q2arb(q):
    return arb(q.numerator) / q.denominator


def work(args):
    bq, nq = args
    ctx.prec = 128
    B.TOL_BITS = 70
    beta, nu = q2arb(bq), q2arb(nq)
    tz = C.find_first_zero(beta, nu)
    tz = C.refine_zero(beta, nu, tz)
    r = max(tz, 1.0) * 2.0 ** -40
    return (str(bq), str(nq), repr(tz - r), repr(tz + r))


if __name__ == '__main__':
    t0 = time.time()
    jobs = [(b, n) for n in NUS for b in BETAS]
    with Pool(int(os.environ.get('NPROC', '8'))) as p:
        res = p.map(work, jobs, chunksize=4)
    out = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'certificates', 'grid_tintervals.csv')
    with open(out, 'w') as f:
        f.write('beta,nu,t_lo,t_hi\n')
        for r in res:
            f.write(','.join(r) + '\n')
    print('wrote', len(res), 'rows in %.0f s' % (time.time() - t0))
