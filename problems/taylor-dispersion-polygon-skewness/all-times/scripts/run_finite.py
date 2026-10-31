#!/usr/bin/env python3
"""run_finite.py nmin nmax [tau0] [procs]: run finite_n.certify for a range of n in parallel; writes a JSON summary."""
import sys, os, json, time
from multiprocessing import Pool
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))


def one(args):
    n, tau0 = args
    import finite_n
    try:
        ok, C = finite_n.certify(n, tau0, verbose=False)
        return dict(n=n, ok=ok, M=C.Mm, Delta0=float(C.Delta0.upper()), rho=float(C.rho.upper()), L1=float(C.L1.upper()),
                    omega_L=float(C.omega_L.upper()), osc_L=float(C.osc_L.upper()), Lip_L=float(C.Lip_L.upper()),
                    delta_L=float(C.delta_L.upper()), U_L=float(C.U_L.upper()), worst=C.grid['worst'], grid=C.grid, short=C.short,
                    mu1=C.mu_f[0], slope=C.grid['slope'], e1=C.grid['e1'], params=finite_n.params(n))
    except Exception as e:
        return dict(n=n, ok=False, error=repr(e))


if __name__ == '__main__':
    nmin, nmax = int(sys.argv[1]), int(sys.argv[2])
    tau0 = float(sys.argv[3]) if len(sys.argv) > 3 else 0.015
    procs = int(sys.argv[4]) if len(sys.argv) > 4 else 3
    t0 = time.time()
    with Pool(procs) as P:
        res = []
        for r in P.imap(one, [(n, tau0) for n in range(nmin, nmax + 1)]):
            res.append(r)
            print({k: r[k] for k in r if k != 'grid'}, flush=True)
    out = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'certificates', f'finite_{nmin}_{nmax}.json')
    os.makedirs(os.path.dirname(out), exist_ok=True)
    json.dump(res, open(out, 'w'), indent=1)
    allok = all(r['ok'] for r in res)
    print(("PASS" if allok else "FAIL") + f": K > 0 on [{tau0}, inf) for {nmin} <= n <= {nmax} ({time.time()-t0:.0f}s)")
    sys.exit(0 if allok else 1)
