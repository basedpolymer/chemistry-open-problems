"""Validation on the 8 five-component cases of Tumbalam Gooty, Agrawal & Tawarmalani (arXiv:2010.12113,
Table of duality gaps): minimise the total reboiler vapor of each of the 6128 configurations (strict model,
multistart ipopt, 2 starts; 6 more starts when within 0.2% of the FTC value) and compare the minimum over
configurations with the published optimum and with the FTC closed form.

Environment: WSL Fedora, Python 3, pyomo 6.10 + ipopt 3.14, multiprocessing (W workers).  Floating point.
Runtime: ~1 h with 4 workers.
Usage: python tg_cases_N5.py outfile.json [workers] [cases]
"""
import json
import multiprocessing as mpc
import sys
import time

from hsmodel import all_configurations, petlyuk_value

CASES = {'A': ([20, 30, 20, 20, 10], [4.1, 3.6, 2.1, 1.42, 1], 402.7),
         'B': ([5, 15, 25, 20, 35], [7.98, 3.99, 3, 1.25, 1], 272.5),
         'C': ([25, 10, 25, 20, 20], [13.72, 3.92, 3.267, 1.21, 1], 260.0),
         'D': ([42.5, 42.5, 5, 5, 5], [3.3275, 3.025, 1.21, 1.1, 1], 896.4),
         'E': ([30, 30, 5, 5, 30], [1.4641, 1.331, 1.21, 1.1, 1], 695.6),
         'F': ([5, 5, 5, 42.5, 42.5], [3.3275, 1.331, 1.21, 1.1, 1], 929.1),
         'G': ([5, 5, 5, 42.5, 42.5], [17.1875, 6.875, 2.75, 1.1, 1], 902.7),
         'H': ([20, 20, 20, 20, 20], [7.5625, 3.025, 1.21, 1.1, 1], 542.0)}
_C = None


def _init():
    global _C
    _C = all_configurations(5)


def work(args):
    from nlp import solve_config
    case, cid = args
    F, al, ref = CASES[case]
    F = [x / 100 for x in F]  # normalise to F = 1 (vapor scales linearly)
    vt, vp = petlyuk_value(al, F, 1.0)
    vp = float(vp)
    best, _ = solve_config(_C[cid], al, F, 1.0, starts=2, seed=cid)
    r = best[0] / vp if best else None
    if r is not None and r < 1.002:
        b2, _ = solve_config(_C[cid], al, F, 1.0, starts=6, seed=cid + 99991, vscale=3.0)
        if b2 and b2[0] / vp < r:
            r = b2[0] / vp
    return case, cid, r


if __name__ == '__main__':
    out = sys.argv[1]
    W = int(sys.argv[2]) if len(sys.argv) > 2 else 4
    cases = list(sys.argv[3]) if len(sys.argv) > 3 else list(CASES)
    _init()
    tasks = [(k, c['id']) for k in cases for c in _C]
    res = {k: {} for k in cases}
    t0 = time.time()
    with mpc.Pool(W, initializer=_init) as pool:
        for n, (k, cid, r) in enumerate(pool.imap_unordered(work, tasks, chunksize=8)):
            res[k][cid] = r
            if n % 2000 == 0:
                print(n, round(time.time() - t0), flush=True)
    summary = {}
    for k in cases:
        F, al, ref = CASES[k]
        vals = {c: v for c, v in res[k].items() if v is not None}
        mn = min(vals.values())
        ftc = [c['id'] for c in _C if c['streams'] and len(c['streams']) == 10 and all(
            v == 'TC' for v in c['hx'].values())][0]
        vt, vp = petlyuk_value(al, [x / 100 for x in F], 1.0)
        summary[k] = {'published_optimum': ref, 'FTC_closed_form_x100': float(vp) * 100,
                      'min_ratio_over_configs': mn, 'min_value_x100': mn * float(vp) * 100,
                      'FTC_ratio': vals.get(ftc), 'n_ties_1e-6': sum(1 for v in vals.values() if v < 1 + 1e-6),
                      'n_failed': len(res[k]) - len(vals)}
        print(k, summary[k], flush=True)
    json.dump({'summary': summary, 'ratios': {k: {str(c): v for c, v in res[k].items()} for k in cases}},
              open(out, 'w'))
