"""Case definitions for run_case.py (Agent C).  Binary columns: ('B', eta) = top total f_p - eta (component p forced
down the bottom section into product p+1), ('T', eta) = top total f_p + eta (component p+1 forced up into product p).
eta_nb_default = Astra slack of every non-binary column (circulation through it); delta = R - max_j H_j."""
CASES = {}

# N = 3, network.md 11.3 data (alpha = (4, 2, 1), z = 1/3, q = 1, delta = 0.05, eta_A = eta_C = 0.003; the
# prefractionator carries half of the circulation, as in ftc3_forced.py)
CASES['n3_base'] = dict(alpha=[4, 2, 1], z=[1, 1, 1], q=1, delta=0.05, eta_nb_default=0.025,
                        eta_bin={(0, 1): ('B', 0.003), (1, 2): ('T', 0.003)}, n=[16, 24, 32, 48, 64, 96, 128])

# N = 4, geometric volatilities, equimolar saturated-liquid feed
CASES['n4_base'] = dict(alpha=[8, 4, 2, 1], z=[1, 1, 1, 1], q=1, delta=0.05, eta_nb_default=0.01,
                        eta_bin={(0, 1): ('B', 0.003), (1, 2): ('B', 0.003), (2, 3): ('T', 0.003)},
                        n=[16, 32, 64, 96])


def _pat(N, pat, eta):
    return {(p, p + 1): (pat[p], eta) for p in range(N - 1)}


# N = 3: all four forcing patterns on the base data
for _p in ['BB', 'BT', 'TB', 'TT']:
    CASES['n3_' + _p] = dict(alpha=[4, 2, 1], z=[1, 1, 1], q=1, delta=0.05, eta_nb_default=0.025,
                             eta_bin=_pat(3, _p, 0.003), n=[16, 32, 64, 128])
# N = 3: feed quality, volatilities, compositions, small delta
for _q in [-2, -0.5, 0, 0.5, 1.5, 3]:
    CASES[f'n3_q{_q}'] = dict(alpha=[4, 2, 1], z=[1, 1, 1], q=_q, delta=0.05, eta_nb_default=0.025,
                              eta_bin=_pat(3, 'BT', 0.003), n=[16, 32, 64, 128])
CASES['n3_close'] = dict(alpha=[1.21, 1.1, 1], z=[1, 1, 1], q=1, delta=0.05, eta_nb_default=0.02,
                         eta_bin=_pat(3, 'BT', 0.003), n=[32, 64, 128, 256, 384], max_dps=600)
CASES['n3_zmid'] = dict(alpha=[4, 2, 1], z=[0.1, 0.8, 0.1], q=1, delta=0.05, eta_nb_default=0.02,
                        eta_bin=_pat(3, 'BT', 0.003), n=[16, 32, 64, 128])
CASES['n3_zlight'] = dict(alpha=[4, 2, 1], z=[0.8, 0.1, 0.1], q=1, delta=0.05, eta_nb_default=0.02,
                          eta_bin=_pat(3, 'BT', 0.003), n=[16, 32, 64, 128])
CASES['n3_ztiny'] = dict(alpha=[4, 2, 1], z=[0.01, 0.01, 0.98], q=1, delta=0.05, eta_nb_default=0.005,
                         eta_bin=_pat(3, 'BT', 0.001), n=[64, 128, 256, 512, 768])
CASES['n3_dsmall'] = dict(alpha=[4, 2, 1], z=[1, 1, 1], q=1, delta=0.002, eta_nb_default=0.001,
                          eta_bin=_pat(3, 'BT', 0.0001), n=[32, 64, 128, 256])
CASES['n3_wide'] = dict(alpha=[100, 10, 1], z=[1, 1, 1], q=1, delta=0.05, eta_nb_default=0.025,
                        eta_bin=_pat(3, 'BT', 0.003), n=[8, 16, 32, 64], trace_depth=1.0)

# N = 4: the eight forcing patterns on the base data
for _p in ['BBB', 'BBT', 'BTB', 'BTT', 'TBB', 'TBT', 'TTB', 'TTT']:
    CASES['n4_' + _p] = dict(alpha=[8, 4, 2, 1], z=[1, 1, 1, 1], q=1, delta=0.05, eta_nb_default=0.01,
                             eta_bin=_pat(4, _p, 0.003), n=[16, 32, 64, 96])
CASES['n4_q0'] = dict(alpha=[6, 3.5, 1.8, 1], z=[0.3, 0.2, 0.15, 0.35], q=0, delta=0.03, eta_nb_default=0.005,
                      eta_bin=_pat(4, 'BBB', 0.002), n=[24, 48, 96, 144])
CASES['n4_q13'] = dict(alpha=[3, 2, 1.4, 1], z=[0.25, 0.3, 0.25, 0.2], q=1.3, delta=0.03, eta_nb_default=0.005,
                       eta_bin=_pat(4, 'BTB', 0.002), n=[32, 64, 128, 192])
# N = 5, 6 (default all-B pattern and a mixed one)
CASES['n5_base'] = dict(alpha=[16, 8, 4, 2, 1], z=[1, 1, 1, 1, 1], q=1, delta=0.05, eta_nb_default=0.005,
                        eta_bin=_pat(5, 'BBBB', 0.002), n=[16, 32, 64, 96])
CASES['n5_mixed'] = dict(alpha=[5, 3.6, 2.4, 1.6, 1], z=[0.15, 0.25, 0.2, 0.3, 0.1], q=0.6, delta=0.03,
                         eta_nb_default=0.003, eta_bin=_pat(5, 'BTBT', 0.002), n=[24, 48, 96, 144])
CASES['n6_base'] = dict(alpha=[32, 16, 8, 4, 2, 1], z=[1, 1, 1, 1, 1, 1], q=1, delta=0.05, eta_nb_default=0.003,
                        eta_bin=_pat(6, 'BBBBB', 0.002), n=[16, 32, 64])
CASES['n6_mixed'] = dict(alpha=[7, 5, 3.4, 2.3, 1.5, 1], z=[0.2, 0.1, 0.25, 0.15, 0.1, 0.2], q=0.8, delta=0.03,
                         eta_nb_default=0.002, eta_bin=_pat(6, 'TBTBT', 0.001), n=[24, 48, 96])


# Agent D's exported specifications (../runs/D_exports/*.json, made with run_D_export.py): finite-stage runs on D's
# own data and parameters.  n chosen from D's rho: largest n with rho^n ~ 1e-15.
def _load_D():
    import os, glob, json, math
    here = os.path.dirname(os.path.abspath(__file__))
    for f in sorted(glob.glob(os.path.join(here, '..', 'runs', 'D_exports', '*.json'))):
        try:
            D = json.load(open(f))
        except Exception:
            continue
        data, par = D['data'], D['params']
        N = int(data['N'])
        pat = par['pattern']
        def pc(s):
            a, b = s.strip('[]').split(',')
            return int(a) - 1, int(b) - 1
        rho = max(float(s['rho']) for s in D['sections'].values())
        nmax = int(math.ceil(15 / -math.log10(rho)))
        nmax = max(32, min(nmax, 600 if N <= 4 else 300))
        name = 'D_' + os.path.basename(f)[:-5]
        CASES[name] = dict(alpha=[str(x) for x in data['alpha']], z=[str(x) for x in data['z']], q=str(data['q']),
                           delta=str(data['delta']),
                           eta_nb={pc(k): str(v) for k, v in par['eta_nb'].items()},
                           eta_bin={(int(p) - 1, int(p)): (pat[int(p) - 1], str(v)) for p, v in par['eta_bin'].items()},
                           n=sorted(set([max(16, nmax // 4), max(24, nmax // 2), nmax])), max_dps=1200)


_load_D()


# red team: equal highest peaks (equal_peaks.py), extreme feed qualities, large forcing
_z3eq = ['0.5258692647380476303989228779107740057576', '0.2370653676309761848005385610446129971212',
         '0.2370653676309761848005385610446129971212']
_z4eq = ['0.4177764105787172124321252210087958284968', '0.226764539753083108704436661477172369543',
         '0.1054590496681996788634381175140318019602', '0.25']
CASES['n3_equal_BT'] = dict(alpha=[4, 2, 1], z=_z3eq, q=1, delta=0.02, eta_nb_default=0.01,
                            eta_bin=_pat(3, 'BT', 0.002), n=[32, 64, 128, 192])
CASES['n3_equal_BB'] = dict(alpha=[4, 2, 1], z=_z3eq, q=1, delta=0.02, eta_nb_default=0.01,
                            eta_bin=_pat(3, 'BB', 0.002), n=[32, 64, 128, 192])
CASES['n4_equal_BBB'] = dict(alpha=[8, 4, 2, 1], z=_z4eq, q=1, delta=0.02, eta_nb_default=0.004,
                             eta_bin=_pat(4, 'BBB', 0.002), n=[32, 64, 128, 192])
CASES['n4_equal_TBT'] = dict(alpha=[8, 4, 2, 1], z=_z4eq, q=1, delta=0.02, eta_nb_default=0.004,
                             eta_bin=_pat(4, 'TBT', 0.002), n=[32, 64, 128, 192])
CASES['n3_q-10'] = dict(alpha=[4, 2, 1], z=[1, 1, 1], q=-10, delta=0.05, eta_nb_default=0.025,
                        eta_bin=_pat(3, 'BT', 0.003), n=[32, 64, 128, 256])
CASES['n3_q10'] = dict(alpha=[4, 2, 1], z=[1, 1, 1], q=10, delta=0.05, eta_nb_default=0.025,
                       eta_bin=_pat(3, 'BT', 0.003), n=[32, 64, 128, 256])
CASES['n4_q-5'] = dict(alpha=[8, 4, 2, 1], z=[1, 1, 1, 1], q=-5, delta=0.05, eta_nb_default=0.01,
                       eta_bin=_pat(4, 'BBT', 0.003), n=[32, 64, 128])
CASES['n4_q6'] = dict(alpha=[8, 4, 2, 1], z=[1, 1, 1, 1], q=6, delta=0.05, eta_nb_default=0.01,
                      eta_bin=_pat(4, 'BBT', 0.003), n=[32, 64, 128])
CASES['n3_bigeta'] = dict(alpha=[4, 2, 1], z=[1, 1, 1], q=1, delta=0.05, eta_nb_default=0.025,
                          eta_bin=_pat(3, 'BT', 0.1), n=[16, 32, 64, 128])
CASES['n4_bigeta'] = dict(alpha=[8, 4, 2, 1], z=[1, 1, 1, 1], q=1, delta=0.05, eta_nb_default=0.01,
                          eta_bin=_pat(4, 'BBT', 0.08), n=[16, 32, 64, 128])
CASES['n4_smallz2'] = dict(alpha=[8, 4, 2, 1], z=[0.3, 0.02, 0.38, 0.3], q=1, delta=0.05, eta_nb_default=0.002,
                           eta_bin=_pat(4, 'BBT', 0.001), n=[32, 64, 128, 192])
CASES['n4_close'] = dict(alpha=[1.6, 1.4, 1.2, 1], z=[1, 1, 1, 1], q=1, delta=0.05, eta_nb_default=0.005,
                         eta_bin=_pat(4, 'BTB', 0.002), n=[64, 128, 256, 384], max_dps=500)
CASES['n7_base'] = dict(alpha=[64, 32, 16, 8, 4, 2, 1], z=[1] * 7, q=1, delta=0.05, eta_nb_default=0.002,
                        eta_bin=_pat(7, 'BBBBBB', 0.002), n=[16, 32])


# red team: OUTSIDE D's parameter region (prefractionator slack too large for the margin delta: kappa_2 < 0)
CASES['n3_invalid'] = dict(alpha=[4, 2, 1], z=[1, 1, 1], q=1, delta=0.002, eta_nb_default=0.05,
                           eta_bin=_pat(3, 'BT', 0.003), n=[32, 64, 128])
CASES['n8_base'] = dict(alpha=[128, 64, 32, 16, 8, 4, 2, 1], z=[1] * 8, q=1, delta=0.05, eta_nb_default=0.001,
                        eta_bin=_pat(8, 'BBBBBBB', 0.002), n=[16])
# red team: extreme data (cheap even for huge n: the chart equations are exact for any n)
CASES['n3_veryclose'] = dict(alpha=[1.05, 1.02, 1], z=[1, 1, 1], q=1, delta=0.05, eta_nb_default=0.01,
                             eta_bin=_pat(3, 'BT', 0.002), n=[500, 1000, 2000, 4000], max_dps=400)
CASES['n3_zB_1e-4'] = dict(alpha=[4, 2, 1], z=['0.4999', '0.0002', '0.4999'], q=1, delta=0.02, eta_nb_default=0.002,
                           eta_bin=_pat(3, 'BT', '1e-5'), n=[64, 128, 256, 512])
CASES['n3_zA_1e-4'] = dict(alpha=[4, 2, 1], z=['0.0001', '0.5', '0.4999'], q=1, delta=0.02, eta_nb_default=0.002,
                           eta_bin=_pat(3, 'BT', '1e-5'), n=[64, 128, 256, 512])
