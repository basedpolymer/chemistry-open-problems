"""Presentation and comparison (NOT part of the proof).

1. Writes certificates/table_grid.tex (Table of the manuscript) from certificates/grid_results.csv.
2. Compares the certified enclosures with the floating-point values published by the source in the
   data files of its arXiv source (arXiv:2102.06751v2): BDoscillate-theta-kappa.csv (eta0 = theta0,
   eta2 = theta2, kappa0, kappa2) and BDoscillate-resonance-transversality.csv (|L(2)|, Re ghat_1(1)),
   for the nine values nu = alpha/(1-alpha), alpha = 0.1..0.9, beta = 0, 0.1, ..., 10.
   The files are read from certificates/sources/ (git-ignored); if absent they are downloaded from
   arXiv (only the two csv files are extracted).
Prints the largest relative deviation per quantity and PASS if it is below 1e-6 (the source's values
come from adaptive floating-point quadrature).

Runtime: a few seconds.  Environment: Python 3 (standard library only).
"""
import csv, io, os, sys, tarfile, urllib.request
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
CERT = os.path.join(HERE, '..', 'certificates')
SRC = os.path.join(CERT, 'sources')
FILES = ['BDoscillate-theta-kappa.csv', 'BDoscillate-resonance-transversality.csv', 'BDoscillate-odd_beta3.csv']


def mid(s):
    """midpoint of an Arb string '[m +/- r]' or a plain number"""
    s = s.strip()
    if s.startswith('['):
        s = s[1:-1]
        return float(s.split('+/-')[0]) if '+/-' in s else float(s)
    return float(s)


def rad(s):
    s = s.strip()
    if s.startswith('[') and '+/-' in s:
        return float(s[1:-1].split('+/-')[1])
    return 0.0


def load_results():
    rows = {}
    with open(os.path.join(CERT, 'grid_results.csv')) as fh:
        rd = csv.DictReader(fh, delimiter=';')
        for r in rd:
            rows[(Fraction(r['beta']), Fraction(r['nu']))] = r
    return rows


def write_table(rows):
    nus = [Fraction(1, 9), Fraction(1, 4), Fraction(3, 7), Fraction(1, 2), Fraction(2, 3), Fraction(1),
           Fraction(3, 2), Fraction(2), Fraction(7, 3), Fraction(4), Fraction(9)]
    betas = [Fraction(0), Fraction(1), Fraction(2), Fraction(5), Fraction(10)]
    out = []
    out.append('\\begin{table}[ht]\n\\centering\\scriptsize\n\\setlength{\\tabcolsep}{3pt}\n')
    out.append('\\begin{tabular}{l' + 'rr' * len(betas) + '}\n\\toprule\n')
    out.append(' & ' + ' & '.join('\\multicolumn{2}{c}{$\\beta=%s$}' % (str(b) if b.denominator == 1 else '%d/%d' % (b.numerator, b.denominator)) for b in betas) + '\\\\\n')
    out.append('$\\nu$ & ' + ' & '.join(['$\\eta_0$ & $\\eta_2$'] * len(betas)) + '\\\\\n\\midrule\n')
    for n in nus:
        cells = []
        for b in betas:
            r = rows[(b, n)]
            e0, e2 = mid(r['eta0']), mid(r['eta2'])
            if e0 >= 1e-3:
                cells.append('$%.4g$ & $%.4f$' % (e0, e2))
            else:
                mant, ex = ('%.3e' % e0).split('e')
                cells.append('$%s\\cdot10^{%d}$ & $%.4f$' % (mant, int(ex), e2))
        ns = str(n) if n.denominator == 1 else '%d/%d' % (n.numerator, n.denominator)
        out.append('$%s$ & ' % ns + ' & '.join(cells) + '\\\\\n')
    out.append('\\bottomrule\n\\end{tabular}\n')
    out.append('\\caption{Certified critical removal rates $\\eta_0$ and direction coefficients $\\eta_2$ '
               '(midpoints of the enclosures, rounded; the full enclosures for all 1111 points are in '
               '\\texttt{certificates/grid\\_results.csv}). Every $\\eta_2$ is certified negative.}\\label{tab:grid}\n'
               '\\end{table}\n')
    with open(os.path.join(CERT, 'table_grid.tex'), 'w') as fh:
        fh.write(''.join(out))
    print('wrote certificates/table_grid.tex')


def get_source():
    if all(os.path.exists(os.path.join(SRC, f)) for f in FILES):
        return True
    try:
        os.makedirs(SRC, exist_ok=True)
        data = urllib.request.urlopen('https://arxiv.org/e-print/2102.06751v2', timeout=60).read()
        with tarfile.open(fileobj=io.BytesIO(data), mode='r:gz') as tf:
            for f in FILES:
                m = tf.getmember(f)
                with open(os.path.join(SRC, f), 'wb') as out:
                    out.write(tf.extractfile(m).read())
        return True
    except Exception as e:
        print('source data not available:', e)
        return False


def main():
    rows = load_results()
    write_table(rows)
    if not get_source():
        print('SKIP comparison with the source')
        return
    alphas = [Fraction(j, 10) for j in range(1, 10)]
    nus = [a / (1 - a) for a in alphas]
    tk = list(csv.DictReader(open(os.path.join(SRC, FILES[0])), delimiter='\t'))
    rt = list(csv.DictReader(open(os.path.join(SRC, FILES[1])), delimiter='\t'))
    worst = {}
    n = 0
    for rec in tk:
        b = Fraction(rec['beta']).limit_denominator(10)
        for j, nu in enumerate(nus, start=1):
            r = rows.get((b, nu))
            if r is None:
                continue
            n += 1
            for (ours, theirs) in [('eta0', 'theta0_%d' % j), ('eta2', 'theta2_%d' % j), ('kappa0', 'kappa0_%d' % j), ('kappa2', 'kappa2_%d' % j)]:
                a, s = mid(r[ours]), float(rec[theirs])
                d = abs(a - s) / max(abs(a), 1e-300)
                if d > worst.get(ours, (0,))[0]:
                    worst[ours] = (d, str(b), str(nu))
    for rec in rt:
        b = Fraction(rec['beta']).limit_denominator(10)
        for j, nu in enumerate(nus, start=1):
            r = rows.get((b, nu))
            if r is None:
                continue
            for (ours, theirs) in [('abs_L2', 'abshatL2%d' % j), ('Re_ghat1_1', 'rehatg11%d' % j)]:
                a, s = mid(r[ours]), float(rec[theirs])
                d = abs(a - s) / max(abs(a), 1e-300)
                if d > worst.get(ours, (0,))[0]:
                    worst[ours] = (d, str(b), str(nu))
    print('compared %d (beta, nu) points with the source data files' % n)
    ok = True
    for k, (d, b, nu) in sorted(worst.items()):
        print('  %-10s max relative deviation %.2e at beta=%s nu=%s' % (k, d, b, nu))
        ok &= d < 1e-6
    print(('PASS' if ok else 'NOTE'), 'certified values agree with the source floating-point values to 1e-6')
    # Table 1 of the source (nu = 0, beta = 3), file BDoscillate-odd_beta3.csv
    oddf = os.path.join(SRC, 'BDoscillate-odd_beta3.csv')
    ours_path = os.path.join(CERT, 'odd_beta3.csv')
    if os.path.exists(oddf) and os.path.exists(ours_path):
        src = list(csv.DictReader(open(oddf), delimiter='\t'))
        ours = list(csv.DictReader(open(ours_path), delimiter=';'))
        print('Table 1 of the source (nu=0, beta=3), %d certified crossings:' % len(ours))
        dev = {'eta0': 0, 'eta2': 0, 'kappa0': 0, 'vartheta0': 0, 'Re_ghat1_1': 0}
        for o, s_ in zip(ours, src):
            for k, sk in [('eta0', 'theta0'), ('eta2', 'theta2'), ('kappa0', 'kappa0'), ('vartheta0', 'vartheta0'), ('Re_ghat1_1', 'rehatg11')]:
                a, b = mid(o[k]), float(s_[sk])
                dev[k] = max(dev[k], abs(a - b) / abs(a))
        for k, d in dev.items():
            print('  %-10s max relative deviation %.2e' % (k, d))
        print('  kappa2: certified value vs source column vs source theta2/4 (rows 1-4):')
        for o, s_ in list(zip(ours, src))[:4]:
            print('    %s  %s  %.10f  %.10f' % (o['j'], o['kappa2'], float(s_['kappa2']), float(s_['theta2']) / 4))
        same = all(abs(float(s_['kappa2']) - float(s_['theta2']) / 4) < 1e-12 for s_ in src)
        print('  NOTE the source column kappa2 equals theta2/(beta+1) = theta2/4 in every row: %s;'
              ' formula (3.23) gives the certified values instead' % same)


if __name__ == '__main__':
    main()
