"""class_pairs.py -- FLOATING-POINT EXPLORATION (not a certificate).

dump mode : read graph6 on stdin, print every gadget class (rooted graph B, root x, class of equal (rho, delta)) with
            rho > 0 as  graph6 x g m rho delta a.
pairs mode: read dumped classes; for every pair of classes of two gadgets with equal rho and equal c = delta - a,
            with a1 > 0 > a2, print the limiting proportion (|a2| g1 + a1 g2)/(|a2| m1 + a1 m2) of the mixture of the two
            gadget types at one hub (Section 5 of the manuscript), and the classes with a = 0 (proportion g/m).
Usage:  for n in 3 4 5 6 7 8; do geng -c -q $n | python3 class_pairs.py dump; done > classes_n8.txt
        python3 class_pairs.py pairs classes_n8.txt
Runtime: about 1 min.
"""
import sys
from collections import defaultdict


def dump():
    from kfnb import scan, parse
    batch = []

    def flush(b):
        n = ord(b[0][0]) - 63
        res = scan(parse(b, n), n, 0.0, 1e-7)
        for row in res:
            gi, x, g, r, d, a = int(row[0]), int(row[1]), int(row[2]), row[3], row[4], row[5]
            if r > 1e-9:
                print('%s %d %d %d %.10f %.10f %.10f' % (b[gi], x, g, n - 1, r, d, a))
    for line in sys.stdin:
        s = line.strip()
        if s:
            batch.append(s)
            if len(batch) >= 20000:
                flush(batch)
                batch = []
    if batch:
        flush(batch)


def pairs(files):
    rows = []
    for fn in files:
        for line in open(fn):
            s, x, g, m, r, d, a = line.split()
            rows.append((s, int(x), int(g), int(m), float(r), float(d), float(a)))
    key = defaultdict(list)
    for row in rows:
        key[(round(row[4], 7), round(row[5] - row[6], 7))].append(row)
    best = []
    for k, L in key.items():
        pos = [t for t in L if t[6] > 1e-7]
        neg = [t for t in L if t[6] < -1e-7]
        zer = [t for t in L if abs(t[6]) <= 1e-7]
        for p in pos:
            for q in neg:
                a1, a2 = p[6], -q[6]
                best.append(((a2 * p[2] + a1 * q[2]) / (a2 * p[3] + a1 * q[3]), p, q))
        for z in zer:
            best.append((z[2] / z[3], z, None))
    best.sort(key=lambda t: -t[0])
    for t in best[:25]:
        print('%.4f' % t[0], t[1], t[2])
    print(len(rows), 'classes;', len(best), 'combinations')


if __name__ == '__main__':
    if sys.argv[1] == 'dump':
        dump()
    else:
        pairs(sys.argv[2:])
