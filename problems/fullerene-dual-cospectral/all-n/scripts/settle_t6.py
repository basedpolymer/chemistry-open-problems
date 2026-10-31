#!/usr/bin/env python3
"""settle_t6.py N LOGDIR -- exact settlement of the T^6 hash collisions of the scan of C_N (dualspec, mask 2).
Prints for every key group the partition into exact T^6 characteristic polynomials, and writes LOGDIR/nN.settled.json.
"""
import collections
import json
import os
import subprocess
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fdc import *


def main():
    N = int(sys.argv[1])
    logdir = sys.argv[2]
    PC = os.path.expanduser('~/fdc/bin/pcselect')
    parent = {}

    def find(x):
        parent.setdefault(x, x)
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    summary = ''
    for line in open(os.path.join(logdir, 'n%d.txt' % N)):
        if line.startswith('C '):
            _, ty, i, j = line.split()
            parent[find(int(i))] = find(int(j))
        elif line.startswith('SUMMARY'):
            summary = line.strip()
    idx = sorted(parent)
    res = dict(n=N, summary=summary, indices=len(idx), groups=[], exact_cospectral_classes=[])
    if idx:
        with tempfile.NamedTemporaryFile('w', suffix='.idx', delete=False) as f:
            f.write(' '.join(map(str, idx)))
            name = f.name
        p1 = subprocess.Popen(['buckygen', '-q', '%dd' % N], stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
        p2 = subprocess.run([PC, name], stdin=p1.stdout, capture_output=True)
        p1.stdout.close()   # pcselect stops after the last selected graph
        p1.kill()
        p1.wait()
        os.unlink(name)
        graphs = dict(zip(idx, read_planar_code(p2.stdout)))
        assert len(graphs) == len(idx)
        groups = collections.defaultdict(list)
        for i in idx:
            groups[find(i)].append(i)
        for g in groups.values():
            byp = collections.defaultdict(list)
            polys = {}
            for i in g:
                rot = graphs[i]
                S = [v for v in range(len(rot)) if len(rot[v]) == 6]
                q = charpoly(rot, S)
                polys[i] = q
                byp[str(q)].append(i)
            res['groups'].append([len(v) for v in byp.values()])
            # first differing coefficient c_k between the first graph of the group and the others
            q0 = polys[g[0]]
            fd = []
            for i in g[1:]:
                d = q0 - polys[i]
                fd.append(None if d.degree() < 0 else q0.degree() - d.degree())
            res.setdefault('first_diff', []).append(fd)
            for k, v in byp.items():
                if len(v) > 1:
                    res['exact_cospectral_classes'].append(dict(indices=v, graph6=[to_graph6(graphs[i]) for i in v]))
    print('n', N, 'indices', len(idx), 'groups', res['groups'], 'first differing c_k', res.get('first_diff'), 'EXACT CLASSES', len(res['exact_cospectral_classes']), flush=True)
    json.dump(res, open(os.path.join(logdir, 'n%d.settled.json' % N), 'w'))


if __name__ == '__main__':
    main()
