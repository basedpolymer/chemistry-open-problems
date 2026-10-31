import sys, time
sys.path.insert(0, '..')
from hubgraph import *
for (p, r, b) in [(1, 2, 2), (2, 3, 3), (3, 4, 4)]:
    t0 = time.time()
    for s in range(p, 3000*p, p):
        c = construct(p, r, b, s)
        if c is not None:
            print(p, r, b, 'first s', s, 'x', c['x'], 'a', c['a'], 'n', c['n'], 'k', c['k'], 'l', c['l'], 't', c['t'], round(time.time()-t0, 2), flush=True)
            break
        if s % (200*p) == 0:
            # diagnostics
            B = b*(b-1); x, Q = q_list(b, s); a = (b-1)*s*Q[1]; n = a*b*r//p; m = n-a*b-1-a
            print('s', s, 'x', x, 'a', a, 'm', m, 'Sig', sum(Q), 'E', float(-delta_closed(a, b)), flush=True)
    else:
        print('none', p, r, b)
