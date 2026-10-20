"""
test_section_sign.py -- numerical test of Agent N's lemma (N1) on MESH steady states.

Agent P (Claude Opus 5.5, effort xhigh).  For every adiabatic CMO section of a random network (randnet.py) at its
steady state (stagenet.py), with u(x; theta) = sum_i alpha_i x_i/(alpha_i - theta) for each liquid on the chain
[top node, stage 1, ..., stage n, bottom node]:
  (I)  Underwood's identity at every cross-section (liquid x_a above, vapor y_b below):
       Psi_s(theta) := V - sum alpha_i w_i/(alpha_i - theta) = L u(x_a) - theta V u(x_b)/S(x_b),
       with w = V y_b - L x_a (the same for every cross-section);
  (N1) along the chain read upward, the sign of u changes at most once: - to + if Psi_s > 0, + to - if
       Psi_s < 0, never if Psi_s = 0.
Tested at the feed roots theta_j and at random theta > 0 (not at a pole).
Run: python test_section_sign.py [ntrials] [seed]  (~1-3 min).  Prints PASS/FAIL; exit code 1 on failure.
"""
import sys, json, time
import numpy as np
from randnet import random_network, random_flows
from refvals import feed_roots

ntr = int(sys.argv[1]) if len(sys.argv) > 1 else 300
seed = int(sys.argv[2]) if len(sys.argv) > 2 else 2024
rng = np.random.default_rng(seed)


def u_of(x, th, a):
    return np.sum(a * x / (a - th), axis=-1)


t0 = time.time()
st = dict(networks=0, solved=0, sections=0, crossections=0, worst_identity=0.0, n1_viol=0, sign_changes=0,
          reverse_sections=0)
for tr in range(ntr):
    N = int(rng.integers(2, 6))
    net, a, z, q = random_network(rng, N)
    f = random_flows(rng, net)
    if f is None:
        continue
    st['networks'] += 1
    x, info = net.solve(f)
    if not info['ok'] or info['res'] > 1e-11:
        continue
    st['solved'] += 1
    y = net.yof(x)
    S = (x * a).sum(axis=1)
    ths = [float(t) for t in feed_roots(a, z, q)] + list(rng.uniform(0.05, 1.5 * a[0], 3))
    for s in net.sections:
        V = f[net.fidx['V:' + s['name']]]; L = f[net.fidx['L:' + s['name']]]
        chain = [s['top']] + list(range(s['first'], s['first'] + s['n'])) + [s['bot']]
        w = V * y[chain[1]] - L * x[chain[0]]
        sg = np.sign(w[np.abs(w) > 1e-10])
        if len(sg) > 1 and np.any(np.diff(sg) > 0):
            st['reverse_sections'] += 1
        st['sections'] += 1
        for th in ths:
            if np.min(np.abs(a - th) / a) < 1e-3:
                continue
            Psi = V - np.sum(a * w / (a - th))
            big = V + np.sum(np.abs(a * w / (a - th)))
            u = u_of(x[chain], th, a)
            for i in range(len(chain) - 1):
                aa, bb = chain[i], chain[i + 1]
                rhs = L * u[i] - th * V * u[i + 1] / S[bb]
                sc = 1 + big + abs(L * u[i]) + abs(th * V * u[i + 1] / S[bb])
                st['worst_identity'] = max(st['worst_identity'], abs(Psi - rhs) / sc)
                st['crossections'] += 1
            # upward reading: bottom node first
            uu = u[::-1]
            tol = 1e-12 * (1 + np.max(np.abs(uu)))
            sig = np.sign(np.where(np.abs(uu) > tol, uu, 0.0))
            nz = sig[sig != 0]
            ch = np.diff(nz)
            nch = int(np.sum(ch != 0))
            st['sign_changes'] += nch
            scale_psi = 1e-9 * (V + np.sum(np.abs(a * w / (a - th))))
            if nch > 1:
                st['n1_viol'] += 1
            elif nch == 1:
                up = ch[ch != 0][0] > 0      # - to + reading upward
                if (up and Psi < -scale_psi) or ((not up) and Psi > scale_psi):
                    st['n1_viol'] += 1
st['seconds'] = round(time.time() - t0, 1)
print(json.dumps(st))
ok1 = st['worst_identity'] < 1e-9
ok2 = st['n1_viol'] == 0
print('(I)  Underwood identity per cross-section :', 'PASS' if ok1 else 'FAIL', st['worst_identity'])
print('(N1) at most one sign change, right sense  :', 'PASS' if ok2 else 'FAIL', st['n1_viol'])
json.dump(st, open('../runs/test_section_sign_%d.json' % seed, 'w'), indent=1)
sys.exit(0 if (ok1 and ok2 and st['solved'] > 0) else 1)
