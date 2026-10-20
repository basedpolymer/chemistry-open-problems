#!/usr/bin/env python3
"""
petlyuk_fixed_flows.py (Agent S) -- numerical exploration (evidence only) for the attainment question, N = 3.
Ternary Petlyuk arrangement (prefractionator + main column with liquid side draw, two thermal couplings), flows FIXED
at the Halvorsen-Skogestad point plus a relative margin delta:
  prefractionator at its preferred split (both feed roots active): top vapor V1T, distillate D1 = z_A + d_B;
  main column D = z_A, S = z_B, B = z_C, reboiler vapor V_B = max(H_1, H_2) (1 + delta);
  prefractionator vapor and liquid increased by delta * V1T at fixed D1.
For n = stages per section (all six sections), solves the MESH/CMO steady state (stagenet.py, continuation in n)
and prints the product purities and R / max_j H_j. Question tested: with fixed flows slightly above the HS point,
do the purities tend to 1 as n -> infinity (route "constructive with slack" / N's compactness route)?
Environment: Python 3, numpy, scipy. Runtime: about a minute.
Usage: python petlyuk_fixed_flows.py alphaA alphaB alphaC zA zB zC q delta
"""
import sys
import numpy as np
from scipy.optimize import brentq
from stagenet import Net
def petlyuk(N, ns, a, z, q, VB, V1, L1, S, D):
    """ns = stages per section (prefrac top, prefrac bottom, main top, mid-upper, mid-lower, main bottom)"""
    net = Net(a, z, q)
    B = 1 - D - S
    LB = VB + B
    Lml = LB - L1 - q
    Lmu = Lml + S
    LT = Lmu + L1
    VT = VB + 1 - q
    Vm = VB - V1
    flows = [VB, V1, L1, S, D, B, LB, Lml, Lmu, LT, VT, Vm, V1 + 1 - q, L1 + q]
    if min(flows) <= 1e-6:
        return None
    # prefractionator
    pt = [net.add_node('stage', f'pt{i}') for i in range(ns[0])]
    pf = net.add_node('stage', 'pf')
    pb = [net.add_node('stage', f'pb{i}') for i in range(ns[1])]
    # main column
    cond = net.add_node('mix', 'cond')
    mt = [net.add_node('stage', f'mt{i}') for i in range(ns[2])]
    J1 = net.add_node('stage', 'J1')
    mu = [net.add_node('stage', f'mu{i}') for i in range(ns[3])]
    Sd = net.add_node('stage', 'S')
    ml = [net.add_node('stage', f'ml{i}') for i in range(ns[4])]
    J2 = net.add_node('stage', 'J2')
    mb = [net.add_node('stage', f'mb{i}') for i in range(ns[5])]
    reb = net.add_node('stage', 'reb')

    def link(chain, L, V):
        for k in range(len(chain) - 1):
            net.add_stream(chain[k], chain[k + 1], 'L', L)
            net.add_stream(chain[k + 1], chain[k], 'V', V)

    net.add_stream('F', pf, 'F', 1.0)
    link([J1] + pt + [pf], L1, V1 + 1 - q)
    link([pf] + pb + [J2], L1 + q, V1)
    link([cond] + mt + [J1], LT, VT)
    link([J1] + mu + [Sd], Lmu, Vm)
    link([Sd] + ml + [J2], Lml, Vm)
    link([J2] + mb + [reb], LB, VB)
    net.add_stream(cond, ('P', 'D'), 'L', D)
    net.add_stream(Sd, ('P', 'S'), 'L', S)
    net.add_stream(reb, ('P', 'B'), 'L', B)
    return net




def hs_point(a, z, q):
    f = lambda t: np.sum(a * z / (a - t)) - (1 - q)
    th = [brentq(f, a[1] + 1e-12, a[0] - 1e-12), brentq(f, a[2] + 1e-12, a[1] - 1e-12)]
    H = [sum(a[p] * z[p] / (t - a[p]) for p in range(j + 1, 3)) for j, t in enumerate(th)]
    # preferred split: V1T = aA zA/(aA - t) + aB dB/(aB - t) for t = th0 and th1
    g = lambda t: a[0] * z[0] / (a[0] - t)
    hB = lambda t: a[1] / (a[1] - t)
    dB = (g(th[1]) - g(th[0])) / (hB(th[0]) - hB(th[1]))
    V1T = g(th[0]) + hB(th[0]) * dB
    return th, H, dB, V1T


def section_layout(n):
    """index lists of the node blocks of petlyuk() with n stages in every section"""
    blocks, k = [], 0
    for size in [n, 1, n, 1, n, 1, n, 1, n, 1, n, 1]:
        blocks.append(list(range(k, k + size)))
        k += size
    return blocks


def interp_guess(Cprev, nprev, n):
    bp, bn = section_layout(nprev), section_layout(n)
    G = np.zeros((sum(len(b) for b in bn), Cprev.shape[1]))
    for b_old, b_new in zip(bp, bn):
        m_old, m_new = len(b_old), len(b_new)
        for i, node in enumerate(b_new):
            pos = 0 if m_new == 1 else i * (m_old - 1) / (m_new - 1)
            lo = int(np.floor(pos)); hi = min(lo + 1, m_old - 1); w = pos - lo
            G[node] = (1 - w) * Cprev[b_old[lo]] + w * Cprev[b_old[hi]]
    return G


def main():
    args = [float(v) for v in sys.argv[1:9]] if len(sys.argv) >= 9 else [4, 2, 1, 1 / 3, 1 / 3, 1 / 3, 1.0, 0.02]
    a = np.array(args[0:3]); z = np.array(args[3:6]); z = z / z.sum(); q = args[6]; delta = args[7]
    extra = float(sys.argv[9]) if len(sys.argv) > 9 else 0.0
    th, H, dB, V1T = hs_point(a, z, q)
    print(f'theta = {th}, H = {H}, max H = {max(H):.6f}, preferred split d_B = {dB:.6f} (z_B = {z[1]:.6f}), '
          f'V1T = {V1T:.6f}')
    D1 = z[0] + dB
    V1T_d = V1T * (1 + delta)
    V1 = V1T_d - (1 - q)          # vapor sent from J2 to the prefractionator bottom
    L1 = V1T_d - D1               # liquid sent from J1 to the prefractionator top
    VB = max(H) * (1 + delta + extra)
    print(f'flows: VB = {VB:.6f}, V1 = {V1:.6f}, L1 = {L1:.6f}, D = {z[0]:.6f}, S = {z[1]:.6f}')
    Cprev, nprev = None, None
    for n in [4, 8, 16, 32, 64, 128]:
        net = petlyuk(3, [n] * 6, a, z, q, VB, V1, L1, z[1], z[0])
        if net is None:
            print('infeasible flows'); return
        C0 = None if Cprev is None else interp_guess(Cprev, nprev, n)
        C = net.solve(maxsweep=3000 if Cprev is None else 300, C0=C0)
        Cprev, nprev = net.C, n
        R, Q, nv = net.duties()
        prods = {st[1][1]: net.outcomp(st[0], net.C[st[0]], st[2]) for st in net.products()}
        print(f'n = {n:3d}: residual {net.res:.1e}  impurities: D {1 - prods["D"][0]:.2e}  S {1 - prods["S"][1]:.2e} '
              f'(A {prods["S"][0]:.2e}, C {prods["S"][2]:.2e})  B {1 - prods["B"][2]:.2e}   R/maxH = {R / max(H):.6f}',
              flush=True)


if __name__ == '__main__':
    main()
