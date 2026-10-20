#!/usr/bin/env python3
"""Section/junction builder on top of stagenet.Net (Agent N).

A network = junction stages + sections.  Section s = chain of n_s equilibrium stages from junction b(s) (bottom
end) to junction t(s) (top end), vapor V_s up, liquid L_s down.  Junction j sends vapor V_s up every section with
b(s) = j and liquid L_s down every section with t(s) = j (a splitter when there are several), receives the liquid of
the sections with b(s) = j and the vapor of the sections with t(s) = j (mixing), may receive the feed, may deliver
liquid or vapor products, and carries the net duty Delta_j = V_out - V_in - Fv.  Total condenser = junction without
vapor outflow; partial reboiler = junction without liquid outflow except the product.
"""
import numpy as np
from stagenet import Net


def build(alpha, z, q, njunc, sections, feed_junc, products, feed_split=None):
    """sections: list of (b, t, nstages, V, L); products: list of (junction, phase 'L'/'V', flow, name).
    feed_split: optional list of (junction, fraction) to split the feed among several junctions."""
    net = Net(alpha, z, q)
    J = [net.add_stage(f"J{j}") for j in range(njunc)]
    sec_stages = []
    for si, (b, t, ns, V, L) in enumerate(sections):
        st = [net.add_stage(f"S{si}.{i}") for i in range(ns)]
        sec_stages.append(st)
        # vapor: junction b -> st[0] -> ... -> st[-1] -> junction t ; liquid reverse
        if ns == 0:
            net.add_vap(J[b], J[t], V)
            net.add_liq(J[t], J[b], L)
            continue
        net.add_vap(J[b], st[0], V)
        for i in range(ns - 1):
            net.add_vap(st[i], st[i + 1], V)
            net.add_liq(st[i + 1], st[i], L)
        net.add_vap(st[-1], J[t], V)
        net.add_liq(J[t], st[-1], L)
        net.add_liq(st[0], J[b], L)
    for (j, ph, f, name) in products:
        if ph == 'L':
            net.add_liq(J[j], ('P', name), f)
        else:
            net.add_vap(J[j], ('P', name), f)
    if feed_split is None:
        net.add_feed(J[feed_junc], 1.0)
    else:
        for j, fr in feed_split:
            net.add_feed(J[j], fr)
    return net, J, sec_stages
