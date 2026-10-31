#!/usr/bin/env python3
"""Write figure.svg: the three kinds of extremal graph (C_9, the double star D_{2,3}, the tree T_{1,2,2}).
Python 3 standard library. The Wiener indices printed in the figure are recomputed here, and the script
prints PASS or FAIL for them."""
import math
import os
import sys
from collections import deque

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "figure.svg")


def wiener(n, edges):
    adj = [[] for _ in range(n)]
    for a, b in edges:
        adj[a].append(b)
        adj[b].append(a)
    tot = 0
    for s in range(n):
        dist = [-1] * n
        dist[s] = 0
        q = deque([s])
        while q:
            u = q.popleft()
            for v in adj[u]:
                if dist[v] < 0:
                    dist[v] = dist[u] + 1
                    q.append(v)
        tot += sum(dist)
    return tot // 2


# C_9 at centre (130,120)
c9 = [(130 + 80 * math.sin(2 * math.pi * i / 9), 120 - 80 * math.cos(2 * math.pi * i / 9)) for i in range(9)]
c9e = [(i, (i + 1) % 9) for i in range(9)]
# double star D_{2,3}: centres 0,1; leaves 2,3 on 0; 4,5,6 on 1
ds = [(340, 120), (420, 120), (290, 75), (290, 165), (480, 60), (495, 120), (480, 180)]
dse = [(0, 1), (0, 2), (0, 3), (1, 4), (1, 5), (1, 6)]
# T_{1,2,2}: root 0; children 1,2,3; leaves 4 | 5,6 | 7,8
tr = [(690, 50), (610, 120), (690, 120), (770, 120), (610, 190), (672, 190), (708, 190), (752, 190), (788, 190)]
tre = [(0, 1), (0, 2), (0, 3), (1, 4), (2, 5), (2, 6), (3, 7), (3, 8)]

checks = [(wiener(9, c9e), 90, "C_9"), (wiener(7, dse), 42, "D_{2,3}"), (wiener(9, tre), 90, "T_{1,2,2}")]
ok = all(a == b for a, b, _ in checks)
for a, b, name in checks:
    print("PASS" if a == b else "FAIL", name, "W =", a)

parts = ['<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 840 270" width="840" height="270" '
         'font-family="Helvetica,Arial,sans-serif">',
         '<rect width="840" height="270" fill="#ffffff"/>']


def draw(pts, edges):
    for a, b in edges:
        parts.append(f'<line x1="{pts[a][0]:.1f}" y1="{pts[a][1]:.1f}" x2="{pts[b][0]:.1f}" y2="{pts[b][1]:.1f}" '
                     'stroke="#333333" stroke-width="2"/>')
    for x, y in pts:
        parts.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="6" fill="#1f5f8b"/>')


draw(c9, c9e)
draw(ds, dse)
draw(tr, tre)
for x, t1, t2 in [(130, "cycle C9 (any d: here d = 4)", "W = 90"),
                  (395, "double star, d = 3, n = 7", "W = 42 = W(C7)"),
                  (690, "tree, d = 4, n = 9", "W = 90 = W(C9)")]:
    parts.append(f'<text x="{x}" y="235" font-size="15" text-anchor="middle" fill="#222222">{t1}</text>')
    parts.append(f'<text x="{x}" y="256" font-size="15" text-anchor="middle" fill="#222222">{t2}</text>')
parts.append("</svg>")
with open(OUT, "w", encoding="utf-8") as f:
    f.write("\n".join(parts) + "\n")
sys.exit(0 if ok else 1)
