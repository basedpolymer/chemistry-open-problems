#!/usr/bin/env python3
"""Write figure.svg: why five reactions stay connected and three two-way sources do not.

Left panel: the polynomials (x+a)(x-1)(x-2), which form one straight segment of rate constants of the
five-reaction network {A->0, A->2A, 2A->3A, 3A->2A, 4A->5A} after dividing by x.
Right panel: the only polynomials with roots 1 and 2 when there are three two-way sources.

Python 3 standard library, exact arithmetic for the claims, floats only for drawing. Prints PASS or
FAIL for the coefficient signs printed in the figure. Runtime under 1 s.
"""
import os
import sys
from fractions import Fraction as Q

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "figure.svg")
W, H = 820, 400
PANELS = ((64, 396), (464, 796))           # x pixel range of each plot
TOP, BOTTOM = 112, 340                     # y pixel range of each plot
XMIN, XMAX, YMIN, YMAX = 0.0, 3.2, -3.0, 3.0


def coeffs(a):
    """Coefficients of (x+a)(x-1)(x-2) = x^3 + (a-3)x^2 + (2-3a)x + 2a, lowest degree first."""
    return [2 * a, 2 - 3 * a, a - 3, Q(1)]


def signs(a):
    return "(" + ", ".join("0" if c == 0 else ("+" if c > 0 else "−") for c in coeffs(a)) + ")"


LEFT = [(Q(1, 2), "series-1", "a = 1/2"), (Q(0), "series-2", "a = 0"), (Q(-3), "series-3", "a = −3")]


def px(panel, x):
    lo, hi = PANELS[panel]
    return lo + (x - XMIN) / (XMAX - XMIN) * (hi - lo)


def py(y):
    return BOTTOM - (y - YMIN) / (YMAX - YMIN) * (BOTTOM - TOP)


def curve(panel, fn, cls):
    pts = []
    n = 320
    for i in range(n + 1):
        x = XMIN + (XMAX - XMIN) * i / n
        pts.append(f"{px(panel, x):.1f},{py(fn(x)):.1f}")
    return f'<polyline class="{cls}" clip-path="url(#clip{panel})" points="{" ".join(pts)}"/>'


def axes(panel):
    lo, hi = PANELS[panel]
    out = [f'<line class="axis" x1="{lo}" y1="{py(0):.1f}" x2="{hi}" y2="{py(0):.1f}"/>']
    for x in (0, 1, 2, 3):
        out.append(f'<line class="axis" x1="{px(panel, x):.1f}" y1="{py(0):.1f}" '
                   f'x2="{px(panel, x):.1f}" y2="{py(0) + 5:.1f}"/>')
        out.append(f'<text class="tick" x="{px(panel, x):.1f}" y="{py(0) + 18:.1f}" '
                   f'text-anchor="middle">{x}</text>')
    out.append(f'<text class="tick" x="{hi}" y="{py(0) + 18:.1f}" text-anchor="end">x</text>')
    return out


def roots(panel):
    return [f'<circle class="root" cx="{px(panel, x):.1f}" cy="{py(0):.1f}" r="5"/>' for x in (1, 2)]


def legend(panel, items):
    lo = PANELS[panel][0]
    out = []
    y = 72
    for cls, label in items:
        out.append(f'<line class="{cls}" x1="{lo}" y1="{y}" x2="{lo + 22}" y2="{y}"/>')
        out.append(f'<text class="label" x="{lo + 30}" y="{y + 4}">{label}</text>')
        y += 19
    return out


def main():
    ok = (signs(Q(1, 2)) == "(+, +, −, +)" and signs(Q(0)) == "(0, +, −, +)"
          and signs(Q(-3)) == "(−, +, −, +)")
    print(("PASS " if ok else "FAIL ") + "sign patterns printed in the left panel")
    ok2 = all(sum(c * Q(x) ** i for i, c in enumerate(coeffs(a))) == 0
              for a, _, _ in LEFT for x in (1, 2))
    print(("PASS " if ok2 else "FAIL ") + "every left-panel curve vanishes at 1 and 2")

    body = [f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" width="{W}" height="{H}" '
            'role="img" aria-labelledby="t d">',
            '<title id="t">Polynomials with roots 1 and 2 in two kinds of network</title>',
            '<desc id="d">Left: with five reactions, the polynomials (x+a)(x-1)(x-2) join a sign pattern '
            'to a different one along one straight line of rate constants. Right: with three two-way '
            'sources, only (x-1)(x-2) and its negative are available, and they meet only at the '
            'identically zero polynomial.</desc>',
            '<style>',
            ':root{--surface:#fcfcfb;--ink:#0b0b0b;--ink2:#52514e;--rule:#b9b8b2;'
            '--s1:#2a78d6;--s2:#eb6834;--s3:#1baf7a}',
            '@media (prefers-color-scheme: dark){:root{--surface:#1a1a19;--ink:#f0efec;--ink2:#c3c2b7;'
            '--rule:#5c5b56;--s1:#3987e5;--s2:#d95926;--s3:#199e70}}',
            'text{font-family:system-ui,-apple-system,"Segoe UI",Helvetica,Arial,sans-serif}',
            '.title{font-size:16px;font-weight:600;fill:var(--ink)}',
            '.label{font-size:13px;fill:var(--ink2)}.tick{font-size:12px;fill:var(--ink2)}',
            '.axis{stroke:var(--rule);stroke-width:1;fill:none}',
            '.series-1,.series-2,.series-3{fill:none;stroke-width:2;stroke-linecap:round;'
            'stroke-linejoin:round}',
            '.series-1{stroke:var(--s1)}.series-2{stroke:var(--s2)}.series-3{stroke:var(--s3)}',
            '.dash{stroke-dasharray:6 5}',
            '.root{fill:var(--ink);stroke:var(--surface);stroke-width:2}',
            '</style>',
            f'<rect width="{W}" height="{H}" fill="var(--surface)"/>']
    for i, (lo, hi) in enumerate(PANELS):
        body.append(f'<clipPath id="clip{i}"><rect x="{lo}" y="{TOP}" width="{hi - lo}" '
                    f'height="{BOTTOM - TOP}"/></clipPath>')

    body.append(f'<text class="title" x="{PANELS[0][0]}" y="30">Five reactions: one straight family</text>')
    body.append(f'<text class="label" x="{PANELS[0][0]}" y="50">(x + a)(x − 1)(x − 2), '
                'coefficient signs change with a</text>')
    body += legend(0, [(cls, f"{name}, signs {signs(a)}") for a, cls, name in LEFT])
    body += axes(0)
    for a, cls, _ in LEFT:
        av = float(a)
        body.append(curve(0, lambda x, av=av: (x + av) * (x - 1) * (x - 2), cls))
    body += roots(0)

    body.append(f'<text class="title" x="{PANELS[1][0]}" y="30">Three two-way sources: two rays</text>')
    body.append(f'<text class="label" x="{PANELS[1][0]}" y="50">only ±(x − 1)(x − 2), '
                'which meet only at f ≡ 0</text>')
    body += legend(1, [("series-1", "signs (+, −, +)"), ("series-2 dash", "signs (−, +, −)")])
    body += axes(1)
    body.append(curve(1, lambda x: (x - 1) * (x - 2), "series-1"))
    body.append(curve(1, lambda x: -(x - 1) * (x - 2), "series-2 dash"))
    body += roots(1)
    body.append(f'<text class="label" x="{PANELS[0][0]}" y="{H - 14}">Dots: the two steady states, '
                'x = 1 and x = 2, shared by every curve.</text>')
    body.append('</svg>')
    with open(OUT, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("\n".join(body) + "\n")
    print("written", os.path.normpath(OUT))
    return 0 if ok and ok2 else 1


if __name__ == "__main__":
    sys.exit(main())
