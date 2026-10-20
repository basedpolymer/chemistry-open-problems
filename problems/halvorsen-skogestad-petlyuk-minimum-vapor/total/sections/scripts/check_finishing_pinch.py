#!/usr/bin/env python3
"""
check_finishing_pinch.py (Agent S) -- exact (sympy) identities behind sections.md §8.6: for the two binary
finishing columns of a ternary Petlyuk arrangement at leading order, the pinch inequality of the binary shooting
argument (Prop. 10 / 10') is equivalent to the vapor exceeding the Underwood minimum computed with the
prefractionator's actual section roots.
  AB column: prefractionator top section with net flows (z_A, d_B, 0), vapor V_P, root r in (b, a):
     V_P = a z_A/(a - r) + b d_B/(b - r).  Main top section V_T, L_T = V_T - z_A (D = z_A); upper-middle section
     V_U = V_T - V_P, L_U = V_U + d_B.  Pinches (binary A/B, relative volatility a/b):
     xT = (V_T/L_T - 1)/(a/b - 1) (B present at zero net flow), xU = ((a/b) V_U/L_U - 1)/(a/b - 1) (A present).
     g(V_T) = (a/b) V_U L_T - V_T L_U  (xT < xU  <=>  g > 0 when L_T, L_U > 0) has the roots
     V_T = a z_A/(a - r) (the larger) and a smaller one; leading coefficient (a - b)/b > 0.
  BC column: prefractionator bottom section with net flows -(0, b_B, z_C), vapor V_PB, root s in (c, b):
     V_PB + b b_B/(b - s) + c z_C/(c - s) = 0.  Main bottom section V_Bo, L_Bo = V_Bo + z_C (B = z_C); lower-middle
     section V_W = V_Bo - V_PB, L_W = L_Bo - L_PB, L_PB = V_PB + b_B + z_C.  g2(V_Bo) = (b/c) V_Bo L_W - V_W L_Bo has
     the roots V_Bo = c z_C/(s - c) (the larger) and a smaller one; leading coefficient (b - c)/c > 0.
Prints PASS/FAIL. Environment: Python 3, sympy. Runtime: seconds.
"""
import sys
import sympy as sp

ok = True
a, b, c, zA, zC, dB, bB, r, s, VT, VBo = sp.symbols('a b c z_A z_C d_B b_B r s V_T V_Bo', positive=True)

VP = a * zA / (a - r) + b * dB / (b - r)
LT = VT - zA; VU = VT - VP; LU = VU + dB
g = (a / b) * VU * LT - VT * LU
Vmin = a * zA / (a - r)
Voth = a * zA / (a - b) - b * dB * (a - r) / ((a - b) * (r - b))
c1 = sp.simplify(g.subs(VT, Vmin)) == 0
c2 = sp.simplify(g.subs(VT, Voth)) == 0
c3 = sp.simplify(sp.Poly(sp.expand(g), VT).all_coeffs()[0] - (a - b) / b) == 0
# Vmin - Voth = a zA [1/(a-r) - 1/(a-b)] + b dB (a-r)/((a-b)(r-b)) > 0 for b < r < a (both terms positive)
c4 = sp.simplify(Vmin - Voth - (a * zA * (1 / (a - r) - 1 / (a - b)) + b * dB * (a - r) / ((a - b) * (r - b)))) == 0
print('AB column: g(a zA/(a-r)) = 0:', c1, '| other root:', c2, '| leading coeff (a-b)/b:', c3, '| difference form:', c4)
ok &= c1 and c2 and c3 and c4

VPB = -(b * bB / (b - s) + c * zC / (c - s))
LPB = VPB + bB + zC
LBo = VBo + zC; VW = VBo - VPB; LW = LBo - LPB
g2 = (b / c) * VBo * LW - VW * LBo
Vmin2 = c * zC / (s - c)
Voth2 = c * zC / (b - c) - b * bB * (s - c) / ((b - s) * (b - c))
d1 = sp.simplify(g2.subs(VBo, Vmin2)) == 0
d2 = sp.simplify(g2.subs(VBo, Voth2)) == 0
d3 = sp.simplify(sp.Poly(sp.expand(g2), VBo).all_coeffs()[0] - (b - c) / c) == 0
d4 = sp.simplify(VW - LW - bB) == 0
d5 = sp.simplify(Vmin2 - Voth2 - (c * zC * (1 / (s - c) - 1 / (b - c)) + b * bB * (s - c) / ((b - s) * (b - c)))) == 0
print('BC column: g2(c zC/(s-c)) = 0:', d1, '| other root:', d2, '| leading coeff (b-c)/c:', d3,
      '| V_W - L_W = b_B:', d4, '| difference form:', d5)
ok &= d1 and d2 and d3 and d4 and d5
print('PASS' if ok else 'FAIL')
sys.exit(0 if ok else 1)
