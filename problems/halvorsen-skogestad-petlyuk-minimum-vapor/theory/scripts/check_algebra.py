"""Symbolic checks (sympy, exact) of the algebraic identities used in the proofs of the paper.
Each claim prints PASS/FAIL; exit code 1 on failure.
Environment: Python 3 + sympy (WSL scientific environment).  Runtime: a few seconds.

 1. Recovery-ordering lemma (feed column of family K), diagonal form: with d = (zA, r zB, r zC, 0),
    g_T(x) = r g_F(x) + (1-r) aA zA/(aA-x) + r aD zD/(x-aD),  g_F(x) = sum_p a_p z_p/(a_p - x).
 2. Recovery-ordering lemma, ratio:  A12 D32 / (A32 D12) = (aA - t3)(t1 - aD) / ((aA - t1)(t3 - aD)).
 3. Family-K lemma: Delta_C(x)/Delta_B(x) = (aB - x)(aB - t)/((x - aC)(t - aC)),
    Delta_p(x) := 1/(a_p - x) - 1/(a_p - t).
 4. X/(X+Y) + Z/(Z+W) - 1 = (XZ - YW)/((X+Y)(Z+W)).
 5. Attainment theorem (FTC), partial fractions for a 4-component column [i,i+3] (3 roots):
    with a_p d_p = V prod_r (a_p - t_r)/prod_{p' != p} (a_p - a_p'),
    V - sum_p a_p d_p/(a_p - x) = V prod_r (x - t_r)/prod_p (x - a_p); and the bottom analogue.
 6. Attainment theorem: at x = 0 the identity gives L_T = V_T prod t_r / prod a_p (top) and
    L_B = V_B prod t_r / prod a_{p+1} (bottom).
"""
import sys
import sympy as sp

ok_all = True


def report(name, expr):
    global ok_all
    ok = sp.simplify(expr) == 0
    print(('PASS' if ok else 'FAIL'), name)
    ok_all &= ok


x, r = sp.symbols('x r')
aA, aB, aC, aD = sp.symbols('aA aB aC aD', positive=True)
zA, zB, zC, zD = sp.symbols('zA zB zC zD', positive=True)
t, t1, t2, t3 = sp.symbols('t t1 t2 t3')

# 1
gT = aA*zA/(aA - x) + r*aB*zB/(aB - x) + r*aC*zC/(aC - x)
gF = aA*zA/(aA - x) + aB*zB/(aB - x) + aC*zC/(aC - x) + aD*zD/(aD - x)
report('1 diagonal form of g_T', gT - (r*gF + (1 - r)*aA*zA/(aA - x) + r*aD*zD/(x - aD)))

# 2
A = lambda y: aA*zA/(aA - y)
D = lambda y: aD*zD/(y - aD)
A12, D12 = A(t1) - A(t2), D(t2) - D(t1)
A32, D32 = A(t2) - A(t3), D(t3) - D(t2)
report('2 ratio A12 D32/(A32 D12)', A12*D32/(A32*D12) - (aA - t3)*(t1 - aD)/((aA - t1)*(t3 - aD)))

# 3
Dl = lambda a, y: 1/(a - y) - 1/(a - t)
report('3 Delta_C/Delta_B', Dl(aC, x)/Dl(aB, x) - (aB - x)*(aB - t)/((x - aC)*(t - aC)))

# 4
X, Y, Z, W = sp.symbols('X Y Z W', positive=True)
report('4 two-fraction identity', X/(X + Y) + Z/(Z + W) - 1 - (X*Z - Y*W)/((X + Y)*(Z + W)))

# 5, 6
V, VB = sp.symbols('V VB')
a = sp.symbols('a0:4', positive=True)          # a0 > a1 > a2 > a3 (components i..i+3)
th = sp.symbols('s0:3')                         # roots s_r in (a_{r+1}, a_r)
top = range(0, 3)                                # top section components i..i+2
ad = {p: V*sp.prod([a[p] - th[q] for q in range(3)])/sp.prod([a[p] - a[pp] for pp in top if pp != p]) for p in top}
lhs = V - sum(ad[p]/(a[p] - x) for p in top)
rhs = V*sp.prod([x - th[q] for q in range(3)])/sp.prod([x - a[p] for p in top])
report('5 top partial fractions', sp.together(lhs - rhs))
bot = range(1, 4)                                # bottom section components i+1..i+3
ab = {p: -VB*sp.prod([a[p] - th[q] for q in range(3)])/sp.prod([a[p] - a[pp] for pp in bot if pp != p]) for p in bot}
lhsB = VB + sum(ab[p]/(a[p] - x) for p in bot)   # V_B - g_B(x), g_B = -sum a b/(a - x)
rhsB = VB*sp.prod([x - th[q] for q in range(3)])/sp.prod([x - a[p] for p in bot])
report('5 bottom partial fractions', sp.together(lhsB - rhsB))
LT = V - sum(ad[p]/a[p] for p in top)
report('6 L_T at x = 0', LT - V*sp.prod(th)/sp.prod([a[p] for p in top]))
LB = VB + sum(ab[p]/a[p] for p in bot)
report('6 L_B at x = 0', LB - VB*sp.prod(th)/sp.prod([a[p] for p in bot]))

print('ALL PASS' if ok_all else 'SOME FAIL')
sys.exit(0 if ok_all else 1)
