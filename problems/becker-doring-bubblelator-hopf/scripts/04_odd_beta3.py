"""The case nu = 0, beta = 3 (alpha = 0, r = 3) of Section 4.3 of the source: the first NZ crossings.

Claims (PASS/FAIL, exit code 1 on failure):
  OB-1  H(t) = Re G_{3,0}(it) has exactly M zeros t_1 < ... < t_M in (0, T_scan], each simple, with
        H'(t_j) of sign (-1)^j, and Im G(it_j) < 0;  T_scan >= sqrt(2/vartheta_NZ).
  OB-2  vartheta_j = -Im G(it_j)/t_j is strictly decreasing for j = 1..NZ, and vartheta_j < vartheta_NZ
        for NZ < j <= M.  Together with |G(it)| <= 2/t (nu = 0, integration by parts) this shows
        that for vartheta >= vartheta_NZ the only purely imaginary characteristic roots are
        +-i t_j at vartheta = vartheta_j (j <= NZ): hypotheses (H1)-(H3) hold at each of the NZ
        crossings, and (H4) at the first one.
  OB-3  Re ghat_1(1) != 0 and eta_2 != 0 at each crossing (enclosures printed); signs reported.
Writes certificates/odd_beta3.csv.

Runtime: about 2 minutes.  Environment: WSL Fedora, Python 3, python-flint 0.9.  Deterministic
(the floating-point search for the zeros is part of the script but the claims only use the
rigorous checks that follow it).
"""
import sys, os, time
from flint import arb, acb, ctx
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import bdhopf as B
import certify_lib as C

ctx.prec = 128
B.TOL_BITS = 90
NZ = int(os.environ.get('NZ', '11'))
FAIL = []


def report(name, ok, msg=''):
    print(('PASS' if ok else 'FAIL'), name, msg)
    if not ok:
        FAIL.append(name)


def Hf(beta, nu, t):
    return float(B.H3_line(arb(t)).mid())


def Hball(t):
    return B.H3_line(t)


def Hpball(t):
    return B.H3_line(t, deriv=1)


def sign_on(a, b, sgn, depth=0):
    """certify sign(H') = sgn on [a, b] by bisection"""
    hp = Hpball(C.iv(a, b))
    if (bool(hp < 0) if sgn < 0 else bool(hp > 0)):
        return True
    if depth > 30:
        return False
    m = 0.5 * (a + b)
    return sign_on(a, m, sgn, depth + 1) and sign_on(m, b, sgn, depth + 1)


def H_nonzero_on(a, b):
    """centered form: H(J) in H(m) + H'(J) [-w/2, w/2]"""
    m = 0.5 * (a + b)
    J = C.iv(a, b)
    hm = Hball(arb(m))
    hp = Hpball(J)
    w2 = arb((b - a) / 2)
    enc = hm + hp * arb(0, 1) * w2
    return bool(enc > 0) or bool(enc < 0)


def main():
    t00 = time.time()
    beta, nu = arb(3), arb(0)
    # ---- floating search for zeros (exploration) ----
    zs = []
    t, h = 0.05, Hf(beta, nu, 0.05)
    step = 0.01
    tstop = None
    while True:
        t2 = t + step
        h2 = Hf(beta, nu, t2)
        if (h > 0) != (h2 > 0):
            a, b = t, t2
            for _ in range(80):
                m = 0.5 * (a + b)
                if (Hf(beta, nu, m) > 0) == (h > 0):
                    a = m
                else:
                    b = m
            zs.append(0.5 * (a + b))
            if len(zs) == NZ:
                th = -float(B.G(beta, nu, acb(0, zs[-1])).imag.mid()) / zs[-1]
                tstop = 1.02 * (2.0 / th) ** 0.5
        t, h = t2, h2
        if tstop is not None and t > tstop:
            break
    # ---- rigorous part ----
    g0 = B.G0(beta, nu)
    Tscan = tstop
    recs = []
    covered = []
    ok1 = True
    for j, z in enumerate(zs, start=1):
        r = max(z, 1.0) * 2.0 ** -45
        tlo, thi = z - r, z + r
        Hlo, Hhi = Hball(arb(tlo)), Hball(arb(thi))
        sgn = -1 if j % 2 == 1 else 1        # expected sign of H'
        c1 = (bool(Hlo > 0) and bool(Hhi < 0)) if sgn < 0 else (bool(Hlo < 0) and bool(Hhi > 0))
        # derivative of constant sign on [tlo-d, thi+d]
        d = 1e-3 * z
        good = sign_on(tlo - d, thi + d, sgn)
        T = C.iv(tlo, thi)
        GT = B.G(beta, nu, acb(0, T))
        c3 = bool(GT.imag < 0)
        ok1 &= (c1 and good and c3)
        covered.append((tlo - d, thi + d))
        recs.append(dict(j=j, T=T, sgn=sgn, vth=-GT.imag / T, c=(c1, good, c3)))
    # scan the gaps: H != 0
    origin_ok = bool(g0 - arb(2.0 ** -20) * B.G0(beta, nu + 1) > 0)
    gaps = []
    left = 2.0 ** -20
    for (a, b) in covered:
        gaps.append((left, a))
        left = b
    gaps.append((left, Tscan))
    nev = 0
    scan_ok = origin_ok
    for (a0, b0) in gaps:
        stack = [(a0 + (b0 - a0) * k / 16, a0 + (b0 - a0) * (k + 1) / 16 if k < 15 else b0, 0) for k in range(16)]
        while stack:
            a, b, dep = stack.pop()
            nev += 1
            if H_nonzero_on(a, b):
                continue
            if dep > 40:
                scan_ok = False
                print('   scan failure on', a, b)
                stack = []
                break
            m = 0.5 * (a + b)
            stack += [(a, m, dep + 1), (m, b, dep + 1)]
    M = len(zs)
    report('OB-1 exactly %d simple zeros of H in (0, %.3f], alternating H\' signs, Im G < 0' % (M, Tscan),
           ok1 and scan_ok, '(%d scan evaluations)' % nev)
    # OB-2
    v = [rr['vth'] for rr in recs]
    dec = all(bool(v[i + 1] < v[i]) for i in range(NZ - 1))
    rest = all(bool(v[i] < v[NZ - 1]) for i in range(NZ, M))
    tail = bool(arb(Tscan) ** 2 * v[NZ - 1] > 2)
    report('OB-2 vartheta_1 > ... > vartheta_%d, later zeros below, T_scan^2 vartheta_%d > 2' % (NZ, NZ),
           dec and rest and tail)
    # OB-3
    rows = []
    ok3 = True
    for rr in recs[:NZ]:
        q = B.hopf_quantities(beta, nu, rr['T'])
        s11 = '+' if q['ghat1_1'].real > 0 else ('-' if q['ghat1_1'].real < 0 else '?')
        se = '+' if q['eta2'] > 0 else ('-' if q['eta2'] < 0 else '?')
        ok3 &= (s11 != '?' and se != '?' and q['eta0'].is_finite())
        rows.append((rr['j'], q['vartheta0'].str(10, radius=True), q['eta0'].str(10, radius=True), q['kappa0'].str(10, radius=True),
                     q['eta2'].str(8, radius=True), q['kappa2'].str(8, radius=True), q['ghat1_1'].real.str(8, radius=True), s11, se))
        print('   crossing %2d: vartheta0=%s eta0=%s kappa0=%s eta2=%s kappa2=%s Re ghat1(1)=%s' % rows[-1][:7])
    report('OB-3 transversality and direction certified at each crossing', ok3)
    with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'certificates', 'odd_beta3.csv'), 'w') as fh:
        fh.write('j;vartheta0;eta0;kappa0;eta2;kappa2;Re_ghat1_1;sign_Re_ghat1_1;sign_eta2\n')
        for r in rows:
            fh.write(';'.join(str(x) for x in r) + '\n')
    print('runtime %.0f s' % (time.time() - t00))
    if FAIL:
        print('FAILED:', FAIL)
        sys.exit(1)
    print('ALL PASS')


if __name__ == '__main__':
    main()
