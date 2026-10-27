"""Numerical checks of explicit physical prefractionator profiles.

Python 3 + mpmath, 70 digits, normally < 30 s. No file writes.
Deterministic N=3..7, q=-0.5,0,0.5,1,1.5; PASS/FAIL and nonzero exit on failure.
These are numerical checks of a proved formula, not interval certificates.
Attainment: GPT-6-Astra (effort xhigh).
"""
import os
for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[name] = "2"
import mpmath as mp
mp.mp.dps = 70


def product(vals):
    return mp.fprod(list(vals))


def bisect(fun, lo, hi):
    # Avoid poles without a scale-sensitive fixed absolute offset.
    pad = (hi - lo) * mp.mpf("1e-65")
    lo, hi = lo + pad, hi - pad
    fl, fh = fun(lo), fun(hi)
    assert fl * fh < 0
    for _ in range(240):
        mid = (lo + hi) / 2
        fm = fun(mid)
        if fm == 0:
            return mid
        if fl * fm > 0:
            lo, fl = mid, fm
        else:
            hi = mid
    return (lo + hi) / 2


def preferred(a, z, q):
    n = len(a)
    U = 1 - q
    fun = lambda t: sum(a[i] * z[i] / (a[i] - t) for i in range(n)) - U
    th = [bisect(fun, a[j+1], a[j]) for j in range(n-1)]
    coeff = [product(a[i] - t for t in th) /
             (a[i] * product(a[i] - a[j] for j in range(n-1) if j != i)) for i in range(n-1)]
    VT = z[0] / coeff[0]
    d = [VT * c for c in coeff] + [mp.mpf(0)]
    b = [z[i] - d[i] for i in range(n)]
    b[0] = mp.mpf(0)
    return th, VT, VT-U, d, b


def construct(a, z, q, eta_ratio=mp.mpf("0.01")):
    n = len(a)
    th, VT0, VB0, d0, b0 = preferred(a, z, q)
    beta = product(a[1:-1]) / product(th[1:-1])
    h = [mp.mpf(0)] * n
    for i in range(1, n-1):
        h[i] = ((a[i]-beta) * product(a[i]-t for t in th[1:-1]) /
                (a[i] * product(a[i]-a[j] for j in range(1,n-1) if j != i)))
    eta = eta_ratio * min(VT0, VB0)
    for _ in range(60):
        d = [d0[i]+eta*h[i] for i in range(n)]
        b = [b0[i]-eta*h[i] for i in range(n)]
        VT, VB = VT0+eta, VB0+eta
        LT, LB = VT-sum(d), VB+sum(b)
        ec = [a[-1]*d[i]/(LT*(a[i]-a[-1])) for i in range(n-1)]
        ec += [1-sum(ec)]
        ea = [a[0]*b[i]/(LB*(a[0]-a[i])) for i in range(1,n)]
        ea = [1-sum(ea)]+ea
        if min(d[:-1]+b[1:]+ec+ea) > 0:
            break
        eta /= 2
    else:
        raise AssertionError("could not find small positive slack")
    pt = lambda t: VT-sum(a[i]*d[i]/(a[i]-t) for i in range(n))
    pb = lambda t: VB+sum(a[i]*b[i]/(a[i]-t) for i in range(n))
    r, u = bisect(pt, a[1], a[0]), bisect(pt, a[-1], a[-2])
    v, ss = bisect(pb, a[1], a[0]), bisect(pb, a[-1], a[-2])
    assert a[0] > r > v > a[1] and a[-2] > u > ss > a[-1]
    chi = [r]+th[1:-1]+[ss]
    weights = [product(a[i]-t for t in chi) /
               (a[i]*product(a[i]-a[j] for j in range(n) if j != i)) for i in range(n)]
    x = [v/sum(weights) for v in weights]
    eu = [u*d[i]/(LT*(a[i]-u)) for i in range(n-1)]+[mp.mpf(0)]
    ev = [mp.mpf(0)]+[v*b[i]/(LB*(v-a[i])) for i in range(1,n)]
    A = x[-1]/ec[-1]
    B = 1-A
    D = x[0]/ea[0]
    C = 1-D
    assert min(A,B,C,D) > 0
    assert max(abs(x[i]-(A*ec[i]+B*eu[i])) for i in range(n)) < mp.mpf("1e-55")
    assert max(abs(x[i]-(C*ev[i]+D*ea[i])) for i in range(n)) < mp.mpf("1e-55")
    top = lambda k: [(A*(a[-1]/u)**k*ec[i]+B*eu[i])/(A*(a[-1]/u)**k+B) for i in range(n)]
    bottom = lambda k: [(C*ev[i]+D*(v/a[0])**k*ea[i])/(C+D*(v/a[0])**k) for i in range(n)]
    return dict(a=a,z=z,q=q,eta=eta,VT=VT,VB=VB,LT=LT,LB=LB,d=d,b=b,
                x=x,top=top,bottom=bottom,r=r,u=u,v=v,s=ss,A=A,B=B,C=C,D=D)


def equilibrium(a, x):
    S = sum(a[i]*x[i] for i in range(len(a)))
    return [a[i]*x[i]/S for i in range(len(a))]


def check_case(a,z,q):
    d = construct(a,z,q)
    n = len(a)
    residual = mp.mpf(0)
    for k in range(65):
        xt, xt1 = d["top"](k), d["top"](k+1)
        xb, xb1 = d["bottom"](k), d["bottom"](k+1)
        assert min(xt+xt1+xb+xb1) > 0
        yt, yb1 = equilibrium(a,xt), equilibrium(a,xb1)
        residual = max(residual,
            max(abs(d["VT"]*yt[i]-d["LT"]*xt1[i]-d["d"][i]) for i in range(n)),
            max(abs(d["VB"]*yb1[i]-d["LB"]*xb[i]+d["b"][i]) for i in range(n)),
            abs(sum(xt)-1),abs(sum(xb)-1))
    feed = [d["LT"]*d["top"](1)[i]+d["VB"]*equilibrium(a,d["bottom"](1))[i]+z[i]
            -d["VT"]*equilibrium(a,d["x"])[i]-d["LB"]*d["x"][i] for i in range(n)]
    residual = max(residual,max(map(abs,feed)))
    assert residual < mp.mpf("1e-50"), residual
    top_bound = d["A"]/d["B"]*(a[-1]/d["u"])**64
    bot_bound = d["D"]/d["C"]*(d["v"]/a[0])**64
    assert d["top"](64)[-1] <= top_bound and d["bottom"](64)[0] <= bot_bound
    print(f"PASS numerical N={n} q={float(q):g} residual={float(residual):.2e} "
          f"top_heavy_64={float(d['top'](64)[-1]):.3e} bottom_light_64={float(d['bottom'](64)[0]):.3e}",
          flush=True)


def main():
    for n in range(3,8):
        a = [mp.mpf(2)**(n-i-1) for i in range(n)]
        raw = [mp.mpf(i+2) for i in range(n)]
        z = [x/sum(raw) for x in raw]
        for q in map(mp.mpf,("-0.5","0","0.5","1","1.5")):
            check_case(a,z,q)
    print("PASS all explicit-profile numerical checks; full FTC matching is not certified")


if __name__ == "__main__":
    try:
        main()
    except Exception:
        print("FAIL explicit-profile numerical check",flush=True)
        raise
