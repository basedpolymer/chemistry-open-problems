"""Theorem 2 return-ratio construction; deterministic numerical checks only.

Python 3 + mpmath, 120 digits; normally <10 s. No files written.
Checks every cross-section and interior composition for 9 constructions.
Attainment: GPT-6-Astra via Codex CLI (effort xhigh).
"""
import os
for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[name] = "2"
import sys
sys.dont_write_bytecode = True
import mpmath as mp
from check_prefractionator import construct, equilibrium
mp.mp.dps = 120


def normalize(x):
    return x/sum(x)


def main():
    a = list(map(mp.mpf,[4,2,1]))
    z = [mp.mpf(1)/3]*3
    dat = construct(a,z,mp.mpf(1),mp.mpf("0.02"))
    av = mp.matrix([a])
    MT = dat["VT"]*mp.diag(a)-mp.matrix(dat["d"])*av
    MB = dat["VB"]*mp.diag(a)+mp.matrix(dat["b"])*av
    T = MT/dat["VT"]
    B = MB/dat["VB"]
    for xi,zeta in (("0.2","0.8"),("0.65","0.5"),("0.9","0.15")):
        xi,zeta = mp.mpf(xi),mp.mpf(zeta)
        fU = mp.matrix([[1-xi,-xi,0]])
        fD = mp.matrix([[0,1-zeta,-zeta]])
        for n in (20,40,80):
            k = n+3
            rt = fU*(T/dat["r"])**n
            rb = fD*(dat["s"]*B**-1)**k
            J = mp.matrix([[1,1,1],list(rt),list(rb)])
            X = mp.lu_solve(J,mp.matrix([1,0,0]))
            top = [normalize(T**j*X) for j in range(n+1)]
            bot = [normalize(B**(-j)*X) for j in range(k+1)]
            minimum = min(min(x) for x in top+bot)
            assert minimum > 0
            residual = mp.mpf(0)
            for j in range(n):
                y = equilibrium(a,top[j])
                residual = max(residual,max(abs(dat["VT"]*y[i]-dat["LT"]*top[j+1][i]-dat["d"][i]) for i in range(3)))
            for j in range(k):
                y = equilibrium(a,bot[j+1])
                residual = max(residual,max(abs(dat["VB"]*y[i]-dat["LB"]*bot[j][i]+dat["b"][i]) for i in range(3)))
            ratio_error = max(abs(top[-1][0]/(top[-1][0]+top[-1][1])-xi),
                              abs(bot[-1][1]/(bot[-1][1]+bot[-1][2])-zeta))
            assert max(residual,ratio_error) < mp.mpf("1e-80")
            print(f"PASS numerical xi={float(xi):g} zeta={float(zeta):g} n={n} k={k} "
                  f"balance={float(residual):.3e} ratio_error={float(ratio_error):.3e} "
                  f"upper_C={float(top[-1][2]):.3e} lower_A={float(bot[-1][0]):.3e}",
                  flush=True)
    print("PASS prescribed return ratios, positive complete profiles; isolated prefractionator only",flush=True)


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print("FAIL "+str(exc),flush=True)
        sys.exit(1)
