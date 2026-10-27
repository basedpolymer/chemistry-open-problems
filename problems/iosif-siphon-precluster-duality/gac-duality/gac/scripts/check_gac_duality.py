"""Exact certificates for agent G's transpose-duality results.

Environment: Python 3.12+, standard library only; repository scripts/crn.py.
Expected runtime: under 10 s. Arithmetic is integer/Fraction throughout.
Prints PASS/FAIL per claim, exits nonzero on any failure.
"""
from fractions import Fraction as F
from itertools import combinations, product
from pathlib import Path
import platform
import sys
import time

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "scripts"))
from crn import Network, rank, nullspace

FAILED = []
COUNT = 0

def check(name, condition):
    global COUNT
    COUNT += 1
    print(("PASS " if condition else "FAIL ") + name)
    if not condition:
        FAILED.append(name)

def criterion(N):
    supports = N.mips_semiflow()
    return all(any(T <= S for T in supports) for S in N.minimal_siphons())

def dot(a, b):
    return sum(x * y for x, y in zip(a, b))

def field(N, x, k):
    rates = [k[j] * monomial(x, a) for j, (a, _) in enumerate(N.reactions)]
    return [dot(row, rates) for row in N.Gamma]

def monomial(x, a):
    out = F(1)
    for xi, ai in zip(x, a):
        out *= xi ** ai
    return out

def incidence(N):
    return [[int(b == c) - int(a == c) for a, b in N.reactions] for c in N.complexes()]

def t_trap(N, J):
    consumed = {i for j in J for i, a in enumerate(N.reactions[j][0]) if a}
    produced = {i for j in J for i, b in enumerate(N.reactions[j][1]) if b}
    return produced <= consumed

def endotactic_violation(N, w):
    # Sweep from high to low; ignore reactions orthogonal to w.
    active = [(a, b) for a, b in N.reactions if dot(w, tuple(y-x for x,y in zip(a,b)))]
    if not active:
        return False
    top = max(dot(w, a) for a, b in active)
    return any(dot(w, a) == top and dot(w, tuple(y-x for x,y in zip(a,b))) > 0 for a,b in active)

def main():
    print(f"Environment: Python {platform.python_version()}, {platform.platform()}; exact Fraction arithmetic; expected runtime <10 s")
    start = time.perf_counter()
    M = Network(2, [((3,0),(2,1)), ((2,1),(3,0))], names=["A","B"])
    D = M.dual()
    check("counterexample dual reactions", D.reactions == [((3,2),(2,3)),((0,1),(1,0))])
    for label, N in [("M", M), ("D", D)]:
        check(label + " exact matrix/rank/nullities", N.Gamma == [[-1,1],[1,-1]] and N.rank()==1 and N.n-N.rank()==1)
        check(label + " H1 H2 H3", N.is_conservative() and N.is_consistent() and N.rows_distinct())
        check(label + " positive equilibrium at all ones", field(N,[F(1),F(1)],[F(1),F(1)])==[0,0])
    check("M reversible, one linkage, two complexes, deficiency 0", M.weakly_reversible() and M.linkage_classes()==1 and len(M.complexes())==2 and M.deficiency()==0)
    check("D not WR, two linkages, four complexes, deficiency 1", not D.weakly_reversible() and D.linkage_classes()==2 and len(D.complexes())==4 and D.deficiency()==1)
    check("D violates endotacticity in rational direction (0,1)", endotactic_violation(D,(0,1)))
    u,v=F(15,8),F(1,8)
    check("extinction initial state shares equilibrium's class", u+v==2)
    check("extinction rational threshold", u**3*v==F(3375,4096) and u**3*v<1)
    check("extinction vector field at initial state", field(D,[u,v],[F(1),F(1)])[1] == -F(721,32768))
    # h(v)=(2-v)^3 v, h'(v)=(2-v)^2(2-4v)>0 on [0,1/8].
    check("rational lower bounds for factors of h' on interval", 2-F(1,8)>0 and 2-4*F(1,8)>0)
    check("dual incidence circulations do not contain positive vector", nullspace(incidence(D),D.r)==[])
    check("M complex balanced at all ones", all(dot(row,[1,1])==0 for row in incidence(M)))
    # Product/source row matching characterizes dual WR when primal WR and H3.
    C = Network(3,[((1,0,0),(0,1,0)),((0,1,0),(0,0,1)),((0,0,1),(1,0,0))])
    W = Network(3,[((1,0,0),(0,1,0)),((1,0,0),(0,0,1)),((0,1,0),(1,0,0)),((0,0,1),(1,0,0))])
    Q = Network(2,[((0,1),(2,0)),((2,0),(0,1))])
    Z = Network(1,[((0,),(1,)),((1,),(0,))])
    H = Network(2,[((2,0),(0,2)),((0,2),(2,0)),((2,0),(1,1)),((1,1),(2,0))])
    check("H ADLS criterion passes, dual criterion fails", criterion(H) and not criterion(H.dual()))
    check("H and dual H satisfy H1 H2 H3", all(N.is_conservative() and N.is_consistent() and N.rows_distinct() for N in [H,H.dual()]))
    check("dual H singleton 4 is critical", H.dual().is_siphon({3}) and H.dual().Gamma[3]==[1,-1])
    L = Network(2,[((1,0),(0,1)),((0,1),(1,0))])
    check("same primal matrix/counts can have different dual counts", L.Gamma==M.Gamma and (len(L.complexes()),L.linkage_classes(),L.deficiency())==(len(M.complexes()),M.linkage_classes(),M.deficiency()) and len(L.dual().complexes())!=len(M.dual().complexes()))
    check("bounded-class property is not unconditional: 0 <-> X", not Z.is_conservative() and Z.is_consistent() and Z.dual().is_conservative())
    networks=[M,D,C,W,Q,Z,Z.dual(),H]
    for label,N in zip(["M","D","C","W","Q","Z","Zstar","H"],networks):
        Nd=N.dual()
        check(label+" transposition/involution", Nd.Gamma==[list(c) for c in zip(*N.Gamma)] and Nd.dual().reactions==N.reactions)
        check(label+" nullity=graph cycle rank+deficiency", N.n-N.rank()==Nd.r-rank(incidence(Nd))+Nd.deficiency())
        check(label+" conservativity-consistency swap", N.is_conservative()==Nd.is_consistent() and N.is_consistent()==Nd.is_conservative())
        if N.weakly_reversible() and N.rows_distinct():
            check(label+" row-permutation WR characterization", Nd.weakly_reversible()==(sorted(map(tuple,N.Ye))==sorted(map(tuple,N.Yp))))
        for size in range(1,N.n+1):
            for Sigma in combinations(range(N.n),size):
                check(label+f" siphon/T-trap {Sigma}", N.is_siphon(Sigma)==t_trap(Nd,Sigma))
                rows=[N.Gamma[i] for i in Sigma]
                restricted=Network(Nd.n,[Nd.reactions[i] for i in Sigma])
                check(label+f" face projection/subnetwork rank {Sigma}", rank(rows)==restricted.rank())
        check(label+" ADLS/T-trap-cycle equivalence", criterion(N)==all(any(T<=frozenset(Sigma) for T in Nd.internal_cycles()) for size in range(1,N.n+1) for Sigma in combinations(range(N.n),size) if t_trap(Nd,Sigma)))
    # A rational facet certificate and Wiuf condition (a), using N: A+B <-> 2A, C <-> D.
    Fnet=Network(4,[((1,1,0,0),(2,0,0,0)),((2,0,0,0),(1,1,0,0)),((0,0,1,0),(0,0,0,1)),((0,0,0,1),(0,0,1,0))])
    Sigma=[0]
    J0=[2,3]
    Fs=Fnet.rank()
    Fi=[Fnet.Gamma[i] for i in Sigma]
    Fb=[[row[j] for j in J0] for row in Fnet.Gamma]
    check("face certificate positive point and boundary point differ in S", rank([list(c) for c in zip(*Fnet.Gamma)]+[[-1,1,0,0]])==Fs)
    check("face is a facet: s=2, projected rank=1", Fs==2 and rank(Fi)==1)
    check("Wiuf (a) rank sum for blind reactions", rank(Fi)+rank(Fb)==Fs)
    FD=Fnet.dual()
    check("Wiuf (a) identical transposed rank sum", rank([[FD.Gamma[j][i] for i in Sigma] for j in range(FD.n)])+rank([FD.Gamma[j] for j in J0])==Fs)
    # Only active embedded class is A <-> 2A; im(A)=R contains a positive vector.
    A=[[1,-1]]
    Atranspose=[[1],[-1]]
    check("Wiuf diamond positive image / dual no positive kernel", rank(A)==1 and nullspace(Atranspose,1)==[])
    # Failure of (a) in a single-linkage CB example: A+B <-> 2A <-> A+C.
    Wfail=Network(3,[((1,1,0),(2,0,0)),((2,0,0),(1,1,0)),((2,0,0),(1,0,1)),((1,0,1),(2,0,0))])
    check("Wiuf (a) can fail in a WR strongly endotactic network", Wfail.weakly_reversible() and Wfail.linkage_classes()==1 and Wfail.is_siphon({0}) and Wfail.rank()==2 and rank([Wfail.Gamma[0]])==1 and all(a[0]>0 for a,b in Wfail.reactions))
    # Complex balance does not give a common norm over the monotone structural class.
    Bnet=Network(3,[((3,0,2),(2,3,0)),((2,3,0),(0,2,3)),((0,2,3),(3,0,2))])
    Balt=Network(3,[((10,0,2),(9,3,0)),((2,1,0),(0,0,3)),((0,10,1),(3,8,0))])
    check("BG example both orientations WR and deficiency zero", all(N.weakly_reversible() and N.deficiency()==0 for N in [Bnet,Bnet.dual()]))
    check("BG example complex balanced at ones", all(dot(row,[1,1,1])==0 for row in incidence(Bnet)))
    check("BG support condition includes a nonpure network", all((Bnet.Ye[i][j]>0)==(Bnet.Gamma[i][j]<0) for i in range(Bnet.n) for j in range(Bnet.r)) and any(Bnet.Ye[i][j]>0 and Bnet.Yp[i][j]>0 for i in range(Bnet.n) for j in range(Bnet.r)))
    check("BG alternate realization has same net matrix and dependence support", Bnet.Gamma==Balt.Gamma and all((Bnet.Ye[i][j]>0)==(Balt.Ye[i][j]>0) for i in range(3) for j in range(3)))
    check("BG alternate kinetics positive equilibrium at ones", field(Balt,[F(1)]*3,[F(1)]*3)==[0,0,0])
    G=[list(c) for c in zip(*Balt.Ye)]
    J=[[sum(Balt.Gamma[i][j]*G[j][k] for j in range(3)) for k in range(3)] for i in range(3)]
    reduced=[[J[i][k]-J[i][2] for k in range(2)] for i in range(2)]
    check("BG alternate restricted Jacobian exact", reduced==[[-15,27],[24,-25]])
    check("BG alternate positive eigenvalue: determinant -273", reduced[0][0]*reduced[1][1]-reduced[0][1]*reduced[1][0]==-273)
    # Prime-size codimension-one certificate: two cyclic orbits in 5 species.
    def shift(v):
        return (v[-1],)+v[:-1]
    seed=(1,0,1,0,0)
    orbit=[]
    for i in range(5):
        orbit.append(seed)
        seed=shift(seed)
    P5=Network(5,[(tuple(scale*x for x in y),tuple(scale*x for x in shift(y))) for scale in [1,2] for y in orbit])
    check("prime example dimensions and complex graph counts", P5.rank()==4 and len(P5.complexes())==10 and P5.linkage_classes()==2 and P5.deficiency()==4)
    check("prime example WR, H3, dual WR and dual deficiency zero", P5.weakly_reversible() and P5.rows_distinct() and P5.dual().weakly_reversible() and P5.dual().deficiency()==0)
    check("prime example every linkage spans full S", all(Network(5,P5.reactions[j:j+5]).rank()==4 for j in [0,5]))
    check("prime example complex balanced at ones", all(dot(row,[1]*10)==0 for row in incidence(P5)))
    check("prime example critical siphon ABC", P5.is_siphon({0,1,2}) and P5.conservation_rays()==[(1,1,1,1,1)] and not criterion(P5))
    check("prime example realizable intermediate face", rank([P5.Gamma[i] for i in [0,1,2]])==3 and sum([0,0,0,2,3])==sum([1]*5))
    check("prime example fails Wiuf a", all(any(a[i]>0 for i in [0,1,2]) for a,b in P5.reactions) and P5.rank()-rank([P5.Gamma[i] for i in [0,1,2]])==1)
    e4=[tuple(int(i==j) for i in range(4)) for j in range(4)]
    C4=Network(4,[(e4[i],e4[(i+1)%4]) for i in range(4)]+[((1,0,1,0),(0,1,0,1)),((0,1,0,1),(1,0,1,0))])
    w=(3,2,1,2)
    heights=[dot(w,a) for a,b in C4.reactions]
    drifts=[dot(w,tuple(y-x for x,y in zip(a,b))) for a,b in C4.reactions]
    check("composite example q=1, both orientations WR and H3", C4.rank()==3 and C4.rows_distinct() and C4.weakly_reversible() and C4.dual().weakly_reversible())
    check("composite example strong-endotacticity violation", any(drifts) and max(heights)==4 and all(drifts[j]==0 for j,h in enumerate(heights) if h==4))
    for p in [2,3,5,7]:
        full_rank=True
        tested=0
        for y in product([0,1],repeat=p):
            if len(set(y))==1:
                continue
            y0=y
            reactions=[]
            for i in range(p):
                reactions.append((y,shift(y)))
                y=shift(y)
            full_rank &= Network(p,reactions).rank()==p-1
            tested+=1
        check(f"prime p={p}: all {tested} nonconstant binary complex orbits span sum-zero space", full_rank and tested==2**p-2)
    # Arbitrarily high dimension at codimension one, sampled n=2,...,10.
    for n in range(2,11):
        e=[tuple(int(i==j) for i in range(n)) for j in range(n)]
        N=Network(n,[(e[j],e[(j+1)%n]) for j in range(n)])
        check(f"cycle family n={n}, q=1, s=n-1 in both orientations", N.rank()==n-1 and N.dual().rank()==n-1 and N.dual().weakly_reversible())
    print(f"Claims: {COUNT}; failures: {len(FAILED)}; runtime: {time.perf_counter()-start:.3f} s")
    return 1 if FAILED else 0

if __name__ == "__main__":
    sys.exit(main())
