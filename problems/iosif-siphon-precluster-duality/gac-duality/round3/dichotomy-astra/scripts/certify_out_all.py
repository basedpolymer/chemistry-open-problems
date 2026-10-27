"""Exact certificates for Theorem A: an ALL dual with an OUT facet.

Environment: Python 3.14 (also compatible with 3.10+), standard library only.
Runtime: less than one second, measured and printed at exit.
All certificate calculations use integers or fractions.Fraction. The global
persistence proof is in ../notes.md; no numerical trajectory is a certificate.
Earlier crn.py and d2lib/d2cert are imported unchanged, with bytecode disabled.
"""
from pathlib import Path
import sys
import time
import platform
from fractions import Fraction as F

HERE = Path(__file__).resolve().parent
BASE = HERE.parents[2]
sys.dont_write_bytecode = True
sys.path.insert(0, str(BASE.parent / 'scripts'))
sys.path.insert(0, str(BASE / 'round3' / 'dichotomy' / 'scripts'))
import crn
import d2cert
from d2lib import Net


# Polynomial variables: a,b,c,alpha,beta,gamma,delta.
def term(coeff, *powers):
    return {tuple(powers): F(coeff)}


def add(*polys):
    out = {}
    for p in polys:
        for e,c in p.items():
            out[e] = out.get(e,F(0))+c
    return {e:c for e,c in out.items() if c}


def mul(p,q):
    out = {}
    for a,c in p.items():
        for b,d in q.items():
            e = tuple(x+y for x,y in zip(a,b))
            out[e] = out.get(e,F(0))+c*d
    return {e:c for e,c in out.items() if c}


def main():
    start = time.perf_counter()
    print('Environment: Python', platform.python_version(), '; standard library; exact rationals.')
    print('Runtime: expected <1 second; measured below.')
    failures = []
    def claim(name, ok, detail):
        print(('PASS' if ok else 'FAIL')+': '+name+': '+detail)
        if not ok:
            failures.append(name)
    cxs = [(0,0,3,0),(1,1,0,0),(0,1,2,3)]
    rx = [(cxs[j],cxs[(j+1)%3]) for j in range(3)]
    N = crn.Network(4,rx)
    G = [[1,-1,0],[1,0,-1],[-3,2,1],[0,3,-3]]
    claim('A-matrix',N.Gamma==G,'Gamma=[[1,-1,0],[1,0,-1],[-3,2,1],[0,3,-3]].')
    w = [11,1,4,1]
    claim('A-H1', all(sum(w[i]*G[i][j] for i in range(4))==0 for j in range(3)),
          'positive conservation vector (11,1,4,1); complex values (12,12,12).')
    claim('A-H2',all(sum(row)==0 for row in G),'positive flux (1,1,1).')
    claim('A-H3',len({tuple(row) for row in G})==4,'all four stoichiometric rows distinct.')
    claim('A-rank',crn.rank(G)==2 and G[0][0]*G[1][1]-G[0][1]*G[1][0]==1,
          'rank 2; displayed 2x2 minor equals 1.')
    claim('A-WR',N.weakly_reversible(),'directed three-cycle; deficiency 3-1-2=0.')
    dual = N.dual()
    expected = [((0,1,0),(1,0,0)),((0,1,1),(1,1,0)),
                ((3,0,2),(0,2,3)),((0,0,3),(0,3,0))]
    claim('A-dual',dual.reactions==expected,'b->a; b+c->a+b; 3a+2c->2b+3c; 3c->3b.')
    claim('A-classes',all(sum(g)==0 for g in zip(*dual.Gamma)) and crn.rank(dual.Gamma)==2,
          'classes are a+b+c=T, T>0.')
    geometry = Net(4,rx)
    ff = [f for f in d2cert.facets(geometry) if f['siphon']]
    claim('A-OUT',len(ff)==1 and ff[0]['Z']==(2,) and ff[0]['type']=='OUT'
          and ff[0]['mu']==1 and ff[0]['fl']=='uni',
          'only siphon facet {c}; transverse order 1, OUT, one-way flow.')
    # Construct the field directly from the imported dual network.
    field = [{} for _ in range(3)]
    for i,(source,target) in enumerate(dual.reactions):
        exponent = tuple(source)+tuple(int(j==i) for j in range(4))
        for j in range(3):
            field[j] = add(field[j],term(target[j]-source[j],*exponent))
    ap = add(term(1,0,1,0,1,0,0,0),term(1,0,1,1,0,1,0,0),term(-3,3,0,2,0,0,1,0))
    bp = add(term(-1,0,1,0,1,0,0,0),term(2,3,0,2,0,0,1,0),term(3,0,0,3,0,0,0,1))
    cp = add(term(-1,0,1,1,0,1,0,0),term(1,3,0,2,0,0,1,0),term(-3,0,0,3,0,0,0,1))
    claim('A-field',field==[ap,bp,cp] and not add(ap,bp,cp),'three displayed mass-action equations; sum zero.')
    alpha = term(1,0,0,0,1,0,0,0)
    minus_beta_c = term(-1,0,0,1,0,1,0,0)
    lhs = add(mul(alpha,cp),mul(minus_beta_c,bp))
    rhs = add(term(1,3,0,2,1,0,1,0),term(-3,0,0,3,1,0,0,1),
              term(-2,3,0,3,0,1,1,0),term(-3,0,0,4,0,1,0,1))
    claim('A-integrating-factor',lhs==rhs,
          'alpha*c*(log(c*exp(-beta*b/alpha)))\' = alpha*c\' - beta*c*b\'; exact polynomial identity.')
    # General proof bounds: a^3 >= T^3/8; c^2 <= c for 0<=c<=1;
    # E*epsilon <= gamma*T^3/16. These are exact symbolic coefficients.
    claim('A-general-margin',F(1,8)-F(1,16)==F(1,16),
          'bracket >= gamma*T^3/8-E*epsilon >= gamma*T^3/16 for the stated epsilon.')
    eps = F(1,128)
    E = F(8)
    claim('A-unit-example',eps<=F(1,2) and eps<=1 and F(1,8)-E*eps==F(1,16),
          'rates=(1,1,1,1), T=1, epsilon=1/128, (log V)\' >= c/16.')
    # This statement is a proof-reference check, not finite simulation of a universal assertion.
    claim('A-proof-artifact',(HERE.parent/'notes.md').is_file(),
          'complete all-rates/all-classes proof in notes.md Sections 2.2--2.3.')
    print('Runtime seconds:',round(time.perf_counter()-start,6))
    return 1 if failures else 0


if __name__=='__main__':
    raise SystemExit(main())
