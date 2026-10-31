/* Exploration: base K_N, a hubs on blocks of size b, N0 hubless base vertices, tree W with m vertices, sum of depths T,
   glued at (mode 0) hub h_1  or (mode 1) a hubless base vertex.
   Formula: Phi = delta - m rho + T = 0 (poles at N, N-1, N+1 as in notes).  Enumerate with P = a-1+m multiple of step.
   usage: x5_search mode bmin bmax amax N0max target */
#include <stdio.h>
#include <stdlib.h>
typedef __int128 i128;
typedef long long ll;
static ll gcdll(ll a, ll b){ while(b){ ll t=a%b; a=b; b=t; } return a<0?-a:a; }
int main(int c, char**v){
  int mode=atoi(v[1]); ll bmin=atoll(v[2]), bmax=atoll(v[3]), amax=atoll(v[4]), N0max=atoll(v[5]); double target=atof(v[6]);
  for(ll b=bmin;b<=bmax;b++){
    ll B=b*(b-1);
    for(ll a=2;a<=amax;a++){
      for(ll N0=(mode==1?1:0);N0<=N0max;N0++){
        ll N=a*b+N0, ap=a+N0;
        ll good=(mode==0?(a-1)*b:a*b);
        double nmaxd=good/target;
        ll mmax=(ll)(nmaxd-(N+a));
        if(mmax<0) break;
        ll g=gcdll(N+1,1+ap), step=(N+1)/g;
        ll P0=a-1; ll first=((P0+step-1)/step)*step;
        for(ll P=first;P<=P0+mmax;P+=step){
          ll m=P-a+1;
          i128 Yh=(i128)a*(b-1), Yph=Yh+1;
          i128 S,polyM;
          if(mode==0){ S=(i128)B*(P+1)*(ap-2)+(i128)m*(b-1)*(a-2); polyM=Yh+B-(N+a-1)+m*(b-1)-m; }
          else { S=(i128)B*((i128)(P+1)*(ap-2)+(i128)m*a); polyM=Yh+B-(N+a-1)-m; }
          i128 EN=Yh*P+(i128)B*(1+ap)*P+S, ENm=Yph*P+S, ENp=(i128)B*(1+ap)*P;
          /* BT = -polyM - EN/N + ENm/(N-1) + ENp/(N+1) */
          i128 frN=EN%N; if(frN<0)frN+=N;
          /* quick filter: need sum of fractions integer */
          i128 D=(i128)N*(N-1)*(N+1);
          i128 num=-(EN%D)*(N-1)%D*(N+1)%D;
          /* careful with overflow: compute residues modulo each factor instead */
          i128 q0=EN/N, r0=EN%N, q1=ENm/(N-1), r1=ENm%(N-1), q2=ENp/(N+1), r2=ENp%(N+1);
          i128 R=-(r0*(i128)(N-1)*(N+1))+r1*(i128)N*(N+1)+r2*(i128)N*(N-1);
          if(R%D!=0) continue;
          i128 BT=-polyM-q0+q1+q2+R/D;
          if(BT%B!=0) continue;
          i128 T=BT/B;
          (void)num;(void)frN;
          if(m==0){ if(T!=0) continue; }
          else if(T<m || T>(i128)m*(m+1)/2) continue;
          ll n=N+a+m;
          printf("mode %d b %lld a %lld N0 %lld m %lld T %lld n %lld good>= %lld prop %.5f\n",mode,b,a,N0,m,(ll)T,n,good,(double)good/n);
          fflush(stdout);
        }
      }
    }
  }
  return 0;
}
