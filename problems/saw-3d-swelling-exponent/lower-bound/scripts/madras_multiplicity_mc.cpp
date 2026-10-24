// HEURISTIC EXPLORATION ONLY (Monte Carlo, floating point): not a certificate.
// Pivot-algorithm sampling of uniform N-step SAWs on Z^3, and for walks whose endpoint
// satisfies ||w_N||_inf <= m, the exact multiplicity K(T w) of Madras' reflection map T
// (Madras, Canad. Math. Bull. 57 (2014), proof of Lemma 2.3): the number of walks theta in
// S_{N,<=m} with T theta = T w.  Also: renewal-depth counts near the top, and validity of
// transverse pivots at the top point.
// Build: g++ -O2 -o mc madras_multiplicity_mc.cpp ; run: ./mc N samples seed
// Runtime ~ minutes single core for N=1000, 1e5 samples.
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <vector>
#include <array>
#include <random>
#include <algorithm>
#include <cmath>
using namespace std;
typedef array<int,3> P3;

struct HSet { // open addressing hash set with generation stamps
  vector<uint64_t> key; vector<uint32_t> gen; uint32_t g=1; size_t mask;
  HSet(size_t cap){ size_t c=1; while(c<cap*4) c<<=1; key.assign(c,0); gen.assign(c,0); mask=c-1; }
  static uint64_t pk(const P3&p){ return ((uint64_t)(p[0]+(1<<20))<<42)|((uint64_t)(p[1]+(1<<20))<<21)|(uint64_t)(p[2]+(1<<20)); }
  void clear(){ g++; }
  bool insert(const P3&p){ uint64_t k=pk(p); size_t h=(k*0x9E3779B97F4A7C15ULL)>>20 & mask;
    while(gen[h]==g){ if(key[h]==k) return false; h=(h+1)&mask; } gen[h]=g; key[h]=k; return true; }
  bool has(const P3&p) const { uint64_t k=pk(p); size_t h=(k*0x9E3779B97F4A7C15ULL)>>20 & mask;
    while(gen[h]==g){ if(key[h]==k) return true; h=(h+1)&mask; } return false; }
};

// 48 lattice symmetries: permutation + signs
struct Sym { int perm[3]; int sgn[3]; };
vector<Sym> syms;
inline P3 applyS(const Sym&s,const P3&v){ P3 r; for(int i=0;i<3;i++) r[i]=s.sgn[i]*v[s.perm[i]]; return r; }

int N; vector<P3> w;
mt19937_64 rng;
HSet *hs;

bool isSAW(const vector<P3>&a){ hs->clear(); for(auto&p:a) if(!hs->insert(p)) return false; return true; }
int ninf(const P3&p){ return max(abs(p[0]),max(abs(p[1]),abs(p[2]))); }

bool pivot_attempt(){
  int k = 1 + rng()%(N-1); const Sym&s = syms[1 + rng()%47];
  hs->clear(); for(int i=0;i<=k;i++) hs->insert(w[i]);
  static vector<P3> nt; nt.resize(N+1);
  for(int i=k+1;i<=N;i++){ P3 d={w[i][0]-w[k][0],w[i][1]-w[k][1],w[i][2]-w[k][2]}; P3 r=applyS(s,d);
    P3 q={w[k][0]+r[0],w[k][1]+r[1],w[k][2]+r[2]}; if(hs->has(q)) return false; nt[i]=q; }
  for(int i=k+1;i<=N;i++) w[i]=nt[i]; return true;
}

struct MData { int span,L,I,H; };
MData madras(const vector<P3>&a){ MData d; d.span=0; for(auto&p:a) d.span=max(d.span,ninf(p));
  d.L=-1; for(int k=1;k<=N;k++) if(ninf(a[k])==d.span){ d.L=k; break; }
  d.I=-1; for(int i=0;i<3;i++) if(abs(a[d.L][i])==d.span){ d.I=i; break; } d.H=a[d.L][d.I]; return d; }
vector<P3> reflectAfter(const vector<P3>&a,int L,int I,int H){ vector<P3> b=a; for(int k=L+1;k<=N;k++) b[k][I]=2*H-a[k][I]; return b; }

int main(int argc,char**argv){
  N = argc>1?atoi(argv[1]):1000; long samples = argc>2?atol(argv[2]):100000; rng.seed(argc>3?atoll(argv[3]):12345);
  int perms[6][3]={{0,1,2},{0,2,1},{1,0,2},{1,2,0},{2,0,1},{2,1,0}};
  for(int p=0;p<6;p++) for(int sg=0;sg<8;sg++){ Sym s; for(int i=0;i<3;i++){ s.perm[i]=perms[p][i]; s.sgn[i]=((sg>>i)&1)?-1:1; } syms.push_back(s); }
  // identity first
  for(size_t i=0;i<syms.size();i++){ bool id=true; for(int j=0;j<3;j++) if(syms[i].perm[j]!=j||syms[i].sgn[j]!=1) id=false; if(id){ swap(syms[0],syms[i]); break; } }
  hs = new HSet(N+10);
  w.resize(N+1); for(int i=0;i<=N;i++) w[i]={i,0,0};
  long acc=0; for(long t=0;t<40L*N;t++) acc+=pivot_attempt();
  const int NM=6; int ms[NM]={2,4,8,16,32,64};
  long nnear[NM]={0}; double sumK[NM]={0}, sumRt[NM]={0}, sumRh[NM]={0}, sumCand[NM]={0}, sumCandRen[NM]={0}, sumCandSAW[NM]={0};
  long histK[NM][8]={{0}}; long pivok[NM][8]={{0}}; long smallspan[NM]={0}; double sumspan[NM]={0};
  double sumE2=0;
  int stride = 10;
  for(long s=0;s<samples;s++){
    for(int t=0;t<stride;t++) pivot_attempt();
    double e2=0; for(int i=0;i<3;i++) e2+= (double)w[N][i]*w[N][i]; sumE2+=e2;
    int en = ninf(w[N]);
    if(en>ms[NM-1]) continue;
    MData d = madras(w);
    vector<P3> psi = reflectAfter(w,d.L,d.I,d.H);
    int b = psi[N][d.I]; int bn = ninf(psi[N]); int sg = (b>0)?1:-1;
    // renewal depth counts (orientation: top at H, sign sg0 = sign(H))
    int s0 = (d.H>0)?1:-1; int Ha = abs(d.H);
    for(int mi=0;mi<NM;mi++){ int m=ms[mi]; if(en>m) continue; nnear[mi]++; if(d.span<=m) smallspan[mi]++; sumspan[mi]+=d.span;
      // renewal depths of tail and head within m of the top
      int Rt=0,Rh=0;
      for(int j=1;j<=m && j<=2*Ha;j++){ int lev=Ha-j; // in oriented coordinate y = s0*x_I
        int tau=-1; for(int k=d.L;k<=N;k++) if(s0*w[k][d.I]==lev){ tau=k; break; }
        if(tau>=0){ bool ok=true; for(int k=tau;k<=N;k++) if(s0*w[k][d.I]>lev){ ok=false; break; } if(ok) Rt++; }
        int sig=-1; for(int k=0;k<=d.L;k++) if(s0*w[k][d.I]==lev){ sig=k; break; }
        if(sig>=0){ bool ok=true; for(int k=sig;k<=d.L;k++) if(s0*w[k][d.I]<lev){ ok=false; break; } if(ok) Rh++; }
      }
      sumRt[mi]+=Rt; sumRh[mi]+=Rh;
      // exact multiplicity
      int K=0, cand=0, candRen=0, candSAW=0;
      int lo = (bn-m+1)/2; if(lo<1) lo=1; int hi=(bn+m)/2;
      for(int Hs=lo; Hs<=hi; Hs++){ int Hv = sg*Hs; if((b - Hv) % 1) continue;
        // endpoint of theta in coord I would be 2Hv - b; need |.|<=m
        if(abs(2*Hv-b)>m) continue; cand++;
        int Lp=-1; for(int k=0;k<=N;k++) if(psi[k][d.I]==Hv){ Lp=k; break; } if(Lp<0) continue;
        // renewal-type necessary condition: psi[Lp..N] on the far side of the plane
        bool ren=true; for(int k=Lp;k<=N;k++) if(sg*psi[k][d.I]<Hs){ ren=false; break; }
        if(!ren) continue; candRen++;
        vector<P3> th = reflectAfter(psi,Lp,d.I,Hv);
        if(!isSAW(th)) continue; candSAW++;
        MData d2 = madras(th); if(d2.L!=Lp||d2.I!=d.I||d2.H!=Hv) continue;
        if(ninf(th[N])>m) continue;
        vector<P3> chk = reflectAfter(th,d2.L,d2.I,d2.H); bool same=true; for(int k=0;k<=N;k++) if(chk[k]!=psi[k]){ same=false; break; }
        if(same) K++; }
      sumK[mi]+=K; sumCand[mi]+=cand; sumCandRen[mi]+=candRen; sumCandSAW[mi]+=candSAW; histK[mi][min(K,7)]++;
      // transverse pivots at the top point: g in D4 acting on the two coordinates != I, about w[L]
      int a1=(d.I+1)%3, a2=(d.I+2)%3; int gi=0;
      for(int sw=0;sw<2;sw++) for(int s1=-1;s1<=1;s1+=2) for(int s2=-1;s2<=1;s2+=2){
        hs->clear(); for(int k=0;k<=d.L;k++) hs->insert(w[k]); bool ok=true;
        for(int k=d.L+1;k<=N && ok;k++){ int u=w[k][a1]-w[d.L][a1], v=w[k][a2]-w[d.L][a2]; int nu = sw? v:u, nv = sw? u:v; nu*=s1; nv*=s2;
          P3 q=w[k]; q[a1]=w[d.L][a1]+nu; q[a2]=w[d.L][a2]+nv; if(hs->has(q)) ok=false; }
        if(ok) pivok[mi][gi]++; gi++; }
    }
  }
  printf("# N=%d samples=%ld stride=%d  mean |w_N|^2 = %.2f\n",N,samples,stride,sumE2/samples);
  printf("# m  n_near  P(near)  meanK  m+1  meanCand  meanCandRenewal  meanCandSAW  meanRtail  meanRhead  histK[0..7+]  transversePivotValidFrac[g0..g7]  fracSpan<=m meanSpan\n");
  for(int mi=0;mi<NM;mi++){ double n=nnear[mi]; if(n==0){ printf("%d 0\n",ms[mi]); continue; }
    printf("%d %ld %.5f %.3f %d %.3f %.3f %.3f %.3f %.3f  ",ms[mi],nnear[mi],n/samples,sumK[mi]/n,ms[mi]+1,sumCand[mi]/n,sumCandRen[mi]/n,sumCandSAW[mi]/n,sumRt[mi]/n,sumRh[mi]/n);
    for(int k=0;k<8;k++) printf("%ld ",histK[mi][k]); printf(" ");
    for(int g=0;g<8;g++) printf("%.3f ",pivok[mi][g]/n); printf(" fracSpan<=m=%.4f meanSpan=%.1f\n",smallspan[mi]/n,sumspan[mi]/n); }
  return 0;
}
