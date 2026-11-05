// Computer-assisted proof of a periodic orbit of the mixed-mechanism network in a compatibility class.
//
// Reduced 6-dim system in z = (x3, x4, x5, x7, x8, x9) with
//   x1 = S - x3 - x4 - x5 - x7 - x8 - x9,  x2 = K - x3 - x4,  x6 = P - x7 - x9   (S, K, P = Stot, Ktot, Ptot).
// All 13 parameters (k1..k10, K, P, S) are intervals read from a parameter file (lines "name lo hi", hex floats),
// so a single run proves the statement for every parameter value in the box.
// Method: interval Krawczyk operator for the Poincare map P on the coordinate section {z_j = c}, with
// C^1 enclosures from the CAPD rigorous Taylor integrator (Lohner algorithm, C1Rect2Set).
// If K(X) is in the interior of X, P has a unique fixed point in X, i.e. there is a periodic orbit through X.
// Optional ("stable"): bound the spectral radius of every matrix in the enclosure [DP(X)] by
// || (S^{-1} [DP] S)^m ||_inf^{1/m} < 1  (m = 2^s by repeated squaring), which proves that all nontrivial Floquet
// multipliers lie in the open unit disc, i.e. the orbit is orbitally asymptotically stable.
//
// Build:  g++ po_proof.cpp `capd-config --cflags --libs` -o po_proof     (CAPD 6.1, g++ 15, WSL Fedora)
// Usage:  po_proof paramfile j c dir z1..z6 rad order mode(stable|existence) label
// Runtime: under one second per orbit.
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>
#include <map>
#include "capd/capdlib.h"
using namespace capd;
using namespace std;

static const char* FIELD =
  "par:k1,k2,k3,k4,k5,k6,k7,k8,k9,k10,K,P,S;var:a,b,c,d,e,f;fun:"
  "k1*(S-a-b-c-d-e-f)*(K-a-b)-(k2+k3)*a,"   // x3' = k1 x1 x2 - (k2+k3) x3
  "k3*a-k4*b,"                              // x4' = k3 x3 - k4 x4
  "k4*b-k5*c*(P-d-f)+k6*d,"                 // x5' = k4 x4 - k5 x5 x6 + k6 x7
  "k5*c*(P-d-f)-(k6+k7)*d,"                 // x7' = k5 x5 x6 - (k6+k7) x7
  "k7*d-k8*(P-d-f)*e+k9*f,"                 // x8' = k7 x7 - k8 x6 x8 + k9 x9
  "k8*(P-d-f)*e-(k9+k10)*f;";               // x9' = k8 x6 x8 - (k9+k10) x9
static const char* PNAMES[13] = {"k1","k2","k3","k4","k5","k6","k7","k8","k9","k10","K","P","S"};

static double vnorm(const IVector& v) { double r = 0; for (int i = 0; i < (int)v.dimension(); ++i) r = std::max(r, abs(v[i]).rightBound()); return r; }
static double mnorm(const IMatrix& A) { double r = 0; for (int i = 0; i < (int)A.numberOfRows(); ++i) { interval s(0.); for (int k = 0; k < (int)A.numberOfColumns(); ++k) s += abs(A[i][k]); r = std::max(r, s.rightBound()); } return r; }

int main(int argc, char* argv[]) {
  if (argc < 15) { cerr << "usage: see header\n"; return 2; }
  cout.precision(17);
  std::map<std::string, interval> prm;
  { ifstream in(argv[1]); string name, lo, hi;
    while (in >> name >> lo >> hi) {
      double l = strtod(lo.c_str(), nullptr), h = strtod(hi.c_str(), nullptr);
      if (!(l <= h)) { cerr << "bad interval for " << name << endl; return 2; }
      prm[name] = interval(l, h); } }
  for (int i = 0; i < 13; ++i) if (!prm.count(PNAMES[i])) { cerr << "missing parameter " << PNAMES[i] << endl; return 2; }
  int j = atoi(argv[2]); double c = atof(argv[3]); int dir = atoi(argv[4]);
  DVector z(6); for (int i = 0; i < 6; ++i) z[i] = atof(argv[5 + i]);
  double rad = atof(argv[11]); int order = atoi(argv[12]); string mode = argv[13]; string label = argv[14];
  const int n = 6, m = 5;
  int idx[5]; { int t = 0; for (int i = 0; i < n; ++i) if (i != j) idx[t++] = i; }
  poincare::CrossingDirection cd = dir > 0 ? poincare::MinusPlus : poincare::PlusMinus;
  bool allok = true;
  try {
    // ---------- non-rigorous refinement (Newton for P(y) = y in the 5 section coordinates)
    DMap dvf(FIELD);
    for (int i = 0; i < 13; ++i) dvf.setParameter(PNAMES[i], prm[PNAMES[i]].mid().leftBound());
    DOdeSolver dsolver(dvf, 20); DCoordinateSection dsec(n, j, c); DPoincareMap dpm(dsolver, dsec, cd);
    z[j] = c;
    for (int it = 0; it < 8; ++it) {
      DMatrix mono(n, n); double rt;
      DVector Pz = dpm(z, mono, rt);
      DMatrix DP = dpm.computeDP(Pz, mono, rt);
      DMatrix A(m, m); DVector F(m);
      for (int r = 0; r < m; ++r) { F[r] = Pz[idx[r]] - z[idx[r]]; for (int s = 0; s < m; ++s) A[r][s] = DP[idx[r]][idx[s]] - (r == s ? 1. : 0.); }
      DVector dz = matrixAlgorithms::gauss(A, F);
      for (int r = 0; r < m; ++r) z[idx[r]] -= dz[r];
      if (it == 7) cout << label << ": Newton residual " << euclNorm(F) << ", return time " << rt << endl;
    }
    // ---------- rigorous part
    IMap ivf(FIELD);
    for (int i = 0; i < 13; ++i) { ivf.setParameter(PNAMES[i], prm[PNAMES[i]]);
      cout << label << ": parameter " << PNAMES[i] << " in " << prm[PNAMES[i]] << endl; }
    IOdeSolver isolver(ivf, order); ICoordinateSection isec(n, j, interval(c)); IPoincareMap ipm(isolver, isec, cd);
    IVector zc(n), X(n);
    for (int i = 0; i < n; ++i) { zc[i] = interval(z[i]); X[i] = zc[i]; }
    zc[j] = interval(c); X[j] = interval(c);
    for (int r = 0; r < m; ++r) { int i = idx[r]; double s = rad * (1.0 + fabs(z[i])); X[i] = zc[i] + interval(-s, s); }
    // image of the centre
    C0HORect2Set s0(zc); interval rt0;
    IVector Pc = ipm(s0, rt0);
    // C^1 enclosure on the box
    C1Rect2Set s1(X); IMatrix mono(n, n); interval rt1;
    IVector PX = ipm(s1, mono, rt1);
    IMatrix DP = ipm.computeDP(PX, mono, rt1);
    IVector Fc(m), Xr(m), cr(m); IMatrix DF(m, m), DPr(m, m);
    for (int r = 0; r < m; ++r) {
      Fc[r] = Pc[idx[r]] - zc[idx[r]]; Xr[r] = X[idx[r]]; cr[r] = zc[idx[r]];
      for (int s = 0; s < m; ++s) { DPr[r][s] = DP[idx[r]][idx[s]]; DF[r][s] = DPr[r][s] - (r == s ? interval(1.) : interval(0.)); }
    }
    DMatrix Cd(m, m);
    { DMatrix mid(m, m); for (int r = 0; r < m; ++r) for (int s = 0; s < m; ++s) mid[r][s] = DF[r][s].mid().leftBound();
      Cd = matrixAlgorithms::gaussInverseMatrix(mid); }
    IMatrix C(m, m); for (int r = 0; r < m; ++r) for (int s = 0; s < m; ++s) C[r][s] = interval(Cd[r][s]);
    IMatrix Id = IMatrix::Identity(m);
    IVector K = cr - C * Fc + (Id - C * DF) * (Xr - cr);
    bool kraw = subsetInterior(K, Xr);
    cout << label << ": return time in " << rt1 << endl;
    cout << label << ": |P(c)-c| <= " << vnorm(Fc) << ", box radius (max) " << vnorm(Xr - cr) << endl;
    cout << label << ": ||I - C DF(X)||_inf <= " << mnorm(Id - C * DF) << endl;
    cout << (kraw ? "PASS " : "FAIL ") << label << ": Krawczyk K(X) in int X => unique periodic orbit through X" << endl;
    allok = allok && kraw;
    // positivity of the concentrations on the box (the orbit passes through X)
    interval x1 = prm["S"] - X[0] - X[1] - X[2] - X[3] - X[4] - X[5];
    interval x2 = prm["K"] - X[0] - X[1];
    interval x6 = prm["P"] - X[3] - X[5];
    bool pos = x1 > 0 && x2 > 0 && x6 > 0; for (int i = 0; i < n; ++i) pos = pos && (X[i] > 0);
    cout << (pos ? "PASS " : "FAIL ") << label << ": all nine concentrations positive on X" << endl;
    allok = allok && pos;
    if (mode == "stable") {
      // diagonal scaling S = diag(1 + |z_i|) and repeated squaring of S^{-1} [DP] S
      IMatrix B(m, m);
      for (int r = 0; r < m; ++r) for (int s = 0; s < m; ++s) B[r][s] = DPr[r][s] * interval(1.0 + fabs(z[idx[s]])) / interval(1.0 + fabs(z[idx[r]]));
      int steps = 0; interval nb;
      for (steps = 1; steps <= 12; ++steps) {
        B = B * B;
        nb = interval(0.);
        for (int r = 0; r < m; ++r) { interval rs(0.); for (int s = 0; s < m; ++s) rs += abs(B[r][s]); nb = (rs.rightBound() > nb.rightBound()) ? rs : nb; }
        if (nb < 1) break;
      }
      bool st = nb < 1;
      cout << label << ": || (S^-1 DP S)^(2^" << steps << ") ||_inf <= " << nb.rightBound() << endl;
      cout << (st ? "PASS " : "FAIL ") << label << ": spectral radius of DP < 1 on X => orbitally asymptotically stable" << endl;
      allok = allok && st;
    }
    cout << label << ": orbit point (section coordinates) X = " << X << endl;
  } catch (exception& e) {
    cout << "FAIL " << label << ": exception " << e.what() << endl;
    return 1;
  }
  return allok ? 0 : 1;
}
