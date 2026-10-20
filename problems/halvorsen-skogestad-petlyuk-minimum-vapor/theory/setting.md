# Setting (Agent T) — definitions used by the theorems

**Version 2 (10:45): the paper (paper.tex, Section 2) is now the reference.** It widens the class to
"two-section column networks" (any number of columns per label, any number of incoming links per node, side
draws, TC or HX on any link, heaters/coolers at any node, interval streams, one feed stage per column, nonnegative
net flows (H1), single feed and saturated-liquid pure products (H2)); regular-column configurations (sec. 2 below)
are a subclass. Agent G's results use the non-degeneracy assumption (N): every column receives a positive flow of
every component of its label.

Version 1 (2026-10-20). Agent C's model (computation/scripts/nlp.py) is a special case of the model below
(it forces V = threshold at interior roots); every lower bound proved here for model (TG) therefore also holds
for Agent C's model.

## 1. Data
Components 1..N, constant relative volatilities a_1 > a_2 > ... > a_N > 0; feed F = 1 with z_p > 0,
sum z_p = 1, liquid fraction q (any real q; physically q in [0,1]). Constant molar overflow, constant pressure,
infinitely many stages in every section. I_r := (a_{r+1}, a_r), r = 1..N-1.
Feed roots: theta_r in I_r, the unique root in I_r of  sum_p a_p z_p/(a_p - theta) = 1 - q.
Peaks: P_r^T = sum_{p<=r} a_p z_p/(a_p - theta_r) (top-referenced),
       H_r   = sum_{p>r} a_p z_p/(theta_r - a_p)  = P_r^T - (1-q)  (= reboiler-referenced peak).
Petlyuk (FTC) value of total reboiler vapor: V_FTC = max_r H_r = max_r P_r^T - (1-q).

## 2. Competitor class: regular-column configurations (Agrawal; Shah-Agrawal 2010; TG = arXiv:2010.12113)
Streams are intervals [i,j]. A configuration is a set S of mixture streams containing [1,N] such that the
split of each [i,j] in S is  distillate [i,k], k = max{k'<j : [i,k'] in S or k'=i},
residue [l,j], l = min{l'>i : [l',j] in S or l'=j}, with l <= k+1, and every stream other than [1,N]
(pure products included) is produced by at least one split. (152 configurations for N = 4, counting the
transfer-stream choices below.) Each stream is the distillate of at most one stream and the residue of at most
one stream.
Transfer streams: a submixture produced by one split is either thermally coupled (TC: two-way vapor/liquid
link) or has a heat exchanger (HX): a partial condenser at the top of the producing section (removes the
reflux L, the stream leaves as saturated vapor of amount D) if it is a distillate, a reboiler at the bottom
(generates the boil-up V, the stream leaves as saturated liquid of amount B) if it is a residue.
A stream produced by two splits is a side draw (TC on both sides, no HX). A pure product produced once is
withdrawn through a total condenser / reboiler; produced twice, as a liquid side draw between the two
sections (vapor of the lower section passes into the upper one). All products are saturated liquid.
Objective: R = total vapor generated in all reboilers.

## 3. Model (TG) of a configuration
Column c (stream [i,j], split [i,k]/[l,j]):  top net flows d_p >= 0 (p=i..k), bottom net (downward) flows
b_p >= 0 (p=l..j), feed net flows f_p = d_p[p<=k] + b_p[p>=l] (sum of the producing sections' flows, or z
for [1,N]); vapors V_T, V_B >= 0 with V_T - V_B = U_c, where U_c = net vapor fed:
   (1-q) for [1,N];  + V_T(parent) for a TC/side-draw distillate link, + D(parent) for an HX distillate link;
   - V_B(parent) for a TC/side-draw residue link, 0 for an HX residue link.
Underwood: for every r in {l-1,...,k}: theta_{c,r} in [a_{r+1}, a_r] with sum_p a_p f_p/(a_p - theta_{c,r}) = U_c,
   g_T(theta) := sum_{p<=k} a_p d_p/(a_p - theta),  g_B(theta) := - sum_{p>=l} a_p b_p/(a_p - theta),
   Y_T >= g_T(theta_{c,r}) for r in {l-1, k},  Y_T = g_T(theta_{c,r}) for l <= r <= k-1,
   Y_T - Y_B = U_c,  0 <= Y_T <= V_T,  0 <= Y_B <= V_B.
Pure product produced twice (side draw): V_T(upper producer's top) ... (see Agent C: VT[t] = VB[u]).

## 4. Underwood flux (the basic tool)
For a section s with vapor V_s and upward net flows w_s (w = d for a top section, w = -b for a bottom section):
   Lam_s(theta) := V_s - sum_p a_p w_{s,p}/(a_p - theta).
Lam is linear in (V, w), so it is conserved at every junction (CMO, no heat exchange), except:
condenser duty Q_c removes Q_c, reboiler duty Q_r adds Q_r, a liquid product of pure p removes
-a_p z_p/(a_p - theta) (i.e. injects a_p z_p/(a_p - theta) if p is light w.r.t. theta), and the main feed
injects Lam_F(theta) = (1-q) - sum_p a_p z_p/(a_p - theta), which is 0 at theta = theta_j.
For column c: Lam_T - Lam_B = Lam_{c,F} := U_c - sum_p a_p f_p/(a_p - theta) on I_j; Lam_T is decreasing,
Lam_B increasing, Lam_{c,F} decreasing on each I_r (nonnegative flows).
TG constraints give, for every column c and every r with both r, r+1 in its stream (implied ones included):
   Lam_T(theta_{c,r}) = Lam_B(theta_{c,r}) >= 0.                                            (U)

## 5. Regions for the split j
L = nodes of streams inside [1,j] and light products; H = streams inside [j+1,N] and heavy products;
M = mixed streams [i,m], i <= j < m (incl. HX nodes on mixed transfer streams). Sections cross only
M -> L ("light exits": top sections of mixed columns with distillate [i,k], k <= j) and H -> M ("heavy exits").
Exact identity at theta = theta_j:
   R - H_j = R(L) + Q_cond(H) + Q_cond(M) + sum_{light exits} Lam_s(theta_j)          (*)
           = R(L) + R(M) + Q_cond(H) + sum_{heavy exits} Lam_s(theta_j).
