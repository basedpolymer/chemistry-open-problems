"""Minimal steady-state solver for CMO equilibrium-stage networks (class A of class.md) and a check of the
stage-level bound  R >= sum_p max(0, phi_p(theta_j))  (Agents S and F).

Agent F, Claude Opus 5.5 (effort xhigh). Environment: Python 3 + numpy (WSL Fedora sci env). Runtime: ~1 min.
Deterministic (fixed seed). Prints PASS/FAIL per test and exits non-zero on failure. Floating point: evidence only;
the bound itself is proved analytically (class.md §5, total-result.md).

Model. Units: 'E' equilibrium stage (outflows: liquid x, vapor y = alpha*x/S(x)), 'C' total condenser (vapor in,
liquid out, same composition), 'R' total vaporizer, 'M' single-phase mixer/splitter. All stream flows are given
(CMO: the flows are design/operating variables); the compositions are computed by successive substitution on
S_u = sum_i alpha_i x_{u,i} (per component the balances are linear once the S_u are fixed).
A stream is (src, dst, phase, flow): src is a unit index or 'F' (feed portion, flow = fraction of the feed),
dst is a unit index or ('P', p) (product stream assigned to component p). For a unit, the outflow of a phase is
the total of its outgoing streams of that phase. Duties follow from the CMO energy balance
Q_u = sum_out v - sum_in v (feed portion f brings vapor (1-q) f).
"""
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla


class Net:
    def __init__(self, alpha, z, q):
        self.alpha = np.asarray(alpha, float)
        self.z = np.asarray(z, float)
        self.q = float(q)
        self.N = len(alpha)
        self.types = []
        self.streams = []

    def unit(self, t):
        self.types.append(t)
        return len(self.types) - 1

    def s(self, src, dst, phase, flow):
        self.streams.append((src, dst, phase, float(flow)))

    # ---------------------------------------------------------------- solve compositions
    def solve(self, iters=2000, tol=1e-13):
        a, N, U = self.alpha, self.N, len(self.types)
        S = np.full(U, a.mean())
        for it in range(iters):
            X = np.zeros((U, N))
            Lout = np.zeros(U); Vout = np.zeros(U)
            for (s_, d, ph, f) in self.streams:
                if isinstance(s_, int):
                    (Lout if ph == 'L' else Vout)[s_] += f
            for i in range(N):
                rows, cols, vals = [], [], []
                b = np.zeros(U)
                for u in range(U):
                    k = a[i] / S[u] if self.types[u] == 'E' else 1.0
                    rows.append(u); cols.append(u); vals.append(Lout[u] + Vout[u] * k)
                for (s_, d, ph, f) in self.streams:
                    if isinstance(d, int):
                        if s_ == 'F':
                            b[d] += f * self.z[i]
                        else:
                            k = a[i] / S[s_] if (self.types[s_] == 'E' and ph == 'V') else 1.0
                            rows.append(d); cols.append(s_); vals.append(-f * k)
                M = sp.csc_matrix((vals, (rows, cols)), shape=(U, U))
                X[:, i] = spla.spsolve(M, b)
            Snew = X @ a / X.sum(axis=1)
            if np.max(np.abs(Snew - S)) < tol:
                S = Snew
                break
            S = 0.5 * S + 0.5 * Snew
        self.X = X / X.sum(axis=1, keepdims=True)
        self.S = self.X @ a
        self.closure = np.max(np.abs(X.sum(axis=1) - 1))
        return it


    def newton(self, iters=100, tol=1e-12):
        """Newton polish on the full balance system (unknowns x_{u,i}; equations: N-1 component balances and
        sum x = 1 per unit). Uses the analytic Jacobian. Call after solve()."""
        a, N, U = self.alpha, self.N, len(self.types)
        x = self.X.copy()
        Lout = np.zeros(U); Vout = np.zeros(U)
        for (s_, d, ph, f) in self.streams:
            if isinstance(s_, int):
                (Lout if ph == 'L' else Vout)[s_] += f
        def yv(u, xu):
            if self.types[u] != 'E':
                return xu, np.eye(N)
            S = a @ xu
            y = a * xu / S
            J = np.diag(a / S) - np.outer(a * xu, a) / S**2
            return y, J
        err = np.inf
        for it in range(iters):
            F = np.zeros((U, N))
            blocks = {}
            def add(r, c, B):
                if (r, c) in blocks:
                    blocks[(r, c)] = blocks[(r, c)] + B
                else:
                    blocks[(r, c)] = B.copy()
            Y = []; JY = []
            for u in range(U):
                y, J = yv(u, x[u]); Y.append(y); JY.append(J)
            for u in range(U):
                F[u] -= x[u] * Lout[u] + Vout[u] * Y[u]
                add(u, u, -(Lout[u] * np.eye(N) + Vout[u] * JY[u]))
            for (s_, d, ph, f) in self.streams:
                if not isinstance(d, int):
                    continue
                if s_ == 'F':
                    F[d] += f * self.z
                elif ph == 'V':
                    F[d] += f * Y[s_]
                    add(d, s_, f * JY[s_])
                else:
                    F[d] += f * x[s_]
                    add(d, s_, f * np.eye(N))
            for u in range(U):
                F[u, N-1] = x[u].sum() - 1
            rows, cols, vals = [], [], []
            for (r, c), B in blocks.items():
                B = B.copy()
                B[N - 1, :] = 1.0 if r == c else 0.0
                rr, cc = np.nonzero(B)
                rows.extend(r * N + rr); cols.extend(c * N + cc); vals.extend(B[rr, cc])
            Jm = sp.csc_matrix((vals, (rows, cols)), shape=(U * N, U * N))
            err = np.max(np.abs(F))
            if err < tol:
                break
            dx = spla.spsolve(Jm, -F.reshape(-1)).reshape(U, N)
            lam = 1.0
            while lam > 1e-4 and np.any(x + lam * dx < -1e-12):
                lam *= 0.5
            x = x + lam * dx
        self.X = x
        self.S = x @ a
        self.closure = err
        return it

    # ---------------------------------------------------------------- stream data
    def comp(self, src, ph):
        if src == 'F':
            return self.z
        x = self.X[src]
        if self.types[src] == 'E' and ph == 'V':
            return self.alpha * x / (self.alpha @ x)
        return x

    def phi(self, e, th):
        src, dst, ph, f = e
        n = f * self.comp(src, ph)
        if src == 'F':
            v = (1 - self.q) * f
        else:
            v = f if ph == 'V' else 0.0
        return v - np.sum(self.alpha * n / (self.alpha - th))

    def duties(self):
        U = len(self.types)
        Q = np.zeros(U)
        for (s_, d, ph, f) in self.streams:
            vin = (1 - self.q) * f if s_ == 'F' else (f if ph == 'V' else 0.0)
            if isinstance(d, int):
                Q[d] -= vin
            if isinstance(s_, int):
                Q[s_] += f if ph == 'V' else 0.0
        return Q

    def feed_roots(self):
        a, z, q = self.alpha, self.z, self.q
        roots = []
        for j in range(self.N - 1):
            lo, hi = a[j + 1], a[j]
            g = lambda t: np.sum(a * z / (a - t)) - (1 - q)
            l, h = lo + 1e-14 * (hi - lo), hi - 1e-14 * (hi - lo)
            for _ in range(200):
                m = 0.5 * (l + h)
                if g(m) > 0:   # g increasing on the interval
                    h = m
                else:
                    l = m
            roots.append(0.5 * (l + h))
        return np.array(roots)

    def check(self, verbose=True):
        """Return (R, list of (j, bound_j, H_j), max node-balance residual, sign-uniformity violations)."""
        a, z, q = self.alpha, self.z, self.q
        Q = self.duties()
        R = np.sum(np.maximum(Q, 0))
        th = self.feed_roots()
        res, viol, out = 0.0, 0, []
        for j, t in enumerate(th):
            U = len(self.types)
            bal = -Q.copy()
            signs = [set() for _ in range(U)]
            for e in self.streams:
                p = self.phi(e, t)
                if isinstance(e[0], int):
                    bal[e[0]] += p
                    if abs(p) > 1e-12:
                        signs[e[0]].add(np.sign(p))
                if isinstance(e[1], int):
                    bal[e[1]] -= p
            res = max(res, np.max(np.abs(bal)))
            viol += sum(1 for s_ in signs if len(s_) > 1)
            bound = sum(max(0.0, self.phi(e, t)) for e in self.streams if not isinstance(e[1], int))
            H = np.sum(a[j + 1:] * z[j + 1:] / (t - a[j + 1:]))
            out.append((j + 1, bound, H))
        if verbose:
            pur = []
            for e in self.streams:
                if not isinstance(e[1], int):
                    p = e[1][1]
                    pur.append(self.comp(e[0], e[2])[p])
            print(f"  R = {R:.6f}; min product purity = {min(pur):.6f}; node-balance residual {res:.1e}; "
                  f"sign-uniformity violations {viol}; closure {self.closure:.1e}")
            for (j, b, H) in out:
                print(f"    j={j}: bound sum max(0,phi_p) = {b:.6f}  (H_j = {H:.6f})  R - bound = {R - b:.3e}"
                      f"  {'PASS' if R >= b - 1e-9 else 'FAIL'}")
        return R, out, res, viol


def column(net, nT, nB, L_T, V_T, feed_frac=1.0, feed_src='F', top='total', bottom='partial',
           pd=0, pb=None):
    """Conventional column: nT top stages, feed stage, nB bottom stages, partial reboiler (an E stage)
    and a total condenser. Returns the unit indices. CMO: L_B = L_T + q f, V_B = V_T - (1-q) f."""
    q = net.q
    f = feed_frac
    L_B = L_T + q * f
    V_B = V_T - (1 - q) * f
    top_st = [net.unit('E') for _ in range(nT)]          # top_st[0] = top stage
    fs = net.unit('E')
    bot_st = [net.unit('E') for _ in range(nB)]           # bot_st[-1] = lowest tray
    reb = net.unit('E')
    cond = net.unit('C')
    D = V_T - L_T
    B = L_B - V_B
    chain = top_st + [fs] + bot_st + [reb]
    # vapor up, liquid down
    for k in range(len(chain) - 1):
        up, dn = chain[k], chain[k + 1]
        Lk = L_T if k < nT else L_B
        Vk = V_T if k < nT else V_B
        net.s(up, dn, 'L', Lk)
        net.s(dn, up, 'V', Vk)
    net.s(chain[0], cond, 'V', V_T)
    net.s(cond, chain[0], 'L', L_T)
    net.s(cond, ('P', pd), 'L', D)
    net.s(feed_src, fs, 'L' if feed_src != 'F' else 'F', f)
    net.s(reb, ('P', pb if pb is not None else net.N - 1), 'L', B)
    return dict(top=top_st, feed=fs, bot=bot_st, reb=reb, cond=cond, D=D, B=B)


def random_network(rng, N=3, U=8, ncyc=10, npath=4, q=None):
    """Random arrangement of E (and some C, R, M) units: random circulations (recycles, bypasses, reverse flows)
    plus random feed->product paths, random phases. Total mass balance holds at every unit by construction;
    duties follow from CMO. Products are whatever leaves; the bound sum_out max(0, phi) needs no purity."""
    alpha = np.sort(rng.uniform(1, 12, N))[::-1]
    alpha[-1] = 1.0
    z = rng.dirichlet(np.ones(N))
    if q is None:
        q = rng.uniform(-0.5, 1.5)
    net = Net(alpha, z, q)
    types = ['E'] * U
    for u in range(U):
        r = rng.uniform()
        types[u] = 'E' if r < 0.7 else ('C' if r < 0.8 else ('R' if r < 0.9 else 'M'))
    types[0] = 'E'
    for t in types:
        net.unit(t)
    inph = {u: None for u in range(U)}   # required inflow phase for C (V), R (L), M (fixed)
    outph = {u: None for u in range(U)}
    for u, t in enumerate(types):
        if t == 'C':
            inph[u], outph[u] = 'V', 'L'
        elif t == 'R':
            inph[u], outph[u] = 'L', 'V'
        elif t == 'M':
            ph = 'L' if rng.uniform() < 0.5 else 'V'
            inph[u], outph[u] = ph, ph
    flows = {}

    def add(a, b, f):
        # phase of the stream a->b must match a's outflow phase and b's inflow phase
        pa, pb = outph[a], inph[b]
        if pa and pb and pa != pb:
            return False
        ph = pa or pb or ('L' if rng.uniform() < 0.5 else 'V')
        flows[(a, b, ph)] = flows.get((a, b, ph), 0.0) + f
        return True

    for _ in range(ncyc):
        k = rng.integers(2, U + 1)
        cyc = [int(c) for c in rng.choice(U, size=k, replace=False)]
        f = rng.exponential(1.0)
        ok = all((outph[cyc[i]] is None or inph[cyc[(i + 1) % k]] is None or
                  outph[cyc[i]] == inph[cyc[(i + 1) % k]]) for i in range(k))
        if not ok:
            continue
        for i in range(k):
            add(cyc[i], cyc[(i + 1) % k], f)
    # feed portions and product paths
    fr = rng.dirichlet(np.ones(npath))
    for m in range(npath):
        k = rng.integers(1, U + 1)
        path = [int(c) for c in rng.choice(U, size=k, replace=False)]
        if not all((outph[path[i]] is None or inph[path[i + 1]] is None or outph[path[i]] == inph[path[i + 1]])
                   for i in range(k - 1)):
            path = [0]
        if types[path[0]] in ('C', 'R', 'M'):
            path = [0] + [p for p in path if p != 0]
            if not all((outph[path[i]] is None or inph[path[i + 1]] is None or
                        outph[path[i]] == inph[path[i + 1]]) for i in range(len(path) - 1)):
                path = [0]
        net.s('F', path[0], 'F', fr[m])
        for i in range(len(path) - 1):
            add(path[i], path[i + 1], fr[m])
        last = path[-1]
        ph = outph[last] or ('L' if rng.uniform() < 0.5 else 'V')
        net.s(last, ('P', int(rng.integers(N))), ph, fr[m])
    for (a, b, ph), f in flows.items():
        net.s(a, b, ph, f)
    return net


def random_tests(n=300, seed=1):
    rng = np.random.default_rng(seed)
    worst, done, fails = np.inf, 0, 0
    for t in range(n):
        net = random_network(rng, N=int(rng.integers(2, 6)), U=int(rng.integers(2, 12)))
        try:
            net.solve(iters=500)
            net.newton()
        except (np.linalg.LinAlgError, RuntimeError):
            continue
        if not np.isfinite(net.closure) or net.closure > 1e-9 or np.any(net.X < -1e-9):
            continue
        R, out, res, viol = net.check(verbose=False)
        if res > 1e-8 or viol:
            fails += 1
            print("  residual/sign problem", res, viol)
            continue
        m = min(R - b for (j, b, H) in out)
        worst = min(worst, m)
        done += 1
        if m < -1e-9:
            fails += 1
            print("  FAIL", t, m)
    print(f"Random networks: {done} solved steady states (of {n}); min over j and networks of R - bound_j = "
          f"{worst:.3e}; failures {fails}  {'PASS' if fails == 0 else 'FAIL'}")
    return fails == 0 and done > 0


def main():
    ok = True
    print("Test 1: binary column, alpha=(2,1), z=(.5,.5), q=1, V=1.6 (> King 1.5), 40+40 stages")
    net = Net([2, 1], [.5, .5], 1.0)
    column(net, 40, 40, L_T=1.1, V_T=1.6, pd=0, pb=1)
    net.solve(); net.newton()
    R, out, res, viol = net.check()
    ok &= res < 1e-9 and viol == 0 and all(R >= b - 1e-9 for (_, b, _) in out)
    print("Test 2: same, V lowered to 1.4 < V_min by continuation, 80+80 stages: impure products; the bound holds and"
          " is tight to roundoff (feed pinch)")
    X0 = None
    for V in np.linspace(1.6, 1.4, 21):
        net = Net([2, 1], [.5, .5], 1.0)
        column(net, 80, 80, L_T=V - 0.5, V_T=V, pd=0, pb=1)
        if X0 is None:
            net.solve()
        else:
            net.X = X0
        net.newton()
        X0 = net.X
    R, out, res, viol = net.check()
    ok &= res < 1e-9 and viol == 0 and all(R >= b - 1e-9 for (_, b, _) in out) and abs(R - out[0][1]) < 1e-6
    print("Test 3: random networks (recycles, bypasses, splitters, total condensers/vaporizers, any q)")
    ok &= random_tests(400, seed=7)
    print("ALL PASS" if ok else "SOME FAIL")
    return ok


if __name__ == '__main__':
    import sys
    sys.exit(0 if main() else 1)
