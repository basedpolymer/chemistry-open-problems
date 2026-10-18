# Independent validated integrator for the parameter-free shooting problem
#     q_ss + q_s/s + g(q) = 0,  q(0) = a,  q_s(0) = 0,   g(q) = exp(-1/q) (q > 0), g = 0 (q <= 0),
# written by the verifier from scratch (Julia + IntervalArithmetic.jl, no code shared with the solver).
#
# Design (deliberately different from the solver's scripts):
#  * Derivatives with respect to a are carried as truncated 3-jets  q(a+d) = Q0 + Q1 d + Q2 d^2 + Q3 d^3
#    (so v = Q1, z = 2 Q2, c = 6 Q3); the variational equations are never written by hand.
#  * Time variable x = log(s^2/4):  q_x = p,  p_x = -e^x g(q).
#  * Taylor method of order N in x with jet-valued coefficients; Lagrange remainder from the N-th
#    coefficient evaluated on a Picard a-priori tube.
#  * Origin: real-variable bootstrap (no complex-analytic estimates), see `initial_state`.
#  * Region q <= 1/40: explicit second-order step with crude bounds on g, g', g'', g''' (valid across q = 0).
#  * Parameter boxes: the box and its exact midpoint are integrated in lockstep and intersected after
#    every step with the mean-value form Q_m(a) in Q_m(amid) + (a - amid)(m+1) Q_{m+1}(box), m = 2,1,0
#    (the jets supply the a-derivatives); plain box integration wraps far too much near q = 0.
#  * Every advance of the time x is checked to be exact in the working floating-point type.

using IntervalArithmetic

const NJ = 4                       # jet length
Jet{T} = NTuple{4,T}

@inline jadd(a::Jet, b::Jet) = (a[1]+b[1], a[2]+b[2], a[3]+b[3], a[4]+b[4])
@inline jsub(a::Jet, b::Jet) = (a[1]-b[1], a[2]-b[2], a[3]-b[3], a[4]-b[4])
@inline jscale(s, a::Jet) = (s*a[1], s*a[2], s*a[3], s*a[4])
@inline jneg(a::Jet) = (-a[1], -a[2], -a[3], -a[4])
@inline jmul(a::Jet, b::Jet) = (a[1]*b[1], a[1]*b[2]+a[2]*b[1], a[1]*b[3]+a[2]*b[2]+a[3]*b[1],
                                a[1]*b[4]+a[2]*b[3]+a[3]*b[2]+a[4]*b[1])
@inline function jinv(a::Jet)
    r1 = inv(a[1])
    r2 = -(a[2]*r1)*r1
    r3 = -(a[2]*r2 + a[3]*r1)*r1
    r4 = -(a[2]*r3 + a[3]*r2 + a[4]*r1)*r1
    (r1, r2, r3, r4)
end
@inline function jexp(s::Jet)
    # k e_k = sum_{j=1}^k j s_j e_{k-j}
    e1 = exp(s[1])
    e2 = s[2]*e1
    e3 = (s[2]*e2 + 2*s[3]*e1)/2
    e4 = (s[2]*e3 + 2*s[3]*e2 + 3*s[4]*e1)/3
    (e1, e2, e3, e4)
end
jzero(::Type{T}) where T = (zero(T), zero(T), zero(T), zero(T))
basetype(::Type{Interval{F}}) where F = F
jhull(a::Jet, b::Jet) = (hull(a[1],b[1]), hull(a[2],b[2]), hull(a[3],b[3]), hull(a[4],b[4]))

# g applied to a jet with q strictly positive
function gjet(q::Jet)
    inf(q[1]) > 0 || error("gjet: q not positive")
    jexp(jneg(jinv(q)))
end

const QFLAT = 1/40        # exact? 1/40 is not dyadic; we only use it through an enclosure
# Bounds for g^{(m)}, m=0..3 on q <= 1/40 (g extended by 0 for q <= 0). Each exp(-1/q) q^{-m}
# increases on (0, 1/m) for m <= 6, so its supremum on (0,1/40] is attained at 1/40.
function flat_bounds(::Type{T}) where T
    q = interval(T, 1)/40
    e = exp(-inv(q)); u = inv(q)
    b0 = sup(e); b1 = sup(e*u^2); b2 = sup(e*(u^4 + 2*u^3)); b3 = sup(e*(u^6 + 6*u^5 + 6*u^4))
    (interval(T, 0, b0), interval(T, 0, b1), interval(T, -b2, b2), interval(T, -b3, b3))
end

# jet of g(q) when the whole q-range lies in q <= 1/40 (possibly negative)
function gjet_flat(q::Jet{T}, fb) where T
    sup(q[1]) <= 1/40 - 1e-12 || error("gjet_flat: q too large")
    g0, g1, g2, g3 = fb
    (g0, g1*q[2], g1*q[3] + g2*q[2]^2/2, g1*q[4] + g2*q[2]*q[3] + g3*q[2]^3/6)
end

struct State{T}
    q::Jet{T}
    p::Jet{T}
end
shull(A::State, B::State) = State(jhull(A.q,B.q), jhull(A.p,B.p))

# vector field (q,p)_x = (p, -e^x g(q))
function field(X, S::State{T}) where T
    State(S.p, jscale(-exp(X), gjet(S.q)))
end
function field_flat(X, S::State{T}, fb) where T
    State(S.p, jscale(-exp(X), gjet_flat(S.q, fb)))
end

# Taylor coefficients (orders 0..N) of the solution through (X, S) (X may be an interval)
function taylor_coeffs(X, S::State{T}, N::Int) where T
    Q = Vector{Jet{T}}(undef, N+1); P = Vector{Jet{T}}(undef, N+1)
    R = Vector{Jet{T}}(undef, N+1); E = Vector{Jet{T}}(undef, N+1); Sn = Vector{Jet{T}}(undef, N+1)
    ex = Vector{T}(undef, N+1)
    eX = exp(X); fact = one(T)
    for i in 0:N
        i > 0 && (fact = fact*i)
        ex[i+1] = eX/fact
    end
    F_ = basetype(T)
    Q[1] = S.q; P[1] = S.p
    inf(S.q[1]) > 0 || error("taylor: q not positive")
    R[1] = jinv(S.q[1:4]); Sn[1] = jneg(R[1]); E[1] = jexp(Sn[1])
    for n in 0:N-1
        if n >= 1
            acc = jzero(T)
            for i in 1:n
                acc = jadd(acc, jmul(Q[i+1], R[n-i+1]))
            end
            R[n+1] = jneg(jmul(R[1], acc)); Sn[n+1] = jneg(R[n+1])
            acc = jzero(T)
            for i in 1:n
                acc = jadd(acc, jscale(i, jmul(Sn[i+1], E[n-i+1])))
            end
            E[n+1] = jscale(inv(interval(F_, n)), acc)
        end
        F = jzero(T)
        for i in 0:n
            F = jadd(F, jscale(ex[i+1], E[n-i+1]))
        end
        Q[n+2] = jscale(inv(interval(F_, n+1)), P[n+1])
        P[n+2] = jscale(-inv(interval(F_, n+1)), F)
    end
    Q, P
end

function inflate(x::Interval{T}, fac, abs_) where T
    m = mid(x); r = radius(x)
    rr = fac*r + abs_*(1 + abs(m))
    interval(T, prevfloat(m - rr), nextfloat(m + rr))
end
inflate(J::Jet, fac, abs_) = map(x -> inflate(x, fac, abs_), J)
inflate(S::State, fac, abs_) = State(inflate(S.q, fac, abs_), inflate(S.p, fac, abs_))
jinterior(a::Jet, b::Jet) = all(isinterior(a[i], b[i]) for i in 1:4)
sinterior(A::State, B::State) = jinterior(A.q, B.q) && jinterior(A.p, B.p)
jfinite(a::Jet) = all(isbounded(a[i]) for i in 1:4)

# Picard a-priori tube: B with S0 + [0,h] F(X,B) inside int(B)
function tube(x, h, S0::State{T}, flat::Bool, fb) where T
    F_ = basetype(T)
    H = interval(F_, 0, sup(h)); X = interval(F_, x) + H
    fld(XX, SS) = flat ? field_flat(XX, SS, fb) : field(XX, SS)
    F0 = fld(interval(F_, x), S0)
    B = inflate(State(jadd(S0.q, jscale(H, F0.q)), jadd(S0.p, jscale(H, F0.p))), 2.0, 1e-13)
    for it in 1:30
        flat || inf(B.q[1]) > 0 || error("tube: q not positive")
        flat && (sup(B.q[1]) < 1/40 - 1e-12 || error("tube: flat bound domain"))
        F = fld(X, B)
        Z = State(jadd(S0.q, jscale(H, F.q)), jadd(S0.p, jscale(H, F.p)))
        if sinterior(Z, B)
            return X, B
        end
        B = inflate(shull(B, Z), 1.5, 1e-13)
    end
    error("tube: Picard inclusion failed")
end

# Taylor step of size h (exact float); returns end enclosure and enclosure over the step
function taylor_step(x, h, S0::State{T}, N::Int) where T
    F_ = basetype(T)
    X, B = tube(x, interval(F_, h), S0, false, nothing)
    Qc, Pc = taylor_coeffs(interval(F_, x), S0, N-1)    # orders 0..N-1 at the start
    Qr, Pr = taylor_coeffs(X, B, N)                     # order N on the tube (remainder)
    hh = interval(F_, h); HH = interval(F_, 0, h)
    function horner(C::Vector, Rn::Jet, t)
        acc = Rn
        for n in N:-1:1
            acc = jadd(jmul_s(acc, t), C[n])
        end
        acc
    end
    # polynomial sum_{n<N} C_n t^n + R_N t^N
    qend = horner(Qc, Qr[N+1], hh); pend = horner(Pc, Pr[N+1], hh)
    qwin = horner(Qc, Qr[N+1], HH); pwin = horner(Pc, Pr[N+1], HH)
    # remainder size (diagnostic for step control)
    rem = maximum(mag(Qr[N+1][i]*hh^N) for i in 1:4)
    State(qend, pend), intersect_state(State(qwin, pwin), B), rem
end
@inline jmul_s(a::Jet, t) = jscale(t, a)
function intersect_state(A::State, B::State)
    State(map(intersect_interval, A.q, B.q), map(intersect_interval, A.p, B.p))
end

# second-order step valid in the flat region q <= 1/40 (including q <= 0)
function flat_step(x, h, S0::State{T}, fb) where T
    F_ = basetype(T)
    X, B = tube(x, interval(F_, h), S0, true, fb)
    F = field_flat(X, B, fb)            # F.p encloses p_x on the tube
    hh = interval(F_, h); HH = interval(F_, 0, h)
    qend = jadd(jadd(S0.q, jscale(hh, S0.p)), jscale(hh^2/2, F.p))
    pend = jadd(S0.p, jscale(hh, F.p))
    qwin = jadd(jadd(S0.q, jscale(HH, S0.p)), jscale(HH^2/2, F.p))
    pwin = jadd(S0.p, jscale(HH, F.p))
    State(qend, pend), intersect_state(State(qwin, pwin), B), zero(T)
end

# Real-variable bootstrap at the origin.  In t = s^2/4 = e^x:  (t q_t)_t = -g(q),  p = t q_t,
# so p(t) = -int_0^t G dtau and Q(t) - A = int_0^t p(tau)/tau dtau (jetwise).  If |G_k| <= M_k on the
# jet ball |Q_k - A_k| <= rho and M_k t0 < rho, a continuity argument keeps the solution in the ball on
# [0, t0] and gives |Q_k - A_k| <= M_k t0, |P_k| <= M_k t0, P_0 in [-M_0 t0, 0].
function initial_state(abox::Interval{T}, x0::Int; rho = 1e-8) where T
    A = (abox, interval(T, 1), interval(T, 0), interval(T, 0))
    r = interval(T, -rho, rho)
    ball = jadd(A, (r, r, r, r))
    G = gjet(ball)
    t0 = interval(T, sup(exp(interval(T, x0))))      # an upper bound for e^{x0}
    Mt = map(m -> sup(interval(T, m)*t0), map(mag, G))  # rounded-up bounds for M_k t0
    all(Mt[i] < rho/2 for i in 1:4) || error("bootstrap failed")
    e = map(m -> interval(T, -m, m), Mt)
    q = jadd(A, e)
    p = (interval(T, -Mt[1], 0), e[2], e[3], e[4])
    State(q, p)
end

# v, z, c and their x-derivatives from a state
vzc(S::State) = (S.q[2], 2*S.q[3], 6*S.q[4])
jkl(S::State) = (S.p[2], 2*S.p[3], 6*S.p[4])

# exact time advance (the state computed by a step is at x+h exactly; refuse rounding)
function advance(x, h)
    xn = x + h
    (xn - x == h && xn - h == x) || error("time step not exactly representable")
    xn
end

# First-order mean-value recentring in the parameter, done with the jets themselves:
# d/da Q_m = (m+1) Q_{m+1}, so Q_m(a) in Q_m(amid) + (a - amid) (m+1) Q_{m+1}(box), for m = 2, 1, 0
# (in this order, each using the already refined next component).  B and C must enclose the box and
# the midpoint solution at the same time (or over the same time interval, pointwise in time).
function recenter(B::State{T}, C::State{T}, d) where T
    q = collect(B.q); p = collect(B.p)
    for i in 3:-1:1               # 1-based jet index i <-> order m = i-1, factor m+1 = i
        q[i] = intersect_interval(q[i], C.q[i] + d*i*q[i+1])
        p[i] = intersect_interval(p[i], C.p[i] + d*i*p[i+1])
    end
    any(isempty_interval, q) && error("recenter: empty intersection")
    any(isempty_interval, p) && error("recenter: empty intersection")
    State(Tuple(q), Tuple(p))
end

"""
    shoot_event(abox; refine, mv)

Integrate the box `abox` (and, if `mv`, its exact midpoint in lockstep for mean-value recentring)
from x0 until the first step end at which v < 0 for every member.  The last step end with v > 0
for every member is xL; the window [xL, xR] is covered by small steps whose step enclosures are
hulled into `W`.  Rigorous content: every member has v(xL) > 0 > v(xR), and W encloses all states
on [xL, xR].
"""
function shoot_event(abox::Interval{T}; N = 16, x0 = -30, hw = 1/512, xmax = 6.0, tol = 1e-15,
                     qflat = 0.024, refine = 0, mv = true) where T
    fb = flat_bounds(T)
    amid = interval(T, mid(abox))
    mv = mv && radius(abox) > 0
    d = abox - amid
    S = initial_state(abox, x0)
    Sc = mv ? initial_state(amid, x0) : S
    x = T(x0)          # time kept in the base number type; every advance is checked exact
    nsteps = 0; nflat = 0
    isflat(S) = sup(S.q[1]) < qflat && sup(S.p[1]) < 0
    function rawstep(x, h, S)
        if isflat(S)
            return flat_step(x, h, S, fb)
        else
            Sn, W, rem = taylor_step(x, h, S, N)
            (rem < tol || h <= 2.0^-14) || error("remainder above tolerance")
            return Sn, W, rem
        end
    end
    function onestep(x, S, Sc, hmax)
        h = hmax
        for k in 1:25
            try
                Sn, W, _ = rawstep(x, h, S)
                if mv
                    Scn, Wc, _ = rawstep(x, h, Sc)
                    Sn = recenter(Sn, Scn, d); W = recenter(W, Wc, d)
                else
                    Scn = Sn
                end
                all(jfinite(J) for J in (Sn.q, Sn.p, W.q, W.p)) || error("non-finite enclosure")
                return h, Sn, Scn, W, isflat(S)
            catch err
                (err isa ErrorException) || rethrow()
            end
            h /= 2
        end
        error("step failed at x = $x")
    end
    # phase 1: up to the window
    while true
        hmax = x < -4 ? 0.5 : 1/32
        h, Sn, Scn, W, fl = onestep(x, S, Sc, hmax)
        if inf(Sn.q[2]) > 0
            x = advance(x, h); S = Sn; Sc = Scn; nsteps += 1; nflat += fl
            x < xmax || error("no event before xmax")
            continue
        end
        break
    end
    # phase 2: window, optionally refined.  Each pass restarts from the last step end with v > 0
    # (for all members) and uses smaller steps; the final window is [xL, xR].
    hcur = hw
    local xL, SL, xR, SR, Wh
    for pass in 0:refine
        xs = x; Ss = S; Scs = Sc
        lastpos = (x, S, Sc); Wsince = nothing
        while true
            h, Sn, Scn, W, fl = onestep(xs, Ss, Scs, hcur)
            xs = advance(xs, h); Ss = Sn; Scs = Scn; nsteps += 1; nflat += fl
            Wsince = Wsince === nothing ? W : shull(Wsince, W)
            if inf(Ss.q[2]) > 0          # still positive for every member: move the left edge
                lastpos = (xs, Ss, Scs); Wsince = nothing
            end
            sup(Ss.q[2]) < 0 && break
            xs - lastpos[1] < 0.5 || error("event window too long")
        end
        xL, SL, _ = lastpos; xR = xs; SR = Ss; Wh = Wsince
        x, S, Sc = lastpos
        hcur = hcur/16
    end
    (xL = xL, xR = xR, SL = SL, SR = SR, W = Wh, steps = nsteps, flat = nflat)
end

# Quantities on the window
function event_quantities(ev)
    W = ev.W
    q = W.q[1]; p = W.p[1]
    v, z, c = vzc(W); j, k, l = jkl(W)
    D = c - k*z/j
    (q = q, p = p, v = v, z = z, c = c, j = j, k = k, D = D,
     vL = ev.SL.q[2], vR = ev.SR.q[2])
end
