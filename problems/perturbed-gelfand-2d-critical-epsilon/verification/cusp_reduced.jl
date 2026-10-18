# Verifier's independent enclosure of the critical value eps0 = max E(a) and of the cusp,
# computed entirely in the parameter-free reduced family (no use of the solver's local certificate).
# Run (from this folder):  julia cusp_reduced.jl        (one thread, about 1 minute; Float64 and 160-bit BigFloat intervals)
#
# 1. Float64 interval bisection locates the sign change of zeta(a) = z(T(a), a).
# 2. 160-bit interval evaluations at the bracket ends a- < a+ re-prove zeta(a-) < 0 < zeta(a+); the cover shows
#    zeta' = D > 0 there, so the unique maximiser a* of E lies in [a-, a+].
# 3. A 160-bit box integration of [a-, a+] encloses eps0 = E(a*), the event time T(a*), and the
#    original-variable cusp (alpha_c, eps_c, lambda_c); one interval Newton step (zeta(a_mid) in 160-bit
#    arithmetic, zeta' = D enclosed on the box) narrows the enclosure of a*:
#       alpha_c = (a* - eps0)/eps0^2,   log(lambda_c/4) = T(a*) - 2 log eps0 - 1/eps0.
include(joinpath(@__DIR__, "jets.jl"))

check(name, ok) = (println(ok ? "PASS " : "FAIL ", name); flush(stdout); ok || exit(1))
full(x) = string("[", round(inf(x), RoundDown, sigdigits = 22), ", ", round(sup(x), RoundUp, sigdigits = 22), "]")

function zeta(a; T = Float64, refine = 6, N = 16, tol = 1e-15, x0 = -30)
    ev = shoot_event(interval(T, a); refine = refine, N = N, tol = tol, x0 = x0)
    event_quantities(ev), ev
end

function main()
    t0 = time()
    # --- 1. bisection with rigorous Float64 point evaluations
    lo, hi = 0.5859375, 0.59375
    Elo, _ = zeta(lo); Ehi, _ = zeta(hi)
    check("zeta(0.5859375) < 0", sup(Elo.z) < 0)
    check("zeta(0.59375) > 0", inf(Ehi.z) > 0)
    for it in 1:60
        m = (lo + hi)/2
        (m == lo || m == hi) && break
        E, _ = zeta(m)
        if sup(E.z) < 0
            lo = m
        elseif inf(E.z) > 0
            hi = m
        else
            break
        end
    end
    println("Float64 bracket: [", lo, ", ", hi, "]  width ", hi - lo); flush(stdout)
    # --- 2. 160-bit endpoints slightly inside/outside
    setprecision(BigFloat, 160)
    # lo and hi already carry rigorous signs; re-check them in 160-bit arithmetic
    aminus = BigFloat(lo); aplus = BigFloat(hi)
    t1 = time()
    Em, _ = zeta(aminus; T = BigFloat, refine = 10, N = 20, tol = 1e-22, x0 = -60)
    Ep, _ = zeta(aplus;  T = BigFloat, refine = 10, N = 20, tol = 1e-22, x0 = -60)
    println("zeta(a-) = ", Em.z); println("zeta(a+) = ", Ep.z)
    println("160-bit point evaluations: ", round(time() - t1, digits = 1), " s"); flush(stdout)
    check("zeta(a-) < 0 (160-bit)", sup(Em.z) < 0)
    check("zeta(a+) > 0 (160-bit)", inf(Ep.z) > 0)
    # --- 3. box [a-, a+]
    box = interval(BigFloat, aminus, aplus)
    ev = shoot_event(box; refine = 10, N = 20, tol = 1e-22, x0 = -60)
    E = event_quantities(ev)
    check("box window: v(xL) > 0 > v(xR), v_x < 0, q > 0, D > 0",
          inf(E.vL) > 0 && sup(E.vR) < 0 && sup(E.j) < 0 && inf(E.q) > 0 && inf(E.D) > 0)
    eps0 = E.q
    Tx = interval(BigFloat, ev.xL, ev.xR)
    # interval Newton step for a*: 0 = zeta(amid) + (a* - amid) D(xi), D(xi) in E.D (box enclosure)
    amid = (aminus + aplus)/2
    Emid, _ = zeta(amid; T = BigFloat, refine = 10, N = 20, tol = 1e-22, x0 = -60)
    println("zeta(amid) = ", Emid.z); flush(stdout)
    astar = intersect_interval(interval(BigFloat, amid) - Emid.z/E.D, box)
    check("Newton enclosure of a* nonempty and inside the bracket", !isempty_interval(astar))
    println("a*    in ", full(astar), "   width ", diam(astar))
    box = astar
    println("eps0  in ", full(eps0), "   width ", diam(eps0))
    println("T(a*) in ", full(Tx))
    println("D(a*) in ", full(E.D))
    e0 = interval(BigFloat, 24210616559524695 // 10^17)
    check("|eps0 - 0.24210616559524695| < 1e-12", sup(abs(eps0 - e0)) < 1//10^12)
    alpha = (box - eps0)/eps0^2
    a0 = interval(BigFloat, 5943243406485813 // 10^15)
    println("alpha_c in ", full(alpha))
    check("|alpha_c - 5.943243406485813| < 1e-10", sup(abs(alpha - a0)) < 1//10^10)
    xorig = Tx - 2*log(eps0) - inv(eps0)
    x0 = interval(BigFloat, -28558378244629573 // 10^17)
    println("log(lambda_c/4) in ", full(xorig))
    check("|log(lambda_c/4) - x0| < 1e-12", sup(abs(xorig - x0)) < 1//10^12)
    lam = 4*exp(xorig)
    println("lambda_c in ", full(lam))
    check("3.006301478866 < lambda_c < 3.006301478873",
          inf(lam) > 3006301478866 // 10^12 && sup(lam) < 3006301478873 // 10^12)
    a_c = box
    check("a* inside (0.5904717279294, 0.5904717279490)",
          inf(a_c) > 5904717279294 // 10^13 && sup(a_c) < 5904717279490 // 10^13)
    println("Runtime seconds: ", round(time() - t0, digits = 1))
end

main()
