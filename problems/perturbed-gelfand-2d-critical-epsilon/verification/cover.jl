# Verifier's independent global cover of the stopping-level curve E(a) on [1/2, 8/5].
# Run (from this folder):  julia -t 6 cover.jl
# Environment: Julia 1.13, IntervalArithmetic.jl 1.0.12; Float64 intervals; about 4 minutes on 6 threads.
#
# For every dyadic box I in [1/2, 1.6015625] (initial width 1/512, adaptive bisection) we certify, for
# the extended problem (g = 0 for q <= 0), a window [xL, xR] with v(xL) > 0 > v(xR) and v_x < 0 on the
# window for every a in I, and one of
#   N:  q > 0 and z < 0 on the window          (E increasing, event physical)
#   C:  q > 0 and D = c - k z / j > 0           (event value of z strictly increasing)
#   P:  z > 0 on the window                     (E decreasing wherever the event is physical)
# and check the pattern N...N C...C P...P with at least one C box.
include(joinpath(@__DIR__, "jets.jl"))

function classify(lo::Float64, hi::Float64)
    ev = shoot_event(interval(lo, hi))
    E = event_quantities(ev)
    ok_win = sup(E.vR) < 0 && inf(E.vL) > 0 && sup(E.j) < 0 && sup(E.p) < 0
    ok_win || return (:none, E)
    if inf(E.q) > 0 && sup(E.z) < 0
        return (:N, E)
    elseif inf(E.z) > 0
        return (:P, E)
    elseif inf(E.q) > 0 && inf(E.D) > 0
        return (:C, E)
    end
    (:none, E)
end

function certify_box(lo, hi, depth)
    kind = :none; E = nothing
    try
        kind, E = classify(lo, hi)
    catch err
        (err isa ErrorException) || rethrow()
    end
    if kind != :none
        return [(lo, hi, kind, E)]
    end
    depth >= 10 && error("FAIL box [$lo, $hi]")
    m = (lo + hi)/2
    vcat(certify_box(lo, m, depth+1), certify_box(m, hi, depth+1))
end

function main()
    t0 = time()
    w = 1/512
    starts = [0.5 + k*w for k in 0:563]           # 564 boxes, last ends at 1.6015625
    @assert starts[end] + w == 1.6015625
    results = Vector{Any}(undef, length(starts))
    Threads.@threads for i in eachindex(starts)
        results[i] = certify_box(starts[i], starts[i] + w, 0)
    end
    leaves = reduce(vcat, results)
    sort!(leaves, by = r -> r[1])
    check(name, ok) = (println(ok ? "PASS " : "FAIL ", name); ok || exit(1))
    check("cover starts at 1/2", leaves[1][1] == 0.5)
    check("cover reaches past 8/5", leaves[end][2] >= 1.6)
    check("no gaps", all(leaves[i][2] == leaves[i+1][1] for i in 1:length(leaves)-1))
    kinds = [r[3] for r in leaves]
    iC = findall(==(:C), kinds)
    check("at least one C box", !isempty(iC))
    c1, c2 = first(iC), last(iC)
    check("pattern N* C+ P*", all(kinds[1:c1-1] .== :N) && all(kinds[c1:c2] .== :C) && all(kinds[c2+1:end] .== :P))
    check("N block nonempty and P block nonempty", c1 > 1 && c2 < length(kinds))
    println("leaves: ", length(leaves), "  N: ", count(==(:N), kinds), "  C: ", length(iC), "  P: ", count(==(:P), kinds))
    println("C region: [", leaves[c1][1], ", ", leaves[c2][2], "]")
    zN = maximum(sup(r[4].z) for r in leaves if r[3] == :N)
    zP = minimum(inf(r[4].z) for r in leaves if r[3] == :P)
    Dmin = minimum(inf(r[4].D) for r in leaves if r[3] == :C)
    qmin = minimum(inf(r[4].q) for r in leaves if r[3] != :P)
    println("max z on N boxes: ", zN, "   min z on P boxes: ", zP, "   min D on C boxes: ", Dmin,
            "   min q on N/C boxes: ", qmin)
    # minimum width of leaves
    println("min leaf width: ", minimum(r[2]-r[1] for r in leaves))
    # point values (8/5 is enclosed in a one-ulp interval; the claims hold for every member)
    for (a, name) in [(interval(Float64, 1//2), "1/2"), (interval(Float64, 3//2), "3/2"), (interval(Float64, 8//5), "8/5")]
        ev = shoot_event(a; refine = 6)
        E = event_quantities(ev)
        println("E($name) in ", E.q, "   z = ", E.z, "   window [", ev.xL, ", ", ev.xR, "]")
        check("window at a = $name", inf(E.vL) > 0 && sup(E.vR) < 0 && sup(E.j) < 0)
        name == "1/2" && check("E(1/2) > 0.238", inf(E.q) > 238 // 1000)
        name == "3/2" && check("E(3/2) > 0.0078 (so a_m > 3/2)", inf(E.q) > 78 // 10000)
        name == "8/5" && check("event level at 8/5 negative (so a_m < 8/5)", sup(E.q) < 0)
    end
    open(joinpath(@__DIR__, "cover-leaves.txt"), "w") do io
        for r in leaves
            E = r[4]
            println(io, r[1], " ", r[2], " ", r[3], " q=", E.q, " z=", E.z, " D=", E.D, " j=", E.j)
        end
    end
    println("Runtime seconds: ", round(time() - t0, digits = 1))
end

main()
