# bk_po_branch.jl -- periodic-orbit branch emanating from the Hopf point (BifurcationKit, orthogonal collocation)
# Usage: julia bk_po_branch.jl      Runtime ~3 min.  Evidence only (no PASS/FAIL): prints K - Kc and the n_M
# half-amplitude along the branch of periodic orbits for M = 15 and M = 16. Supercritical <=> the branch leaves
# towards K > Kc with amp^2 ~ (K - Kc); subcritical <=> towards K < Kc.  The first printed points are the
# branch-switching predictor; read the monotone part.
using BifurcationKit, LinearAlgebra
include_string(Main, join(readlines(joinpath(@__DIR__, "bk_check.jl"))[11:25], "\n"))   # reuse bd_reduced from bk_check.jl
for (M, Kc) in [(15, 212.895216331), (16, 20.77381284)]
    K0 = 0.999 * Kc
    y0 = fill(1 + K0, M - 1)
    prob = BifurcationProblem(bd_reduced, y0, (K = K0,), (@optic _.K))
    opts = ContinuationPar(p_min = 0.99*Kc, p_max = 1.01 * Kc, ds = 2e-5 * Kc, dsmin = 1e-9 * Kc, dsmax = 1e-4 * Kc,
                           detect_bifurcation = 3, nev = M - 1, n_inversion = 8, max_bisection_steps = 40,
                           tol_stability = 1e-10, newton_options = NewtonPar(tol = 1e-10 * (1 + Kc)^2))
    br = continuation(prob, PALC(), opts; verbosity = 0)
    idx = findfirst(sp -> sp.type == :hopf, br.specialpoint)
    optspo = ContinuationPar(opts; p_min = 0.98*Kc, p_max = 1.02*Kc, ds = 1e-6*Kc, dsmax = 1e-4*Kc, max_steps = 40,
                             detect_bifurcation = 0, newton_options = NewtonPar(tol = 1e-8 * (1 + Kc)^2, max_iterations = 25))
    brpo = continuation(br, idx, optspo, Collocation(30, 3); δp = 1e-6*Kc, verbosity = 0,
                        record_from_solution = (x, p; k...) -> begin
                            xx = get_periodic_orbit(p.prob, x, p.p)
                            (max = maximum(xx[end, :]), min = minimum(xx[end, :]), period = getperiod(p.prob, x, p.p))
                        end)
    println("M=$M Kc=$(br.specialpoint[idx].param)")
    for k in 1:5:length(brpo)
        println("   K-Kc = $(brpo.param[k] - br.specialpoint[idx].param)  amp(n_M) = $((brpo.branch[k].max - brpo.branch[k].min)/2)  stable=$(brpo.stable[k])")
    end
end
