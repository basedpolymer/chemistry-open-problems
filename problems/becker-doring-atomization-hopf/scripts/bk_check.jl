# bk_check.jl -- independent check of the Hopf normal form with BifurcationKit v0.8 (Julia 1.13).
# Usage: julia bk_check.jl M1 Kc1 [M2 Kc2 ...]      (Kc = approximate critical K, e.g. from the sweep files)
# Runtime: ~1-2 min (compilation) + a few seconds per point.
# The mass-conserving system is reduced to y = (n_2..n_M) by n_1 = m(K) - sum_{l>=2} l n_l, with
# m(K) = (1+K) M (M+1)/2 the mass of the constant equilibrium, so y = (1+K) is an equilibrium for every K.
# Derivatives (Jacobian, second and third differentials) come from ForwardDiff inside BifurcationKit; the
# Hopf point is located by BifurcationKit's own continuation/bisection and its normal form
#     dz/dt = (i omega + a dK) z + b z |z|^2
# is computed by get_normal_form. Nothing is shared with the Python code. Prints RESULT lines.
using BifurcationKit, LinearAlgebra

function bd_reduced!(dy, y, p)
    K = p.K; M = length(y) + 1; N = M - 1
    A = 1 + K
    m = A * M * (M + 1) / 2
    n1 = m - sum((l + 1) * y[l] for l in 1:N)        # y[l] = n_{l+1}
    n(l) = l == 1 ? n1 : y[l-1]
    J(l) = n(l) * n1 - n(l + 1)                      # flux J_l, l = 1..N
    for l in 2:N
        dy[l-1] = J(l - 1) - J(l)
    end
    dy[M-1] = J(N) - K * n(M)
    return dy
end
bd_reduced(y, p) = bd_reduced!(similar(y, promote_type(eltype(y), typeof(p.K))), y, p)

function check(M, Kc)
    K0 = 0.999 * Kc
    y0 = fill(1 + K0, M - 1)
    prob = BifurcationProblem(bd_reduced, y0, (K = K0,), (@optic _.K))
    opts = ContinuationPar(p_min = K0, p_max = 1.001 * Kc, ds = 2e-5 * Kc, dsmin = 1e-9 * Kc, dsmax = 1e-4 * Kc,
                           detect_bifurcation = 3, nev = M - 1, n_inversion = 8, max_bisection_steps = 40,
                           tol_stability = 1e-10, newton_options = NewtonPar(tol = 1e-10 * (1 + Kc)^2))
    br = continuation(prob, PALC(), opts; verbosity = 0)
    idx = findfirst(sp -> sp.type == :hopf, br.specialpoint)
    if idx === nothing
        println("RESULT M=$M NO_HOPF_DETECTED")
        return
    end
    hp = get_normal_form(br, idx; verbose = false)
    println("RESULT M=$M Kc=$(br.specialpoint[idx].param) omega=$(hp.ω) real(b)=$(real(hp.nf.b)) type=$(hp.type)")
end

for i in 1:2:length(ARGS)
    check(parse(Int, ARGS[i]), parse(Float64, ARGS[i+1]))
end
