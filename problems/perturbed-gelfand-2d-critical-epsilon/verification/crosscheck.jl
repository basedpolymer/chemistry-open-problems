# Verifier's cross-check between the two independent computations (about 10 s, one thread):
# the folds found by NON-RIGOROUS direct shooting on the original equation (direct_original.log, values of
# a = eps + eps^2 alpha printed to 8 decimals) must satisfy E(a) = eps for the reduced stopping-level
# function, computed here with the validated integrator (Float64 intervals), with z < 0 at the first
# fold and z > 0 at the second.  This checks the fold equivalence lambda_alpha = 0 <=> eps = E(a).
include(joinpath(@__DIR__, "jets.jl"))
for (a, eps) in [(0.05390218,0.05),(1.39642155,0.05),(0.19885668,0.15),(1.11213391,0.15),(0.43624794,0.23),(0.76177042,0.23),(0.57534657,0.242),(0.60574651,0.242)]
    E = event_quantities(shoot_event(interval(a); refine=4, xmax=80.0))
    println("a = ", a, "  E(a) in ", E.q, "   (direct-shooting eps = ", eps, ")   z = ", E.z)
end
