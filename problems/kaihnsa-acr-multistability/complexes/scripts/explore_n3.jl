# Exploration (floating point, not a certificate): search 6-complex bimolecular networks on 3 species
# (full-dimensional, one linkage class) for ACR in X together with two exponentially stable positive steady states.
# Inverse approach: choose steady states on the plane x = c, K = span of their monomial vectors,
# V = A*W with W a basis of K^perp, A such that each column A*w_i lies in the cone of reaction vectors of reactant i.
# Usage: julia explore_n3.jl [samples_per_case] [seed]   (runtime: about 80 s for 1000 samples per case)
# Not a certificate: prints candidate hits only (no PASS/FAIL). Seed 3 with 1000 samples found the first hit
# (complexes Z, 2Z, Y+Z, 2Y, X, X+Z; weight (2,1,1)); see report.md, Section 9.
using LinearAlgebra, Random

function combinations(v, k)
    out = Vector{Vector{Int}}()
    n = length(v)
    function rec(start, cur)
        if length(cur) == k
            push!(out, copy(cur)); return
        end
        for i in start:n
            push!(cur, v[i]); rec(i+1, cur); pop!(cur)
        end
    end
    rec(1, Int[])
    out
end

const SP = ["X","Y","Z"]
function allcomplexes()
    cs = Vector{Vector{Int}}()
    for a in 0:2, b in 0:2, c in 0:2
        if a+b+c <= 2
            push!(cs, [a,b,c])
        end
    end
    cs
end
cname(c) = begin
    parts = String[]
    for k in 1:3
        c[k] == 1 && push!(parts, SP[k])
        c[k] == 2 && push!(parts, "2"*SP[k])
    end
    isempty(parts) ? "0" : join(parts, "+")
end

# valid facet normals of cone(G) (columns of G) in R^3
function facetnormals(G)
    H = Vector{Vector{Float64}}()
    m = size(G,2)
    for a in 1:m, b in a+1:m
        h = cross(G[:,a], G[:,b])
        norm(h) < 1e-12 && continue
        s = [dot(h, G[:,k]) for k in 1:m]
        if all(s .>= -1e-12)
            push!(H, h/norm(h))
        elseif all(s .<= 1e-12)
            push!(H, -h/norm(h))
        end
    end
    H
end
incone(H, v; tol=1e-10) = all(dot(h, v) >= -tol*norm(v) for h in H)

monos(Cm, x) = [prod(x .^ Cm[:,i]) for i in 1:size(Cm,2)]

function hurwitz3(J)
    a1 = -tr(J)
    a2 = J[1,1]*J[2,2]-J[1,2]*J[2,1] + J[1,1]*J[3,3]-J[1,3]*J[3,1] + J[2,2]*J[3,3]-J[2,3]*J[3,2]
    a3 = -det(J)
    (a1 > 0 && a3 > 0 && a1*a2 > a3, (a1, a2, a3))
end

function sample_cone(G)
    m = size(G,2)
    w = zeros(m)
    # random sparse positive weights
    k = rand(1:m)
    idx = randperm(m)[1:k]
    w[idx] .= exp.(2 .* randn(k))
    G*w
end

function run(nsamp, seed)
    Random.seed!(seed)
    cs = allcomplexes()
    hits = 0
    stats = zeros(Int, 5)
    for S6 in combinations(1:length(cs), 6)
        C = hcat(cs[S6]...)
        rank(C[:,2:end] .- C[:,1]) == 3 || continue
        for nonr in 0:6
            R = nonr == 0 ? collect(1:6) : setdiff(1:6, [nonr])
            Cr = C[:,R]
            rank(Cr[:,2:end] .- Cr[:,1]) == 3 || continue
            r = length(R)
            # SF pairs (differ by e_X) among reactants
            sfpairs = [(i,j) for i in 1:r for j in 1:r if Cr[:,j] - Cr[:,i] == [1,0,0]]
            if r == 5 && isempty(sfpairs)
                continue
            end
            Gs = [Float64.(hcat([C[:,j] - C[:,R[i]] for j in 1:6 if j != R[i]]...)) for i in 1:r]
            Hs = [facetnormals(G) for G in Gs]
            npts = r == 6 ? 3 : 2
            for s in 1:nsamp
                c = exp(2*(rand()-0.5)*2)
                P = [[c, exp(4*(rand()-0.5)), exp(4*(rand()-0.5))] for _ in 1:npts]
                K = hcat([monos(Cr, p) for p in P]...)
                rank(K) == npts || continue
                W = nullspace(K')'   # 3 x r
                size(W,1) == 3 || continue
                # choose basis B of 3 reactant columns
                B = randperm(r)[1:3]
                abs(det(W[:,B])) < 1e-8 && continue
                Vb = hcat([sample_cone(Gs[i]) for i in B]...)
                A = Vb / W[:,B]
                V = A*W
                ok = all(incone(Hs[i], V[:,i]) for i in 1:r)
                stats[1] += 1
                ok || continue
                stats[2] += 1
                # stability at the chosen points (+ third root search for r = 5 later)
                nst = 0
                for p in P
                    J = V*Diagonal(monos(Cr,p))*Cr'*Diagonal(1 ./ p)
                    ev = eigvals(J)
                    nJ = opnorm(J)
                    if maximum(real.(ev)) < -1e-4*nJ && minimum(abs.(ev)) > 1e-4*nJ
                        nst += 1
                    end
                end
                nst >= 1 && (stats[3] += 1)
                if nst >= 2
                    hits += 1
                    println("HIT r=$r complexes=", join(cname.(eachcol(C)), ", "), " reactants=", join(cname.(eachcol(Cr)), ", "))
                    println("  c=$c points=$P")
                    println("  V=$V")
                    flush(stdout)
                end
            end
        end
    end
    println("total hits: $hits  stats(tried,conefeasible,>=1stable) = $stats")
end

nsamp = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 2000
seed = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 1
run(nsamp, seed)
