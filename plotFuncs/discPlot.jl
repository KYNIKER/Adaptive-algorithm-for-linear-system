# Based on the paper JuliaReach: a Toolbox for Set-Based Reachability
using Plots, LazySets, LinearAlgebra, BenchmarkTools, FastExpm, Profile, PProf


include("../helperfunctions.jl")
#include("../models.jl")
include("../models/heat/heat_load.jl")
include("../models/motor/motor_load.jl")
include("../models/motor/motor_load.jl")
include("../models/building/building_load.jl")
include("../models/PDE/pde_load.jl")
include("../models/ISS/iss_load.jl")
include("../models/beam/beam_load.jl")
include("../models/FOM/fom_load.jl")
include("../models/MNA1/mna1_load.jl")
include("../models/MNA5/mna5_load.jl")
#include("CegarFunctions.jl")
#include("../CegarInhomogenous.jl")
include("plotHelper.jl")

function getDiscOne(A, B, initialTimeStep, X0::Zonotope{N,Vector{N},Matrix{N}}, U::Zonotope, iterations) where {N}
    stepsBeforeReduce = 4
    maxOrder = 50
    reduceToOrder = 1
    finalReduce = maxOrder # We reduce the final P to this

    ANorm = norm(A, Inf)
    XNorm = norm(X0, Inf)::Float64
    XG = copy(genmat(X0))
    XDim, p = size(XG)
    XC = copy(X0.center)
    m = initialTimeStep / 2^iterations
    #m = initialTimeStep / 2^(ceil(Integer, log2(initialTimeStep)) + ceil(Integer, -log2(10.0^(-Digits))) - 1)   #Calculate the smallest number larger than 10^-Digits obtained by repeatedly dividing initialTimeStep by 2.
    R = Zonotope[]
    changedTimeStep = true
    #println("m: ", m)
    #discritezationDict = Dict()#Dict{Float64, Tuple{Zonotope{N,Vector{N},Matrix{N}}, Matrix{Float64}}}()
    phiDict = Dict{Float64,Matrix{Float64}}()
    discritezationDict = Dict{Float64,Zonotope{N,Vector{N},Matrix{N}}}()
    inputDiscritezationDict = Dict()#Dict{Float64, Zonotope{N,Vector{N},Matrix{N}}}()

    k = size(U.generators, 2)

    U = concretize(B * U)

    let ϕ::Matrix{Float64} = fastExpm(A .* m; threshold=eps(Float64), nonzero_tol=eps(Float64))
        tempM = similar(ϕ)
        d = m
        dia::Matrix{Float64} = diagm(ones(XDim))

        if !(ρ(U.center, U) < norm(U.center)) #Origin is *not* in input
            println("Origin is not in input...")
            #println(d)
            û = copy(U.center)
            invA = inv(Matrix(A))
            Ut = Zonotope(U.center - û, genmat(U))
            dU = overapproximate(d * Ut, Zonotope)
            P = minkowski_sum(dU, E_ψ(Ut, d, A))
            #println(typeof(P))
            #P = PCA_reduce(P)

            P̂ = invA * (ϕ - dia) * û
            lt = concretize(minkowski_sum(convert(Zonotope, ϕ * X0), dU))
            rt = concretize(minkowski_sum(E_ψ(Ut, d, A), E⁺(X0, d, A)))
            PZ = Zonotope(P̂, zeros(Float64, size(U.center, 1), 1))
            f = concretize(minkowski_sum(lt, rt))
            disc = overapproximate(CH(X0, minkowski_sum(f, PZ)), Zonotope)

            while d < initialTimeStep
                #inputDiscritezationDict[d] = P
                P = minkowski_sum(P, linear_map(ϕ, P))
                if LazySets.order(P) > maxOrder
                    #println("REDUCE!")
                    P = reduce_order(P, reduceToOrder)
                end

                #=P̂ = invA * (ϕ - dia) * û
                lt = concretize(minkowski_sum(convert(Zonotope, ϕ * X0), box_approximation_symmetric(d * Ut)))
                rt = concretize(minkowski_sum(E_ψ(Ut, d, A), E⁺(X0, d, A)))
                PZ = Zonotope(P̂, zeros(Float64, size(U.center, 1), 1))
                f = concretize(minkowski_sum(lt, rt))
                disc = overapproximate(CH(X0, minkowski_sum(f, PZ)), Zonotope)=#
                discritezationDict[d] = copy(disc)
                phiDict[d] = copy(ϕ)
                disc = overapproximate(CH(disc, linear_map(ϕ, disc)), Zonotope)

                mul!(tempM, ϕ, ϕ)
                copy!(ϕ, tempM)
                d = d * 2
            end
            #=P̂ = invA * (ϕ - dia) * û

            lt = concretize(minkowski_sum(convert(Zonotope, ϕ * X0), box_approximation_symmetric(d * Ut)))
            rt = concretize(minkowski_sum(E_ψ(Ut, d, A), E⁺(X0, d, A)))
            PZ = Zonotope(P̂, zeros(Float64, size(U.center, 1), 1))
            f = concretize(minkowski_sum(lt, rt))

            disc = overapproximate(CH(X0, minkowski_sum(f, PZ)), Zonotope)=#
            discritezationDict[d] = copy(disc)
            phiDict[d] = copy(ϕ)
            #dU = box_approximation_symmetric(initialTimeStep * Ut)
            #P = minkowski_sum(dU, E_ψ(Ut, initialTimeStep, A))
            if LazySets.order(P) > finalReduce
                #println("REDUCE!")
                P = reduce_order(P, finalReduce)
            end
            inputDiscritezationDict[initialTimeStep] = P
        else
            dU = overapproximate(d * U, Zonotope)
            P = minkowski_sum(dU, E_ψ(U, d, A))
            lt = concretize(minkowski_sum(convert(Zonotope, ϕ * X0), dU))
            rt = concretize(minkowski_sum(E_ψ(U, d, A), E⁺(X0, d, A)))
            f = concretize(minkowski_sum(lt, rt))
            disc = overapproximate(CH(X0, f), Zonotope)
            while d < initialTimeStep
                #inputDiscritezationDict[d] = P
                P = minkowski_sum(P, linear_map(ϕ, P))
                if LazySets.order(P) > maxOrder
                    #println("REDUCE!")
                    P = reduce_order(P, reduceToOrder)
                end


                phiDict[d] = copy(ϕ)
                discritezationDict[d] = copy(disc)
                disc = overapproximate(CH(disc, linear_map(ϕ, disc)), Zonotope)
                mul!(tempM, ϕ, ϕ)
                copy!(ϕ, tempM)
                d = d * 2
            end
            # lt = concretize(minkowski_sum(convert(Zonotope, ϕ * X0), box_approximation_symmetric(d * U)))
            # rt = concretize(minkowski_sum(E_ψ(U, d, A), E⁺(X0, d, A)))
            # disc = overapproximate( CH(X0, concretize(minkowski_sum(lt, rt))), Zonotope)
            discritezationDict[d] = copy(disc)
            phiDict[d] = copy(ϕ)
            if LazySets.order(P) > finalReduce
                #println("REDUCE!")
                P = reduce_order(P, finalReduce)
            end
            inputDiscritezationDict[initialTimeStep] = P
        end
    end

    println("Dicts done!")
    finalHomogeneous = discritezationDict[initialTimeStep]
    finalInput = inputDiscritezationDict[initialTimeStep]

    return MinkowskiSum(finalHomogeneous, finalInput)
end


initialTimeStep = 0.002 * 2^6
levels = 5
dimToPlot1 = 1
dimToPlot2 = 2

A, B, ballβ, P₁, T, constraint, dimToPlot = load_building()


finalZonotope = getDiscOne(A, B, initialTimeStep, P₁, ballβ, 0) # Start with no extra disc

#corners2 = Vector(undef, size(boxes2, 1))

H = box_approximation(finalZonotope)
H_proj = LazySets.project(H, [dimToPlot1, dimToPlot2]) # Only get the dimension we care about
tope = vertices_list(H_proj)

# Switch 3rd and 4th index for Plotting
tope[3], tope[4] = tope[4], tope[3]

cordsX = getindex.(tope, 1)
cordsY = getindex.(tope, 2)

# println(cordsX)
# println(cordsY)

maxValx = maximum(cordsX)
minValx = minimum(cordsX)
maxValy = maximum(cordsY)
minValy = minimum(cordsY)

p = plot(dpi=300, thickness_scaling=1, xlabel="Value 1", ylabel="Value 2")


res = Shape(cordsX, cordsY)



palette = Plots.palette(:tab10)

color = palette[1]

plot!(p, res, color=:pink, alpha=0.2, label="δ = $(initialTimeStep)")

for level in 1:5
    finalZonotope = getDiscOne(A, B, initialTimeStep, P₁, ballβ, level)
    H = box_approximation(finalZonotope)
    H_proj = LazySets.project(H, [dimToPlot1, dimToPlot2]) # Only get the dimension we care about
    tope = vertices_list(H_proj)

    tope[3], tope[4] = tope[4], tope[3]
    cordsX = getindex.(tope, 1)
    cordsY = getindex.(tope, 2)

    res2 = Shape(cordsX, cordsY)
    timeval = initialTimeStep / (2^level)

    color = palette[level]
    plot!(p, res2, color=color, alpha=0.1 * level, label="δ = $(timeval)")
    println("Tope: ", tope)
end

savefig(p, "discExample2.png")
plot(p)
