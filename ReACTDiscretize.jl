export ReACTDiscretizeBase
using LinearAlgebra, LazySets, ReachabilityAnalysis
include("helperfunctions.jl")
isinvertible(x) = applicable(inv, x) && isone(inv(Matrix(x)) * x)

function ReACTDiscretize(A, B, X0::Zonotope{N,Vector{N},Matrix{N}}, U::Zonotope, δ⁻, δ⁺, alg::ReachabilityAnalysis.Exponentiation.AbstractExpAlg=ReachabilityAnalysis.Exponentiation.BaseExp, maxOrder::Int=10, reduceOrder::Int=10) where {N}
    phiDict = Dict{Float64,Matrix{Float64}}()
    discritezationDict = Dict{Float64,Zonotope{N,Vector{N},Matrix{N}}}()
    inputDiscritezationDict = Dict{Float64,Zonotope{N,Vector{N},Matrix{N}}}()

    U = linear_map(B, U) #overapproximate(concretize(B * U), Zonotope)

    let ϕ::Matrix{Float64} = ReachabilityAnalysis.Exponentiation._exp(A, δ⁻, alg)
        tempM = similar(ϕ)
        d = δ⁻
        isInvA = isinvertible(A)

        A_abs = ReachabilityAnalysis.Exponentiation.elementwise_abs(A)
        Φcache = sum(A) == abs(sum(A)) ? ϕ : nothing
        P2A_abs = ReachabilityAnalysis.Exponentiation.Φ₂(A_abs, δ⁻, alg, isInvA, nothing)


        dU = overapproximate(d * U, Zonotope)
        E_ψ = convert(Zonotope, symmetric_interval_hull(linear_map(P2A_abs, symmetric_interval_hull(linear_map(A, U)))))
        P = minkowski_sum(dU, E_ψ)
        E⁺ = convert(Zonotope, symmetric_interval_hull(linear_map(P2A_abs, symmetric_interval_hull(linear_map(A * A, X0)))))
        lt = minkowski_sum(linear_map(ϕ, X0), dU)
        rt = minkowski_sum(E_ψ, E⁺)
        f = minkowski_sum(lt, rt)
        disc = overapproximate(CH(X0, f), Zonotope)
        Φ₁ = ReachabilityAnalysis.Exponentiation.Φ₁(A, d, alg, isInvA, nothing)
        #P = linear_map(Φ₁, U)
        #println(disc)
        while d < δ⁺
            discritezationDict[d] = copy(disc)
            inputDiscritezationDict[d] = copy(P)
            if maxOrder > 0
                if LazySets.order(P) > maxOrder
                    P = reduce_order(P, reduceOrder)
                end
                if LazySets.order(disc) > maxOrder
                    disc = reduce_order(disc, reduceOrder)
                end
            end
            phiDict[d] = copy(ϕ)
            #tΦ₁ = ReachabilityAnalysis.Exponentiation.Φ₁(A, d, alg, isInvA, Φcache)
            disc = overapproximate(CH(disc, minkowski_sum(P, linear_map(ϕ, disc))), Zonotope)
            P = minkowski_sum(P, linear_map(ϕ, P))
            mul!(tempM, ϕ, ϕ)
            copy!(ϕ, tempM)
            d = d * 2
        end
        if maxOrder > 0
            if LazySets.order(disc) > maxOrder
                disc = reduce_order(disc, reduceOrder)
            end
        end
        discritezationDict[d] = copy(disc)
        phiDict[d] = copy(ϕ)
        inputDiscritezationDict[d] = P

    end
    return discritezationDict, inputDiscritezationDict, phiDict
end


function ReACTDiscretize(A, B, X0::Zonotope{N,Vector{N},Matrix{N}}, U::Nothing, δ⁻, δ⁺, alg::ReachabilityAnalysis.Exponentiation.AbstractExpAlg=ReachabilityAnalysis.Exponentiation.BaseExp, maxOrder::Int=5, reduceOrder::Int=5) where {N}
    XDim, _ = size(genmat(X0))
    phiDict = Dict{Float64,Matrix{Float64}}()
    discritezationDict = Dict{Float64,Zonotope{N,Vector{N},Matrix{N}}}()



    let ϕ::Matrix{Float64} = ReachabilityAnalysis.Exponentiation._exp(A, δ⁻, alg)
        tempM = similar(ϕ)
        d = δ⁻
        dia::Matrix{Float64} = diagm(ones(XDim))
        isInvA = isinvertible(A)
        Φ = ReachabilityAnalysis.Exponentiation._exp(A, δ⁻, alg)
        A_abs = ReachabilityAnalysis.Exponentiation.elementwise_abs(A)
        Φcache = sum(A) == abs(sum(A)) ? Φ : nothing
        P2A_abs = ReachabilityAnalysis.Exponentiation.Φ₂(A_abs, δ⁻, alg, isInvA, Φcache)
        P1A = ReachabilityAnalysis.Exponentiation.Φ₁(A, δ⁻, alg, isInvA, Φcache)

        û = B
        P̂ = P1A * û
        PZ = Zonotope(P̂, zeros(Float64, XDim, 1))
        E⁺ = convert(Zonotope, symmetric_interval_hull(linear_map(P2A_abs, symmetric_interval_hull(linear_map(A * A, X0)))))

        f = minkowski_sum(linear_map(ϕ, X0), E⁺)
        disc = overapproximate(CH(X0, minkowski_sum(f, PZ)), Zonotope)
        while d < δ⁺
            discritezationDict[d] = disc
            if maxOrder > 0
                if LazySets.order(disc) > maxOrder
                    disc = reduce_order(disc, reduceOrder)
                end
            end
            phiDict[d] = copy(ϕ)
            disc = overapproximate(CH(disc, linear_map(ϕ, disc)), Zonotope)
            mul!(tempM, ϕ, ϕ)
            copy!(ϕ, tempM)
            d = d * 2
        end
        if maxOrder > 0
            if LazySets.order(disc) > maxOrder
                disc = reduce_order(disc, reduceOrder)
            end
        end
        discritezationDict[d] = copy(disc)
        phiDict[d] = copy(ϕ)
    end
    return discritezationDict, Nothing, phiDict
end