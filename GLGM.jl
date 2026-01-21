using LazySets, LinearAlgebra, FastExpm, Polyhedra, CDDLib
include("helperfunctions.jl")

function GLGM(A, B, initialTimeStep, interval, X0::Zonotope{N,Vector{N},Matrix{N}}, U::Zonotope, constraint) where {N}
    XG = copy(genmat(X0))
    XDim, p = size(XG)
    time::Float64 = minimum(interval)
    endtime::Float64 = maximum(interval)
    U = concretize(B * U)
    ϕ::Matrix{Float64} = fastExpm(A .* initialTimeStep; threshold=eps(Float64), nonzero_tol=eps(Float64))
    V::Zonotope{N,Vector{N},Matrix{N}} = U
    S::Zonotope{N,Vector{N},Matrix{N}} = X0
    dia::Matrix{Float64} = diagm(ones(XDim))
    if !(zeros(size(U.center)) ∈ U)
        û = copy(U.center)
        invA = inv(Matrix(A))
        Ut = Zonotope(U.center - û, genmat(U))
        dU = overapproximate(initialTimeStep * Ut, Zonotope)
        P = minkowski_sum(dU, E_ψ(Ut, initialTimeStep, A))

        P̂ = invA * (ϕ - I) * û
        lt = minkowski_sum(convert(Zonotope, ϕ * X0), dU)
        rt = minkowski_sum(E_ψ(Ut, initialTimeStep, A), E⁺(X0, initialTimeStep, A))
        PZ = Zonotope(P̂, zeros(Float64, size(U.center, 1), 1))
        f = minkowski_sum(lt, rt)
        disc = overapproximate(CH(X0, minkowski_sum(f, PZ)), Zonotope)
        V = copy(P)
        S = copy(disc)
    else

        dU = overapproximate(initialTimeStep * U, Zonotope)
        P = minkowski_sum(dU, E_ψ(U, initialTimeStep, A))
        lt = concretize(minkowski_sum(convert(Zonotope, ϕ * X0), dU))
        rt = concretize(minkowski_sum(E_ψ(U, initialTimeStep, A), E⁺(X0, initialTimeStep, A)))
        f = concretize(minkowski_sum(lt, rt))
        disc = overapproximate(CH(X0, f), Zonotope)
        V = P
        S = disc
    end
    R = minkowski_sum(S, V)
    #=println(order(R))
    println(V == U)
    println(S == X0)=#
    #=while time <= endtime
        R = minkowski_sum(S, V)
        if order(R) > 50
            R = reduce_order(R, 10)
        end
        #R = reduce_order(R, 5)
        if !mapreduce(c -> ρ(c.a, R) <= c.b, &, constraint)
            println(time)
            println(map(c -> ρ(c.a, R), constraint))
            println(map(c -> c.b, constraint))
            return false
        end
        S = linear_map(ϕ, S)
        V = minkowski_sum(V, linear_map(ϕ, V))
        if order(V) > 50
            println("reduce")
            V = reduce_order(V, 10)
        end
        time = time + initialTimeStep
    end=#
    while time <= endtime
        if order(R) > 5
            R = reduce_order(R, 5)
        end
        #R = reduce_order(R, 5)
        if !mapreduce(c -> ρ(c.a, R) <= c.b, &, constraint)
            println(time)
            println(map(c -> ρ(c.a, R), constraint))
            println(map(c -> c.b, constraint))
            return false
        end
        R = minkowski_sum(V, linear_map(ϕ, R))
        time = time + initialTimeStep
    end
    return true
end