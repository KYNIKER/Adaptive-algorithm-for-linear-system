using LazySets, LinearAlgebra

function rectangleFromHBox(corners)
    Shape(getindex.(corners, 1), getindex.(corners, 2))
end

function rectangleFromHBox!(res::AbstractVector{Shape}, cornerss, timestepsize, dim)

    for i in 1:size(cornerss, 1)#ændrer til getindex -> max min 
        tope = getindex(cornerss, i)
        dimCoords = getindex.(tope, dim)
        maxcor = maximum(dimCoords)
        mincor = minimum(dimCoords)
        res[i] = Shape([timestepsize*(i-1), timestepsize*i, timestepsize*i, timestepsize*(i-1)], [mincor, mincor, maxcor, maxcor])
    end
    return res
end

function reachsets(A, timestepsize, interval, X₀, μ)
    T = maximum(interval)
    N = floor(Int, T/timestepsize)
    
    ANorm = norm(A, Inf)
    α = (exp(ANorm*timestepsize)-1-timestepsize*ANorm)/norm(X₀, Inf)
    β = (exp(ANorm*timestepsize)-1)*μ/ANorm

    ϕ = exp(A*timestepsize)

    ϕp = (I+ϕ)/2
    ϕm = (I-ϕ)/2
    gens = hcat(ϕp*X₀.generators,ϕm*X₀.center, ϕm*X₀.generators)

    R₁ = minkowski_sum(Zonotope(ϕp*X₀.center, gens), Zonotope(zeros(2), (α+β)*I(2)))
    #res = [R₁]
    R = [R₁]

    ballβ = Zonotope(zeros(2), β*I(2))

    for i in 2:N
        push!(R, minkowski_sum(linear_map(ϕ, R[i-1]), ballβ))
        #res[i] = R[i]
    end
    return R
end