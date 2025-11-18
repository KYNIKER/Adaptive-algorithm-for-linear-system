using LazySets, LinearAlgebra, Plots, FastExpm

using Random # For strategy 4

function rectangleFromHBox(res::AbstractVector{Shape}, corners)
    for i in 1:size(corners,1)
        tope = getindex(corners, i)
        res[i] = Shape(getindex.(tope, 1), getindex.(tope, 2))
    end
    return res
end

function rectangleFromHBox!(res::AbstractVector{Shape}, cornerss, timestepsize, dim, timescale::Float64)
    deltt = timestepsize
    currentTime = 0
    for i in 1:size(cornerss, 1)#ændrer til getindex -> max min 
        tope = getindex(cornerss, i)
        dimCoords = getindex.(tope, dim)
        maxcor = maximum(dimCoords)
        mincor = minimum(dimCoords)
        res[i] = Shape([currentTime, currentTime + deltt, currentTime + deltt, currentTime], [mincor, mincor, maxcor, maxcor])
        #res[i] = Shape([deltt*(i-1), deltt*i, deltt*i, deltt*(i-1)], [mincor, mincor, maxcor, maxcor])
        currentTime = currentTime + deltt
        deltt = deltt+timescale*timestepsize
    end
    return res
end

function rectangleFromHBox2Dims(res::AbstractVector{Shape}, corners, dim1, dim2)
    for i in 1:size(corners, 1)
        tope = getindex(corners, i)
        # Switch 3rd and 4th index
        tope[3], tope[4] = tope[4], tope[3]

        res[i] = Shape(getindex.(tope, dim1), getindex.(tope, dim2))
    end
end

function rectangleFromHBoxWithTimestepArray(res::AbstractVector{Shape}, cornerss, timesteps::Vector{Float64}, startTime, dim::Int64)

    currentTime = startTime
    for i in 1:size(cornerss, 1)
        timestep::Float64 = timesteps[i]::Float64

        tope::Vector{Vector{Float64}} = getindex(cornerss, i)::Vector{Vector{Float64}}
        dimCoords::Vector{Float64} = getindex.(tope, dim)
        maxcor::Float64 = maximum(dimCoords)::Float64
        mincor::Float64 = minimum(dimCoords)::Float64
        res[i] = Shape([currentTime, currentTime + timestep, currentTime + timestep, currentTime], [mincor, mincor, maxcor, maxcor])
        #res[i] = Shape([deltt*(i-1), deltt*i, deltt*i, deltt*(i-1)], [mincor, mincor, maxcor, maxcor])
        currentTime::Float64 = currentTime + timestep::Float64
    end
    return res
end

function rectangleFromHBox!(res::AbstractVector{Shape}, cornerss, timestepsize, dim, timesteps)

    for i in 1:size(cornerss, 1)#ændrer til getindex -> max min 
        tope = getindex(cornerss, i)
        dimCoords = getindex.(tope, dim)
        maxcor = maximum(dimCoords)
        mincor = minimum(dimCoords)
        curtime = timesteps[i]
        res[i] = Shape([timestepsize*(curtime-1), timestepsize*curtime, timestepsize*curtime, timestepsize*(curtime-1)], [mincor, mincor, maxcor, maxcor])
    end
    return res
end

function reachsets(A, timestepsize, interval, X₀, μ, timescale)
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
    i = 2
    deltt=timestepsize
    time = timestepsize # We have already done one timestep
    #while i < N
    while time < T
        deltt = deltt+timescale*timestepsize
        #α = (exp(ANorm*timestepsize)-1-timestepsize*ANorm)/norm(X₀, Inf)
        β = (exp(ANorm*deltt)-1)*μ/ANorm

        ϕ = exp(A*deltt)

        ballβ = Zonotope(zeros(dim(X₀)), β*I(2))
        push!(R, minkowski_sum(linear_map(ϕ, R[i-1]), ballβ))
        #res[i] = R[i]



        #N = floor(Int, T/deltt) # Re-estimate amount of timesteps
        time = time + deltt
        i = i + 1
    end
    return R
end

function initialStep(A::Matrix{Float64}, ANorm::Float64, timestep::Float64, X₀::Zonotope{Float64, Vector{Float64}, Matrix{Float64}}, μ::Float64)
    
    # Splitting this up so the compiler doesn't get confused
    step1::Float64 = exp(ANorm * timestep)::Float64 - 1 - timestep * ANorm
    step2::Float64 = manualNormInf(X₀)
    α::Float64 = step1 / step2
    #α = (exp(ANorm*timestep)-1-timestep*ANorm)/norm(X₀, Inf)
    β::Float64 = (exp(ANorm*timestep)::Float64-1)*μ/ANorm

    ϕ = exp(A*timestep)::Matrix{Float64}

    ϕp::Matrix{Float64} = (I+ϕ)/2
    ϕm::Matrix{Float64} = (I-ϕ)/2
    gens = hcat(ϕp*X₀.generators,ϕm*X₀.center, ϕm*X₀.generators)::Matrix{Float64}

    R::Zonotope{Float64, Vector{Float64}, Matrix{Float64}} = minkowski_sum(Zonotope(ϕp*X₀.center, gens), Zonotope(zeros(dim(X₀)), (α+β)*I(dim(X₀))))

    ballβ = Zonotope(zeros(dim(X₀)), diagm(β*ones(dim(X₀))))

    return (R, ballβ, ϕ)
end

function manualNormInf(Z::Zonotope{Float64, Vector{Float64}, Matrix{Float64}})
    dimensions = dim(Z)::Int64
    maxVal::Float64 = Z.center[1] + sum(Z.generators[1])
    for i in 2:dimensions
        maxVal = max(maxVal, Z.center[i] + sum(Z.generators[i]))
    end

    return maxVal
end


function initialStepNoInput(A::Matrix{Float64}, Anorm::Float64, timestep::Float64, X₀::Zonotope{Float64, Vector{Float64}, Matrix{Float64}})
    
    # Splitting this up so the compiler doesn't get confused
    step1 = exp(Anorm * timestep)::Float64 - 1 - timestep * Anorm
    step2::Float64 = manualNormInf(X₀)
    α = step1 / step2
    
    #α::Float64 = (exp(ANorm*timestep)-1-timestep*ANorm)/manualNormInf(X₀)::Float64

    ϕ = exp(A*timestep)::Matrix{Float64}

    ϕp::Matrix{Float64} = (I+ϕ)/2
    ϕm::Matrix{Float64} = (I-ϕ)/2

    gens = hcat(ϕp*X₀.generators,ϕm*X₀.center, ϕm*X₀.generators)::Matrix{Float64}

    R::Zonotope{Float64, Vector{Float64}, Matrix{Float64}} = minkowski_sum(Zonotope(ϕp*X₀.center, gens), Zonotope(zeros(dim(X₀)), α*I(dim(X₀))))

    return R, ϕ
end


function reachSetsForTimesteps(A, timesteps, interval, X₀, μ)
    ANorm = norm(A, Inf)

    timestep = timesteps[1] # Use first element

    R₁, ballβ, ϕ = initialStep(A, ANorm, timestep, X₀, μ)

    
    R = [R₁]

    # We only further utilize these if we have an input
    Ω = []
    S = ballβ
    V = ballβ
    if (μ != 0)
        push!(Ω, minkowski_sum(R₁, ballβ))
        V = linear_map(ϕ, V)
    end

    totalLength = length(timesteps)
    currentTime = timestep

    for i in 2:totalLength
        timestep = timesteps[i]
        if timestep == timesteps[i-1] # Check if timestep is the same
            prevR = R[i-1] # Reuse previous result

            newR = linear_map(ϕ, prevR)
            push!(R, newR)

            # If we have an input
            if (μ != 0)
                S = minkowski_sum(S, V)
                V = linear_map(ϕ, V)

                push!(Ω, minkowski_sum(newR, S))
            end

        else # Recalculate values for new timestep
            newR, ballβ, ϕ = initialStep(A, ANorm, timestep, X₀, μ)
            
            # We move R to the current time
            ϕCurrentTime = exp(A*currentTime)
            newR = Zonotope(ϕCurrentTime * newR.center, ϕCurrentTime * newR.generators)

            push!(R, newR)

            if (μ != 0)
                # We have to retroactivily calculate the new S
                steps = floor(Int, currentTime / timestep)
                S = ballβ
                V = ballβ
                for j in 1:steps
                    S = minkowski_sum(S, V)
                    V = linear_map(ϕ, V)
                end
                
                push!(Ω, minkowski_sum(newR, S))
            end
        end
        currentTime = currentTime + timestep
    end
    # If no input just return R
    if (μ == 0)
        return R
    end
    # Else return Ω
    return Ω
end

# Move a zonotope forward in time, assumes no input
function forwardTimeNoInput(A, R, currentTime) 
    ϕCurrentTime = exp(A*currentTime)
    return linear_map(ϕCurrentTime, R)
end

# We reuse the previous S, and do not recompute it
function forwardTimeReUse(A, R, currentTime, ϕ, ballβ, prevS) 
    newR = linear_map(exp(A*currentTime), R)
    #newR = forwardTimeNoInput(A, R, currentTime)
    V = linear_map(exp(A*currentTime), ballβ)

    S = minkowski_sum(prevS, V)
    V = linear_map(ϕ, V)

    newΩ = minkowski_sum(newR, S)

    #println("Set S, V: $S, $V with steps $steps")
    return (newR, S, V, newΩ)
end

function forwardTime(A :: Matrix{Float64}, R :: Zonotope{Float64, Vector{Float64}, Matrix{Float64}}, currentTime :: Float64, timestep :: Float64, ϕ :: Matrix{Float64}, prevTime :: Float64, prevS :: Zonotope{Float64, Vector{Float64}, Matrix{Float64}}, prevV :: Zonotope{Float64, Vector{Float64}, Matrix{Float64}})
    newR :: Zonotope{Float64, Vector{Float64}, Matrix{Float64}} = linear_map(exp(A*currentTime), R)
    #newR :: Zonotope = forwardTimeNoInput(A, R, currentTime)
    steps = floor(Int, (currentTime - prevTime) / timestep)
    S = prevS 
    V = prevV 
    for j in 1:steps
        S :: Zonotope{Float64, Vector{Float64}, Matrix{Float64}} = minkowski_sum(S, V) # this one
        V :: Zonotope{Float64, Vector{Float64}, Matrix{Float64}} = linear_map(ϕ, V) # this one
    end
    newΩ :: Zonotope{Float64, Vector{Float64}, Matrix{Float64}} = minkowski_sum(newR, S)

    #println("Set S, V: $S, $V with steps $steps")
    return (newR, S, V, newΩ)
end

function reachsetslog(A, timestepsize, interval, X₀, μ)
    T = maximum(interval)
    N = floor(Int, T/timestepsize)
    
    bits = ceil(Int, log2(N))
    println(N)
    println(bits)
    ϕs = Vector{AbstractMatrix}(undef, bits)

    ballβs = Vector{Zonotope}(undef, bits)



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
    ts = [1]

    ballβ = Zonotope(zeros(2), β*I(2))

    ϕs[1] = ϕ
    ballβs[1] = ballβ

    for i in 2:bits
        ϕs[i] = ϕs[i-1]*ϕs[i-1]
        ballβs[i] = minkowski_sum(ballβs[i-1], ballβs[i-1])
    end

    for i in 2:bits
        if ((N >> (i)) & 1) == 1
            push!(R, minkowski_sum(linear_map(ϕs[i], R[end]), ballβs[i]))
            push!(ts, ts[end]+(2^(i)))
        end
            #res[i] = R[i]
    end
    return R, ts
end

function approximateArea(zonotope::Zonotope)
    zonotope |> box_approximation |> area
end

function approximateArea(zonotope::AbstractVector{Zonotope{Float64, Vector{Float64}, Matrix{Float64}}})
    area(reduce(∪, map(box_approximation, zonotope)))
end

function approximateArea(zonotope::Vector{Any})
    area(reduce(∪, map(box_approximation, zonotope)))
end

#intersects function is based upon "Set operations and order reductions for constrained zonotopes" - Vignesh Raghuraman, Justin P. Koeln
#Note that this only returns true if they intersect and does not depend on whether the intersection between them is empty.
function intersectss(c::Vector{Float64}, G::Matrix{Float64}, h::Vector{Float64}, f::Float64, tempM :: Matrix{Float64}) :: Bool #Kan måske speedes op ved at give et ekstra argument tempM som bliver brugt til at udregne højresiden så at vi ikke altid assigner memory til den operation.
    #c::Vector{Float64} = zonotope.center
    #G::Matrix{Float64} = genmat(zonotope)
    #println(h, reshape(halfspace.a, (1,:)))
    #l::Float64 = abs(f - dot(h, c))
    #r::Float64 = sum(abs.(G .* h))
    #mul!(tempM, G, h)
    return abs(f - dot(h, c)) <= sum(abs.(G .* h))#mapreduce(abs, +, G .* h)#sum(abs.(G * h))
end

#=function intersects(zonotope::Zonotope, halfspace::HalfSpace)
    h, f = tosimplehrep(halfspace)
    l = abs(only(f) - dot(h', zonotope.center))
    r = sum(abs.(zonotope.generators .* h'))
    return l <= r
end=#
function intersects(zonotope::Zonotope, halfspace::LazySets.HalfSpace) :: Bool
    c::Vector{Float64} = zonotope.center
    G::Matrix{Float64} = genmat(zonotope)
    h::Vector{Float64} = halfspace.a
    f::Float64 = halfspace.b 
    #d = c + h
    #d::Vector{Float64} = c.*h
    l::Float64 = abs(f - dot(c, h))#sum(c.* h)
    #l = abs(l)# < 0 ? -1 * l : l
    #t::Matrix{Float64} = h'
    #mul!(tempM, h', G)
    r::Float64 = abs_sum(h, G)
    return l <= r#false
end

function intersects(halfspace::HalfSpace, zonotope::Zonotope)
    intersects(zonotope, halfspace)
end

function intersects(z1::Zonotope, z2::Zonotope) #Should be possible to compute this from above paper by viewing the zonotopes as being constrained but with no constrains and whether their intersecting zonotope has empty 
    isdisjoint(z1, z2)
end

#=function intersectss(c::Vector{Float64}, G::Matrix{Float64}, h::Vector{Float64}, f::Float64) :: Bool
    #c::Vector{Float64} = zonotope.center
    #G::Matrix{Float64} = genmat(zonotope)
    #println(h, reshape(halfspace.a, (1,:)))
    l::Float64 = abs(f - dot(h, c))
    r::Float64 = sum(abs.(G .* h))
    return l <= r
end=#

function smallestDeterminant(A, timestepsizes)
    trace = tr(A)
    return timestepsizes[argmin(map(exp, timestepsizes * trace))]
end


function linear_map_zonotope_nD(M::AbstractMatrix, c::Vector{Float64}, G::Matrix{Float64})
    #c = M * LazySets.center(Z)
    #gi = M * genmat(Z)
    c = M * c
    G = M * G
    return Zonotope(c, G)
end
# https://github.com/JuliaReach/ReachabilityAnalysis.jl/blob/master/src/Discretization/Exponentiation.jl l:341
function Φ₂(A, δ, invertible = false)
    let A = abs.(A)
        if invertible
            Aδ = A .* δ
            Φ = fastExpm(Aδ)#; threshold=eps(Float64), nonzero_tol=eps(Float64))
            n = size(A, 1)
            N = eltype(A)
            In = Matrix(one(N) * I, n, n)
            B = Φ - In - Aδ
            Ainv = inv(Matrix(A))
            Ainvsqr = Ainv^2
            return Ainvsqr * B
        else
            n = checksquare(A)
            B = _P_3n(A, δ, n)
            P = fastExpm(B; threshold=eps(Float64), nonzero_tol=eps(Float64))
            return _P₂_blk(P, n)
        end 
    end
end

@inline function _P_3n(A::AbstractMatrix{N}, δ, n) where {N}
    return [Matrix(A * δ) Matrix(δ * I, n, n) zeros(n, n);
            zeros(n, 2 * n) Matrix(δ * I, n, n);
            zeros(n, 3 * n)]::Matrix{N}
end

@inline _P₂_blk(P, n) = P[1:n, (2 * n + 1):(3 * n)]

function E⁺(X, δ, A)
    return convert(Zonotope, symmetric_interval_hull(Φ₂(A, δ, isinvertible(A)) * symmetric_interval_hull(A * A * X)))
end

function E_ψ(U, δ, A)
    return convert(Zonotope, symmetric_interval_hull(Φ₂(A, δ, isinvertible(A)) * symmetric_interval_hull(A * U)))
end