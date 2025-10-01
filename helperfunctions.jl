using LazySets, LinearAlgebra

function rectangleFromHBox(res::AbstractVector{Shape}, corners)
    for i in 1:size(corners,1)
        tope = getindex(corners, i)
        res[i] = Shape(getindex.(tope, 1), getindex.(tope, 2))
    end
    return res
end

function rectangleFromHBox!(res::AbstractVector{Shape}, cornerss, timestepsize, dim, timescale)
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

function rectangleFromHbox!(res::AbstractVector{Shape}, corners, dim1, dim2)
    for i in 1:size(corners, 1)
        tope = getindex(corners, i)
        res[i] = Shape(getindex.(tope, dim1), getindex.(tope, dim2))
    end
end

function rectangleFromHBoxWithTimestepArray(res::AbstractVector{Shape}, cornerss, timesteps, startTime, dim)

    currentTime = startTime
    for i in 1:size(cornerss, 1)#
        timestep = timesteps[i]

        tope = getindex(cornerss, i)
        dimCoords = getindex.(tope, dim)
        maxcor = maximum(dimCoords)
        mincor = minimum(dimCoords)
        res[i] = Shape([currentTime, currentTime + timestep, currentTime + timestep, currentTime], [mincor, mincor, maxcor, maxcor])
        #res[i] = Shape([deltt*(i-1), deltt*i, deltt*i, deltt*(i-1)], [mincor, mincor, maxcor, maxcor])
        currentTime = currentTime + timestep
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

function initialStep(A, ANorm, timestep, X₀, μ)
    α = (exp(ANorm*timestep)-1-timestep*ANorm)/norm(X₀, Inf)
    β = (exp(ANorm*timestep)-1)*μ/ANorm

    ϕ = exp(A*timestep)

    ϕp = (I+ϕ)/2
    ϕm = (I-ϕ)/2
    gens = hcat(ϕp*X₀.generators,ϕm*X₀.center, ϕm*X₀.generators)


    R = minkowski_sum(Zonotope(ϕp*X₀.center, gens), Zonotope(zeros(dim(X₀)), (α+β)*I(dim(X₀))))

    ballβ = Zonotope(zeros(dim(X₀)), β*I(dim(X₀)))

    return (R, ballβ, ϕ)
end

function reachSetsForTimesteps(A, timesteps, interval, X₀, μ)
    ANorm = norm(A, Inf)

    timestep = timesteps[1] # Use first element

    R₁, ballβ, ϕ = initialStep(A, ANorm, timestep, X₀, μ)

    R = [R₁]

    totalLength = length(timesteps)
    currentTime = timestep

    for i in 2:totalLength
        timestep = timesteps[i]
        if timestep == timesteps[i-1] # Check if timestep is the same
            prevR = R[i-1] # Reuse previous result
            
            # If no input we can avoid the minkowski_sum
            # and thus not add more generators to the zonotope
            if (μ == 0)
                push!(R, linear_map(ϕ, prevR))
            else 
                push!(R, minkowski_sum(linear_map(ϕ, prevR), ballβ))
            end

        else # Recalculate values for new timestep
            newR, ballβ, ϕ = initialStep(A, ANorm, timestep, X₀, μ)
            
            # We move R to the current time
            ϕCurrentTime = exp(A*currentTime)
            newR = Zonotope(ϕCurrentTime * newR.center, ϕCurrentTime * newR.generators)

            push!(R, newR)
        end
        currentTime = currentTime + timestep

    end
    return R
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