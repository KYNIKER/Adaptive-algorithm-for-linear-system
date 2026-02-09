using LazySets, LinearAlgebra

function ReACT(A, B, initialTimeStep, interval, X0::Zonotope{N,Vector{N},Matrix{N}}, U::Zonotope, constraint, Digits::Integer, STRATEGY::Integer) where {N}
    XDim, _ = size(genmat(X0))
    m = initialTimeStep / 2^(ceil(Integer, log2(initialTimeStep)) + ceil(Integer, -log2(10.0^(-Digits))) - 1)   #Calculate the smallest number larger than 10^-Digits obtained by repeatedly dividing initialTimeStep by 2.
    changedTimeStep = true
    elems = (ceil(Integer, log2(initialTimeStep)) + ceil(Integer, -log2(10.0^(-Digits))) - 1)
    phiDict = Dict{Float64,Matrix{Float64}}()
    sizehint!(phiDict, elems)
    discritezationDict = Dict{Float64,Zonotope{N,Vector{N},Matrix{N}}}()
    sizehint!(discritezationDict, elems)
    inputDiscritezationDict = Dict{Float64,Zonotope{N,Vector{N},Matrix{N}}}()
    sizehint!(inputDiscritezationDict, elems)

    constraintProjVectors = map(x -> x.a, constraint)
    constraintProjBounds = ρ.(constraintProjVectors, constraint)



    time::Float64 = minimum(interval)
    endtime::Float64 = maximum(interval)

    currentTimeStep = copy(initialTimeStep)

    attemptsRecorder = Integer[]

    U = concretize(B * U)

    let ϕ::Matrix{Float64} = exp(A .* m)
        tempM = similar(ϕ)
        d = m
        dia::Matrix{Float64} = diagm(ones(XDim))

        if !(zeros(XDim) ∈ U) #Origin is *not* in input
            û = copy(U.center)
            invA = inv(Matrix(A))
            Ut = Zonotope(U.center - û, genmat(U))
            dU = overapproximate(d * Ut, Zonotope)
            P = minkowski_sum(dU, E_ψ(Ut, d, A))


            P̂ = invA * (ϕ - dia) * û
            lt = minkowski_sum(convert(Zonotope, ϕ * X0), dU)
            rt = minkowski_sum(E_ψ(Ut, d, A), E⁺(X0, d, A))
            PZ = Zonotope(P̂, zeros(Float64, size(U.center, 1), 1))
            f = minkowski_sum(lt, rt)
            disc = overapproximate(CH(X0, minkowski_sum(f, PZ)), Zonotope)
            lt = missing
            rt = missing
            PZ = missing
            f = missing

            while d < initialTimeStep
                inputDiscritezationDict[d] = P
                P = minkowski_sum(P, linear_map(ϕ, P))
                discritezationDict[d] = disc
                phiDict[d] = copy(ϕ)
                disc = overapproximate(CH(disc, linear_map(ϕ, disc)), Zonotope)

                mul!(tempM, ϕ, ϕ)
                copy!(ϕ, tempM)
                d = d * 2
            end

            discritezationDict[d] = copy(disc)
            phiDict[d] = copy(ϕ)
            inputDiscritezationDict[initialTimeStep] = P
        else
            dU = overapproximate(d * U, Zonotope)
            P = minkowski_sum(dU, E_ψ(U, d, A))
            lt = concretize(minkowski_sum(convert(Zonotope, ϕ * X0), dU))
            rt = concretize(minkowski_sum(E_ψ(U, d, A), E⁺(X0, d, A)))
            f = concretize(minkowski_sum(lt, rt))
            disc = overapproximate(CH(X0, f), Zonotope)
            while d < initialTimeStep
                inputDiscritezationDict[d] = P
                P = minkowski_sum(P, linear_map(ϕ, P))


                phiDict[d] = copy(ϕ)
                discritezationDict[d] = copy(disc)
                disc = overapproximate(CH(disc, linear_map(ϕ, disc)), Zonotope)
                mul!(tempM, ϕ, ϕ)
                copy!(ϕ, tempM)
                d = d * 2
            end
            discritezationDict[d] = copy(disc)
            phiDict[d] = copy(ϕ)
            inputDiscritezationDict[initialTimeStep] = P
        end
    end

    V::Zonotope{N,Vector{N},Matrix{N}} = copy(inputDiscritezationDict[initialTimeStep])
    Sρ = zeros(Float64, length(constraint))
    newR::Zonotope{N,Vector{N},Matrix{N}} = discritezationDict[initialTimeStep]
    i = 1


    Φ::Matrix{Float64} = diagm(ones(Float64, size(A, 2)))
    tempM = similar(Φ)
    ϕt = similar(Φ)
    newRR = copy(newR)


    while time < endtime

        attempts = 1
        approveFlag = false

        while !approveFlag
            if currentTimeStep < m
                return false
            end

            if changedTimeStep
                newR = discritezationDict[currentTimeStep]
                V = copy(inputDiscritezationDict[currentTimeStep])
                ϕt = phiDict[currentTimeStep]
                newRR = linear_map(Φ, newR)
                V = linear_map(Φ, V)
            else
                newRR = linear_map(ϕt, newRR)
                V = linear_map(ϕt, V)
            end

            changedTimeStep = false
            hom = map(x -> ρ(x, newRR), constraintProjVectors)
            inhom = map(x -> ρ(x, V), constraintProjVectors)

            if reduce(&, <=(Sρ + hom + inhom, constraintProjBounds))
                approveFlag = true
                Sρ += map(x -> ρ(x, V), constraintProjVectors) + map(x -> dot(V.center, x), constraintProjVectors)
                mul!(tempM, Φ, ϕt)
                copy!(Φ, tempM)
            else
                newR = copy(newR)
                currentTimeStep = currentTimeStep / 2
                changedTimeStep = true
                attempts = attempts + 1
            end
        end


        push!(attemptsRecorder, attempts)
        i = i + 1
        time = time + currentTimeStep

        # Reset / apply strategy
        # Only do this if the current timestep is less than the initial
        if STRATEGY == 0
            # Only reduce
        elseif STRATEGY == 1
            # always try double
            if currentTimeStep < initialTimeStep
                currentTimeStep = currentTimeStep * 2
                changedTimeStep = true
            end
        elseif STRATEGY == 2
            # If attemptsrecorder past 4 are successes, double timestep
            if currentTimeStep < initialTimeStep
                lowest = min(4, i - 1)
                window = @view attemptsRecorder[i-lowest:i-1]
                if all(window .== 1)
                    currentTimeStep = currentTimeStep * 2
                    changedTimeStep = true
                end
            end
        end
    end

    return true
end

function ReACT(A, B, initialTimeStep, interval, X0::Zonotope{N,Vector{N},Matrix{N}}, U::Nothing, constraint, Digits::Integer, STRATEGY::Integer) where {N}
    XG = copy(genmat(X0))
    XDim, _ = size(XG)
    m = initialTimeStep / 2^(ceil(Integer, log2(initialTimeStep)) + ceil(Integer, -log2(10.0^(-Digits))) - 1)   #Calculate the smallest number larger than 10^-Digits obtained by repeatedly dividing initialTimeStep by 2.
    changedTimeStep = true

    elems = (ceil(Integer, log2(initialTimeStep)) + ceil(Integer, -log2(10.0^(-Digits))) - 1)
    phiDict = Dict{Float64,Matrix{Float64}}()
    sizehint!(phiDict, elems)
    discritezationDict = Dict{Float64,Zonotope{N,Vector{N},Matrix{N}}}()
    sizehint!(discritezationDict, elems)

    constraintProjVectors = map(x -> x.a, constraint)
    constraintProjBounds = ρ.(constraintProjVectors, constraint)



    time::Float64 = minimum(interval)
    endtime::Float64 = maximum(interval)

    currentTimeStep = copy(initialTimeStep)

    attemptsRecorder = Integer[]


    let ϕ::Matrix{Float64} = exp(A .* m)
        tempM = similar(ϕ)
        d = m
        dia::Matrix{Float64} = diagm(ones(XDim))
        û = B
        P̂ = (A \ (ϕ - dia)) * û
        PZ = Zonotope(P̂, zeros(Float64, XDim, 1))
        f = minkowski_sum(convert(Zonotope, ϕ * X0), E⁺(X0, d, A))
        disc = overapproximate(CH(X0, minkowski_sum(f, PZ)), Zonotope)
        PZ = missing
        f = missing

        while d < initialTimeStep
            discritezationDict[d] = disc
            phiDict[d] = copy(ϕ)
            disc = overapproximate(CH(disc, linear_map(ϕ, disc)), Zonotope)

            mul!(tempM, ϕ, ϕ)
            copy!(ϕ, tempM)
            d = d * 2
        end

        discritezationDict[d] = copy(disc)
        phiDict[d] = copy(ϕ)
    end

    newR::Zonotope{N,Vector{N},Matrix{N}} = discritezationDict[initialTimeStep]
    i = 1


    Φ::Matrix{Float64} = diagm(ones(Float64, size(A, 2)))
    tempM = similar(Φ)
    ϕt = similar(Φ)
    newRR = copy(newR)
    initialϕ = phiDict[initialTimeStep]

    while time < endtime

        attempts = 1
        approveFlag = false

        while !approveFlag
            if currentTimeStep < m

                return false
            end

            if changedTimeStep
                newR = discritezationDict[currentTimeStep]
                ϕt = phiDict[currentTimeStep]
                newRR = linear_map(Φ, newR)
            else
                newRR = linear_map(ϕt, newRR)
            end

            changedTimeStep = false

            hom = map(x -> ρ(x, newRR), constraintProjVectors)
            if reduce(&, <=(hom, constraintProjBounds))
                approveFlag = true
                mul!(tempM, Φ, ϕt)
                copy!(Φ, tempM)
            else
                newR = copy(newR)
                currentTimeStep = currentTimeStep / 2
                changedTimeStep = true
                attempts = attempts + 1
            end
        end

        push!(attemptsRecorder, attempts)
        i = i + 1
        time = time + currentTimeStep

        # Reset / apply strategy
        # Only do this if the current timestep is less than the initial
        if STRATEGY == 0
            # Only reduce
        elseif STRATEGY == 1
            # always try double
            if currentTimeStep < initialTimeStep
                currentTimeStep = currentTimeStep * 2
                changedTimeStep = true
            end
        elseif STRATEGY == 2
            # If attemptsrecorder past 4 are successes, double timestep
            if currentTimeStep < initialTimeStep
                lowest = min(4, i - 1)
                window = @view attemptsRecorder[i-lowest:i-1]
                if all(window .== 1)
                    currentTimeStep = currentTimeStep * 2
                    changedTimeStep = true
                end
            end
        end
    end

    return true
end