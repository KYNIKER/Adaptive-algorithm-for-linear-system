using LazySets, LinearAlgebra
include("ReACTDiscretize.jl")

function ReACT(A, B, initialTimeStep, interval, X0::Zonotope{N,Vector{N},Matrix{N}}, U::Zonotope, constraint, δ⁻, STRATEGY::Integer, alg::ReachabilityAnalysis.Exponentiation.AbstractExpAlg=ReachabilityAnalysis.Exponentiation.BaseExp, maxOrder::Int=5, reduceOrder::Int=5) where {N}

    m = copy(δ⁻)
    changedTimeStep = true
    #elems = (ceil(Integer, log2(initialTimeStep)) + ceil(Integer, -log2(10.0^(-Digits))) - 1)
    phiDict = Dict{Float64,Matrix{Float64}}()
    #sizehint!(phiDict, elems)
    discritezationDict = Dict{Float64,Zonotope{N,Vector{N},Matrix{N}}}()
    #sizehint!(discritezationDict, elems)
    inputDiscritezationDict = Dict{Float64,Zonotope{N,Vector{N},Matrix{N}}}()
    #sizehint!(inputDiscritezationDict, elems)

    constraintProjVectors = map(x -> x.a, constraint)
    constraintProjBounds = ρ.(constraintProjVectors, constraint)

    discritezationDict, inputDiscritezationDict, phiDict = ReACTDiscretize(A, B, X0, U, m, initialTimeStep, alg, maxOrder, reduceOrder)

    time::Float64 = minimum(interval)
    endtime::Float64 = maximum(interval)

    currentTimeStep = copy(initialTimeStep)

    attemptsRecorder = Integer[]

    V::Zonotope{N,Vector{N},Matrix{N}} = copy(inputDiscritezationDict[initialTimeStep])
    Sρ = zeros(Float64, length(constraint))
    newR::Zonotope{N,Vector{N},Matrix{N}} = discritezationDict[initialTimeStep]
    i = 1


    Φ::Matrix{Float64} = diagm(ones(Float64, size(A, 2)))
    tempM = similar(Φ)
    ϕt = similar(Φ)
    newRR = copy(newR)

    # Main simulation loop

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
            #hom = map(x -> ρ(x, newRR), constraintProjVectors)
            inhom = map(x -> ρ(x, V), constraintProjVectors)

            if all((inputAccumulated + inputCurrent + ρ(x, newRR)) <= y for (inputAccumulated, inputCurrent, x, y) in zip(Sρ, inhom, constraintProjVectors, constraintProjBounds))
                #if reduce(&, <=(Sρ + hom + inhom, constraintProjBounds))
                approveFlag = true
                Sρ += inhom
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

#  React with no input
function ReACT(A, B, initialTimeStep, interval, X0::Zonotope{N,Vector{N},Matrix{N}}, U::Nothing, constraint, δ⁻, STRATEGY::Integer, alg::ReachabilityAnalysis.Exponentiation.AbstractExpAlg=ReachabilityAnalysis.Exponentiation.BaseExp, maxOrder::Int=5, reduceOrder::Int=5) where {N}
    XG = copy(genmat(X0))
    XDim, _ = size(XG)
    m = copy(δ⁻)
    changedTimeStep = true

    #elems = (ceil(Integer, log2(initialTimeStep)) + ceil(Integer, -log2(10.0^(-Digits))) - 1)
    phiDict = Dict{Float64,Matrix{Float64}}()
    #sizehint!(phiDict, elems)
    discritezationDict = Dict{Float64,Zonotope{N,Vector{N},Matrix{N}}}()
    #sizehint!(discritezationDict, elems)

    discritezationDict, _, phiDict = ReACTDiscretize(A, B, X0, U, m, initialTimeStep, alg, maxOrder, reduceOrder)

    constraintProjVectors = map(x -> x.a, constraint)
    constraintProjBounds = ρ.(constraintProjVectors, constraint)

    time::Float64 = minimum(interval)
    endtime::Float64 = maximum(interval)

    currentTimeStep = copy(initialTimeStep)

    attemptsRecorder = Integer[]

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

            #hom = map(x -> ρ(x, newRR), constraintProjVectors)
            if all(ρ(x, newRR) <= y for (x, y) in zip(constraintProjVectors, constraintProjBounds))
                #if reduce(&, <=(hom, constraintProjBounds))
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

# This is solely used for plotting, not benchmarking
function PlotReACT(A, B, initialTimeStep, interval, X0::Zonotope{N,Vector{N},Matrix{N}}, U::Zonotope, constraint, δ⁻, STRATEGY::Integer, alg::ReachabilityAnalysis.Exponentiation.AbstractExpAlg=ReachabilityAnalysis.Exponentiation.BaseExp, maxOrder::Int=5, reduceOrder::Int=5) where {N}
    m = copy(δ⁻)
    # Override for 1 timestep system for plotting
    #m = initialTimeStep / 2^(ceil(Integer, log2(initialTimeStep)) + ceil(Integer, -log2(10.0^(-Digits))) - 1)   #Calculate the smallest number larger than 10^-Digits obtained by repeatedly dividing initialTimeStep by 2.
    changedTimeStep = true
    #elems = (ceil(Integer, log2(initialTimeStep)) + ceil(Integer, -log2(10.0^(-Digits))) - 1)
    phiDict = Dict{Float64,Matrix{Float64}}()
    #sizehint!(phiDict, elems)
    discritezationDict = Dict{Float64,Zonotope{N,Vector{N},Matrix{N}}}()
    #sizehint!(discritezationDict, elems)
    inputDiscritezationDict = Dict{Float64,Zonotope{N,Vector{N},Matrix{N}}}()
    #sizehint!(inputDiscritezationDict, elems)

    constraintProjVectors = map(x -> x.a, constraint)
    constraintProjBounds = map(x -> x.b, constraint)#ρ.(constraintProjVectors, constraint)#


    discritezationDict, inputDiscritezationDict, phiDict = ReACTDiscretize(A, B, X0, U, m, initialTimeStep, alg, maxOrder, reduceOrder)

    time::Float64 = minimum(interval)
    endtime::Float64 = maximum(interval)

    currentTimeStep = copy(initialTimeStep)

    attemptsRecorder = Integer[]
    V::Zonotope{N,Vector{N},Matrix{N}} = copy(inputDiscritezationDict[initialTimeStep])
    Sρ = zeros(Float64, length(constraint))
    newR::Zonotope{N,Vector{N},Matrix{N}} = discritezationDict[initialTimeStep]
    i = 1

    reachSet = Vector{Zonotope{N,Vector{N},Matrix{N}}}()
    timeStepRecorder = Vector{Float64}()

    #println(discritezationDict[m].center)

    Φ::Matrix{Float64} = diagm(ones(Float64, size(A, 2)))
    tempM = similar(Φ)
    ϕt = similar(Φ)
    newRR = copy(newR)
    #lastV = Zonotope(zero(discritezationDict[m].center), [zero(discritezationDict[m].center)])
    Ub = overapproximate(linear_map(B, U), Zonotope)
    while time < endtime

        attempts = 1
        approveFlag = false

        while !approveFlag
            if currentTimeStep < m
                #currentTimeStep = m
                #println(constraintProjBounds)
                return reachSet, timeStepRecorder, attemptsRecorder
            end

            if changedTimeStep
                newR = discritezationDict[currentTimeStep]
                #V = copy(inputDiscritezationDict[currentTimeStep])
                ϕt = phiDict[currentTimeStep]
                newRR = linear_map(Φ, discritezationDict[currentTimeStep])
                #V = linear_map(Φ, V)
            else
                newRR = linear_map(ϕt, newRR)
                #V = linear_map(ϕt, V)
            end
            V = linear_map(ReachabilityAnalysis.Exponentiation.Φ₁(A, time, alg, false, nothing), Ub)
            changedTimeStep = false
            #hom = map(x -> ρ(x, newRR), constraintProjVectors)
            #println(hom)
            inhom = map(x -> ρ(x, V), constraintProjVectors)

            if all((input + ρ(x, newRR)) < y for (input, x, y) in zip(Sρ, constraintProjVectors, constraintProjBounds))
                approveFlag = true
                push!(reachSet, minkowski_sum(newRR, V)) #minkowski_sum(newRR, V)) # Add homogeneous and input to reachSet
                Sρ = copy(inhom)

                mul!(tempM, Φ, ϕt)
                copy!(Φ, tempM)
            else
                #=if attempts == 1
                    println(Sρ, " ", hom, " ", inhom)
                end=#

                newR = copy(newR)
                currentTimeStep = currentTimeStep / 2
                changedTimeStep = true
                attempts = attempts + 1
            end
        end


        push!(attemptsRecorder, attempts)
        push!(timeStepRecorder, currentTimeStep)

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
    return reachSet, timeStepRecorder, attemptsRecorder
end

# This shows off our approach with spaceex discretization. Also only for plotting
function PlotReACTIndividualDisc(A, B, initialTimeStep, interval, X0::Zonotope{N,Vector{N},Matrix{N}}, U::Zonotope, constraint, δ⁻, STRATEGY::Integer, alg::ReachabilityAnalysis.Exponentiation.AbstractExpAlg=ReachabilityAnalysis.Exponentiation.BaseExp, maxOrder::Int=5, reduceOrder::Int=5) where {N}
    # Here we use typical spaceex discretization for each timestep size

    m = copy(δ⁻)

    #m = initialTimeStep / 2^(ceil(Integer, log2(initialTimeStep)) + ceil(Integer, -log2(10.0^(-Digits))) - 1)   #Calculate the smallest number larger than 10^-Digits obtained by repeatedly dividing initialTimeStep by 2.
    changedTimeStep = true
    #elems = (ceil(Integer, log2(initialTimeStep)) + ceil(Integer, -log2(10.0^(-Digits))) - 1)
    phiDict = Dict{Float64,Matrix{Float64}}()
    #sizehint!(phiDict, elems)
    discritezationDict = Dict{Float64,Zonotope{N,Vector{N},Matrix{N}}}()
    #sizehint!(discritezationDict, elems)
    inputDiscritezationDict = Dict{Float64,Zonotope{N,Vector{N},Matrix{N}}}()
    #sizehint!(inputDiscritezationDict, elems)

    constraintProjVectors = map(x -> x.a, constraint)
    constraintProjBounds = map(x -> x.b, constraint)#ρ.(constraintProjVectors, constraint)#


    phiDict = Dict{Float64,Matrix{Float64}}()
    discritezationDict = Dict{Float64,Zonotope{N,Vector{N},Matrix{N}}}()
    inputDiscritezationDict = Dict{Float64,Zonotope{N,Vector{N},Matrix{N}}}()

    tempU = linear_map(B, U) #overapproximate(concretize(B * U), Zonotope)
    d = m
    while d <= initialTimeStep
        let ϕ::Matrix{Float64} = ReachabilityAnalysis.Exponentiation._exp(A, d, alg)
            tempM = similar(ϕ)

            isInvA = isinvertible(A)

            A_abs = ReachabilityAnalysis.Exponentiation.elementwise_abs(A)
            Φcache = sum(A) == abs(sum(A)) ? ϕ : nothing
            P2A_abs = ReachabilityAnalysis.Exponentiation.Φ₂(A_abs, d, alg, isInvA, nothing)


            dU = overapproximate(d * tempU, Zonotope)
            E_ψ = convert(Zonotope, symmetric_interval_hull(linear_map(P2A_abs, symmetric_interval_hull(linear_map(A, tempU)))))
            #P = minkowski_sum(dU, E_ψ)
            E⁺ = convert(Zonotope, symmetric_interval_hull(linear_map(P2A_abs, symmetric_interval_hull(linear_map(A * A, X0)))))
            lt = minkowski_sum(linear_map(ϕ, X0), dU)
            rt = minkowski_sum(E_ψ, E⁺)
            f = minkowski_sum(lt, rt)
            disc = overapproximate(CH(X0, f), Zonotope)
            Φ₁ = ReachabilityAnalysis.Exponentiation.Φ₁(A, d, alg, isInvA, nothing)
            P = linear_map(Φ₁, tempU)
            discritezationDict[d] = copy(disc)
            inputDiscritezationDict[d] = copy(P)
            phiDict[d] = copy(ϕ)
            d *= 2
        end
    end

    #    discritezationDict, inputDiscritezationDict, phiDict = ReACTDiscretize(A, B, X0, U, m, initialTimeStep, alg, maxOrder, reduceOrder)

    time::Float64 = minimum(interval)
    endtime::Float64 = maximum(interval)

    currentTimeStep = copy(initialTimeStep)

    attemptsRecorder = Integer[]

    V::Zonotope{N,Vector{N},Matrix{N}} = copy(inputDiscritezationDict[initialTimeStep])
    Sρ = zeros(Float64, length(constraint))
    newR::Zonotope{N,Vector{N},Matrix{N}} = discritezationDict[initialTimeStep]
    i = 1

    reachSet = Vector{Zonotope{N,Vector{N},Matrix{N}}}()
    timeStepRecorder = Vector{Float64}()

    #println(discritezationDict[m].center)

    Φ::Matrix{Float64} = diagm(ones(Float64, size(A, 2)))
    tempM = similar(Φ)
    ϕt = similar(Φ)
    newRR = copy(newR)
    #lastV = Zonotope(zero(discritezationDict[m].center), [zero(discritezationDict[m].center)])
    Ub = overapproximate(linear_map(B, U), Zonotope)
    while time < endtime

        attempts = 1
        approveFlag = false

        while !approveFlag
            if currentTimeStep < m
                #currentTimeStep = m
                #println(constraintProjBounds)
                return reachSet, timeStepRecorder, attemptsRecorder
            end

            if changedTimeStep
                newR = discritezationDict[currentTimeStep]
                #V = copy(inputDiscritezationDict[currentTimeStep])
                ϕt = phiDict[currentTimeStep]
                newRR = linear_map(Φ, discritezationDict[currentTimeStep])
                #V = linear_map(Φ, V)
            else
                newRR = linear_map(ϕt, newRR)
                #V = linear_map(ϕt, V)
            end
            V = linear_map(ReachabilityAnalysis.Exponentiation.Φ₁(A, time, alg, false, nothing), Ub)
            changedTimeStep = false
            #hom = map(x -> ρ(x, newRR), constraintProjVectors)
            #println(hom)
            inhom = map(x -> ρ(x, V), constraintProjVectors)

            if all((input + ρ(x, newRR)) < y for (input, x, y) in zip(Sρ, constraintProjVectors, constraintProjBounds))
                approveFlag = true
                push!(reachSet, minkowski_sum(newRR, V)) #minkowski_sum(newRR, V)) # Add homogeneous and input to reachSet
                Sρ = copy(inhom)

                mul!(tempM, Φ, ϕt)
                copy!(Φ, tempM)
            else
                #=if attempts == 1
                    println(Sρ, " ", hom, " ", inhom)
                end=#

                newR = copy(newR)
                currentTimeStep = currentTimeStep / 2
                changedTimeStep = true
                attempts = attempts + 1
            end
        end


        push!(attemptsRecorder, attempts)
        push!(timeStepRecorder, currentTimeStep)

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
    return reachSet, timeStepRecorder, attemptsRecorder
end
