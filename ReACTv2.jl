using LazySets, LinearAlgebra

function ReACTWithSupport(A, B, initialTimeStep, interval, X0::Zonotope{N,Vector{N},Matrix{N}}, U::Zonotope, constraint, δ⁻, STRATEGY::Integer, alg::ReachabilityAnalysis.Exponentiation.AbstractExpAlg=ReachabilityAnalysis.Exponentiation.BaseExp, maxOrder::Int=5, reduceOrder::Int=5) where {N}
    m = δ⁻
    changedTimeStep = true
    elems = (ceil(Integer, log2(initialTimeStep)) + ceil(Integer, -log2(m)) - 1)
    phiDict = Dict{Float64,Matrix{Float64}}()
    sizehint!(phiDict, elems)
    discritezationDict = Dict{Float64,Zonotope{N,Vector{N},Matrix{N}}}()
    sizehint!(discritezationDict, elems)
    inputDiscritezationDict = Dict{Float64,Zonotope{N,Vector{N},Matrix{N}}}()
    sizehint!(inputDiscritezationDict, elems)

    constraintProjVectors = map(x -> x.a, constraint)
    oldConstraintProjVectors = copy(constraintProjVectors)
    constraintProjBounds = map(x -> x.b, constraint)



    time::Float64 = minimum(interval)
    endtime::Float64 = maximum(interval)

    currentTimeStep = copy(initialTimeStep)

    attemptsRecorder = Integer[]


    discritezationDict, inputDiscritezationDict, phiDict = ReACTDiscretize(A, B, X0, U, m, initialTimeStep, alg, maxOrder, reduceOrder)

    for key in keys(phiDict)
        phiDict[key] = permutedims(phiDict[key])
    end

    Sρ = zeros(Float64, length(constraint))
    newRR = discritezationDict[initialTimeStep]
    i = 1

    ϕt::Matrix{Float64} = diagm(ones(Float64, size(A, 2)))


    while time < endtime

        attempts = 1
        approveFlag = false

        while !approveFlag
            if currentTimeStep < m
                return false
            end

            if changedTimeStep
                newRR = discritezationDict[currentTimeStep]
                ϕt = phiDict[currentTimeStep]
            end
            constraintProjVectors = map(x -> ϕt * x, oldConstraintProjVectors)

            changedTimeStep = false
            if all((input + ρ(x, newRR)) <= y for (input, x, y) in zip(Sρ, constraintProjVectors, constraintProjBounds))
                Sρ += map(x -> ρ(x, inputDiscritezationDict[currentTimeStep]), oldConstraintProjVectors)
                approveFlag = true
                oldConstraintProjVectors = constraintProjVectors
            else
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

# No input version (U::Nothing)
function ReACTWithSupport(A, B, initialTimeStep, interval, X0::Zonotope{N,Vector{N},Matrix{N}}, U::Nothing, constraint, δ⁻, STRATEGY::Integer, alg::ReachabilityAnalysis.Exponentiation.AbstractExpAlg=ReachabilityAnalysis.Exponentiation.BaseExp, maxOrder::Int=5, reduceOrder::Int=5) where {N}
    m = copy(δ⁻)
    changedTimeStep = true
    phiDict = Dict{Float64,Matrix{Float64}}()
    discritezationDict = Dict{Float64,Zonotope{N,Vector{N},Matrix{N}}}()


    constraintProjVectors = map(x -> x.a, constraint)
    oldConstraintProjVectors = copy(constraintProjVectors)
    constraintProjBounds = map(x -> x.b, constraint)



    time::Float64 = minimum(interval)
    endtime::Float64 = maximum(interval)

    currentTimeStep = copy(initialTimeStep)

    attemptsRecorder = Integer[]


    discritezationDict, _, phiDict = ReACTDiscretize(A, B, X0, U, m, initialTimeStep, alg, maxOrder, reduceOrder)

    for key in keys(phiDict)
        phiDict[key] = permutedims(phiDict[key])
    end

    newR::Zonotope{N,Vector{N},Matrix{N}} = discritezationDict[initialTimeStep]
    i = 1

    ϕt::Matrix{Float64} = diagm(ones(Float64, size(A, 2)))
    newRR = copy(newR)


    while time < endtime
        attempts = 1
        approveFlag = false
        while !approveFlag
            if currentTimeStep < m
                return false
            end

            if changedTimeStep
                newRR = discritezationDict[currentTimeStep]
                ϕt = phiDict[currentTimeStep]
            end
            constraintProjVectors = map(x -> ϕt * x, oldConstraintProjVectors)
            changedTimeStep = false

            if all(ρ(x, newRR) <= y for (x, y) in zip(constraintProjVectors, constraintProjBounds))
                approveFlag = true
                oldConstraintProjVectors = constraintProjVectors
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