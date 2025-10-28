using LazySets, LinearAlgebra

include("helperfunctions.jl")


"""
    cegarInputSystem(A, initialTimeStep, interval, X0::Zonotope, U::Zonotope, constraint::HalfSpace, Digits :: Integer)
Based on "Efficient Computation of Reachable Sets of Linear Time-Invariant Systems with Inputs" - A. Girard, C. Le Guernic, and O. Maler
# System description
Given a LTI system: x' = Ax + Bu(t)
"""
function cegarInputSystem(A, initialTimeStep, interval, X0::Zonotope, U::Zonotope, constraint::HalfSpace, Digits :: Integer)
    ANorm = norm(A, Inf)
    initialTimeStepSize = initialTimeStep
    m = initialTimeStep / 2^(log2(initialTimeStep)+ceil(-log2(10.0^(-Digits))))   #Calculate the smallest number larger than 10^-Digits obtained by repeatedly dividing initialTimeStep by 2.
    R = []
    Ω = []
    changedTimeStep = true

    previouslyCalculatedDict = Dict()
    previousInputValuesDict = Dict()

    time = minimum(interval)
    endtime = maximum(interval)

    currentTimeStep = initialTimeStep
    timeStepRecorder = []
    attemptsRecorder = []

    S = Vector{Zonotope}()
    sizehint!(S, ceil(endtime / m))
    push!(S, U)
    let ϕ = exp(A*m)
        for id in 2:(ceil(endtime / m))
            push!(S,ϕ * S[end])
        end
    end

    newR = nothing
    ϕ = nothing
    i = 1


    while time < endtime
        attempts = 1
        approveFlag = false
        prevTime = 0.

        while !approveFlag
            if currentTimeStep == 0
                throw(AssertionError("Error model fails at time $time, constraint is not satisfied"))
            end

            if changedTimeStep
                if haskey(previouslyCalculatedDict, currentTimeStep)
                    newR, ϕ = previouslyCalculatedDict[currentTimeStep]
                    prevTime = previousInputValuesDict[currentTimeStep]
                else 
                    newR, ϕ = initialStepNoInput(A, ANorm, currentTimeStep, X0)
                    previouslyCalculatedDict[currentTimeStep] = (newR, ϕ)
                    previousInputValuesDict[currentTimeStep] = time
                end
                # So we have calculated the discritezation with the current timestep and then we teleport it.
                # We should probably store the teleportation matrix and just perform the teleportation without calling the function.
                newR = forwardTimeNoInput(A, newR, time)
                changedTimeStep = false
            else 
                newR = linear_map(ϕ, R[i - 1])
            end

            msum = minkowski_sum(newR, S[ceil(Integer, (time + currentTimeStep) / m)])
            if !intersects(constraint, msum)
                approveFlag = true
            else 
                currentTimeStep = currentTimeStep / 2
                changedTimeStep = true
                attempts = attempts + 1
            end
        end

        push!(R, newR)
    
        push!(timeStepRecorder, currentTimeStep)
        push!(attemptsRecorder, attempts)
        i = i + 1
        time = round(time + currentTimeStep, digits=Digits)

        # Reset / apply strategy
        if currentTimeStep == initialTimeStepSize
            currentTimeStep = initialTimeStep
            changedTimeStep = true
        end
    end

    return (R, timeStepRecorder, attemptsRecorder)
end