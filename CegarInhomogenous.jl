using LazySets, LinearAlgebra, Printf

include("helperfunctions.jl")
include("reductionMethods.jl")

"""
    cegarInputSystem(A, initialTimeStep, interval, X0::Zonotope, U::Zonotope, constraint::HalfSpace, Digits :: Integer)
Based on "Efficient Computation of Reachable Sets of Linear Time-Invariant Systems with Inputs" - A. Girard, C. Le Guernic, and O. Maler
# System description
Given a LTI system: x' = Ax + Bu(t)
"""
function cegarInputSystem(A, initialTimeStep, interval, X0::Zonotope, U::Zonotope, constraint::HalfSpace, Digits :: Integer)
    h::Vector{Float64} = constraint.a
    f::Float64 = constraint.b 
    ANorm = norm(A, Inf)
    initialTimeStepSize = initialTimeStep
    m = initialTimeStep / 2^(ceil(Integer, log2(initialTimeStep)) + ceil(Integer, -log2(10.0^(-Digits))) - 1)   #Calculate the smallest number larger than 10^-Digits obtained by repeatedly dividing initialTimeStep by 2.
    R = []
    changedTimeStep = true

    discritezationDict = Dict()
    inputDiscritezationDict = Dict()

    time = minimum(interval)
    endtime = maximum(interval)

    currentTimeStep = initialTimeStep
    timeStepRecorder = Float64[]
    attemptsRecorder = Integer[]

    #println("m:$m")
    k = size(U.generators, 2)

    let ϕ = exp(A * m)
        d = m
        P = PCA_reduce(minkowski_sum(U, linear_map(ϕ, U)))
        while d < initialTimeStep
            inputDiscritezationDict[d] = P
            P = minkowski_sum(P, linear_map(ϕ, P))
            ϕ = ϕ * ϕ
            d = d * 2
        end
        inputDiscritezationDict[initialTimeStep] = PCA_reduce(P)
    end


    V = Vector{Zonotope}()
    sizehint!(V, ceil(Integer, endtime / initialTimeStep) + 1)
    push!(V, inputDiscritezationDict[initialTimeStep])
    let ϕ = exp(A * initialTimeStep)
        for id in 2:(ceil(Integer, endtime / initialTimeStep) + 1)
            push!(V, Zonotope(ϕ * V[id - 1].center, ϕ * V[id - 1].generators))
        end
    end


    S = Vector{Zonotope}()
    sizehint!(S, size(V, 1))
    push!(S, V[1])
    for id in 2:size(V, 1)
        push!(S, PCA_reduce(minkowski_sum(S[id - 1], V[id])))
        #push!(S, minkowski_sum(S[id - 1], V[id]))

    end

    newR = nothing
    i = 1
    #=for i in eachindex(S)
        @printf "volume of S[%d]: %.12f\n" i area(S[i])
        #println("volume of V[",i,"]: ", area(V[i]))
    end=#

    while time < endtime
        attempts = 1
        approveFlag = false

        while !approveFlag
            if currentTimeStep < m
                println(timeStepRecorder)
                return (R, timeStepRecorder, attemptsRecorder)
                throw(AssertionError("Error model fails at time $time, constraint is not satisfied after $attempts attempts"))
            end

            if changedTimeStep
                if haskey(discritezationDict, currentTimeStep)
                    newR, ϕt = discritezationDict[currentTimeStep]
                else 
                    newR, ϕt = initialStepNoInput(A, ANorm, currentTimeStep, X0)
                    discritezationDict[currentTimeStep] = (newR, ϕt)
                end

                ϕT = exp(A * time)
                
                newR = linear_map(ϕT, newR)
                changedTimeStep = false
            else 
                _, ϕt = discritezationDict[currentTimeStep]
                newR = linear_map(ϕt, R[i - 1])
            end

            msum::Zonotope = minkowski_sum(newR, S[ceil(Integer, (time + currentTimeStep) / initialTimeStep)])
            if !intersectss(msum.center, genmat(msum), h, f)
                approveFlag = true
                newR = msum
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
        if currentTimeStep < initialTimeStepSize
            currentTimeStep = currentTimeStep * 2
            changedTimeStep = true
        end
    end

    return (R, timeStepRecorder, attemptsRecorder)
end