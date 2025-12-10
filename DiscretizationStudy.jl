# Based on the paper JuliaReach: a Toolbox for Set-Based Reachability
using LazySets, LinearAlgebra, Printf, FastExpm#, Polyhedra#, CDDLib

include("helperfunctions.jl")
include("models.jl")
include("models/heat/heat_load.jl")
include("models/motor/motor_load.jl")
include("models/motor/motor_load.jl")
include("models/building/building_load.jl")
include("models/PDE/pde_load.jl")
#include("CegarFunctions.jl")
include("CegarInhomogenous.jl")


function cegarInputSystemNoOutputOldDisc(A, B, initialTimeStep, interval, X0::Zonotope{N,Vector{N},Matrix{N}}, U::Zonotope, constraint, Digits :: Integer, STRATEGY :: Integer) where {N}
    maxOrder = 50
    reduceToOrder = 1
    finalReduce = maxOrder # We reduce the final P to this

    XG = copy(genmat(X0))
    XDim, p = size(XG)
    m = initialTimeStep / 2^(ceil(Integer, log2(initialTimeStep)) + ceil(Integer, -log2(10.0^(-Digits))) - 1)   #Calculate the smallest number larger than 10^-Digits obtained by repeatedly dividing initialTimeStep by 2.
    changedTimeStep = true

    elems = (ceil(Integer, log2(initialTimeStep)) + ceil(Integer, -log2(10.0^(-Digits))) - 1)
    phiDict = Dict{Float64,  Matrix{Float64}}()
    sizehint!(phiDict, elems)
    discritezationDict = Dict{Float64, Zonotope{N,Vector{N},Matrix{N}}}()
    sizehint!(discritezationDict, elems)
    inputDiscritezationDict = Dict{Float64, Zonotope{N,Vector{N},Matrix{N}}}()
    sizehint!(inputDiscritezationDict, elems)

    time::Float64 = minimum(interval)
    endtime::Float64 = maximum(interval)

    currentTimeStep = copy(initialTimeStep)

    attemptsRecorder = Integer[]

    U = concretize(B * U)

    d = m
    dia :: Matrix{Float64} = diagm(ones(XDim))
    #ϕ :: Matrix{Float64} = nothing
    
    if !(ρ(U.center, U) < norm(U.center)) #Origin is *not* in input
        println("Origin is not in input...")
        # Calculate homogenous part
        invA = inv(Matrix(A))
        û = copy(U.center)
        while d < initialTimeStep
            ϕ :: Matrix{Float64} = fastExpm(A .* d; threshold=eps(Float64), nonzero_tol=eps(Float64))
    
            û = copy(U.center)
            invA = inv(Matrix(A))
            Ut = Zonotope(U.center - û, genmat(U))
            dU = overapproximate(d * Ut, Zonotope)
            #P = minkowski_sum(dU, E_ψ(Ut, d, A))

            P̂ = invA * (ϕ - dia) * û
            lt = concretize(minkowski_sum(convert(Zonotope, ϕ * X0), dU))
            rt = concretize(minkowski_sum(E_ψ(Ut, d, A), E⁺(X0, d, A)))
            PZ = Zonotope(P̂, zeros(Float64, size(U.center, 1), 1))
            f = concretize(minkowski_sum(lt, rt))
            disc = overapproximate(CH(X0, minkowski_sum(f, PZ)), Zonotope)
            

            discritezationDict[d] = copy(disc)
            phiDict[d] = copy(ϕ)
            d = d * 2
        end

        # Calculate for initual time step
        ϕ = fastExpm(A .* d; threshold=eps(Float64), nonzero_tol=eps(Float64))
        û = copy(U.center)
        invA = inv(Matrix(A))
        Ut = Zonotope(U.center - û, genmat(U))
        dU = overapproximate(d * Ut, Zonotope)
        P = minkowski_sum(dU, E_ψ(Ut, d, A))

        P̂ = invA * (ϕ - dia) * û
        lt = concretize(minkowski_sum(convert(Zonotope, ϕ * X0), dU))
        rt = concretize(minkowski_sum(E_ψ(Ut, d, A), E⁺(X0, d, A)))
        PZ = Zonotope(P̂, zeros(Float64, size(U.center, 1), 1))
        f = concretize(minkowski_sum(lt, rt))
        disc = overapproximate(CH(X0, minkowski_sum(f, PZ)), Zonotope)
        

        discritezationDict[d] = copy(disc)
        phiDict[d] = copy(ϕ)
        inputDiscritezationDict[initialTimeStep] = P
    else
        while d < initialTimeStep
            ϕ = fastExpm(A .* d; threshold=eps(Float64), nonzero_tol=eps(Float64))
            dU = overapproximate(d * U, Zonotope)
            lt = concretize(minkowski_sum(convert(Zonotope, ϕ * X0), dU))
            rt = concretize(minkowski_sum(E_ψ(U, d, A), E⁺(X0, d, A)))
            f = concretize(minkowski_sum(lt, rt))   
            disc = overapproximate(CH(X0, f), Zonotope)

            discritezationDict[d] = copy(disc)
            phiDict[d] = copy(ϕ)
            d = d * 2
        end

        ϕ = fastExpm(A .* d; threshold=eps(Float64), nonzero_tol=eps(Float64))
        dU = overapproximate(d * U, Zonotope)
        P = minkowski_sum(dU, E_ψ(U, d, A))
        lt = concretize(minkowski_sum(convert(Zonotope, ϕ * X0), dU))
        rt = concretize(minkowski_sum(E_ψ(U, d, A), E⁺(X0, d, A)))
        f = concretize(minkowski_sum(lt, rt))
        disc = overapproximate(CH(X0, f), Zonotope)

        discritezationDict[d] = copy(disc)
        phiDict[d] = copy(ϕ)
        inputDiscritezationDict[initialTimeStep] = P
    end
    

    println("Dicts done!")

    V :: Zonotope{N,Vector{N},Matrix{N}} = copy(inputDiscritezationDict[initialTimeStep])

    S ::Zonotope{N,Vector{N},Matrix{N}} = copy(inputDiscritezationDict[initialTimeStep])
    newR :: Zonotope{N,Vector{N},Matrix{N}} = discritezationDict[initialTimeStep]
    i = 1
    inputStepsCounter = 0


    Φ :: Matrix{Float64} = diagm(ones(Float64, size(A, 2)))
    tempM = similar(Φ)
    ϕt = similar(Φ)
    newRR = copy(newR)
    initialϕ = phiDict[initialTimeStep]

    while time < endtime

        attempts = 1
        
        approveFlag = false

        if (time + currentTimeStep) > inputStepsCounter * initialTimeStep
            inputStepsCounter += 1
            V = concretize(linear_map(initialϕ, V))
            S = concretize(minkowski_sum(S, V))
            if LazySets.order(S) > maxOrder
                S = reduce_order(S, reduceToOrder)
            end
        end

        while !approveFlag
            if currentTimeStep < m
                println("Error model fails at time $time, constraint is not satisfied after $attempts attempts")
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
            
            nsum :: Zonotope = minkowski_sum(newRR, S)
            if mapreduce(c -> !intersects(nsum, c), &, constraint) 
                approveFlag = true
                mul!(tempM, Φ, ϕt)
                copy!(Φ, tempM) #Her kan man vente med at udregne Phi * phit indtil at man ved hvor mange gange at man vil gange phit på, så at man kan lave et mere effektivt kald på matmul(Phi, phit, phit, ...)
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
                lowest = min(4, i-1)
                window = @view attemptsRecorder[i-lowest:i-1]   
                if all(window .== 1)
                    currentTimeStep = currentTimeStep * 2
                    changedTimeStep = true
                end
            end
        end 
    end

    println("Successful run")

    return true
end

function cegarInputSystemOldDisc(A, B, initialTimeStep, interval, X0::Zonotope{N,Vector{N},Matrix{N}}, U::Zonotope, constraint, Digits :: Integer, STRATEGY :: Integer) where {N}
    maxOrder = 50
    reduceToOrder = 1

    XG = copy(genmat(X0))
    XDim, p = size(XG)
    XC = copy(X0.center)
    m = initialTimeStep / 2^(ceil(Integer, log2(initialTimeStep)) + ceil(Integer, -log2(10.0^(-Digits))) - 1)   #Calculate the smallest number larger than 10^-Digits obtained by repeatedly dividing initialTimeStep by 2.
    R = Zonotope[]
    changedTimeStep = true
    #println("m: ", m)
    #discritezationDict = Dict()#Dict{Float64, Tuple{Zonotope{N,Vector{N},Matrix{N}}, Matrix{Float64}}}()
    phiDict = Dict{Float64,  Matrix{Float64}}()
    discritezationDict = Dict{Float64, Zonotope{N,Vector{N},Matrix{N}}}()
    inputDiscritezationDict = Dict()#Dict{Float64, Zonotope{N,Vector{N},Matrix{N}}}()

    time::Float64 = minimum(interval)
    endtime::Float64 = maximum(interval)

    currentTimeStep = copy(initialTimeStep)
    timeStepRecorder = Float64[]
    attemptsRecorder = Integer[]


    U = concretize(B * U)
    d = m
    dia :: Matrix{Float64} = diagm(ones(XDim))
    
    if !(ρ(U.center, U) < norm(U.center)) #Origin is *not* in input
        println("Origin is not in input...")
        # Calculate homogenous part
        invA = inv(Matrix(A))
        û = copy(U.center)
        while d < initialTimeStep
            ϕ :: Matrix{Float64} = fastExpm(A .* d; threshold=eps(Float64), nonzero_tol=eps(Float64))
    
            û = copy(U.center)
            invA = inv(Matrix(A))
            Ut = Zonotope(U.center - û, genmat(U))
            dU = overapproximate(d * Ut, Zonotope)
            #P = minkowski_sum(dU, E_ψ(Ut, d, A))

            P̂ = invA * (ϕ - dia) * û
            lt = concretize(minkowski_sum(convert(Zonotope, ϕ * X0), dU))
            rt = concretize(minkowski_sum(E_ψ(Ut, d, A), E⁺(X0, d, A)))
            PZ = Zonotope(P̂, zeros(Float64, size(U.center, 1), 1))
            f = concretize(minkowski_sum(lt, rt))
            disc = overapproximate(CH(X0, minkowski_sum(f, PZ)), Zonotope)
            

            discritezationDict[d] = copy(disc)
            phiDict[d] = copy(ϕ)
            d = d * 2
        end

        # Calculate for initual time step
        ϕ = fastExpm(A .* d; threshold=eps(Float64), nonzero_tol=eps(Float64))
        û = copy(U.center)
        invA = inv(Matrix(A))
        Ut = Zonotope(U.center - û, genmat(U))
        dU = overapproximate(d * Ut, Zonotope)
        P = minkowski_sum(dU, E_ψ(Ut, d, A))

        P̂ = invA * (ϕ - dia) * û
        lt = concretize(minkowski_sum(convert(Zonotope, ϕ * X0), dU))
        rt = concretize(minkowski_sum(E_ψ(Ut, d, A), E⁺(X0, d, A)))
        PZ = Zonotope(P̂, zeros(Float64, size(U.center, 1), 1))
        f = concretize(minkowski_sum(lt, rt))
        disc = overapproximate(CH(X0, minkowski_sum(f, PZ)), Zonotope)
        

        discritezationDict[d] = copy(disc)
        phiDict[d] = copy(ϕ)
        inputDiscritezationDict[initialTimeStep] = P
    else
        while d < initialTimeStep
            ϕ = fastExpm(A .* d; threshold=eps(Float64), nonzero_tol=eps(Float64))
            dU = overapproximate(d * U, Zonotope)
            lt = concretize(minkowski_sum(convert(Zonotope, ϕ * X0), dU))
            rt = concretize(minkowski_sum(E_ψ(U, d, A), E⁺(X0, d, A)))
            f = concretize(minkowski_sum(lt, rt))   
            disc = overapproximate(CH(X0, f), Zonotope)

            discritezationDict[d] = copy(disc)
            phiDict[d] = copy(ϕ)
            d = d * 2
        end

        ϕ = fastExpm(A .* d; threshold=eps(Float64), nonzero_tol=eps(Float64))
        dU = overapproximate(d * U, Zonotope)
        P = minkowski_sum(dU, E_ψ(U, d, A))
        lt = concretize(minkowski_sum(convert(Zonotope, ϕ * X0), dU))
        rt = concretize(minkowski_sum(E_ψ(U, d, A), E⁺(X0, d, A)))
        f = concretize(minkowski_sum(lt, rt))
        disc = overapproximate(CH(X0, f), Zonotope)

        discritezationDict[d] = copy(disc)
        phiDict[d] = copy(ϕ)
        inputDiscritezationDict[initialTimeStep] = P
    end
    

    println("Dicts done!")

    STEPS = ceil(Integer, endtime / initialTimeStep) + 1
    V :: Zonotope{N,Vector{N},Matrix{N}} = inputDiscritezationDict[initialTimeStep]

    S ::Zonotope{N,Vector{N},Matrix{N}} = inputDiscritezationDict[initialTimeStep]
    newR :: Zonotope{N,Vector{N},Matrix{N}} = discritezationDict[initialTimeStep]
    i = 1
    inputStepsCounter = 0


    Φ :: Matrix{Float64} = diagm(ones(Float64, size(A, 2)))
    tempM = similar(Φ)
    ϕt = similar(Φ)
    tempXG = Matrix{Float64}(undef, 1, 150)
    RG = similar(genmat(newR))
    RC = similar(newR.center)
    newRR = copy(newR)
    einmal :: Matrix{Float64} = diagm(ones(Float64, size(A, 2)))
    initialϕ :: Matrix{Float64} = phiDict[initialTimeStep]

    tubes = []

    while time < endtime

        attempts = 1
        
        approveFlag = false

        if (time + currentTimeStep) > inputStepsCounter * initialTimeStep
            inputStepsCounter += 1
            V = linear_map(initialϕ, V)            
            S = minkowski_sum(S, V)
            if LazySets.order(S) > maxOrder 
                S = reduce_order(S, reduceToOrder)
            end
        end

        while !approveFlag
            if currentTimeStep < m
                println("Error model fails at time $time, constraint is not satisfied after $attempts attempts")
                return (tubes, timeStepRecorder, attemptsRecorder)
            end

            if changedTimeStep
                newR= discritezationDict[currentTimeStep]
                ϕt = phiDict[currentTimeStep] 
                newRR = linear_map(Φ, newR)
            else 
                newRR = linear_map(ϕt, newRR)#linear_map!(newRR, ϕt, R[i-1]) #smallStep(newRR, ϕt, RC, RG)
            end

            changedTimeStep = false
            
            nsum :: Zonotope = minkowski_sum(newRR, S)
            if mapreduce(c -> !intersects(nsum, c), &, constraint) #intersectss(msum.center, genmat(msum), h, f, tempXG)
                approveFlag = true
                push!(R, copy(newRR))
                push!(tubes, copy(nsum))
                mul!(tempM, Φ, ϕt)
                copy!(Φ, tempM) #Her kan man vente med at udregne Phi * phit indtil at man ved hvor mange gange at man vil gange phit på, så at man kan lave et mere effektivt kald på matmul(Phi, phit, phit, ...)
            else 
                if i > 1
                    R[i - 1] = copy(R[i - 1])
                end
                currentTimeStep = currentTimeStep / 2 
                changedTimeStep = true
                attempts = attempts + 1
            end
        end

        push!(timeStepRecorder, currentTimeStep)
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
                lowest = min(4, i-1)
                window = @view attemptsRecorder[i-lowest:i-1]   
                if all(window .== 1)
                    currentTimeStep = currentTimeStep * 2
                    changedTimeStep = true
                end
            end
        end 
    end

    return (tubes, timeStepRecorder, attemptsRecorder)
end