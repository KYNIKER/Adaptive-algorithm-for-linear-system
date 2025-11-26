using LazySets, LinearAlgebra, Printf, FastExpm#, Polyhedra#, CDDLib

include("helperfunctions.jl")
include("reductionMethods.jl")

"""
    cegarInputSystem(A, initialTimeStep, interval, X0::Zonotope, U::Zonotope, constraint::HalfSpace, Digits :: Integer)
Based on "Efficient Computation of Reachable Sets of Linear Time-Invariant Systems with Inputs" - A. Girard, C. Le Guernic, and O. Maler
# System description
Given a LTI system: x' = Ax + Bu(t)
"""
function cegarInputSystem(A, initialTimeStep, interval, X0::Zonotope{N,Vector{N},Matrix{N}}, U::Zonotope, constraint, Digits :: Integer) where {N}
    stepsBeforeReduce = 4
    
    ANorm = norm(A, Inf)
    XNorm = norm(X0, Inf)::Float64
    XG = copy(genmat(X0))
    XDim, p = size(XG)
    XC = copy(X0.center)
    m = initialTimeStep / 2^(ceil(Integer, log2(initialTimeStep)) + ceil(Integer, -log2(10.0^(-Digits))) - 1)   #Calculate the smallest number larger than 10^-Digits obtained by repeatedly dividing initialTimeStep by 2.
    R = Zonotope[]
    changedTimeStep = true
    #println("m: ", m)
    discritezationDict = Dict()#Dict{Float64, Tuple{Zonotope{N,Vector{N},Matrix{N}}, Matrix{Float64}}}()
    inputDiscritezationDict = Dict()#Dict{Float64, Zonotope{N,Vector{N},Matrix{N}}}()

    time::Float64 = minimum(interval)
    endtime::Float64 = maximum(interval)

    currentTimeStep = copy(initialTimeStep)
    timeStepRecorder = Float64[]
    attemptsRecorder = Integer[]

    k = size(U.generators, 2)

    

    let ϕ :: Matrix{Float64} = fastExpm(A .* m; threshold=eps(Float64), nonzero_tol=eps(Float64))
        tempM = similar(ϕ)
        d = m
        dia :: Matrix{Float64} = diagm(ones(XDim))
        i = 1
        
        if !(zeros(size(U.center)) ∈ U)
            û = copy(U.center)
            #U.center = U.center * 0
            invA = inv(Matrix(A))
            #P̂ = invA * (ϕ - dia) * û
            Ut = Zonotope(U.center - û, genmat(U))

            dU = box_approximation_symmetric(d * Ut)
            #println(typeof(dU))
            #println(typeof(E_ψ(U, d, A)))
            P = minkowski_sum(dU, E_ψ(Ut, d, A))#minkowski_sum(U, linear_map(ϕ, U))
            
            #disc :: Zonotope{N,Vector{N},Matrix{N}} = copy(X0)
            while d < initialTimeStep
                #println("d: ", d)
                inputDiscritezationDict[d] = P
                P = minkowski_sum(P, linear_map(ϕ, P))
                if i % stepsBeforeReduce == 0
                    P = PCA_reduce(P)
                end
                i += 1
                #=α = (exp(ANorm * d) - 1 - d * ANorm) / XNorm
                ϕp = (dia+ϕ)/2
                ϕm = (dia-ϕ)/2
                gens::Matrix{Float64} = hcat(ϕp * XG, ϕm * XC, ϕm * XG, α*dia) #hcat(ϕp * XG, ϕm * XC, ϕm * XG)

                disc = Zonotope(ϕp*XC, gens)#minkowski_sum(Zonotope(ϕp*XC, gens), Zonotope(zeros(XDim), α*dia))
                =#
                #disc = minkowski_sum(convex_hull(X0, ϕ * minkowski_sum(X0, box_approximation_symmetric(d * U))),minkowski_sum(E_ψ(U, d, A), E⁺(X0, d, A)))
                P̂ = invA * (ϕ - dia) * û
                #println(P̂)
                lt = concretize(minkowski_sum(convert(Zonotope, ϕ * X0), box_approximation_symmetric(d * Ut)))
                rt = concretize(minkowski_sum(E_ψ(Ut, d, A), E⁺(X0, d, A)))

                #println(typeof(lt), " vs ", typeof(rt))
                PZ = Zonotope(P̂, zeros(Float64, size(U.center, 1), 1))
                f = concretize(minkowski_sum(lt, rt))
                
                disc = overapproximate(CH(X0, minkowski_sum(f, PZ)), Zonotope)
                discritezationDict[d] = (copy(disc), copy(ϕ))#(copy(Zonotope(disc.center - PZ.center, genmat(disc))), copy(ϕ))

                mul!(tempM, ϕ , ϕ)
                copy!(ϕ, tempM)
                d = d * 2
            end

            #=α = (exp(ANorm * d) - 1 - d * ANorm) / XNorm
            ϕp = (dia+ϕ)/2
            ϕm = (dia-ϕ)/2
            gens = hcat(ϕp * XG, ϕm * XC, ϕm * XG, α*dia) #hcat(ϕp * XG, ϕm * XC, ϕm * XG)

            disc = Zonotope(ϕp*XC, gens)#minkowski_sum(Zonotope(ϕp*XC, gens), Zonotope(zeros(XDim), α*dia))
            =#
            P̂ = invA * (ϕ - dia) * û
            #println(P̂)
            lt = concretize(minkowski_sum(convert(Zonotope, ϕ * X0), box_approximation_symmetric(d * Ut)))
            rt = concretize(minkowski_sum(E_ψ(Ut, d, A), E⁺(X0, d, A)))

            #println(typeof(lt), " vs ", typeof(rt))
            PZ = Zonotope(P̂, zeros(Float64, size(U.center, 1), 1))
            f = concretize(minkowski_sum(lt, rt))
            
            disc = overapproximate(CH(X0, minkowski_sum(f, PZ)), Zonotope)
            discritezationDict[d] = (copy(disc), copy(ϕ))

            #P = minkowski_sum(U, linear_map(ϕ, U))
            #P = minkowski_sum(dU, E_ψ(U, d, A))
            inputDiscritezationDict[initialTimeStep] = P
        else
        
            dU = box_approximation_symmetric(d * U)
            #println(typeof(dU))
            #println(typeof(E_ψ(U, d, A)))
            P = minkowski_sum(dU, E_ψ(U, d, A))#minkowski_sum(U, linear_map(ϕ, U))
            #disc :: Zonotope{N,Vector{N},Matrix{N}} = copy(X0)
            while d < initialTimeStep
                println("d: ", d)
                inputDiscritezationDict[d] = P
                P = minkowski_sum(P, linear_map(ϕ, P))
                if i % stepsBeforeReduce == 0
                    P = PCA_reduce(P)
                end
                i += 1
                #=α = (exp(ANorm * d) - 1 - d * ANorm) / XNorm
                ϕp = (dia+ϕ)/2
                ϕm = (dia-ϕ)/2
                gens::Matrix{Float64} = hcat(ϕp * XG, ϕm * XC, ϕm * XG, α*dia) #hcat(ϕp * XG, ϕm * XC, ϕm * XG)

                disc = Zonotope(ϕp*XC, gens)#minkowski_sum(Zonotope(ϕp*XC, gens), Zonotope(zeros(XDim), α*dia))
                =#
                #disc = minkowski_sum(convex_hull(X0, ϕ * minkowski_sum(X0, box_approximation_symmetric(d * U))),minkowski_sum(E_ψ(U, d, A), E⁺(X0, d, A)))
                lt = concretize(minkowski_sum(convert(Zonotope, ϕ * X0), box_approximation_symmetric(d * U)))
                rt = concretize(minkowski_sum(E_ψ(U, d, A), E⁺(X0, d, A)))


                #println(typeof(lt), " vs ", typeof(rt))
                f = concretize(minkowski_sum(lt, rt))
                disc = overapproximate(CH(X0, f), Zonotope)
                discritezationDict[d] = (copy(disc), copy(ϕ))

                mul!(tempM, ϕ , ϕ)
                copy!(ϕ, tempM)
                d = d * 2
            end

            #=α = (exp(ANorm * d) - 1 - d * ANorm) / XNorm
            ϕp = (dia+ϕ)/2
            ϕm = (dia-ϕ)/2
            gens = hcat(ϕp * XG, ϕm * XC, ϕm * XG, α*dia) #hcat(ϕp * XG, ϕm * XC, ϕm * XG)

            disc = Zonotope(ϕp*XC, gens)#minkowski_sum(Zonotope(ϕp*XC, gens), Zonotope(zeros(XDim), α*dia))
            =#
            lt = concretize(minkowski_sum(convert(Zonotope, ϕ * X0), box_approximation_symmetric(d * U)))
            rt = concretize(minkowski_sum(E_ψ(U, d, A), E⁺(X0, d, A)))


            #println(typeof(lt), " vs ", typeof(rt))
            f = concretize(minkowski_sum(lt, rt))
            disc = overapproximate(CH(X0, f), Zonotope)
            discritezationDict[d] = (copy(disc), copy(ϕ))
            #P = minkowski_sum(U, linear_map(ϕ, U))
            #P = minkowski_sum(dU, E_ψ(U, d, A))
            inputDiscritezationDict[initialTimeStep] = P
        end
    end

    println("Dicts done!")

    STEPS = ceil(Integer, endtime / initialTimeStep) + 1
    V :: Zonotope{N,Vector{N},Matrix{N}} = inputDiscritezationDict[initialTimeStep]

    S ::Zonotope{N,Vector{N},Matrix{N}} = inputDiscritezationDict[initialTimeStep]
    newR :: Zonotope{N,Vector{N},Matrix{N}}, _ = discritezationDict[initialTimeStep]
    i = 1


    Φ :: Matrix{Float64} = diagm(ones(Float64, size(A, 2)))
    tempM = similar(Φ)
    ϕt = similar(Φ)
    tempXG = Matrix{Float64}(undef, 1, 150)
    RG = similar(genmat(newR))
    RC = similar(newR.center)
    #ϕT :: Matrix{Float64} = Φ
    newRR = copy(newR)
    einmal :: Matrix{Float64} = diagm(ones(Float64, size(A, 2)))
    (_, initialϕ :: Matrix{Float64}) = discritezationDict[initialTimeStep]

    tubes = []

    while time < endtime

        attempts = 1
        
        approveFlag = false

        if ceil(Integer, (time + currentTimeStep) / initialTimeStep) > ceil(Integer, time / initialTimeStep)
            V = linear_map(initialϕ, V)            
            S = minkowski_sum(S, V)
            if i % stepsBeforeReduce == 0
                S = PCA_reduce(S)
            end
        end

        while !approveFlag
            if currentTimeStep < m
                println("Error model fails at time $time, constraint is not satisfied after $attempts attempts")
                return (tubes, timeStepRecorder, attemptsRecorder)
            end

            if changedTimeStep
                newR, ϕt = discritezationDict[currentTimeStep]
                newRR = linear_map(Φ, newR)
                #msum = newRR
            else 
                #newR, ϕt = discritezationDict[currentTimeStep]
                #newRR = R[i - 1] #linear_map!(newRR, einmal, R[i-1])#newRR = copy(R[i - 1])
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

        #push!(R, msum)
    
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



function cegarInputSystemNoOutput(A, initialTimeStep, interval, X0::Zonotope{N,Vector{N},Matrix{N}}, U::Zonotope, constraint, Digits :: Integer, STRATEGY :: Integer) where {N}
    stepsBeforeReduce = 5
    XG = copy(genmat(X0))
    XDim, p = size(XG)
    m = initialTimeStep / 2^(ceil(Integer, log2(initialTimeStep)) + ceil(Integer, -log2(10.0^(-Digits))) - 1)   #Calculate the smallest number larger than 10^-Digits obtained by repeatedly dividing initialTimeStep by 2.
    changedTimeStep = true
    discritezationDict = Dict()#Dict{Float64, Tuple{Zonotope{N,Vector{N},Matrix{N}}, Matrix{Float64}}}()
    inputDiscritezationDict = Dict()#Dict{Float64, Zonotope{N,Vector{N},Matrix{N}}}()

    time::Float64 = minimum(interval)
    endtime::Float64 = maximum(interval)

    currentTimeStep = copy(initialTimeStep)

    attemptsRecorder = Integer[]

    let ϕ :: Matrix{Float64} = fastExpm(A .* m; threshold=eps(Float64), nonzero_tol=eps(Float64))
        tempM = similar(ϕ)
        d = m
        dia :: Matrix{Float64} = diagm(ones(XDim))
        i = 1
        
        if !(zeros(size(U.center)) ∈ U) #Origin is *not* in input
            û = copy(U.center)
            invA = inv(Matrix(A))
            Ut = Zonotope(U.center - û, genmat(U))
            dU = box_approximation_symmetric(d * Ut)
            P = minkowski_sum(dU, E_ψ(Ut, d, A))

            while d < initialTimeStep
                inputDiscritezationDict[d] = P
                P = minkowski_sum(P, linear_map(ϕ, P))
                if i % stepsBeforeReduce == 0
                    P = PCA_reduce(P)
                end
                i += 1
                P̂ = invA * (ϕ - dia) * û
                lt = concretize(minkowski_sum(convert(Zonotope, ϕ * X0), box_approximation_symmetric(d * Ut)))
                rt = concretize(minkowski_sum(E_ψ(Ut, d, A), E⁺(X0, d, A)))
                PZ = Zonotope(P̂, zeros(Float64, size(U.center, 1), 1))
                f = concretize(minkowski_sum(lt, rt))
                disc = overapproximate(CH(X0, minkowski_sum(f, PZ)), Zonotope)
                discritezationDict[d] = (copy(disc), copy(ϕ))

                mul!(tempM, ϕ , ϕ)
                copy!(ϕ, tempM)
                d = d * 2
            end
            P̂ = invA * (ϕ - dia) * û

            lt = concretize(minkowski_sum(convert(Zonotope, ϕ * X0), box_approximation_symmetric(d * Ut)))
            rt = concretize(minkowski_sum(E_ψ(Ut, d, A), E⁺(X0, d, A)))
            PZ = Zonotope(P̂, zeros(Float64, size(U.center, 1), 1))
            f = concretize(minkowski_sum(lt, rt))
            
            disc = overapproximate(CH(X0, minkowski_sum(f, PZ)), Zonotope)
            discritezationDict[d] = (copy(disc), copy(ϕ))
            inputDiscritezationDict[initialTimeStep] = P
        else
            dU = box_approximation_symmetric(d * U)
            P = minkowski_sum(dU, E_ψ(U, d, A))
            while d < initialTimeStep
                inputDiscritezationDict[d] = P
                P = minkowski_sum(P, linear_map(ϕ, P))
                if i % stepsBeforeReduce == 0
                    P = PCA_reduce(P)
                end
                i += 1
                lt = concretize(minkowski_sum(convert(Zonotope, ϕ * X0), box_approximation_symmetric(d * U)))
                rt = concretize(minkowski_sum(E_ψ(U, d, A), E⁺(X0, d, A)))
                f = concretize(minkowski_sum(lt, rt))
                disc = overapproximate(CH(X0, f), Zonotope)
                discritezationDict[d] = (copy(disc), copy(ϕ))

                mul!(tempM, ϕ , ϕ)
                copy!(ϕ, tempM)
                d = d * 2
            end

            lt = concretize(minkowski_sum(convert(Zonotope, ϕ * X0), box_approximation_symmetric(d * U)))
            rt = concretize(minkowski_sum(E_ψ(U, d, A), E⁺(X0, d, A)))
            disc = overapproximate( CH(X0, concretize(minkowski_sum(lt, rt))), Zonotope)
            discritezationDict[d] = (copy(disc), copy(ϕ))
            inputDiscritezationDict[initialTimeStep] = P
        end
    end

    println("Dicts done!")

    #STEPS = ceil(Integer, endtime / initialTimeStep) + 1
    V :: Zonotope{N,Vector{N},Matrix{N}} = inputDiscritezationDict[initialTimeStep]

    S ::Zonotope{N,Vector{N},Matrix{N}} = inputDiscritezationDict[initialTimeStep]
    newR :: Zonotope{N,Vector{N},Matrix{N}}, _ = discritezationDict[initialTimeStep]
    i = 1


    Φ :: Matrix{Float64} = diagm(ones(Float64, size(A, 2)))
    tempM = similar(Φ)
    ϕt = similar(Φ)
    newRR = copy(newR)
    (_, initialϕ) = discritezationDict[initialTimeStep]

    while time < endtime

        attempts = 1
        
        approveFlag = false

        if ceil(Integer, (time + currentTimeStep) / initialTimeStep) > ceil(Integer, time / initialTimeStep)
            V = linear_map(initialϕ, V)            
            S = minkowski_sum(S, V)
            if i % stepsBeforeReduce == 0
                S = PCA_reduce(S)
            end
        end

        while !approveFlag
            if currentTimeStep < m
                println("Error model fails at time $time, constraint is not satisfied after $attempts attempts")
                return false
            end

            if changedTimeStep
                newR, ϕt = discritezationDict[currentTimeStep]
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
                currentTimeStep = currentTimeStep / 2 
                changedTimeStep = true
                attempts = attempts + 1
            end
        end

        #push!(R, msum)
    
        #push!(timeStepRecorder, currentTimeStep)
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



function OneTimeStepSystem(A, timestep::Float64, interval, X0::Zonotope{N,Vector{N},Matrix{N}}, U::Zonotope, constraint, Digits :: Integer) where {N}
    ANorm = norm(A, Inf)
    XNorm = norm(X0, Inf)::Float64
    XG = copy(genmat(X0))
    XDim, p = size(XG)
    XC = copy(X0.center)
    R = Zonotope[]
    m = timestep / 2^(ceil(Integer, log2(timestep)) + ceil(Integer, -log2(10.0^(-Digits))) - 1)   #Calculate the smallest number larger than 10^-Digits obtained by repeatedly dividing initialTimeStep by 2.
    #println("m: ", m)

    time::Float64 = minimum(interval)
    endtime::Float64 = maximum(interval)

    #finalDisc :: Zonotope{N,Vector{N},Matrix{N}} = nothing
    disc :: Zonotope{N,Vector{N},Matrix{N}} = copy(X0)
    initialϕ= nothing

    P = nothing

    # Calculate based on the smallest number, for some digits
    let ϕ :: Matrix{Float64} = fastExpm(A .* m; threshold=eps(Float64), nonzero_tol=eps(Float64))
        tempM = similar(ϕ)
        d = m
        P = minkowski_sum(U, linear_map(ϕ, U))
        # disc :: Zonotope{N,Vector{N},Matrix{N}} = copy(X0)
        dia :: Matrix{Float64} = diagm(ones(XDim))

        while d < timestep
            println("d: ", d)
            P = PCA_reduce(minkowski_sum(P, linear_map(ϕ, P)))
            mul!(tempM, ϕ , ϕ)
            copy!(ϕ, tempM)
            d = d * 2
        end

        α = (exp(ANorm * d) - 1 - d * ANorm) / XNorm
        ϕp = (dia+ϕ)/2
        ϕm = (dia-ϕ)/2
        gens = hcat(ϕp * XG, ϕm * XC, ϕm * XG, α*dia) #hcat(ϕp * XG, ϕm * XC, ϕm * XG)

        disc = Zonotope(ϕp*XC, gens)#minkowski_sum(Zonotope(ϕp*XC, gens), Zonotope(zeros(XDim), α*dia))
        initialϕ = ϕ
        #P = minkowski_sum(U, linear_map(ϕ, U))
        P = PCA_reduce(P)
    end

    println("Dicts done!")

    #STEPS = ceil(Integer, endtime / initialTimeStep) + 1
    V :: Zonotope{N,Vector{N},Matrix{N}} = P

    S ::Zonotope{N,Vector{N},Matrix{N}} = P
    newR :: Zonotope{N,Vector{N},Matrix{N}} = disc
    #i = 1


    Φ :: Matrix{Float64} = diagm(ones(Float64, size(A, 2)))
    tempM = similar(Φ)
    #ϕt = similar(Φ)
    #tempXG = Matrix{Float64}(undef, 1, 150)
    # RG = similar(genmat(newR))
    # RC = similar(newR.center)
    #ϕT :: Matrix{Float64} = Φ
    #newRR = copy(newR)
    einmal :: Matrix{Float64} = diagm(ones(Float64, size(A, 2)))
    #(_, initialϕ) = discritezationDict[initialTimeStep]


    tubes = []

    # Loop
    while time < endtime
        V = linear_map(initialϕ, V)            
        S = PCA_reduce(minkowski_sum(S, V)) :: Zonotope{Float64,Vector{Float64},Matrix{Float64}}
        

        newR = linear_map(initialϕ, newR)#linear_map!(newRR, ϕt, R[i-1]) #smallStep(newRR, ϕt, RC, RG)
            
            
        nsum :: Zonotope = minkowski_sum(newR, S)
        if mapreduce(c -> !intersects(nsum, c), &, constraint) #intersectss(msum.center, genmat(msum), h, f, tempXG)
            push!(R, copy(newR))
            push!(tubes, copy(nsum))
            mul!(tempM, Φ, initialϕ)
            copy!(Φ, tempM) #Her kan man vente med at udregne Phi * phit indtil at man ved hvor mange gange at man vil gange phit på, så at man kan lave et mere effektivt kald på matmul(Phi, phit, phit, ...)
        else # If we hit a constraint, we fail
            println("Error model fails at time $time, constraint is not satisfied with timestep $timestep")
            amountOfElements = length(tubes)
            timeStepRecorder = [timestep for i in 1:amountOfElements]
            attemptsRecorder = [1 for i in 1:amountOfElements]
            return (tubes, timeStepRecorder, attemptsRecorder)
        end
        
        # Update counters
        #i = i + 1
        time = time + timestep
    end

    amountOfElements = length(tubes)
    timeStepRecorder = [timestep for i in 1:amountOfElements]
    attemptsRecorder = [1 for i in 1:amountOfElements]
    return (tubes, timeStepRecorder, attemptsRecorder)
end

#=function cegarInputSystemPreAlloc(A, initialTimeStep, interval, X0::Zonotope{N,Vector{N},Matrix{N}}, U::Zonotope, constraint::HalfSpace, Digits :: Integer) where {N}
    h::Vector{Float64} = constraint.a
    f::Float64 = constraint.b 
    ANorm = norm(A, Inf)
    XNorm = norm(X0, Inf)
    XG = genmat(X0)
    XDim, p = size(XG)
    XC = X0.center
    initialTimeStepSize = initialTimeStep
    m = initialTimeStep / 2^(ceil(Integer, log2(initialTimeStep)) + ceil(Integer, -log2(10.0^(-Digits))) - 1)   #Calculate the smallest number larger than 10^-Digits obtained by repeatedly dividing initialTimeStep by 2.
    R = Zonotope[]
    changedTimeStep = true

    discritezationDict = Dict{Float64, Tuple{Zonotope, Matrix{Float64}}}()
    inputDiscritezationDict = Dict{Float64, Zonotope}()

    time = minimum(interval)
    endtime = maximum(interval)

    currentTimeStep = initialTimeStep
    timeStepRecorder = Float64[]
    attemptsRecorder = Integer[]

    #println("m:$m")
    k = size(U.generators, 2)

    let ϕ :: Matrix{Float64} = exp(A * m)
        tempM = similar(ϕ)
        d = m
        P = PCA_reduce(minkowski_sum(U, linear_map(ϕ, U)))
        #disc :: Zonotope = X0
        dia :: Matrix{Float64} = diagm(ones(XDim))
        while d < initialTimeStep
            inputDiscritezationDict[d] = P
            P = minkowski_sum(P, linear_map(ϕ, P))
            
            α = (exp(ANorm * d) - 1 - d * ANorm) / XNorm
            ϕp = (dia+ϕ)/2
            ϕm = (dia-ϕ)/2
            gens = hcat(ϕp * XG, ϕm * XC, ϕm * XG, α*dia) #hcat(ϕp * XG, ϕm * XC, ϕm * XG)

            disc = Zonotope(ϕp*XC, gens)#minkowski_sum(Zonotope(ϕp*XC, gens), Zonotope(zeros(XDim), α*dia))
            discritezationDict[d] = (disc, ϕ)

            mul!(tempM, ϕ , ϕ)
            copy!(ϕ, tempM)
            d = d * 2
        end
        α = (exp(ANorm * d) - 1 - d * ANorm) / XNorm
        ϕp = (dia+ϕ)/2
        ϕm = (dia-ϕ)/2
        gens = hcat(ϕp * XG, ϕm * XC, ϕm * XG, α*dia) 

            disc = Zonotope(ϕp*XC, gens)
            discritezationDict[d] = (disc, ϕ)
        inputDiscritezationDict[initialTimeStep] = PCA_reduce(P)
    end
    STEPS = ceil(Integer, endtime / initialTimeStep) + 1

    VP = size(genmat(inputDiscritezationDict[initialTimeStep]), 2)
    V = Vector{Zonotope{N,Vector{N},Matrix{N}}}(undef, STEPS)
    V[1] = inputDiscritezationDict[initialTimeStep]

    @inbounds for i in 2:STEPS
        V[i] = Zonotope(Vector{N}(undef, XDim), Matrix{N}(undef, XDim, VP) )
    end

    let (_, ϕ) = discritezationDict[initialTimeStep]
        @inbounds for id in 2:STEPS
            linear_map!(V[id], ϕ, V[id - 1])
        end
    end


    S = Vector{Zonotope}()
    sizehint!(S, size(V, 1))
    push!(S, V[1])
    for id in 2:size(V, 1)
        t = PCA_reduce(minkowski_sum(S[id - 1], V[id]))
        push!(S, t)
    end

    newR :: Zonotope, _ = discritezationDict[initialTimeStep]
    i = 1

    Φ :: Matrix{Float64} = diagm(ones(Float64, size(A, 2)))
    tempM = similar(Φ)
    tempXG = Matrix{Float64}(undef, 1, 150)
    RG = similar(genmat(newR))
    RC = similar(newR.center)
    #ϕT :: Matrix{Float64} = Φ

    while time < endtime
        #=if i% 100 == 0
            println("Time at step $i: $time")
        end=#
        attempts = 1
        approveFlag = false

        while !approveFlag
            if currentTimeStep < m
                println("Error model fails at time $time, constraint is not satisfied after $attempts attempts")
                return (R, timeStepRecorder, attemptsRecorder)
            end

            if changedTimeStep

                newR, ϕt = discritezationDict[currentTimeStep]

                #ϕT = Φ * ϕt
                newR = linear_map(Φ, newR)
                
                changedTimeStep = false
            else 
                _, ϕt = discritezationDict[currentTimeStep]
                newR = R[i - 1]
            end

            #newR = linear_map(ϕt, newR)
            RC = newR.center
            RG = genmat(newR)
            msum::Zonotope = linear_map_zonotope_nD(ϕt, RC, RG) #minkowski_sum(newR, S[ceil(Integer, (time + currentTimeStep) / initialTimeStep)])
            if !intersects(msum, constraint)#intersectss(msum.center, genmat(msum), h, f, tempXG)
                approveFlag = true
                newR = msum
                mul!(tempM, Φ, ϕt)
                copy!(Φ, tempM) #Her kan man vente med at udregne Phi * phit indtil at man ved hvor mange gange at man vil gange phit på, så at man kan lave et mere effektivt kald på matmul(Phi, phit, phit, ...)
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
end=#

function smallStep(newR, ϕt, RC, RG) :: Zonotope
    copy!(RC, newR.center)
    mul!(RG, ϕt, genmat(newR))
    return Zonotope( ϕt * RC, RG)#
end