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
    XNorm = norm(X0, Inf)
    XDim = dim(X0)
    XG = genmat(X0)
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

            #=
                    α = (exp(ANorm*timestep)-1-timestep*ANorm)/norm(X₀, Inf)
                    ϕ = exp(A*timestep)

                    ϕp = (I+ϕ)/2
                    ϕm = (I-ϕ)/2
                    gens = hcat(ϕp*X₀.generators,ϕm*X₀.center, ϕm*X₀.generators)


                    R :: Zonotope = minkowski_sum(Zonotope(ϕp*X₀.center, gens), Zonotope(zeros(dim(X₀)), α*I(dim(X₀))))

                    return (Zonotope(R.center, R.generators), ϕ)
            =#

        end
        α = (exp(ANorm * d) - 1 - d * ANorm) / XNorm
        ϕp = (dia+ϕ)/2
        ϕm = (dia-ϕ)/2
        gens = hcat(ϕp * XG, ϕm * XC, ϕm * XG, α*dia) #hcat(ϕp * XG, ϕm * XC, ϕm * XG)

            disc = Zonotope(ϕp*XC, gens)#minkowski_sum(Zonotope(ϕp*XC, gens), Zonotope(zeros(XDim), α*dia))
            discritezationDict[d] = (disc, ϕ)
        inputDiscritezationDict[initialTimeStep] = PCA_reduce(P)
    end


    V = Vector{Zonotope}()
    sizehint!(V, ceil(Integer, endtime / initialTimeStep) + 1)
    push!(V, inputDiscritezationDict[initialTimeStep])
    let (_, ϕ) = discritezationDict[initialTimeStep]
        for id in 2:(ceil(Integer, endtime / initialTimeStep) + 1)
            #push!(V, linear_map(ϕ, V[id - 1]))
            t :: Zonotope = linear_map(ϕ, V[id - 1])#Zonotope(ϕ * V[id - 1].center, ϕ * V[id - 1].generators)
            push!(V, t)
        end
    end


    S = Vector{Zonotope}()
    sizehint!(S, size(V, 1))
    push!(S, V[1])
    for id in 2:size(V, 1)
        t = PCA_reduce(minkowski_sum(S[id - 1], V[id]))
        push!(S, t)
        #push!(S, minkowski_sum(S[id - 1], V[id]))

    end

    newR :: Zonotope, _ = discritezationDict[initialTimeStep]
    i = 1
    #=for i in eachindex(S)
        @printf "volume of S[%d]: %.12f\n" i area(S[i])
        #println("volume of V[",i,"]: ", area(V[i]))
    end=#

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
                #println(R, timeStepRecorder)
                println("Error model fails at time $time, constraint is not satisfied after $attempts attempts")
                return (R, timeStepRecorder, attemptsRecorder)
                #throw(AssertionError("Error model fails at time $time, constraint is not satisfied after $attempts attempts"))
            end

            if changedTimeStep
                #=if haskey(discritezationDict, currentTimeStep)
                    newR, ϕt = discritezationDict[currentTimeStep]
                else 
                    newR, ϕt = initialStepNoInput(A, ANorm, currentTimeStep, X0)
                    discritezationDict[currentTimeStep] = (newR, ϕt)
                end=#
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
end