# Based on the paper JuliaReach: a Toolbox for Set-Based Reachability
using Plots, LazySets, LinearAlgebra
include("helperfunctions.jl")
include("models.jl")
include("models/heat/heat_load.jl")
include("models/motor/motor_load.jl")
include("models/motor/motor_load.jl")
include("models/building/building_load.jl")
include("models/PDE/pde_load.jl")
#include("CegarFunctions.jl")
include("CegarInhomogenous.jl")


function cegarInputSystemOldDisc(A, initialTimeStep, interval, X0::Zonotope{N,Vector{N},Matrix{N}}, U::Zonotope, constraint, Digits :: Integer) where {N}
    stepsBeforeReduce = 5
    
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
    SVDict = Dict()#Dict{Float64, Tuple{Float64, Zonotope{N, Vector{N}, Matrix{N}}, Zonotope{N, Vector{N}, Matrix{N}}}}

    time::Float64 = minimum(interval)
    endtime::Float64 = maximum(interval)

    currentTimeStep = copy(initialTimeStep)
    timeStepRecorder = Float64[]
    attemptsRecorder = Integer[]

    k = size(U.generators, 2)

    
    #ϕ :: Matrix{Float64} = fastExpm(A .* initialTimeStep; threshold=eps(Float64), nonzero_tol=eps(Float64))

    dU = box_approximation_symmetric(initialTimeStep * U)
        #println(typeof(dU))
        #println(typeof(E_ψ(U, d, A)))
    P = minkowski_sum(dU, E_ψ(U, initialTimeStep, A))#minkowski_sum(U, linear_map(ϕ, U))


    let ϕ :: Matrix{Float64} = fastExpm(A .* m; threshold=eps(Float64), nonzero_tol=eps(Float64))
        tempM = similar(ϕ)
        d = m
        #println(typeof(dU))
        #println(typeof(E_ψ(U, d, A)))
        #disc :: Zonotope{N,Vector{N},Matrix{N}} = copy(X0)
        dia :: Matrix{Float64} = diagm(ones(XDim))
        i = 1
        while d < initialTimeStep
            println("d: ", d)

            i += 1
            #=α = (exp(ANorm * d) - 1 - d * ANorm) / XNorm
            ϕp = (dia+ϕ)/2
            ϕm = (dia-ϕ)/2
            gens::Matrix{Float64} = hcat(ϕp * XG, ϕm * XC, ϕm * XG, α*dia) #hcat(ϕp * XG, ϕm * XC, ϕm * XG)

            disc = Zonotope(ϕp*XC, gens)#minkowski_sum(Zonotope(ϕp*XC, gens), Zonotope(zeros(XDim), α*dia))
            =#
            #disc = minkowski_sum(convex_hull(X0, ϕ * minkowski_sum(X0, box_approximation_symmetric(d * U))),minkowski_sum(E_ψ(U, d, A), E⁺(X0, d, A)))
            lt = concretize(minkowski_sum(convert(Zonotope, ϕ * X0), box_approximation_symmetric(d * U)))
            println("lt")
            rt = concretize(minkowski_sum(E_ψ(U, d, A), E⁺(X0, d, A)))
            println("rt")

            #println(typeof(lt), " vs ", typeof(rt))
            disc = overapproximate(CH(X0, concretize(minkowski_sum(lt, rt))), Zonotope)
            println("disc: ", d)
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
        println("lt")
        rt = concretize(minkowski_sum(E_ψ(U, d, A), E⁺(X0, d, A)))
        println("rt")

        #println(typeof(lt), " vs ", typeof(rt))
        disc = overapproximate( CH(X0, concretize(minkowski_sum(lt, rt))), Zonotope)
        println("disc: ", d)
        discritezationDict[d] = (copy(disc), copy(ϕ))
        #P = minkowski_sum(U, linear_map(ϕ, U))
        #P = minkowski_sum(dU, E_ψ(U, d, A))
    end

    println("Dicts done!")

    STEPS = ceil(Integer, endtime / initialTimeStep) + 1
    V :: Zonotope{N,Vector{N},Matrix{N}} = P

    S ::Zonotope{N,Vector{N},Matrix{N}} = P
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
                
                # push!(tubes, minkowski_sum(newRR, S))
                # push!(timeStepRecorder, currentTimeStep)
                # push!(attemptsRecorder, attempts)
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
            if mapreduce(c -> !intersects(nsum, c), :(&&), constraint) #intersectss(msum.center, genmat(msum), h, f, tempXG)
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










# Running the stuff








const μ = 0.01
const STRATEGY = 0

initialTimeStep = 0.016 * 0.5 ^ 2
#strategy = 1
Digits = 5
reuse = true
plotConstraint = true
input = true

if input
    A,  ballβ, P₁, T, constraint, dimToPlot = load_heat_input()
    A,  ballβ, P₁, T, constraint, dimToPlot = loadCosWave()
    println("A invertible?", isinvertible(A))
    #=ANorm = norm(A, Inf)
    m = initialTimeStep / 2^(ceil(Integer, log2(initialTimeStep)) + ceil(Integer, -log2(10.0^(-Digits))) - 1)
    β = (exp(ANorm*(m))-1)*norm(ballβ)/ANorm
    ballβ = Zonotope(zeros(dim(ballβ)), β*I(dim(ballβ)))=#
else
    A, P₁, constraint, T, dimToPlot = loadHeat01()

    A, P₁, constraint, T, dimToPlot = loadCosWave()
    constraint = [constraint]

    ANorm = norm(A, Inf)
    m = initialTimeStep / 2^(ceil(Integer, log2(initialTimeStep)) + ceil(Integer, -log2(10.0^(-Digits))) - 1)
    β = (exp(ANorm*(m))-1)*μ/ANorm
    #println("original area ballβ: ", area(Zonotope(zeros(dim(P₁)), ((exp(ANorm*(initialTimeStep))-1)*μ/ANorm)*I(dim(P₁)))))
    ballβ = Zonotope(zeros(LazySets.dim(P₁)), β*I(LazySets.dim(P₁)))
end

###
#@profview boxes2, timesteps, attemptsRecorder = reachSetsCegarInput(A, initialTimeStep, T, P₁, constraint, μ, strategy, digits, reuse)

#@time boxes2, timesteps, attemptsRecorder = reachSetsCegarInput(A, initialTimeStep, T, P₁, constraint, μ, strategy, digits, reuse)
println("initialTimeStep: ", initialTimeStep)
@time boxes2, timesteps, attemptsRecorder = cegarInputSystemOldDisc(A, initialTimeStep, T, P₁, ballβ, constraint, Digits)
#@time boxes2, timesteps, attemptsRecorder = reachSetsCegar(A, initialTimeStep, T, P₁, constraint, strategy, digits)
#@profview boxes2, timesteps, attemptsRecorder = cegarInputSystem(A, initialTimeStep, T, P₁, ballβ, constraint, Digits)



corners2 = Vector(undef, size(boxes2, 1))



@time begin
    for i in 1:(size(boxes2, 1))
        H = box_approximation(boxes2[i])
        #println(H)
        #H.radius = abs(H.radius)
        H_proj = LazySets.project(H, [dimToPlot]) # Only get the dimension we care about
        corners2[i] = vertices_list(H_proj)
    end
end

p = nothing
if plotConstraint
    cornersToSearch = [value[1] for box in corners2 for value in box]
    maxVal = maximum(cornersToSearch)
    minVal = minimum(cornersToSearch)
    # Take min and max with constraint
    constraintValAdjusted = constraint[1].b * 1.1
    maxVal = max(maxVal, constraintValAdjusted) 
    minVal = min(minVal, constraintValAdjusted)

    p = plot(dpi=300, thickness_scaling=1, ylims=(minVal, maxVal))
else
    p = plot(dpi=300, thickness_scaling=1)
end


shapes2 = Vector{Shape}(undef, size(boxes2, 1))

@time rectangleFromHBoxWithTimestepArray(shapes2, corners2, timesteps, minimum(T), 1)

for i in eachindex(shapes2)
    plot!(p, shapes2[i], vars=(1,0), c=:green, alpha=:0.2, lab="")
end
##

# We only use the endtime of the first (but worse) model
Tstart = minimum(T)
Tend = sum(timesteps) + Tstart

T = [Tstart, Tend]





@time boxes2, timesteps, attemptsRecorder = cegarInputSystem(A, initialTimeStep, T, P₁, ballβ, constraint, Digits)
#@time boxes2, timesteps, attemptsRecorder = reachSetsCegar(A, initialTimeStep, T, P₁, constraint, strategy, digits)
#@profview boxes2, timesteps, attemptsRecorder = cegarInputSystem(A, initialTimeStep, T, P₁, ballβ, constraint, Digits)

corners2 = Vector(undef, size(boxes2, 1))



@time begin
    for i in 1:(size(boxes2, 1))
        H = box_approximation(boxes2[i])
        #println(H)
        #H.radius = abs(H.radius)
        H_proj = LazySets.project(H, [dimToPlot]) # Only get the dimension we care about
        corners2[i] = vertices_list(H_proj)
    end
end

# p = nothing
# if plotConstraint
#     cornersToSearch = [value[1] for box in corners2 for value in box]
#     maxVal = maximum(cornersToSearch)
#     minVal = minimum(cornersToSearch)
#     # Take min and max with constraint
#     constraintValAdjusted = constraint[1].b * 1.1
#     maxVal = max(maxVal, constraintValAdjusted) 
#     minVal = min(minVal, constraintValAdjusted)

#     p = plot(dpi=300, thickness_scaling=1, ylims=(minVal, maxVal))
# else
#     p = plot(dpi=300, thickness_scaling=1)
# end


shapes2 = Vector{Shape}(undef, size(boxes2, 1))

@time rectangleFromHBoxWithTimestepArray(shapes2, corners2, timesteps, minimum(T), 1)

for i in eachindex(shapes2)
    plot!(p, shapes2[i], vars=(1,0), c=:yellow, alpha=:0.2, lab="")
end
##

#println(shapes2)

if plotConstraint
    if 0 < constraint[1].b
        plot!(LazySets.HalfSpace([0.0, -1.0], -constraint[1].b), lab="constraint", c=:purple)
    else
        plot!(LazySets.HalfSpace([0.0, 1.0], constraint[1].b), lab="constraint", c=:purple)
    end
end



plot(p)

