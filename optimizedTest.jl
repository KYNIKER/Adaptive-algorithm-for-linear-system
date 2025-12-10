# Based on the paper JuliaReach: a Toolbox for Set-Based Reachability
using Plots, LazySets, LinearAlgebra, BenchmarkTools, FastExpm, Profile, PProf
include("helperfunctions.jl")
include("models.jl")
include("models/heat/heat_load.jl")
include("models/motor/motor_load.jl")
include("models/motor/motor_load.jl")
include("models/building/building_load.jl")
include("models/PDE/pde_load.jl")
include("models/ISS/iss_load.jl")
include("models/beam/beam_load.jl")
include("models/FOM/fom_load.jl")
include("models/MNA1/mna1_load.jl")
include("models/MNA5/mna5_load.jl")
#include("CegarFunctions.jl")
include("CegarInhomogenous.jl")

const μ = 0.001
const STRATEGY = 0

initialTimeStep = 0.001
#strategy = 1
Digits = 5
reuse = true
plotConstraint = true
input = true
plotOutput = false

if input
    A, B, ballβ, P₁, T, constraint, dimToPlot = load_beam()
    #t1 = fastExpm(A .* 1; threshold=eps(Float64), nonzero_tol=eps(Float64))
    #t2 = fastExpm(A .* 0.25; threshold=eps(Float64), nonzero_tol=eps(Float64))
    #println(unique(t1 - t2 * t2 * t2 * t2))
    println("A invertible?", isinvertible(A))
    #=ANorm = norm(A, Inf)
    m = initialTimeStep / 2^(ceil(Integer, log2(initialTimeStep)) + ceil(Integer, -log2(10.0^(-Digits))) - 1)
    β = (exp(ANorm*(m))-1)*norm(ballβ)/ANorm
    ballβ = Zonotope(zeros(dim(ballβ)), β*I(dim(ballβ)))=#
else
    A, P₁, constraint, T, dimToPlot = loadHeat01()

    ANorm = norm(A, Inf)
    m = initialTimeStep / 2^(ceil(Integer, log2(initialTimeStep)) + ceil(Integer, -log2(10.0^(-Digits))) - 1)
    β = (exp(ANorm*(m))-1)*μ/ANorm
    #println("original area ballβ: ", area(Zonotope(zeros(dim(P₁)), ((exp(ANorm*(initialTimeStep))-1)*μ/ANorm)*I(dim(P₁)))))
    ballβ = Zonotope(zeros(dim(P₁)), β*I(dim(P₁)))
end

constraint = isa(constraint, Array) ? constraint : [constraint]

###
#@profview boxes2, timesteps, attemptsRecorder = reachSetsCegarInput(A, initialTimeStep, T, P₁, constraint, μ, strategy, digits, reuse)

#@time boxes2, timesteps, attemptsRecorder = reachSetsCegarInput(A, initialTimeStep, T, P₁, constraint, μ, strategy, digits, reuse)
println("initialTimeStep: ", initialTimeStep)
#boxes2, timesteps, attemptsRecorder = cegarInputSystem(A, initialTimeStep, T, P₁, ballβ, constraint, Digits)
#@time boxes2, timesteps, attemptsRecorder = cegarInputSystem(A, B, initialTimeStep, T, P₁, ballβ, constraint, Digits)

#sent = OneTimeStepSystem(A, B, initialTimeStep, T, P₁, ballβ, constraint, Digits, STRATEGY)
res = cegarInputSystemNoOutput(A, B, initialTimeStep, T, P₁, ballβ, constraint, Digits, STRATEGY)

#=Profile.Allocs.clear()
Profile.Allocs.@profile sample_rate=1.0 cegarInputSystemNoOutput(A, B, initialTimeStep, T, P₁, ballβ, constraint, Digits, STRATEGY)

PProf.Allocs.pprof(from_c=false)=#

#@profview cegarInputSystemNoOutput(A, B, initialTimeStep, T, P₁, ballβ, constraint, Digits, STRATEGY)
# cegarInputSystem(A, initialTimeStep, T, P₁, ballβ, constraint, Digits)
#@time boxes2, timesteps, attemptsRecorder = reachSetsCegar(A, initialTimeStep, T, P₁, constraint, strategy, digits)

if plotOutput
    boxes2, timesteps, attemptsRecorder = cegarInputSystem(A, B, initialTimeStep, T, P₁, ballβ, constraint, Digits)
    corners2 = Vector(undef, size(boxes2, 1))

    begin
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

    rectangleFromHBoxWithTimestepArray(shapes2, corners2, timesteps, minimum(T), 1)

    for i in eachindex(shapes2)
        attemptCount = attemptsRecorder[i]
        if attemptCount == 1 # First attempt
            plot!(p, shapes2[i], vars=(1,0), c=:green, alpha=:0.2, lab="")
        elseif attemptCount == 2 # Second
            plot!(p, shapes2[i], vars=(1,0), c=:yellow, alpha=:0.2, lab="") 
        elseif attemptCount == 3 # Third
            plot!(p, shapes2[i], vars=(1,0), c=:red, alpha=:0.2, lab="")
        else # More than three attempts
            plot!(p, shapes2[i], vars=(1,0), c=:red, alpha=:0.2, lab="")
        end
    end
    ##

if plotConstraint
    if 0 < constraint[1].b
        plot!(LazySets.HalfSpace([0.0, -1.0], -constraint[1].b), lab="constraint", c=:purple)
    else
        plot!(LazySets.HalfSpace([0.0, 1.0], constraint[1].b), lab="constraint", c=:purple)
    end
end

    println("Amount of attempts: ", sum(attemptsRecorder))
    println("Unique timesteps ", unique(timesteps))

    #println("List of attempts: ", attemptsRecorder)

    plot(p)
end