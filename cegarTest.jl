# Based on the paper JuliaReach: a Toolbox for Set-Based Reachability
using Plots, LazySets, LinearAlgebra
include("helperfunctions.jl")
include("models.jl")
#include("CegarFunctions.jl")
include("CegarInhomogenous.jl")

const UseCrane = false # Crane usually has T = 15
const μ = 0.00

Tstart = 0.0
T = 4.
initialTimeStep = 0.02
strategy = 2
dimToPlot = 2
Digits = 4
reuse = false
plotConstraint = true

# const tΔ = 0.01
# const r = 1.2

A, P₁, constraint, T, dimToPlot = loadHeat01()
#A, P₁ = loadCrane()
if UseCrane
    A, P₁, constraint, T, dimToPlot = loadCrane()
end

p = nothing

if UseCrane
    p = plot(dpi=300, thickness_scaling=1)
else
    p = plot(dpi=300, thickness_scaling=1, ylims=(-1.85, 2.5))
end
ANorm = norm(A, Inf)
m = initialTimeStep / 2^(ceil(Integer, log2(initialTimeStep)) + ceil(Integer, -log2(10.0^(-Digits))) - 1)
β = (exp(ANorm*(m))-1)*μ/ANorm
#println("original area ballβ: ", area(Zonotope(zeros(dim(P₁)), ((exp(ANorm*(initialTimeStep))-1)*μ/ANorm)*I(dim(P₁)))))
ballβ = Zonotope(zeros(dim(P₁)), β*I(dim(P₁)))
###
#ProfileView.Profile.init()
#using ProfileCanvas
#ProfileCanvas.@profview boxes2, timesteps, attemptsRecorder = cegarInputSystem(A, initialTimeStep, [Tstart, T], P₁, ballβ, constraint, 1)
#@profview boxes2, timesteps, attemptsRecorder = cegarInputSystem(A, initialTimeStep, [Tstart, T], P₁, ballβ, constraint, 2)
t = @elapsed boxes2, timesteps, attemptsRecorder = cegarInputSystem(A, initialTimeStep, T, P₁, ballβ, constraint, Digits)
println(t)

@profview boxes2, timesteps, attemptsRecorder = cegarInputSystem(A, initialTimeStep, T, P₁, ballβ, constraint, Digits)

#@time boxes2, timesteps, attemptsRecorder = reachSetsCegar(A, initialTimeStep, [Tstart, T], P₁, constraint, strategy, 2)


corners2::Vector{Vector{Vector{Float64}}} = Vector(undef, size(boxes2, 1))

@time begin
    for i in 1:(size(boxes2, 1))
        #corners2[i] = vertices_list(box_approximation(boxes2[i]))
        H = box_approximation(boxes2[i])
        H_proj = project(H, [dimToPlot]) # Only get the dimension we care about
        corners2[i] = vertices_list(H_proj)
    end
end

#=if !UseCrane
    plot!(constraint, lab="constraint", c=:purple)
end=#
if plotConstraint
    cornersToSearch = [value[1] for box in corners2 for value in box]
    maxVal = maximum(cornersToSearch)
    minVal = minimum(cornersToSearch)
    # Take min and max with constraint
    constraintValAdjusted = constraint.b * 1.1
    maxVal = max(maxVal, constraintValAdjusted) 
    minVal = min(minVal, constraintValAdjusted)

    p = plot(dpi=300, thickness_scaling=1, ylims=(minVal, maxVal))
else
    p = plot(dpi=300, thickness_scaling=1)
end
shapes2 = Vector{Shape}(undef, size(boxes2, 1))
#@time rectangleFromHBox!(shapes2, corners2, 2*tΔ, 2)7



#@time rectangleFromHBoxWithTimestepArray(shapes2, corners2, timesteps, Tstart, dimToPlot)
@time rectangleFromHBoxWithTimestepArray(shapes2, corners2, timesteps, minimum(T), 1)
#@time rectangleFromHBox2Dims(shapes2, corners2, 1, 2)



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


println("Amount of attempts: ", sum(attemptsRecorder))
println("Unique timesteps ", length(unique(timesteps)))

println("List of attempts: ", attemptsRecorder)

plot(p)