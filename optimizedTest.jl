# Based on the paper JuliaReach: a Toolbox for Set-Based Reachability
using Plots, LazySets, LinearAlgebra
include("helperfunctions.jl")
include("models.jl")
include("CegarFunctions.jl")

const μ = 0.01

initialTimeStep = 0.4
strategy = 2
digits = 4
reuse = true
plotConstraint = true

A, P₁, constraint, T, dimToPlot = loadHeat01()

T = [0, 8]

###
#@profview boxes2, timesteps, attemptsRecorder = reachSetsCegarInput(A, initialTimeStep, T, P₁, constraint, μ, strategy, digits, reuse)

#@time boxes2, timesteps, attemptsRecorder = reachSetsCegarInput(A, initialTimeStep, T, P₁, constraint, μ, strategy, digits, reuse)
@time boxes2, timesteps, attemptsRecorder = reachSetsCegar(A, initialTimeStep, T, P₁, constraint, strategy, digits)
#@profview boxes2, timesteps, attemptsRecorder = reachSetsCegar(A, initialTimeStep, T, P₁, constraint, strategy, digits)


corners2 = Vector(undef, size(boxes2, 1))

@time begin
    for i in 1:(size(boxes2, 1))
        H = box_approximation(boxes2[i])
        H_proj = project(H, [dimToPlot]) # Only get the dimension we care about
        corners2[i] = vertices_list(H_proj)
    end
end

p = nothing
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

@time rectangleFromHBoxWithTimestepArray(shapes2, corners2, timesteps, minimum(T), 1)

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
    if 0 < constraint.b
        plot!(HalfSpace([0.0, -1.0], -constraint.b), lab="constraint", c=:purple)
    else
        plot!(HalfSpace([0.0, 1.0], constraint.b), lab="constraint", c=:purple)
    end
end

println("Amount of attempts: ", sum(attemptsRecorder))
println("Unique timesteps ", length(unique(timesteps)))

#println("List of attempts: ", attemptsRecorder)

plot(p)