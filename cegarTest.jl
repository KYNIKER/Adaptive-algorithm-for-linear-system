# Based on the paper JuliaReach: a Toolbox for Set-Based Reachability
using Plots, LazySets, LinearAlgebra
include("helperfunctions.jl")
include("models.jl")

const UseCrane = false # Crane usually has T = 15
const μ = 0.01

Tstart = 0
T = 4
initialTimeStep = 0.4
strategy = 3
dimToPlot = 2

if !UseCrane
    constraint = HalfSpace([0., 1.], -1.8)
else
    constraint = HalfSpace([0., 0., 0., 0., 0., 1.], -1.57)
end

# const tΔ = 0.01
# const r = 1.2


A, P₁ = loadCosWave()
#A, P₁ = loadCrane()
if UseCrane
    A, P₁ = loadCrane()
end

p = nothing

if UseCrane
    p = plot(dpi=300, thickness_scaling=1)
else
    p = plot(dpi=300, thickness_scaling=1, ylims=(-1.75, 2.5))
end


###
@time boxes2, timesteps, attemptsRecorder = reachSetsCegar(A, initialTimeStep, [Tstart, T], P₁, constraint, strategy)


corners2 = Vector(undef, size(boxes2, 1))

@time begin
    for i in 1:(size(boxes2, 1))
        corners2[i] = vertices_list(box_approximation(boxes2[i]))
    end
end

shapes2 = Vector{Shape}(undef, size(boxes2, 1))
#@time rectangleFromHBox!(shapes2, corners2, 2*tΔ, 2)7



@time rectangleFromHBoxWithTimestepArray(shapes2, corners2, timesteps, Tstart, dimToPlot)

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

if !UseCrane
    plot!(constraint, lab="constraint", c=:purple)
end

println("Amount of attempts: ", sum(attemptsRecorder))
println("Unique timesteps ", length(unique(timesteps)))

println("List of attempts: ", attemptsRecorder)

plot(p)