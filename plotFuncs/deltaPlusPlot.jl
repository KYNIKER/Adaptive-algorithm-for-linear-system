# Based on the paper JuliaReach: a Toolbox for Set-Based Reachability
using Plots, LazySets, LinearAlgebra, BenchmarkTools, FastExpm, Profile, PProf


include("../helperfunctions.jl")
include("../models.jl")
include("../models/heat/heat_load.jl")
include("../models/motor/motor_load.jl")
include("../models/motor/motor_load.jl")
include("../models/building/building_load.jl")
include("../models/PDE/pde_load.jl")
include("../models/ISS/iss_load.jl")
include("../models/beam/beam_load.jl")
include("../models/FOM/fom_load.jl")
include("../models/MNA1/mna1_load.jl")
include("../models/MNA5/mna5_load.jl")
#include("CegarFunctions.jl")
include("../CegarInhomogenous.jl")
include("plotHelper.jl")


# We load a simple coswave
#const μ = 0.
#const STRATEGY = 2

initialTimeStep = 0.5
Digits = 5


A, B, U, P₁, T, constraint, dimToPlot = load_building()



#T = [0, 5]
STRATEGY = 1

constraint = isa(constraint, Array) ? constraint : [constraint]

boxes1, timesteps1, attemptsRecorder1 = cegarInputSystem(A, B, initialTimeStep, T, P₁, U, constraint, Digits, STRATEGY)
shapes1, maxVal1, minVal1 = getShapes(boxes1, timesteps1)

initialTimeStep = 2


boxes2, timesteps2, attemptsRecorder2 = cegarInputSystem(A, B, initialTimeStep, T, P₁, U, constraint, Digits, STRATEGY)
shapes2, maxVal2, minVal2 = getShapes(boxes2, timesteps2)

initialTimeStep = 8

boxes3, timesteps3, attemptsRecorder3 = cegarInputSystem(A, B, initialTimeStep, T, P₁, U, constraint, Digits, STRATEGY)
shapes3, maxVal3, minVal3 = getShapes(boxes3, timesteps3)

println("Finished simulations")

constraintValAdjusted = constraint[1].b * 1.1
maxVal = max(maxVal1, maxVal2, maxVal3, constraintValAdjusted) 
minVal = min(minVal1, minVal2, minVal3, constraintValAdjusted)


p = plot(dpi=300, thickness_scaling=1, ylims=(minVal, maxVal), xlims=(0, maximum(T)), xlabel="Time", ylabel="Value")


for i in eachindex(shapes1)
    plot!(p, shapes1[i], vars=(1,0), c=:forestgreen, alpha=:0.2,
        label = i == 1 ? "δ⁺ = 0.5" : "")
end

for i in eachindex(shapes2)
    plot!(p, shapes2[i], vars=(1,0), c=:blue, alpha=:0.2,
        label = i == 1 ? "δ⁺ = 2" : "")
end

for i in eachindex(shapes3)
    plot!(p, shapes3[i], vars=(1,0), c=:red, alpha=:0.2,
        label = i == 1 ? "δ⁺ = 8" : "")
end

println("Amount of steps strategy 1: ", length(shapes1))
println("Amount of steps strategy 2: ", length(shapes2))
println("Amount of steps strategy 3: ", length(shapes3))

# Plot constraint

#plot!(LazySets.HalfSpace([0.0, -1.0], -constraint[1].b), lab="Unsafe Region", c=:black)
plot!(LazySets.HalfSpace([0.0, -1.0], -constraint[1].b), lab="Unsafe Region", c=:black)

# println("Amount of attempts strategy 1: ", sum(attemptsRecorder2))
# println("Amount of attempts strategy 2: ", sum(attemptsRecorder3))

# println("Timesteps strategy 1: ", timesteps2)
# println("Timesteps strategy 2: ", timesteps3)


# ts = range(0, T; length=100)
# # Y = reduce(hcat, f.(ts))
# # X = Y[1:n, :]   
# x0 = P₁.center
# println(x0)


savefig(p, "deltaPlusPlot.pdf")
plot(p)