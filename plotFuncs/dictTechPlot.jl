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
include("../DiscretizationStudy.jl")
include("plotHelper.jl")


# We load a simple coswave
#const μ = 0.
STRATEGY = 2

initialTimeStep = 0.016
Digits = 5

#A, P₁, constraint, T, dimToPlot= loadCosWave()
A, B, U, P₁, T, constraint, dimToPlot = load_building()

# T = [0, 2]

constraint = isa(constraint, Array) ? constraint : [constraint]

boxes1, timesteps1, attemptsRecorder1 = cegarInputSystem(A, B, initialTimeStep, T, P₁, U, constraint, Digits, STRATEGY)
shapes1, maxVal1, minVal1 = getShapes(boxes1, timesteps1)

boxes2, timesteps2, attemptsRecorder2 = cegarInputSystemOldDisc(A, B, initialTimeStep, T, P₁, U, constraint, Digits, STRATEGY)
shapes2, maxVal2, minVal2 = getShapes(boxes2, timesteps2)

# boxes2, timesteps2, attemptsRecorder2 = OneTimeStepSystem(A, B, initialTimeStep, T, P₁, U, constraint, Digits, STRATEGY)
# shapes2, maxVal2, minVal2 = getShapes(boxes2, timesteps2)


println("Finished simulations")

constraintValAdjusted = constraint[1].b * 1.1
maxVal = max(maxVal1, maxVal2, constraintValAdjusted) 
minVal = min(minVal1, minVal2, constraintValAdjusted)

p = plot(dpi=300, thickness_scaling=1, ylims=(minVal, maxVal), xlims=(0, maximum(T)), xlabel="Time", ylabel="Value")


for i in eachindex(shapes1)
    plot!(p, shapes1[i], vars=(1,0), c=:black, alpha= i == 1 ? 0.6 : 0.2,
        label = i == 1 ? "ReACT with Dict Tech" : "")
end

for i in eachindex(shapes2)
    plot!(p, shapes2[i], vars=(1,0), c=:forestgreen, alpha= i == 1 ? 0.6 : 0.2,
        label = i == 1 ? "ReACT without Dict Tech" : "")
end


# Plot constraint

#plot!(LazySets.HalfSpace([0.0, 1.0], constraint[1].b), lab="Unsafe Region", c=:black)



plot(p)