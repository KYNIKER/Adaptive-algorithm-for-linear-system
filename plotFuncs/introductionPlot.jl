# Based on the paper JuliaReach: a Toolbox for Set-Based Reachability
using Plots, LazySets, LinearAlgebra, BenchmarkTools, Profile, PProf


include("../helperfunctions.jl")
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
include("../ReACT.jl")
include("plotHelper.jl")


# We load a simple coswave
STRATEGY = 2

initialTimeStep = 0.8
Digits = 3



A = [0. 1.;
    -1. 0.]
P₁ = Zonotope([0., 1.5], [[0.0; 0.05]])
constraint = LazySets.HalfSpace([0., 1.], -1.7)
T = [0, 3]
dimToPlot = 2
# No input
B = diagm([0.0, 0.0])
U = LazySets.Zonotope([0.0, 0.0], [[0.0, 0.0]])
#U::Zonotope = BallInf([0.0], 0.0)


constraint = isa(constraint, Array) ? constraint : [constraint]

boxes1, timesteps1, attemptsRecorder1 = PlotReACT(A, B, initialTimeStep, T, P₁, U, constraint, Digits, STRATEGY)
#println(boxes1)
shapes1, maxVal1, minVal1 = getShapes(boxes1, timesteps1)

#println("Finished simulations")

# Get new values
U = Zonotope(zeros(dim(P₁)), [zeros(dim(P₁))])
initialTimeStep = 0.1

println("Starting second simulation with timestep size: ", initialTimeStep)

boxes2, timesteps2, attemptsRecorder2 = PlotReACT(A, B, initialTimeStep, T, P₁, U, constraint, -1, STRATEGY)
shapes2, maxVal2, minVal2 = getShapes(boxes2, timesteps2)


println("Finished simulations")

constraintValAdjusted = constraint[1].b * 1.1
maxVal = max(maxVal2, constraintValAdjusted)
minVal = min(minVal2, constraintValAdjusted)


p = plot(dpi=1200, thickness_scaling=1, guidefontsize=25,
    xguidefont=font(25, "Times"),
    yguidefont=font(25, "Times"),
    xtick=([0, 1], ["0", "T"]),
    ylims=(minVal, maxVal), xlims=(0, maximum(T)), xlabel="Time", ylabel="Value")

#=
for i in eachindex(shapes1)
    if i == 1
        plot!(p, shapes1[i], vars=(1, 0), c=:forestgreen, alpha=:0.2,
            label="ReACT")
    else
        plot!(p, shapes1[i], vars=(1, 0), c=:forestgreen, alpha=:0.2,
            label="")
    end
end
=#
for i in eachindex(shapes2)
    if i == 1
        plot!(p, shapes2[i], vars=(1, 0), c=:blue, alpha=:0.2,
            label="Fixed Timestep")
    else
        plot!(p, shapes2[i], vars=(1, 0), c=:blue, alpha=:0.2,
            label="")
    end
end

# Plot real coswave
#ω = sqrt(2.5)
#t = 0:0.01:10
#x1 = 1.5 .* cos.(ω .* t)

#plot!(p, t, x1, label="Cos(t)", c=:red)

# Plot constraint

plot!(LazySets.HalfSpace([0.0, 1.0], constraint[1].b), lab="Unsafe Region", c=:black)



plot(p)