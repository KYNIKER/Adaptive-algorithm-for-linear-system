# Based on the paper JuliaReach: a Toolbox for Set-Based Reachability
using Plots, LazySets, LinearAlgebra, BenchmarkTools, FastExpm, Profile, PProf, ReachabilityAnalysis


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

load_func = load_building
A, B, ballβ, P₁, T, constraint, dimToPlot = load_func()


# ReACT
initialTimeStep = 0.5
Digits = 3
STRATEGY = 1

boxes1, timesteps1, attemptsRecorder1 = cegarInputSystem(A, B, initialTimeStep, T, P₁, ballβ, constraint, Digits, STRATEGY)
shapes1, maxVal1, minVal1 = getShapes(boxes1, timesteps1)


# GLGM
δ = 0.002
tVal = maximum(T)
n = size(A, 1)

sys = @system(x' = Ax + Bu, x ∈ Universe(n), u ∈ ballβ)
alg = GLGM06(;δ=δ, approx_model=Forward())
prob = InitialValueProblem(sys, P₁)
sol = solve(prob, alg; T=tVal) # Running the actual time
#res = mapreduce(c -> ρ(c.a, sol) <= c.b, &, constraint) # Check if hits constraint

GLGMShapes = []
fp = flowpipe(sol)
for (i, Ω) in enumerate(fp)
    currT = i * δ
    proj = project(Ω, dimToPlot)
    max = ρ([1.0], proj)
    min = -ρ([-1.0], proj)
    newShape = Shape([currT, currT + δ, currT + δ, currT], [min, min, max, max])
    push!(GLGMShapes, newShape)
end

# Plotting 
#p = plot(dpi=300, thickness_scaling=1, ylims=(minVal, maxVal), xlims=(0, maximum(T)), xlabel="Time", ylabel="Value")
#p = plot(dpi=300, thickness_scaling=1, xlims=(0, maximum(T)), xlabel="Time", ylabel="Value")

constraintValAdjusted = constraint[1].b * 1.1
maxVal = max(maxVal1, constraintValAdjusted) 
minVal = min(minVal1, constraintValAdjusted)

p = plot(dpi=300, thickness_scaling=1, ylims=(minVal, maxVal), xlims=(0, maximum(T)), xlabel="Time", ylabel="Value")

for i in eachindex(shapes1)
    plot!(p, shapes1[i], vars=(1,0), c=:forestgreen, alpha=:0.2,
        label = i == 1 ? "ReACT" : "")
end
for i in eachindex(GLGMShapes)
    plot!(p, GLGMShapes[i], vars=(1,0), c=:blue, alpha=:1, 
    label = i == 1 ? "GLGM06" : "")
end

plot!(LazySets.HalfSpace([0.0, -1.0], -constraint[1].b), lab="Unsafe Region", c=:black)

savefig(p, "GlGMvsReACT.pdf")
plot(p)


# println(sol)
# # Plotting
# plot(sol, vars=(1, 2), xlab="x", ylab="v", lw=0.5, color=:blue)
# plot!(constraint, lw=0.5, color=:red)
