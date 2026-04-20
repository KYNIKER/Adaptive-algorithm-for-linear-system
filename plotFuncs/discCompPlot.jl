# Based on the paper JuliaReach: a Toolbox for Set-Based Reachability
using Plots, LazySets, LinearAlgebra, BenchmarkTools, Profile, PProf, Plots.PlotMeasures, LaTeXStrings
gr()

include("../helperfunctions.jl")
include("../models/heat/heat_load.jl")
include("../models/motor/motor_load.jl")
include("../models/building/building_load.jl")
include("../models/PDE/pde_load.jl")
include("../models/ISS/iss_load.jl")
include("../models/beam/beam_load.jl")
include("../models/MNA1/mna1_load.jl")
include("../models/MNA5/mna5_load.jl")
#include("CegarFunctions.jl")
include("../ReACTv2.jl")
include("plotHelper.jl")


const ∞ = Inf

Digits = 1e-3
initialTimeStep = (2.0)^10 * 1e-3
#initialTimeStep = 0.25
palette = Plots.palette(:fes10) #Plots.palette(:cyclic_mrybm_35_75_c68_n256_s25, 2)
c1 = palette[9]
c2 = palette[6]
LazySets.Comparison.set_tolerance(Float64)
LazySets.Comparison.set_ztol(Float64, 1e-10)
load_func = load_heat_input
A, B, ballβ, P₁, T, constraint, dimToPlot = load_func()
name = "heatDiscComp"

p = plot(dpi=1200, thickness_scaling=1, guidefontsize=25, minorgrid=true, #ϵ=dm / 4,
    legendfont=font(12, "Times"),
    #legendcolumn=-1,
    legend_position=:bottomright,
    tickfont=font(8, "Times"),
    xguidefont=font(12, "Times"),
    yguidefont=font(12, "Times"),
    xtick=([0, maximum(T)], [L"0", L"T"]),
    bottom_margin=2mm,
    left_margin=5mm,
    right_margin=5mm,
    top_margin=2mm,
    xlims=(0, maximum(T)),
    xlabel=L"Time", ylabel=L"x_{%$dimToPlot}")

STRATEGY = 2

constraint = isa(constraint, Array) ? constraint : [constraint]
constraintValAdjusted = constraint[1].b * 1.2

maxVal = constraintValAdjusted
minVal = constraintValAdjusted

#T = 0.5




# boxes1, timesteps1, attemptsRecorder1 = PlotReACT(A, B, initialTimeStep, T, P₁, U, constraint, Digits, STRATEGY)
# shapes1, maxVal1, minVal1 = plotProjectedFlowpipe(boxes1, timesteps1, 0, dimToPlot; approx=true)
# maxVal = max(maxVal1, maxVal)
# minVal = min(minVal1, minVal)
# println(minVal, " ", maxVal)

# for i in eachindex(shapes1)
#     plot!(p, shapes1[i], vars=(1, 0), c=palette[2],
#         label=i == 1 ? L"\delta^{+} / \delta^- = %$initialTimeStep / %$dm" : "")
# end


boxes1, timesteps1 = PlotReACTWithSupport(A, B, initialTimeStep, T, P₁, ballβ, constraint, Digits, [constraint[1].a, -constraint[1].a], STRATEGY)

shapes1, maxVal1, minVal1 = plotSupportFlowpipe(boxes1, timesteps1, 1, 2)
maxVal = max(maxVal1, maxVal)
minVal = min(minVal1, minVal)
println(minVal, " ", maxVal)

for i in eachindex(shapes1)
    plot!(p, shapes1[i], vars=(1, 0), c=c2, la=0.0, alpha=0.7, lw=0.0,
        label=i == 1 ? L"Alg.\: 2: \delta^{+} / \delta^- = %$initialTimeStep / %$Digits" : "")
end

println("Max timestep: $(maximum(timesteps1))")
println("Min timestep: $(minimum(timesteps1))")
initialTimeStep = (2.0)^1 * 1e-3
boxes2, timesteps2 = PlotReACTWithSupport(A, B, initialTimeStep, T, P₁, ballβ, constraint, 1e-3, [constraint[1].a, -constraint[1].a], STRATEGY; naive=true)
shapes2, maxVal2, minVal2 = plotSupportFlowpipe(boxes2, timesteps2, 1, 2)
maxVal = max(maxVal2, maxVal)
minVal = min(minVal2, minVal)
println(minVal, " ", maxVal)
println("$(length(timesteps1)) $(length(timesteps2))")


for i in eachindex(shapes2)
    plot!(p, shapes2[i], vars=(1, 0), c=c1, la=0.0, alpha=0.7, lw=0.01,
        label=i == 1 ? L"Alg.\: 1: \delta^{+} / \delta^- = %$initialTimeStep / %$Digits" : "")
end


#=
initialTimeStep = 2

boxes2, timesteps2, attemptsRecorder2 = PlotReACT(A, B, initialTimeStep, T, P₁, U, constraint, Digits, STRATEGY)
shapes2, maxVal2, minVal2 = getShapes(boxes2, timesteps2)

maxVal = max(maxVal2, maxVal)
minVal = min(minVal2, minVal)


for i in eachindex(shapes2)
    plot!(p, shapes2[i], vars=(1, 0), c=:blue, alpha=:0.2,
        label=i == 1 ? "δ⁺ = " * string(initialTimeStep) : "")
end


initialTimeStep = 8

boxes3, timesteps3, attemptsRecorder3 = PlotReACT(A, B, initialTimeStep, T, P₁, U, constraint, Digits, STRATEGY)
shapes3, maxVal3, minVal3 = getShapes(boxes3, timesteps3)

maxVal = max(maxVal3, maxVal)
minVal = min(minVal3, minVal)

for i in eachindex(shapes3)
    plot!(p, shapes3[i], vars=(1, 0), c=:red, alpha=:0.2,
        label=i == 1 ? "δ⁺ = " * string(initialTimeStep) : "")
end
=#
println("Finished simulations")

ylims!((minVal, maxVal))
yticks!([minVal, 0, constraint[1].b], [string(round(minVal; sigdigits=2)), "0.0", string(constraint[1].b)])





#println("Amount of steps strategy 1: ", length(shapes1))
#println("Amount of steps strategy 2: ", length(shapes2))
#println("Amount of steps strategy 3: ", length(shapes3))

# Plot constraint

#plot!(LazySets.HalfSpace([0.0, -1.0], -constraint[1].b), lab="Unsafe Region", c=:black)

# println("Amount of attempts strategy 1: ", sum(attemptsRecorder2))
# println("Amount of attempts strategy 2: ", sum(attemptsRecorder3))

# println("Timesteps strategy 1: ", timesteps2)
# println("Timesteps strategy 2: ", timesteps3)


# ts = range(0, T; length=100)
# # Y = reduce(hcat, f.(ts))
# # X = Y[1:n, :]   
# x0 = P₁.center
# println(x0)

plot!(LazySets.HalfSpace([0.0, -1.0], -constraint[1].b), lab="Unsafe region", c=:black, fillstyle=:/)
xlims!(0, maximum(T))
lens!(p, [9.25, 9.65], [0.096, 0.102], inset=(1, bbox(0.07, 0.7, 0.29, 0.25)), lc=:black, xtick=[], ytick=[], tickfont=font(20, "Times"), subplot=2)

savefig(p, "plots/" * name * ".pdf")
plot(p)