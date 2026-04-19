using Plots, LazySets, LinearAlgebra, BenchmarkTools, Profile, Plots.PlotMeasures, LaTeXStrings


include("../helperfunctions.jl")
include("../models/heat/heat_load.jl")
include("../models/motor/motor_load.jl")
include("../models/building/building_load.jl")
include("../models/PDE/pde_load.jl")
include("../models/ISS/iss_load.jl")
include("../models/beam/beam_load.jl")
include("../models/MNA1/mna1_load.jl")
include("../ReACT.jl")
include("plotHelper.jl")


const ∞ = Inf

Digits = 3
dm = 10.0^-Digits
initialTimeStep = dm * 2^10
palette = Plots.palette(:fes10)
c1 = palette[9]
c2 = palette[6]
LazySets.Comparison.set_tolerance(Float64)
LazySets.Comparison.set_ztol(Float64, 1e-10)
A, B, U, P₁, T, constraint, dimToPlot = load_heat_input()
name = "heatDiscComp"

p = plot(dpi=1200, thickness_scaling=1, guidefontsize=25, minorgrid=true, #ϵ=dm / 4,
    legendfont=font(12, "Times"),
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

STRATEGY = 1

constraint = isa(constraint, Array) ? constraint : [constraint]
constraintValAdjusted = constraint[1].b * 1.2

maxVal = constraintValAdjusted
minVal = constraintValAdjusted

initialTimeStep = initialTimeStep / 2^5 # We reduce the max initial time step as this is simply unnecessary
boxes1, timesteps1, attemptsRecorder1 = PlotReACTIndividualDisc(A, B, initialTimeStep, T, P₁, U, constraint, dm, STRATEGY)
shapes1, maxVal1, minVal1 = plotProjectedFlowpipe(boxes1, timesteps1, 0, dimToPlot; approx=true)
maxVal = max(maxVal1, maxVal)
minVal = min(minVal1, minVal)
println(minVal, " ", maxVal)

for i in eachindex(shapes1)
    plot!(p, shapes1[i], vars=(1, 0), c=c2, la=0.0, alpha=1.0, lw=0.0,
        label=i == 1 ? L"Alg.\: 1" : "")
end

println("Max timestep: $(maximum(timesteps1))")
println("Min timestep: $(minimum(timesteps1))")

initialTimeStep = dm * 2^10
boxes2, timesteps2, attemptsRecorder2 = PlotReACT(A, B, initialTimeStep, T, P₁, U, constraint, dm, STRATEGY)
shapes2, maxVal2, minVal2 = plotProjectedFlowpipe(boxes2, timesteps2, 0, dimToPlot; approx=true)
maxVal = max(maxVal2, maxVal)
minVal = min(minVal2, minVal)
println(minVal, " ", maxVal)


for i in eachindex(shapes2)
    plot!(p, shapes2[i], vars=(1, 0), c=c1, la=0.1, alpha=1.0, lw=0.05,
        label=i == 1 ? L"Alg.\: 2" : "")
end


println("Finished simulations")

ylims!((minVal, maxVal))
yticks!([minVal, 0, constraint[1].b], [string(round(minVal; sigdigits=2)), "0.0", string(constraint[1].b)])

plot!(LazySets.HalfSpace([0.0, -1.0], -constraint[1].b), lab=L"\mathcal{X}_\bot", c=:black, fillstyle=:/)
xlims!(0, maximum(T))

savefig(p, "plots/" * name * ".pdf")
plot(p)