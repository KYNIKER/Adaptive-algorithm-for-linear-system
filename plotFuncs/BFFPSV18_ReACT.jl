using Plots, LazySets, LinearAlgebra, BenchmarkTools, Profile, PProf, ReachabilityAnalysis, LaTeXStrings, Plots.PlotMeasures

include("../helperfunctions.jl")
include("../models/heat/heat_load.jl")
include("../models/motor/motor_load.jl")
include("../models/building/building_load.jl")
include("../models/PDE/pde_load.jl")
include("../models/ISS/iss_load.jl")
include("../models/beam/beam_load.jl")
include("../models/MNA1/mna1_load.jl")
include("plotHelper.jl")
include("../ReACT.jl")

name = "BFFPSV18vsReACTBuilding"
load_func = load_heat_input
A, B, ballβ, P₁, T, constraint, dimToPlot = load_func()

palette = Plots.palette(:fes10)
c1 = palette[9]
c2 = palette[6]
LazySets.Comparison.set_tolerance(Float64)
LazySets.Comparison.set_ztol(Float64, 1e-10)
alp = 0.7
# ReACT
Digits = 1e-4
initialTimeStep = (2.0)^14 * 1e-4
STRATEGY = 1

boxes1, timesteps1, attemptsRecorder1 = PlotReACT(A, B, initialTimeStep, T, P₁, ballβ, constraint, Digits, STRATEGY)
shapes1, maxVal1, minVal1 = plotProjectedFlowpipe(boxes1, timesteps1, 0, dimToPlot; approx=true)
# BFFPSV18
δ = 0.00
tVal = maximum(T)
n = size(A, 1)

sys = @system(x' = Ax + Bu, x ∈ Universe(n), u ∈ ballβ)
prob = InitialValueProblem(sys, P₁)

sol = solve(prob; T=tVal,
    alg=BFFPSV18(δ=1e-3, vars=[133], partition=[i:i for i in 1:200]))

solution_proj = LazySets.project(sol, [dimToPlot])
p = plot(dpi=1200, thickness_scaling=1, guidefontsize=25, minorgrid=false,
    legendfont=font(12, "Times"),
    tickfont=font(8, "Times"),
    xguidefont=font(12, "Times"),
    yguidefont=font(12, "Times"),
    xtick=([0, tVal], [L"0", L"T"]),
    bottom_margin=2mm,
    left_margin=5mm,
    right_margin=5mm,
    top_margin=2mm,
    xlims=(0, tVal), xlabel=L"Time", ylabel=L"x_{%$dimToPlot}")


constraintValAdjusted = constraint[1].b * 1.1
maxVal = max(maxVal1, constraintValAdjusted)
minVal = min(minVal1, constraintValAdjusted)


ylims!((minVal, maxVal))
yticks!([minVal, 0, constraint[1].b], [string(round(minVal; sigdigits=2)), "0.0", string(constraint[1].b)])


plot!(p, flowpipe(solution_proj)[end], vars=(0, dimToPlot), color=c2, c=c2, la=0.0, alpha=1.0, lw=0.0, lab=L"BFFPSV")
plot!(p, solution_proj, vars=(0, dimToPlot), color=c2, c=c2, la=0.0, alpha=1.0, lw=0.0)
for i in eachindex(shapes1)
    plot!(p, shapes1[i], color=c1, c=c1, la=0.1, alpha=1.0, lw=0.01,
        label=i == 1 ? L"Alg.\: 3: \delta^{+} / \delta^- = %$initialTimeStep / %$Digits" : "")
end

plot!(LazySets.HalfSpace([0.0, -1.0], -constraint[1].b), lab="Unsafe Region", c=:black, fillstyle=:/)
xlims!((0, tVal))

savefig(p, "plots/" * name * "Plot.pdf")
plot(p)
