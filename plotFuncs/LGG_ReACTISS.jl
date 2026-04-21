using Plots, LazySets, LinearAlgebra, BenchmarkTools, ReachabilityAnalysis, LaTeXStrings, Plots.PlotMeasures

include("../models/heat/heat_load.jl")
include("../models/motor/motor_load.jl")
include("../models/building/building_load.jl")
include("../models/PDE/pde_load.jl")
include("../models/ISS/iss_load.jl")
include("../models/beam/beam_load.jl")
include("../models/MNA1/mna1_load.jl")
include("plotHelper.jl")
include("../ReACT.jl")

name = "LGGvsReACTISSSupport"
load_func = load_iss
A, B, ballβ, P₁, T, constraint, dimToPlot = load_func()
palette = Plots.palette(:fes10)
c1 = palette[9]
c2 = palette[6]
alp = 0.7
# ReACT
Digits = 6e-4
initialTimeStep = (2.0)^5 * Digits
STRATEGY = 2


boxes1, timesteps1 = PlotReACT(A, B, initialTimeStep, T, P₁, ballβ, constraint, Digits, [constraint[1].a, constraint[2].a], STRATEGY)

shapes1, maxVal1, minVal1 = plotSupportFlowpipe(boxes1, timesteps1, 1, 2)
# LGG
δ = 0.00
tVal = maximum(T)
n = size(A, 1)

sys = @system(x' = Ax + Bu, x ∈ Universe(n), u ∈ ballβ)
prob = InitialValueProblem(sys, P₁)

sol = solve(prob; T=tVal,
    alg=LGG09(δ=6e-4, template=CustomDirections([constraint[1].a, constraint[2].a]), approx_model=Forward()))
println("Ran models")

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
    xlims=(0, tVal), xlabel=L"Time", ylabel=L"y")



constraintValAdjusted = constraint[1].b * 1.1
maxVal = constraintValAdjusted
minVal = -constraintValAdjusted

ylims!((minVal, maxVal))
yticks!([-constraint[1].b, 0, constraint[1].b], [string(-constraint[1].b), "0.0", string(constraint[1].b)])


for i in eachindex(shapes1)
    plot!(p, shapes1[i], color=c1, c=c1, la=0.0, alpha=0.7, lw=0.05,
        label=i == 1 ? L"Alg.\: 3: \delta^{+} / \delta^- = %$initialTimeStep / %$Digits" : "")
end
for (i, rp) in pairs(flowpipe(sol))
    el1 = ρ(constraint[1].a, rp)
    el2 = -ρ(constraint[2].a, rp)

    plot!(p, Shape([Digits * (i - 1), Digits * (i), Digits * (i), Digits * (i - 1)], [el1, el1, el2, el2]), color=c2, c=c2, la=0.0, alpha=1.0, lw=0.05,
        label=i == 1 ? L"LGG" : "")
end

plot!(LazySets.HalfSpace([0.0, -1.0], -constraint[1].b), lab="Unsafe Region", c=:black, fillstyle=:/)
plot!(LazySets.HalfSpace([0.0, 1.0], -constraint[1].b), c=:black, fillstyle=:/)
xlims!((0, tVal))
println(minVal, " ", maxVal)
savefig(p, "plots/" * name * "Plot.pdf")
savefig(p, "plots/" * name * "Plot.png")
display(p)
