using Plots, LazySets, LinearAlgebra, ReachabilityAnalysis, LaTeXStrings, Plots.PlotMeasures
gr()

include("../models/heat/heat_load.jl")
include("../models/motor/motor_load.jl")
include("../models/building/building_load.jl")
include("../models/PDE/pde_load.jl")
include("../models/ISS/iss_load.jl")
include("../models/beam/beam_load.jl")
include("../models/MNA1/mna1_load.jl")
include("plotHelper.jl")
include("../ReACT.jl")

name = "LGGvsReACTBuilding"
load_func = load_building
A, B, ballβ, P₁, T, constraint, dimToPlot = load_func()
palette = Plots.palette(:fes10)
c1 = palette[9]
c2 = palette[6]
alp = 0.7
# ReACT
Digits = 2e-3
initialTimeStep = (2.0)^9 * 2e-3
STRATEGY = 1

boxes1, timesteps1 = PlotReACT(A, B, initialTimeStep, T, P₁, ballβ, constraint, Digits, [constraint[1].a, -constraint[1].a], STRATEGY)

shapes1, maxVal1, minVal1 = plotSupportFlowpipe(boxes1, timesteps1, 1, 2)
tVal = maximum(T)
n = size(A, 1)

println(maxVal1, minVal1)

# LGG

sys = @system(x' = Ax + Bu, x ∈ Universe(n), u ∈ ballβ)
alg = LGG09(δ=2e-3, template=CustomDirections([sparsevec([25], [1.0], 48)]), approx_model=Forward())
prob = InitialValueProblem(sys, P₁)

sol = solve(prob; T=tVal, alg=LGG09(; δ=0.002, vars=(25), n=48))
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

for i in eachindex(shapes1)
    plot!(p, shapes1[i], color=c1, c=c1, la=0.0, alpha=0.7, lw=0.05,
        label=i == 1 ? L"Alg.\: 3: \delta^{+} / \delta^- = %$initialTimeStep / %$Digits" : "")
end


plot!(p, flowpipe(solution_proj)[1], vars=(0, dimToPlot), color=c2, c=c2, la=0.0, alpha=1.0, lw=0.05, lab=L"LGG")
plot!(p, solution_proj, vars=(0, dimToPlot), color=c2, c=c2, la=0.0, alpha=1.0, lw=0.0)

plot!(LazySets.HalfSpace([0.0, -1.0], -constraint[1].b), lab="Unsafe Region", c=:black, fillstyle=:/)
xlims!((0, tVal))
lens!(p, [0.0, 0.6], [0.004, 0.0065], inset=(1, bbox(0.1, 0.7, 0.23, 0.23)), lc=:black, xtick=[], ytick=[], tickfont=font(20, "Times"), subplot=2)

savefig(p, "plots/" * name * "Plot.pdf")
savefig(p, "plots/" * name * "Plot.png")
display(p)
