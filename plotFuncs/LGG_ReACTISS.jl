# Based on the paper JuliaReach: a Toolbox for Set-Based Reachability
using Plots, LazySets, LinearAlgebra, BenchmarkTools, Profile, PProf, ReachabilityAnalysis, LaTeXStrings, Plots.PlotMeasures, Polyhedra, CDDLib


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
#include("../CegarInhomogenous.jl")
include("plotHelper.jl")
include("../ReACT.jl")

name = "LGGvsReACTISSSupport"
load_func = load_iss
A, B, ballβ, P₁, T, constraint, dimToPlot = load_func()
#T = [0.0, 5.0]
palette = Plots.palette(:fes10)
c1 = palette[9]
c2 = palette[6]
LazySets.Comparison.set_tolerance(Float64)
LazySets.Comparison.set_ztol(Float64, 1e-10)
alp = 0.7
# ReACT
Digits = 6e-4 #2e-3 #
initialTimeStep = (2.0)^6 * Digits #(2.0)^10 * Digits #
STRATEGY = 2

boxes1, timesteps1, attemptsRecorder1 = PlotReACT(A, B, initialTimeStep, T, P₁, ballβ, constraint, Digits, STRATEGY)
#pop!(shapes1)

#shapes1, maxVal1, minVal1 = plotSupportFlowpipe(boxes1, timesteps1, constraint[1].a, -constraint[1].a)
# LGG
δ = 0.00
tVal = maximum(T)
n = size(A, 1)

sys = @system(x' = Ax + Bu, x ∈ Universe(n), u ∈ ballβ)
#alg = LGG09(δ=2e-3, template=CustomDirections([sparsevec([25], [1.0], 48)]), approx_model=Forward())
prob = InitialValueProblem(sys, P₁)

sol = solve(prob; T=tVal,
    #alg=LGG09(; δ=2e-3, template=CustomDirections([constraint[1].a]), approx_model=Forward())) #solve(prob, alg; T=20.0) # Running the actual time
    alg=LGG09(δ=6e-4, template=CustomDirections([constraint[1].a, constraint[2].a]), approx_model=Forward()))
println("Ran models")
#=
sol = solve(prob; T=tVal,
    alg=BFFPSV18(δ=2e-3, vars=[25], partition=[i:i for i in 1:48]))
=#
#sol.options[:plot_vars] = [0, 25]
#solution_proj = LazySets.project(sol, [dimToPlot])
#res = mapreduce(c -> ρ(c.a, sol) <= c.b, &, constraint) # Check if hits constraint
p = plot(dpi=1200, thickness_scaling=1, guidefontsize=25, minorgrid=false,
    legendfont=font(12, "Times"),
    #legendcolumn=-1,
    #legend_position=:outertop,
    tickfont=font(8, "Times"),
    xguidefont=font(12, "Times"),
    yguidefont=font(12, "Times"),
    xtick=([0, tVal], [L"0", L"T"]),
    #ytick=([], []),
    bottom_margin=2mm,
    left_margin=5mm,
    right_margin=5mm,
    top_margin=2mm,
    #ylims=(minVal1, maxVal1), 
    xlims=(0, tVal), xlabel=L"Time", ylabel=L"y")



constraintValAdjusted = constraint[1].b * 1.1
maxVal = constraintValAdjusted #max(maxVal1, constraintValAdjusted)
minVal = -constraintValAdjusted #min(minVal1, -constraintValAdjusted)

#p = plot(dpi=300, thickness_scaling=1, xlims=(0, maximum(T)), xlabel="Time", ylabel="Value")

ylims!((minVal, maxVal))
yticks!([-constraint[1].b, 0, constraint[1].b], [string(-constraint[1].b), "0.0", string(constraint[1].b)])
#=
for i in eachindex(shapes1)
    plot!(p, shapes1[i], color=c1, c=c1, la=0.9, alpha=1.0, lw=0.05, label=i == 1 ? L"Alg.\: 3: \delta^{+} / \delta^- = %$initialTimeStep / %$Digits" : "")
end
=#
#plot!(p, flowpipe(solution_proj)[1], vars=(0, dimToPlot), color=c2, c=c2, la=0.0, alpha=1.0, lw=0.0, lab=L"LGG")
#LGGprojection = [Shape([Digits * (i - 1), Digits * (i), Digits * (i), Digits * (i - 1)], [elem[1], elem[1], elem[2], elem[2]]) for (i, elem) in pairs([map(x -> [ρ(constraint[1].a, x), ρ(constraint[2].a, x)], flowpipe(sol))])]
#plot!(p, sol, vars=(0, 1), color=c2, c=c2, la=0.0, alpha=0.7, lw=0.0)
#println(flowpipe(sol)[end])

for (i, rp) in pairs(flowpipe(sol))
    el1 = ρ(constraint[1].a, rp)
    el2 = -ρ(constraint[2].a, rp)

    plot!(p, Shape([Digits * (i - 1), Digits * (i), Digits * (i), Digits * (i - 1)], [el1, el1, el2, el2]), color=c2, c=c2, la=0.0, alpha=0.7, lw=0.05,
        label=i == 1 ? L"LGG" : "")
end
println("plotted lgg")
time = 0.0
for (i, rp) in pairs(boxes1)
    t = timesteps1[i]
    el1 = ρ(constraint[1].a, rp)
    el2 = -ρ(constraint[2].a, rp)
    plot!(p, Shape([time, time + t, time + t, time], [el1, el1, el2, el2]), color=c1, c=c1, la=0.0, alpha=0.7, lw=0.05, label=i == 1 ? L"Alg.\: 3: \delta^{+} / \delta^- = %$initialTimeStep / %$Digits" : "")
    global time += t
end
println("plotted ours")
plot!(LazySets.HalfSpace([0.0, -1.0], -constraint[1].b), lab="Unsafe Region", c=:black, fillstyle=:/)
plot!(LazySets.HalfSpace([0.0, 1.0], -constraint[1].b), c=:black, fillstyle=:/)
xlims!((0, tVal))
println(minVal, " ", maxVal)
savefig(p, "plots/" * name * "Plot.pdf")
plot(p)
