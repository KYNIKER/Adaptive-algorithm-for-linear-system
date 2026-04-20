using Plots, LazySets, LinearAlgebra, Plots.PlotMeasures, LaTeXStrings

include("../models/heat/heat_load.jl")
include("../models/motor/motor_load.jl")
include("../models/building/building_load.jl")
include("../models/PDE/pde_load.jl")
include("../models/ISS/iss_load.jl")
include("../models/beam/beam_load.jl")
include("../models/MNA1/mna1_load.jl")
include("../ReACT.jl")
include("plotHelper.jl")


# We load a simple coswave
STRATEGY = 1

initialTimeStep = 0.8
Digits = 1e-1

palette = Plots.palette(:fes10)
alp = 0.7

A = [0. 1.;
    -2.5 0.]
P₁ = Zonotope([0., 1.5], [[0.0; 0.05]])
constraint = LazySets.HalfSpace(-[0., 1.], 1.65)
T = [0, 8]
dimToPlot = 1
# No input
B = diagm([1.0, 1.0])
U = Zonotope(zeros(LazySets.dim(P₁)), zeros(LazySets.dim(P₁), 0))


constraint = isa(constraint, Array) ? constraint : [constraint]

boxes1, timesteps1 = PlotReACT(A, B, (2.0)^2 * 1e-1, T, P₁, U, constraint, 1e-1, [[0., 1.], [0., -1.]], STRATEGY)
shapes1, maxVal1, minVal1 = plotSupportFlowpipe(boxes1, timesteps1, 1, 2)

# Get new values



println("Starting second simulation with timestep size: ", initialTimeStep)

boxes2, timesteps2 = PlotReACT(A, B, 1e-1, T, P₁, U, constraint, 1e-1, [[0., 1.], [0., -1.]], STRATEGY)
shapes2, maxVal2, minVal2 = plotSupportFlowpipe(boxes2, timesteps2, 1, 2)

println("Finished simulations")

constraintValAdjusted = -constraint[1].b * 1.2
maxVal = max(maxVal1, maxVal2, constraintValAdjusted)
minVal = min(minVal1, minVal2, constraintValAdjusted)


p = plot(dpi=1200, thickness_scaling=1, guidefontsize=25, minorgrid=true,
    legendfont=font(12, "Times"),
    #legendcolumn=-1,
    legend_position=:topright,
    tickfont=font(8, "Times"),
    xguidefont=font(12, "Times"),
    yguidefont=font(12, "Times"),
    xtick=([0, 8], [L"0", L"T"]),
    ytick=([], []),
    bottom_margin=2mm,
    left_margin=5mm,
    right_margin=5mm,
    top_margin=2mm,
    ylims=(minVal, maxVal), xlims=(0, maximum(T)), xlabel=L"Time", ylabel=L"x")



for i in eachindex(shapes1)
    if i == 1
        plot!(p, shapes1[i], vars=(1, 0), c=palette[9], alpha=0.7, lw=0.05,
            label="Our approach")
    else
        plot!(p, shapes1[i], vars=(1, 0), c=palette[9], alpha=0.7, lw=0.05,
            label="")
    end
end
for i in eachindex(shapes2)
    if i == 1
        plot!(p, shapes2[i], vars=(1, 0), c=palette[6], alpha=1.0, la=0.0, lw=0.05, #fa=0.0,
            label="Fixed step")
    else
        plot!(p, shapes2[i], vars=(1, 0), c=palette[6], alpha=1.0, la=0.0, lw=0.05, #fa=0.0,
            label="")
    end
end

plot!(LazySets.HalfSpace(-constraint[1].a, -constraint[1].b), lab="Unsafe region", c=:black, alpha=1.0, fillstyle=:/)


savefig(p, "plots/" * "Introduction.pdf")
savefig(p, "plots/" * "Introduction.png")
display(p)