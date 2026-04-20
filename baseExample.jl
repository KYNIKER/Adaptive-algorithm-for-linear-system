using Plots, LazySets, LinearAlgebra, BenchmarkTools, Profile
include("helperfunctions.jl")
include("plotFuncs/plotHelper.jl")
include("models/heat/heat_load.jl")
include("models/motor/motor_load.jl")
include("models/building/building_load.jl")
include("models/PDE/pde_load.jl")
include("models/ISS/iss_load.jl")
include("models/beam/beam_load.jl")
include("models/MNA1/mna1_load.jl")
include("ReACT.jl")

# Strategy 1 always tries increasing timestep
# Strategy 2 attempts to double timestep if no failure in past 4
STRATEGY = 2

δ⁻ = 0.002
initialTimeStep = δ⁻ * 2^9
plotConstraint = true
plotOutput = true


A, B, ballβ, P₁, T, constraint, dimToPlot = load_building()

constraint = isa(constraint, Array) ? constraint : [constraint]

###

res = ReACT(A, B, initialTimeStep, T, P₁, ballβ, constraint, δ⁻, STRATEGY)
println("Successful simulation?: ", res)

if plotOutput
    boxes, timesteps = PlotReACT(A, B, initialTimeStep, T, P₁, ballβ, constraint, δ⁻, [constraint[1].a, -constraint[1].a], STRATEGY)


    shapes, maxVal, minVal = plotSupportFlowpipe(boxes, timesteps, 1, 2)

    p = nothing
    if plotConstraint
        # Take min and max with constraint
        constraintValAdjusted = constraint[1].b * 1.1
        maxVal = max(maxVal, constraintValAdjusted)
        minVal = min(minVal, constraintValAdjusted)
    end

    p = plot(dpi=300, thickness_scaling=1, ylims=(minVal, maxVal))
    palette = Plots.palette(:fes10)
    c1 = palette[9]

    for i in eachindex(shapes)
        plot!(p, shapes[i], color=c1, c=c1, la=0.1, alpha=0.7, lw=0.05,
            label=i == 1 ? "Alg. 3: δ⁺ / δ⁻ = $initialTimeStep / $δ⁻" : "")
    end

    if plotConstraint
        if 0 < constraint[1].b
            plot!(LazySets.HalfSpace([0.0, -1.0], -constraint[1].b), lab="Unsafe Region", c=:black, fillstyle=:/)
        else
            plot!(LazySets.HalfSpace([0.0, 1.0], constraint[1].b), lab="Unsafe Region", c=:black, fillstyle=:/)
        end
    end

    println("Amount of timesteps: ", length(timesteps))
    println("Unique timesteps ", unique(timesteps))

    plot(p)
end