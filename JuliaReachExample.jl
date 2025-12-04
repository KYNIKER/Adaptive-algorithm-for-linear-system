using ReachabilityAnalysis, Plots

include("models.jl")
include("models/heat/heat_load.jl")
include("models/motor/motor_load.jl")
include("models/motor/motor_load.jl")
include("models/building/building_load.jl")
include("models/PDE/pde_load.jl")
include("models/ISS/iss_load.jl")
include("models/beam/beam_load.jl")
include("models/FOM/fom_load.jl")
include("models/MNA1/mna1_load.jl")
include("models/MNA5/mna5_load.jl")

# Timestep size
δ = 0.002

# System description
A, B, ballβ, P₁, time, constraint, dimToPlot = load_building()
t = maximum(time)
n = size(A, 1)
println(n)
sys = @system(x' = Ax + Bu, x ∈ Universe(n), u ∈ ballβ)

prob = InitialValueProblem(sys, P₁)


# flowpipe computation
alg = GLGM06(;δ=δ, approx_model=Forward())
@time sol = solve(prob, alg; T=t)
@time res = mapreduce(c -> ρ(c.a, sol) <= c.b, &, constraint)
println(res)

#LazySets.set_ztol(Float64, 1e-8)
#fig = plot(sol; vars=(0, dimToPlot), linecolor=:blue, color=:blue,
#           alpha=0.8, lw=1.0, xlab="t", ylab=string(dimToPlot))
#plot(sol, vars=(0, 1), xlab="t", ylab="x(t)")