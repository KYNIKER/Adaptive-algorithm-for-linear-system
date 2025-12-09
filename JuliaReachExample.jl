using ReachabilityAnalysis, Plots, LazySets

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
δ = 0.001

# System description
A, B, ballβ, P₁, time, constraint, dimToPlot = load_beam()
t = 14.0 #maximum(time)
n = size(A, 1)
println(n)
sys = @system(x' = Ax + Bu, x ∈ Universe(n), u ∈ ballβ)

prob = InitialValueProblem(sys, P₁)


# flowpipe computation
#alg = GLGM06(;δ=δ, approx_model=Forward())
partition = [i:i for i in 1:348]
alg = BFFPSV18(;δ=δ, vars=sparsevec([89], [1.0], 348), dim=348, partition=partition)
@time sol = solve(prob, alg; T=t, property=constraint, mode="check")
@time res = mapreduce(c -> ρ([1], sol) <= c.b, &, constraint)
println(res)

LazySets.set_ztol(Float64, 1e-5)
fig = plot(sol; vars=(0, 1), linecolor=:blue, color=:blue,
           alpha=0.8, lw=1.0, xlab="t", ylab=string(dimToPlot))
#plot(sol, vars=(0, 1), xlab="t", ylab="x(t)")