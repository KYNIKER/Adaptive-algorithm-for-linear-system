# Based on the paper JuliaReach: a Toolbox for Set-Based Reachability
using Plots, LazySets, LinearAlgebra, BenchmarkTools, FastExpm, Profile, PProf
include("helperfunctions.jl")
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
#include("CegarFunctions.jl")
include("CegarInhomogenous.jl")

const μ = 0.001
const STRATEGY = 1

initialTimeStep = 0.0005
#strategy = 1
Digits = 5
reuse = true
plotConstraint = true
input = true
plotOutput = false

