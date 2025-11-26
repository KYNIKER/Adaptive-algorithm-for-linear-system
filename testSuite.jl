using Plots, LazySets, LinearAlgebra, BenchmarkTools
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


names = ["beam", "building", "fom", "heat", "iss", "motor", "pde"]
fileLoads = [load_beam, load_building, load_fom, load_heat_input, load_iss, load_motor, load_pde]
timestepList = [0.0016 * 2^10]


const STRATEGY = 0

#strategy = 1
Digits = 8


for (name, initialTimeStep, load_func) in zip(names, timestepList, fileLoads)

    println("Running benchmark for: ", name)

    # Actual run
    GC.gc()# Force garbage collection
    A, ballβ, P₁, T, constraint, dimToPlot = load_func()

    println("T: ", T)

    cegarInputSystem(A, initialTimeStep, T, P₁, ballβ, constraint, Digits)

    #@btime boxes2, timesteps, attemptsRecorder = cegarInputSystem(A, initialTimeStep, T, P₁, ballβ, constraint, Digits)


end




#GC.gc()# Force garbage collection


# Benchmark: @btime