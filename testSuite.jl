using Plots, LazySets, LinearAlgebra, BenchmarkTools, CSV, DataFrames
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


#names = ["beam", "building", "fom", "heat", "iss", "motor", "pde"]
#fileLoads = [load_beam, load_building, load_fom, load_heat_input, load_iss, load_motor, load_pde]


# read the docs https://juliaci.github.io/BenchmarkTools.jl/stable/manual/
function runBenchmark(name, initialTimeStep, Digits, load_func, STRATEGY, maxOrder = 50, reduceOrder = 1)
    println("Running benchmark for: ", name)
    BenchmarkTools.DEFAULT_PARAMETERS.seconds = 3600
    BenchmarkTools.DEFAULT_PARAMETERS.samples = 10

    # Actual run
    GC.gc()# Force garbage collection
    A, B, ballβ, P₁, T, constraint, _ = load_func() # load

    b = @benchmarkable _ = cegarInputSystemNoOutput($A, $B, $initialTimeStep, $T, $P₁, $ballβ, $constraint, $Digits, $STRATEGY, $maxOrder, $reduceOrder)

    #tune!(b; verbose=true) # Tune to find the optimal samples/evals
    println("tuned!")
    y = run(b; verbose=true)

    # Convert time to seconds from nanoseconds
    timeList = []
    for timeVal in y.times
        push!(timeList, timeVal / 1e9)
    end

    # Get timesteps
    timesteps = cegarInputSystemOnlyTiming(A, B, initialTimeStep, T, P₁, ballβ, constraint, Digits, STRATEGY, maxOrder, reduceOrder)

    isSuccess = sum(timesteps, dims=1) >= T# Check if we reach the end
    uniqueTimesteps = unique(timesteps)
    # Write to csv file
    df = DataFrame(strategy = STRATEGY, initialTimeStep = initialTimeStep, Digits = Digits, maxOrder = maxOrder, reduceOrder = reduceOrder, avgTime = mean(timeList), medianTime = median(timeList), success = isSuccess, memory = y.memory, allocs = y.allocs, timesteps = [uniqueTimesteps])

    filename = "results/" * name * "Results" * ".csv"
    if isfile(filename)# Check if file exists
        open(filename, "a") do File
            CSV.write(File, df, delim = ";", append=true)
        end
    else
        open(filename, "w") do File
            CSV.write(File, df, delim = ";",writeheader = true)
        end
    end

    println("Completed run for: ", name)
end

names = ["building", "heat", "motor"]
fileLoads = [load_building, load_heat_input, load_motor]


for (name, modelFunc) in zip(names, fileLoads)
    #name = name * "_NoReduceHomogenous"
    for maxOrder in [20, 35, 50]
        for reduceOrder in [1, 2, 5]
            for initialTimeStep in [5, 8, 10, 20]
                runBenchmark(name, initialTimeStep, 5, modelFunc, 1, maxOrder, reduceOrder)
            end
        end
    end
end
# #runBenchmark("FOM", 0.5, 5, load_fom, 2)
# modelname = "mna5"
# model = load_mna5
# dig = 1
# #runBenchmark(modelname, 1.0, dig, model, 2)
# #GC.gc()
# runBenchmark(modelname, 1.0, dig, model, 1)
# GC.gc()
#runBenchmark(modelname, 0.0025, dig, model, 0)
#GC.gc()