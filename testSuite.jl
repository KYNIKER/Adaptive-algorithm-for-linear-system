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
function runBenchmark(name, initialTimeStep, Digits, load_func, STRATEGY)
    println("Running benchmark for: ", name)
    #BenchmarkTools.DEFAULT_PARAMETERS.seconds = 3600
    BenchmarkTools.DEFAULT_PARAMETERS.samples = 3

    # Actual run
    GC.gc()# Force garbage collection
    A, B, ballβ, P₁, T, constraint, _ = load_func() # load

    b = @benchmarkable _ = cegarInputSystemNoOutput($A, $B, $initialTimeStep, $T, $P₁, $ballβ, $constraint, $Digits, $STRATEGY)

    #tune!(b; verbose=true) # Tune to find the optimal samples/evals
    #println("tuned!")
    y = run(b; verbose=true)
    println("Run completed.")
    # Convert time to seconds from nanoseconds
    timeList = []
    for timeVal in y.times
        push!(timeList, timeVal / 1e9)
    end

    # Get timesteps
    timesteps = cegarInputSystemOnlyTiming(A, B, initialTimeStep, T, P₁, ballβ, constraint, Digits, STRATEGY)

    isSuccess = sum(timesteps, dims=1) >= T# Check if we reach the end
    uniqueTimesteps = unique(timesteps)
    # Write to csv file
    df = DataFrame(strategy=STRATEGY, initialTimeStep=initialTimeStep, Digits=Digits, avgTime=mean(timeList), medianTime=median(timeList), success=isSuccess, memory=y.memory, allocs=y.allocs, timesteps=[uniqueTimesteps])

    filename = "results/" * name * "ResultsNewMethod" * ".csv"
    if isfile(filename)# Check if file exists
        open(filename, "a") do File
            CSV.write(File, df, delim=";", append=true)
        end
    else
        open(filename, "w") do File
            CSV.write(File, df, delim=";", writeheader=true)
        end
    end

    println("Completed run for: ", name)
end

#runBenchmark("FOM", 0.5, 5, load_fom, 2)
#=modelname = "heat"
model = load_heat_input
dig = 1=#
runBenchmark("mna5", 1., 1, load_mna5, 1)
GC.gc()#=
runBenchmark("motor", 1.0, 3, load_motor, 1)
GC.gc()
runBenchmark("building", 1.0, 3, load_building, 1)
GC.gc()

runBenchmark("heatInput", 1.0, 4, load_heat_input, 1)
GC.gc()
runBenchmark("ISS", 0.05, 3, load_iss, 1)
GC.gc()
runBenchmark("beam", 0.0025, 4, load_beam, 1)
GC.gc()
runBenchmark("mna1", 1.0, 4, load_mna1, 1)
GC.gc()
runBenchmark("mna5", 1.0, 1, load_mna5, 1)
GC.gc()
#runBenchmark(modelname, 0.0025, dig, model, 0)
#GC.gc()=#