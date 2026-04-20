using Plots, LazySets, LinearAlgebra, BenchmarkTools, CSV, DataFrames, Expokit, ReachabilityAnalysis
include("models/heat/heat_load.jl")
include("models/motor/motor_load.jl")
include("models/building/building_load.jl")
include("models/PDE/pde_load.jl")
include("models/ISS/iss_load.jl")
include("models/beam/beam_load.jl")
include("models/MNA1/mna1_load.jl")
include("ReACT.jl")


# read the docs https://juliaci.github.io/BenchmarkTools.jl/stable/manual/
function runBenchmark(name, initialTimeStep, δ⁻, load_func, STRATEGY)
    LazySets.load_expokit()
    println("Running benchmark for: ", name)
    #BenchmarkTools.DEFAULT_PARAMETERS.seconds = 3600
    BenchmarkTools.DEFAULT_PARAMETERS.samples = 10

    # Actual run
    GC.gc()# Force garbage collection
    A, B, ballβ, P₁, T, constraint, _ = load_func() # load

    b = @benchmarkable _ = ReACT($A, $B, $initialTimeStep, $T, $P₁, $ballβ, $constraint, $δ⁻, $STRATEGY)

    y = run(b; verbose=true)
    println("Run completed.")

    # Convert time to seconds from nanoseconds
    timeList = []
    for timeVal in y.times
        push!(timeList, timeVal / 1e9)
    end


    isSuccess = ReACT(A, B, initialTimeStep, T, P₁, ballβ, constraint, δ⁻, STRATEGY) # Check if we reach the end

    # Write to csv file
    df = DataFrame(strategy=STRATEGY, initialTimeStep=initialTimeStep, δ⁻=δ⁻, avgTime=mean(timeList), medianTime=median(timeList), success=isSuccess, memory=y.memory, allocs=y.allocs)

    filename = "results/ReportReACTv5_" * name * "Results" * ".csv"
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


runBenchmark("ISS", (2.0)^6 * 6e-4, 6e-4, load_iss, 1)
GC.gc()
runBenchmark("beam", (2.0)^5 * 5e-5, 5e-5, load_beam, 1)
GC.gc()

runBenchmark("motor", (2.0)^4 * 1e-3, 1e-3, load_motor, 1)
GC.gc()

runBenchmark("pde", (2.0)^13 * 3e-4, 3e-4, load_pde, 1)
GC.gc()

runBenchmark("building", (2.0)^10 * 2e-3, 2e-3, load_building, 1)
GC.gc()

runBenchmark("heatInput", (2.0)^11 * 1e-3, 1e-3, load_heat_input, 1)
GC.gc()

runBenchmark("mna1", (2.0)^12 * 4e-4, 4e-4, load_mna1, 1)
GC.gc()
