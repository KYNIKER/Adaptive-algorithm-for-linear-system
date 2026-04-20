using LazySets, LinearAlgebra, BenchmarkTools, CSV, DataFrames, Expokit, ReachabilityAnalysis
include("models/heat/heat_load.jl")
include("models/motor/motor_load.jl")
include("models/building/building_load.jl")
include("models/PDE/pde_load.jl")
include("models/ISS/iss_load.jl")
include("models/beam/beam_load.jl")
include("models/MNA1/mna1_load.jl")
include("ReACT.jl")

strategy = 2

# read the docs https://juliaci.github.io/BenchmarkTools.jl/stable/manual/
function runBenchmark(name, initialTimeStep, δ⁻, load_func, STRATEGY)
    LazySets.load_expokit()
    println("Running benchmark for: ", name)
    #BenchmarkTools.DEFAULT_PARAMETERS.seconds = 3600
    BenchmarkTools.DEFAULT_PARAMETERS.samples = 10

    # Actual run
    GC.gc()# Force garbage collection
    A, B, ballβ, P₁, T, constraint, _ = load_func() # load

    local tsteps = ReACT(A, B, initialTimeStep, T, P₁, ballβ, constraint, δ⁻, STRATEGY) # Check if we reach the end
    b = @benchmarkable _ = ReACT($A, $B, $initialTimeStep, $T, $P₁, $ballβ, $constraint, $δ⁻, $STRATEGY)

    y = run(b; verbose=true)
    println("Run completed.")

    # Convert time to seconds from nanoseconds
    timeList = []
    for timeVal in y.times
        push!(timeList, timeVal / 1e9)
    end


    isSuccess = tsteps > 0
    # Write to csv file
    df = DataFrame(strategy=STRATEGY, initialTimeStep=initialTimeStep, δ⁻=δ⁻, avgTime=mean(timeList), medianTime=median(timeList), success=isSuccess, memory=y.memory, allocs=y.allocs, steps=tsteps)

    filename = "results/ReACT_" * name * "Results" * ".csv"
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


runBenchmark("ISS", (2.0)^5 * 6e-4, 6e-4, load_iss, strategy)
GC.gc()

runBenchmark("beam", (2.0)^5 * 5e-5, 5e-5, load_beam, strategy)
GC.gc()

runBenchmark("motor", (2.0)^3 * 1e-3, 1e-3, load_motor, strategy)
GC.gc()

runBenchmark("pde", (2.0)^10 * 3e-4, 3e-4, load_pde, strategy)
GC.gc()

runBenchmark("building", (2.0)^9 * 2e-3, 2e-3, load_building, strategy)
GC.gc()

runBenchmark("heatInput", (2.0)^10 * 1e-3, 1e-3, load_heat_input, strategy)
GC.gc()

runBenchmark("mna1", (2.0)^11 * 4e-4, 4e-4, load_mna1, strategy)
GC.gc()
