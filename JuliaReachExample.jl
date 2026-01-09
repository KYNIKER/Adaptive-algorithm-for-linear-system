using ReachabilityAnalysis, Plots, BenchmarkTools, CSV, DataFrames

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
#δ = 0.002

# # System description
# A, B, ballβ, P₁, time, constraint, dimToPlot = load_building()
# t = maximum(time)
# n = size(A, 1)
# println(n)
# sys = @system(x' = Ax + Bu, x ∈ Universe(n), u ∈ ballβ)

# prob = InitialValueProblem(sys, P₁)


# # flowpipe computation
# alg = GLGM06(;δ=δ, approx_model=Forward())
# @time sol = solve(prob, alg; T=t) # Running the actual time
# @time res = mapreduce(c -> ρ(c.a, sol) <= c.b, &, constraint) # Check if hits constraint
# println(res)



function RunCode(prob, alg, t, constraint)
    sol = solve(prob, alg; T=t) # Running the actual time
    res = mapreduce(c -> ρ(c.a, sol) <= c.b, &, constraint) # Check if hits constraint
    return res
end

function doGLGMJuliaTest(load_func, name)
    println("Running benchmark for: ", name, "...")
    BenchmarkTools.DEFAULT_PARAMETERS.seconds = 3600
    BenchmarkTools.DEFAULT_PARAMETERS.samples = 3
    δ = 0.0006
    GC.gc() # Force garbage collection
    A, B, ballβ, P₁, time, constraint, _ = load_func()
    t = maximum(time)
    n = size(A, 1)

    sys = @system(x' = Ax + Bu, x ∈ Universe(n), u ∈ ballβ)
    alg = GLGM06(; δ=δ, approx_model=Forward())
    prob = InitialValueProblem(sys, P₁)


    b = @benchmarkable _ = RunCode($prob, $alg, $t, $constraint)
    tune!(b) # Tune to find the optimal samples/evals
    y = run(b)

    # Check if it works
    println("Procceeding to data collection...")
    res = RunCode(prob, alg, t, constraint)

    # Convert time to seconds from nanoseconds
    timeList = []
    for timeVal in y.times
        push!(timeList, timeVal / 1e9)
    end

    df = DataFrame(name=name, timestepsize=δ, avgTime=mean(timeList), medianTime=median(timeList), success=res)
    filename = "results/" * name * "JuliaResults" * ".csv"
    if isfile(filename)# Check if file exists
        open(filename, "a") do File
            CSV.write(File, df, delim=";", append=true)
        end
    else
        open(filename, "w") do File
            CSV.write(File, df, delim=";", writeheader=true)
        end
    end

    println("Finished running benchmark for: ", name, "!")
end

doGLGMJuliaTest(load_iss, "iss")
