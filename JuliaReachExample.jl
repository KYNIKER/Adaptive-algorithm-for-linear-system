using ReachabilityAnalysis, Plots, LazySets, BenchmarkTools, CSV, DataFrames

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

function RunCodeLGG09(prob, alg, t, constraint)
    sol = solve(prob, alg; T=t) # Running the actual time
    res = mapreduce(c -> ρ(c.a, sol) <= c.b, &, constraint) # Check if hits constraint
    return res
end

# read the docs https://juliaci.github.io/BenchmarkTools.jl/stable/manual/
function doLGG09JuliaTest(load_func, name, δ=0.0006)
    name = "LGG09_" * name
    println("Running benchmark for: ", name, "...")
    BenchmarkTools.DEFAULT_PARAMETERS.seconds = 3600
    BenchmarkTools.DEFAULT_PARAMETERS.samples = 3

    GC.gc() # Force garbage collection
    A, B, ballβ, P₁, time, constraint, _ = load_func()
    t = maximum(time)
    n = size(A, 1)


    # Fetch constraint dims
    template = CustomDirections(map(x -> x.a, constraint))
    # idxs = map(x -> x.a.nzind, constraint) #constraint[1].a
    # idxs = vcat(idxs...)

    sys = @system(x' = Ax + Bu, x ∈ Universe(n), u ∈ ballβ)
    alg = LGG09(δ=δ, template=template, approx_model=Forward())
    prob = InitialValueProblem(sys, P₁)

    b = @benchmarkable _ = RunCodeLGG09($prob, $alg, $t, $constraint)
    #tune!(b) # Tune to find the optimal samples/evals
    y = run(b)

    # Check if it works
    println("Procceeding to data collection...")
    res = RunCodeLGG09(prob, alg, t, constraint)

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

function RunCodeBFFPSV18(prob, alg, t, constraint)
    sol = solve(prob, alg; T=t) # Running the actual time
    res = mapreduce(c -> ρ(c.a.nzind, sol) <= c.b, &, constraint) # Check if hits constraint
    return res
end

# read the docs https://juliaci.github.io/BenchmarkTools.jl/stable/manual/
function doBFFPSV18JuliaTest(load_func, name, δ=0.0006)
    name = "BFFPSV18_" * name
    println("Running benchmark for: ", name, "...")
    BenchmarkTools.DEFAULT_PARAMETERS.seconds = 3600
    BenchmarkTools.DEFAULT_PARAMETERS.samples = 3

    GC.gc() # Force garbage collection
    A, B, ballβ, P₁, time, constraint, _ = load_func()
    t = maximum(time)
    n = size(A, 1)

    # Fetch constraint dims
    idxs = map(x -> x.a.nzind, constraint) #constraint[1].a
    idxs = vcat(idxs...)

    amountOfConstraints = length(idxs)

    println("Contraint idxs: ", idxs)


    sys = @system(x' = Ax + Bu, x ∈ Universe(n), u ∈ ballβ)
    partition = [i:i for i in 1:n] 
    alg = BFFPSV18(;δ=δ, vars=sparsevec(idxs, ones(amountOfConstraints), n), dim=n, partition=partition)

    prob = InitialValueProblem(sys, P₁)


    b = @benchmarkable _ = RunCodeBFFPSV18($prob, $alg, $t, $constraint)
    y = run(b)

    # Check if it works
    println("Procceeding to data collection...")
    res = RunCodeBFFPSV18(prob, alg, t, constraint)

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

# Timestep size
#δ = 0.0001


#doBFFPSV18JuliaTest(load_beam, "beam", 0.0001)
#doBFFPSV18JuliaTest(load_heat_input, "heat")
doLGG09JuliaTest(load_heat_input, "heat")

# # System description
# A, B, ballβ, P₁, time, constraint, dimToPlot = load_beam()
# t = 20.0 #maximum(time)
# n = size(A, 1)
# println(n)
# sys = @system(x' = Ax + Bu, x ∈ Universe(n), u ∈ ballβ)

# prob = InitialValueProblem(sys, P₁)


# # flowpipe computation
# alg = GLGM06(;δ=δ, approx_model=Forward())
# #partition = [i:i for i in 1:348]
# #alg = BFFPSV18(;δ=δ, vars=sparsevec([89], [1.0], 348), dim=348, partition=partition)
# @time sol = solve(prob, alg; T=t, property=constraint, mode="check")
# @time res = mapreduce(c -> ρ([1], sol) <= c.b, &, constraint)
# println(res)

# LazySets.set_ztol(Float64, 1e-5)
# fig = plot(sol; vars=(0, 1), linecolor=:blue, color=:blue,
#            alpha=0.8, lw=1.0, xlab="t", ylab=string(dimToPlot))
# #plot(sol, vars=(0, 1), xlab="t", ylab="x(t)")