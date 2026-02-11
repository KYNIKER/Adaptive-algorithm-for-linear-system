using ReachabilityAnalysis, Plots, LazySets, BenchmarkTools, CSV, DataFrames, MathematicalPredicates

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



algCheckDict = Dict(
    "beam" => BFFPSV18(δ=5e-5, vars=[89], partition=[i:i for i in 1:348]),
    "building" => BFFPSV18(δ=2e-3, vars=[25], partition=[i:i for i in 1:48]),
    "heat" => BFFPSV18(δ=1e-3, vars=[133], partition=[i:i for i in 1:200]),

    #"iss" => BFFPSV18(δ = 6e-4, vars = 136:270, partition = vcat([[i] for i in 1:135], [136:270]), sparse = true, lazy_initial_set = false)
    "iss" => BFFPSV18(δ=6e-4, vars=136:270, partition=vcat([i:i for i in 1:135], [136:270])),
    #"iss" => BFFPSV18(; δ = 6e-4, vars = 136:270, partition = [1:135, 136:270])

    "motor" => BFFPSV18(δ=1e-3, vars=[1, 5], partition=[i:i for i in 1:8]),
    "mna1" => BFFPSV18(δ=4e-4, vars=[1], partition=[i:i for i in 1:578]),
    "mna5" => BFFPSV18(δ=3e-1, vars=[1, 2], partition=[i:i for i in 1:10913]),
    "pde" => BFFPSV18(δ=3e-4, vars=1:84, partition=[i:i for i in 1:84])
)

function RunCodeBFFPSV18(prob, alg, t, constraint, isSparse)
    sol = solve(prob, alg; T=t) # Running the actual time
    #println(sol.F)
    if isSparse # Usually true
        println("sparse")
        return mapreduce(c -> ρ(c.a.nzind, sol) <= c.b, &, constraint) # Check if hits constraint
    end
    # Typical array in case of low dimensions
    # Manually fetch indexes for each constraint
    return mapreduce(c -> ρ([0.0, 1.0], sol.F) <= c.b, &, constraint) # Check if hits constraint
end

# read the docs https://juliaci.github.io/BenchmarkTools.jl/stable/manual/
function doBFFPSV18JuliaTest(load_func, name)
    namePrint = "BFFPSV18_" * name
    println("Running benchmark for: ", namePrint, "...")
    BenchmarkTools.DEFAULT_PARAMETERS.seconds = 3600
    BenchmarkTools.DEFAULT_PARAMETERS.samples = 3

    GC.gc() # Force garbage collection
    A, B, ballβ, P₁, time, constraint, _ = load_func()
    t = maximum(time)
    n = size(A, 1)

    _name = lowercase(name)
    if haskey(algCheckDict, _name)
        algCheck = algCheckDict[_name]
    else
        println("No arguments found for ", name)
        return
    end

    sys = @system(x' = Ax + Bu, x ∈ Universe(n), u ∈ ballβ)

    prob = InitialValueProblem(sys, P₁)

    isSparse = constraint[1].a isa SparseVector
    b = @benchmarkable _ = RunCodeBFFPSV18($prob, $algCheck, $t, $constraint, $isSparse)
    y = run(b)

    # Check if it works
    println("Procceeding to data collection...")
    res = RunCodeBFFPSV18(prob, algCheck, t, constraint, isSparse)
    # Convert time to seconds from nanoseconds
    timeList = []
    for timeVal in y.times
        push!(timeList, timeVal / 1e9)
    end

    df = DataFrame(name=name, avgTime=mean(timeList), medianTime=median(timeList), success=res)
    filename = "results/" * namePrint * "JuliaResults" * ".csv"
    if isfile(filename)# Check if file exists
        open(filename, "a") do File
            CSV.write(File, df, delim=";", append=true)
        end
    else
        open(filename, "w") do File
            CSV.write(File, df, delim=";", writeheader=true)
        end
    end

    println("Finished running benchmark for: ", namePrint, "!")
end

names = ["beam", "building", "heat", "iss", "motor", "pde", "mna1", "mna5"]
loadFuncs = [load_beam, load_building, load_heat_input, load_iss, load_motor, load_pde, load_mna1, load_mna5]

for (name, load_func) in zip(names, loadFuncs)
    #doBFFPSV18JuliaTest(load_func, name)
end

#doBFFPSV18JuliaTest(load_building, "building")
doBFFPSV18JuliaTest(load_motor, "motor")
#doBFFPSV18JuliaTest(load_iss, "iss")


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