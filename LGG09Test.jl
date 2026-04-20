using ReachabilityAnalysis, Plots, LazySets, BenchmarkTools, CSV, DataFrames, MathematicalPredicates

include("models/heat/heat_load.jl")
include("models/motor/motor_load.jl")
include("models/building/building_load.jl")
include("models/PDE/pde_load.jl")
include("models/ISS/iss_load.jl")
include("models/beam/beam_load.jl")
include("models/MNA1/mna1_load.jl")

_, _, _, _, _, ISSproperty, _ = load_iss()
_, _, _, _, _, pdeProperty, _ = load_pde()

algCheckDict = Dict(
    "beam" => LGG09(δ=5e-5, template=CustomDirections([sparsevec([89], [1.0], 348)]), approx_model=Forward()),
    "building" => LGG09(δ=2e-3, template=CustomDirections([sparsevec([25], [1.0], 48)]), approx_model=Forward()),
    "heat" => LGG09(δ=1e-3, template=CustomDirections([sparsevec([133], [1.0], 200)]), approx_model=Forward()),
    "iss" => LGG09(δ=6e-4, template=CustomDirections([ISSproperty[1].a]), approx_model=Forward()),
    "motor" => LGG09(δ=1e-3, template=CustomDirections([sparsevec([1], [1.0], 8), sparsevec([5], [1.0], 8)]), approx_model=Forward()),
    "mna1" => LGG09(δ=4e-4, template=CustomDirections([sparsevec([1], [1.0], 578)]), approx_model=Forward()),
    "pde" => LGG09(δ=3e-4, template=CustomDirections([pdeProperty[1].a]), approx_model=Forward())
)


function RunCodeLGG09(prob, alg, t, constraint, isSparse=true)
    sol = solve(prob, alg; T=t) # Running the actual time
    return mapreduce(c -> ρ(c.a, sol) <= c.b, &, constraint) # Check if hits constraint
    if isSparse # Usually true
        println("sparse")
        return mapreduce(c -> ρ(c.a, sol) <= c.b, &, constraint) # Check if hits constraint
    end
    # Typical array in case of low dimensions
    # Manually fetch indexes for each constraint
    return mapreduce(c -> ρ(findall(c.a .== 1), sol.F) <= c.b, &, constraint) # Check if hits constraint
end




# read the docs https://juliaci.github.io/BenchmarkTools.jl/stable/manual/
function doLGG09JuliaTest(load_func, name)
    println("Running benchmark for: ", name, "...")
    BenchmarkTools.DEFAULT_PARAMETERS.seconds = 3600
    BenchmarkTools.DEFAULT_PARAMETERS.samples = 10

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
    if _name == "mna1"
        sys = @system(x' = Ax + B, x ∈ Universe(n))
    else
        sys = @system(x' = Ax + Bu, x ∈ Universe(n), u ∈ ballβ)
    end

    prob = InitialValueProblem(sys, P₁)
    isSparse = constraint[1].a isa SparseVector

    b = @benchmarkable _ = RunCodeLGG09($prob, $algCheck, $t, $constraint, $isSparse)
    y = run(b)

    # Check if it works
    println("Procceeding to data collection...")
    res = RunCodeLGG09(prob, algCheck, t, constraint, isSparse)

    # Convert time to seconds from nanoseconds
    timeList = []
    for timeVal in y.times
        push!(timeList, timeVal / 1e9)
    end

    df = DataFrame(name=name, avgTime=mean(timeList), medianTime=median(timeList), success=res)
    filename = "results/LGG09_" * name * "Results" * ".csv"
    if isfile(filename) # Check if file exists
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

names = ["beam", "building", "heat", "iss", "motor", "pde", "mna1"]
loadFuncs = [load_beam, load_building, load_heat_input, load_iss, load_motor, load_pde, load_mna1]

for (name, load_func) in zip(names, loadFuncs)
    doLGG09JuliaTest(load_func, name)
end