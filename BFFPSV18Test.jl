using ReachabilityAnalysis, Plots, LazySets, BenchmarkTools, CSV, DataFrames, MathematicalPredicates

include("models/heat/heat_load.jl")
include("models/motor/motor_load.jl")
include("models/building/building_load.jl")
include("models/PDE/pde_load.jl")
include("models/ISS/iss_load.jl")
include("models/beam/beam_load.jl")
include("models/MNA1/mna1_load.jl")


algCheckDict = Dict(
    "beam" => BFFPSV18(δ=5e-5, vars=[89], partition=[i:i for i in 1:348]),
    "building" => BFFPSV18(δ=2e-3, vars=[25], partition=[i:i for i in 1:48]),
    "heat" => BFFPSV18(δ=1e-3, vars=[133], partition=[i:i for i in 1:200]),
    "iss" => BFFPSV18(δ=6e-4, vars=136:270, partition=vcat([i:i for i in 1:135], [136:270])),
    "motor" => BFFPSV18(δ=1e-3, vars=[1, 5], partition=[i:i for i in 1:8]),
    "mna1" => BFFPSV18(δ=4e-4, vars=[1], partition=[i:i for i in 1:578]),
    "pde" => BFFPSV18(δ=3e-4, vars=1:84, partition=[i:i for i in 1:84])
)

function RunCodeISS(prob, alg, t, constraint)
    sol = solve(prob, alg; T=t) # Running the actual time
    mapping = constraint[1].a[136:270]
    return ρ(mapping, sol.F) <= constraint[1].b
end

function RunCodePDE(prob, alg, t, constraint)
    sol = solve(prob, alg; T=t) # Running the actual time
    mapping = constraint[1].a[1:84]
    return ρ(mapping, sol.F) <= constraint[1].b
end

function RunCodeBFFPSV18(prob, alg, t, constraint)
    sol = solve(prob, alg; T=t) # Running the actual time
    constraintAmount = length(constraint)
    flag = true
    for i in 1:constraintAmount
        oneHotEncoding = zeros(constraintAmount)
        oneHotEncoding[i] = 1.0
        flag &= ρ(oneHotEncoding, sol.F) <= constraint[i].b
    end
    return flag
end

# read the docs https://juliaci.github.io/BenchmarkTools.jl/stable/manual/
function doBFFPSV18JuliaTest(load_func, name)
    println("Running benchmark for: ", name, "...")
    BenchmarkTools.DEFAULT_PARAMETERS.seconds = 3600
    BenchmarkTools.DEFAULT_PARAMETERS.samples = 1  # reduced from 10 for faster evaluation

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
        b = @view B[diagind(B)]
        sys = @system(x' = Ax + b, x ∈ Universe(n))
    else
        sys = @system(x' = Ax + Bu, x ∈ Universe(n), u ∈ ballβ)
    end

    prob = InitialValueProblem(sys, P₁)

    b = @benchmarkable _ = RunCodeBFFPSV18($prob, $algCheck, $t, $constraint)
    if _name == "iss"
        b = @benchmarkable _ = RunCodeISS($prob, $algCheck, $t, $constraint)
    elseif _name == "pde"
        b = @benchmarkable _ = RunCodePDE($prob, $algCheck, $t, $constraint)
    end
    y = run(b)

    # Check if it works
    println("Procceeding to data collection...")
    if _name == "iss"
        res = RunCodeISS(prob, algCheck, t, constraint)
    elseif _name == "pde"
        res = RunCodePDE(prob, algCheck, t, constraint)
    else
        res = RunCodeBFFPSV18(prob, algCheck, t, constraint)
    end
    # Convert time to seconds from nanoseconds
    timeList = []
    for timeVal in y.times
        push!(timeList, timeVal / 1e9)
    end

    df = DataFrame(name=name, avgTime=mean(timeList), medianTime=median(timeList), success=res)
    filename = "results/BFFPSV_" * name * "JuliaResults" * ".csv"
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

names = ["beam", "building", "heat", "iss", "motor", "pde", "mna1"]
loadFuncs = [load_beam, load_building, load_heat_input, load_iss, load_motor, load_pde, load_mna1]

for (name, load_func) in zip(names, loadFuncs)
    doBFFPSV18JuliaTest(load_func, name)
end

