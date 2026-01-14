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
δplus = 0.01
DIGITS = [2, 3, 4, 5]


# read the docs https://juliaci.github.io/BenchmarkTools.jl/stable/manual/
function runAreaBenchmark(name, load_func)
    println("Running area benchmark for: ", name)


    # Actual run
    GC.gc()# Force garbage collection
    A, B, ballβ, P₁, _, _, _ = load_func() # load
    Pbounds = []
    Dbounds = []
    for digit in DIGITS
        P, disc = discing(A, B, P₁, ballβ, δplus, digit)
        push!(Pbounds, mapreduce(x -> abs.(x), +, eachcol(genmat(P))))
        push!(Dbounds, mapreduce(x -> abs.(x), +, eachcol(genmat(disc))))
    end

    PboundsRelative = map(x -> x .\ Pbounds[1], Pbounds)
    DboundsRelative = map(x -> x .\ Dbounds[1], Dbounds)

    # Write to csv file
    maxRelHom = map(x -> maximum(x), PboundsRelative)
    df = DataFrame(initialTimeStep=δplus, Digits=DIGITS, maxRelativeHomogenous=maxRelHom, maxRelativeInhomogenous=map(x -> maximum(x), DboundsRelative), minRelativeHomogenous=map(x -> minimum(x), PboundsRelative), minRelativeInhomogenous=map(x -> minimum(x), DboundsRelative))

    filename = "area/" * name * "Area" * ".csv"
    if isfile(filename)# Check if file exists
        open(filename, "a") do File
            CSV.write(File, df, delim=";", append=true; bufsize=500000000)
        end
    else
        open(filename, "w") do File
            CSV.write(File, df, delim=";", writeheader=true; bufsize=500000000)
        end
    end

    println("Completed run for: ", name)
end


#runBenchmark(modelname, 1.0, dig, model, 2)
#GC.gc()
runAreaBenchmark("building", load_building)
GC.gc()
runAreaBenchmark("beam", load_beam)
GC.gc()
runAreaBenchmark("heat", load_heat_input)
GC.gc()
runAreaBenchmark("motor", load_motor)
GC.gc()
runAreaBenchmark("pde", load_pde)
GC.gc()
runAreaBenchmark("ISS", load_iss)
GC.gc()
runAreaBenchmark("fom", load_fom)
GC.gc()
#=runAreaBenchmark("mna1", load_mna1)
GC.gc()
runAreaBenchmark("mna5", load_mna5)
GC.gc()=#
#runBenchmark(modelname, 0.0025, dig, model, 0)
#GC.gc()