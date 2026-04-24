# This file runs our testSuite. And also our plots if enabled

using Pkg
Pkg.activate(@__DIR__) # Select current repository
Pkg.instantiate() # Installs all packages from Project.toml

const RUN_PLOTS = true

println("Running Benchmarks")
include("testSuite.jl")
include("BFFPSV18Test.jl")
include("LGG09Test.jl")
println("Finished Benchmarks")

if RUN_PLOTS
    println("Starting plots")
    include("plotFuncs/BFFPSV18_ReACT.jl")
    include("plotFuncs/discCompPlot.jl")
    include("plotFuncs/introductionPlot.jl")
    include("plotFuncs/LGG_ReACT.jl")
    include("plotFuncs/LGG_ReACTISS.jl")
    include("plotFuncs/timeDiscretizationPlot.jl")
    println("Finished plots")
end