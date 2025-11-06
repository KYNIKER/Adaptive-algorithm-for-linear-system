using Plots, LazySets, LinearAlgebra
include("helperfunctions.jl")
include("models.jl")
include("CegarFunctions.jl")

A = [0.0 1.0; -1.0 0.0]
ϕ = [1.0 0.0; 0.0 1.0]
@profview forwardTime(A, Zonotope([0,0], diagm([1, 1])), 10, 1, ϕ,0, Zonotope([0,0], diagm([1, 1])), Zonotope([0,0],diagm([1, 1])))