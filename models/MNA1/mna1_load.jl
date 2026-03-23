using ReachabilityAnalysis, LinearAlgebra

include("./mna1_model.jl")
include("./mna1_specifications.jl")

function load_mna1()
    A, b = mna1_model()
    X0, time_horizon, constraint = mna1_specification()
    T = [0, time_horizon]
    dimToPlot = 1
    X0 = convert(Zonotope, X0)
    X0 = Zonotope(Vector(X0.center), Matrix(X0.generators))

    U = Zonotope(ones(length(b)), [zeros(length(b))])

    return A, diagm(b), U, X0, T, [constraint], dimToPlot
end