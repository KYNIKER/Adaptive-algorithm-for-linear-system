using ReachabilityAnalysis

include("./mna5_model.jl")
include("./mna5_specifications.jl")

function load_mna5()
    A, b = mna5_model()
    X0, time_horizon, constraint = mna5_specification()
    T = [0, time_horizon]
    dimToPlot = 1
    X0 = convert(Zonotope, X0)
    X0 = Zonotope(Vector(X0.center), Matrix(X0.generators))

    return A, b, nothing, X0, T, constraint, dimToPlot
end