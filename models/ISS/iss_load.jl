using ReachabilityAnalysis

include("./iss_model.jl")
include("./iss_specifications.jl")

function load_iss()
    A, B, U = iss_model()
    X0, time_horizon, constraint = iss_specification()
    T = [0, time_horizon]
    dimToPlot = 182
    X0 = convert(Zonotope, X0)
    X0 = Zonotope(Vector(X0.center), Matrix(X0.generators))
    InputZonotope::Zonotope = U

    return A, B, InputZonotope, X0, T, constraint, dimToPlot
end