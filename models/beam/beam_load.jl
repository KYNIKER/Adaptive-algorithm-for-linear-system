using  ReachabilityAnalysis

include("./beam_model.jl")
include("./beam_specifications.jl")

function load_beam()
    A, B, U = beam_model()
    X0, time_horizon, constraint = beam_specification()
    T = [0, time_horizon]
    dimToPlot = 89
    X0 = convert(Zonotope, X0)
    X0 = Zonotope(Vector(X0.center), Matrix(X0.generators))
    InputZonotope :: Zonotope = box_approximation(B*U) 

    return A, InputZonotope, X0, T, [constraint], dimToPlot
end