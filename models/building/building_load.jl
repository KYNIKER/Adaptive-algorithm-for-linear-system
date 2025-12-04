using ReachabilityAnalysis

include("./building_model.jl")
include("./building_specifications.jl")

function load_building()
    A, B, U = building_model()
    X0, time_horizon, constraint = building_specification()
    T = [0, time_horizon]
    dimToPlot = 25
    X0 = convert(Zonotope, X0)
    X0 = Zonotope(Vector(X0.center), Matrix(X0.generators))
    InputZonotope :: Zonotope = U#box_approximation(B*U) 

    return A, B, InputZonotope, X0, T, [constraint], dimToPlot
end