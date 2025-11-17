using  ReachabilityAnalysis

include("./heat_model.jl")
include("./heat_specifications.jl")

function load_heat_input()
    A, B, U = heat_model()
    X0, time_horizon, constraint = heat_specification()
    T = [0, 2]
    dimToPlot = 133
    X0 = convert(Zonotope, X0)
    X0 = Zonotope(Vector(X0.center), Matrix(X0.generators))
    
    InputZonotope :: Zonotope = B*U

    return A, InputZonotope, X0, T, [constraint], dimToPlot
end