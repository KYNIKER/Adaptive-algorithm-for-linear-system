using  ReachabilityAnalysis

include("./fom_model.jl")
include("./fom_specifications.jl")

function load_fom()
    A, B, U = fom_model()
    X0, time_horizon, constraint = fom_specification()
    T = [0, time_horizon]
    dimToPlot = 1
    X0 = convert(Zonotope, X0)
    X0 = Zonotope(Vector(X0.center), Matrix(X0.generators))
    InputZonotope :: Zonotope = U#box_approximation(B*U) 

    return A, B, InputZonotope, X0, T, [constraint], dimToPlot
end