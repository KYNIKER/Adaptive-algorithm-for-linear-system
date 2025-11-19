using  ReachabilityAnalysis

include("./pde_model.jl")
include("./pde_specifications.jl")

function load_pde()
    A, B, U = pde_model()
    X0, time_horizon, constraint = pde_specification()
    T = [0, time_horizon]
    dimToPlot = 1
    X0 = convert(Zonotope, X0)
    X0 = Zonotope(Vector(X0.center), Matrix(X0.generators))
    InputZonotope :: Zonotope = box_approximation(B*U) 

    return A, InputZonotope, X0, T, [constraint], dimToPlot
end 