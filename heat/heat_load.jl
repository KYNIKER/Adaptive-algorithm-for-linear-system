using ReachabilityBenchmarks, MathematicalSystems, Reachability, Plots

include(@current_path "heat_model.jl")
include(@current_path "heat_specifications.jl")

function load_heat_input()
    A, B, U = heat_model()
    X0, time_horizon, constraint = heat_specification()
    T = [0, time_horizon]
    dimToPlot = 133
    X0 = convert(Zonotope, X0)
    X0 = Zonotope(Vector(X0.center), Matrix(X0.generators))
    
    InputZonotope = B*Box(U) # No way det her virker

    return A, InputZonotope, X0, T, constraint, dimToPlot
end