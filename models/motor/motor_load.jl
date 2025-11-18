using  ReachabilityAnalysis

include("./motor_model.jl")
include("./motor_specifications.jl")

function load_motor()
    A, B, U = motor_model()
    X0, time_horizon, constraints = motor_specification()
    T = [0, time_horizon]
    dimToPlot = 0#{0, 5}
    X0 = convert(Zonotope, X0)
    X0 = Zonotope(Vector(X0.center), Matrix(X0.generators))
    InputZonotope :: Zonotope = box_approximation(B*U) 

    return A, InputZonotope, X0, T, constraints, dimToPlot
end