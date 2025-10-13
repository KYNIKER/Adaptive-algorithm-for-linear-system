using Plots, LazySets, LinearAlgebra
using MAT
# https://juliaio.github.io/MAT.jl/stable/methods/




function loadCosWave()
    A = [0. 1.; 
           -2.5 0.]
    P = Zonotope([0., 1.5], [[0.0; 0.05]])
    return (A, P)
end

function loadBuilding()
    # TODO This currently does not work 
    
    filedict = matread("MAT_Files/building.mat")
    A = filedict["A"] # This is a sparse matrix
    
    # Input matrix
    inputMatrix = filedict["B"]
    # input domain
    U = BallInf([0.5], 0.3)

    X0 = Hyperrectangle(; low=[fill(0.0002, 10); zeros(14); -0.0001; zeros(23)], high=[fill(0.00025, 10); zeros(14); 0.0001; zeros(23)])
    X0 = convert(Zonotope, X0)

    return (A, X0)
end

# Recommended T = 15
function loadCrane() # From https://github.com/JuliaReach/ReachabilityModels.jl/blob/master/src/models/crane/
    A = [0.0 1.0 0.0 0.0 0.0 0.0;
        -0.417533 -3.1931759963 39.24 0.0 -14.825331 11.123344;
        0.0 0.0 0.0 1.0 0.0 0.0;
        0.0417533 0.31931759963 -4.905 0.0 1.4825331 -1.1123344;
        0.0638407957 -0.32473339016573 0.0 0.0 -3.7332068901 -0.7007592976;
        0.0853437452 -0.72366802635628 0.0 0.0 -5.9714023436 -2.2736115136]

    X0 = Hyperrectangle([2.5, 0.0, 0.0, 0.0, 0.0, 0.0],
      [2.5, 0.0, 0.2, 0.1, 0.0, 0.0])

    X0 = convert(Zonotope, X0)
    return (A, X0)
end

# T = 15
# https://github.com/JuliaReach/ReachabilityModels.jl/blob/master/src/models/ellipse
function loadEclipse()
    A = [3.0 -9.0; 4.0 -3.0]
    
    X0₁ = Singleton([1.0, 0.0]) ⊕ BallInf(zeros(2), 0.01)
    X0₂ = Singleton([1.5, 0.0]) ⊕ BallInf(zeros(2), 0.02)

    α = 0.8
    Z0 = α * X0₁ + (1 - α) * X0₂;

    convert(Zonotope, Z0)

    return (A, Z0)
end

# T = 5
# https://github.com/JuliaReach/ReachabilityModels.jl/blob/master/src/models/five_dim_sys
function fiveDimSys()
    D = [-1.0 -4.0 0.0 0.0 0.0;
     4.0 -1.0 0.0 0.0 0.0;
     0.0 0.0 -3.0 1.0 0.0;
     0.0 0.0 -1.0 -3.0 0.0;
     0.0 0.0 0.0 0.0 -2.0]
    P = [0.6 -0.1 0.1 0.7 -0.2;
        -0.5 0.7 -0.1 -0.8 0.0;
        0.9 -0.5 0.3 -0.6 0.1;
        0.5 -0.7 0.5 0.6 0.3;
        0.8 0.7 0.6 -0.3 0.2]
    A = P * D * inv(P)

    X0 = BallInf([1.0, 0.0, 0.0, 0.0, 0.0], 0.1)

    Z0 = convert(Zonotope, X0)

    return (A, Z0)
end



