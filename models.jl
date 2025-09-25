using Plots, LazySets, LinearAlgebra
using MAT
# https://juliaio.github.io/MAT.jl/stable/methods/


function loadTest() # Maybe redo in the future with just 1 parameter
    A = [-1. 0.; 
         0. -1.] #
    P = Zonotope([0., 1.5], [0.0 0.0; 0.0 0.5])
    return (A, P)
end

function loadCosWave(tΔ)
    A = [cos(tΔ) sin(tΔ); 
        -sin(tΔ) cos(tΔ)]
    P = Zonotope([0., 0.0], [0.0 0.0; 0.0 1])
    return (A, P)
end

function loadBuilding()
    
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

