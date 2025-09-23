

function loadTest() # Maybe redo in the future with just 1 parameter
    A = [-1. 0.; 
         0. -1.] #
    P = Zonotope([0., 1.5], [0.0 0.0; 0.0 0.5])
    return A, P
end

function loadCosWave(tΔ)
    A = [cos(tΔ) sin(tΔ); 
        -sin(tΔ) cos(tΔ)]
    P = Zonotope([0., 0.0], [0.0 0.0; 0.0 1])
    return A, P
end

