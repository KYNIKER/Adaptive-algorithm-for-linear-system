using LazySets, LinearAlgebra

function PCA_reduce(Z, k) # Reduce to k generators

    c = Z.center # get center
    G = genmat(Z) # Get generators

    U, S, V = svd(G) # Apply singular value decomposition (https://www.geeksforgeeks.org/machine-learning/singular-value-decomposition-svd/)

    # Note that U finds the eigenvectors of the covariance matrix
    # S ranks importance, and is ordered equivalently to eigenvalues
    k = min(k, length(S))   # Ensure k is not greater than the amount of dimensions

    G_red = U[:, 1:k] * Diagonal(S[1:k]) # Get top k

    return Zonotope(c, G_red)
end

function box_reduce(Z::Zonotope, k::Number) # Reduce to k generators

    if k == 1 # If we only want one generator
        hyperRect = interval_hull(Z)
        return Zonotope(Z.center, [radius_hyperrectangle(hyperRect)])
    end

    c::Vector{Float64} = Z.center # get center
    G::Matrix{Float64} = genmat(Z) # Get generators

    amountOfGens = size(G, 2)
    dim = size(G, 1)
    k = min(k, amountOfGens) # Ensure k is not greater than the amount of generators

    com = map(x->norm(x, 1) - norm(x, Inf), eachcol(G))
    G = G[:,sortperm(com)]
    #G = sortslices(G, dims = 2, by = x->norm(x, 1) - norm(x, Inf))


    # We split the generators into two groups. Where G1 will be overapproximated
    G1 = G[:, 1:amountOfGens-(k-(dim-1))]
    G2 = G[:, amountOfGens-(k-dim):amountOfGens]

    Z1::Zonotope = Zonotope(zeros(dim), G1)
    Z2 = Zonotope(c, G2)

    Z1HyperRect = interval_hull(Z1) # Note that Z1 will always return something with dims 2. 
    Z1red = Zonotope(zeros(dim), [radius_hyperrectangle(Z1HyperRect)])

    Zred = minkowski_sum(Z1red, Z2)

    return Zred

end