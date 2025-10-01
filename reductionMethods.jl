using LazySets, LinearAlgebra

function PCA_reduce(Z, k) # Reduce to k generators

    c = center(Z) # get center
    G = genmat(Z) # Get generators

    U, S, V = svd(G) # Apply singular value decomposition (https://www.geeksforgeeks.org/machine-learning/singular-value-decomposition-svd/)

    # Note that U finds the eigenvectors of the covariance matrix
    # S ranks importance, and is ordered equivalently to eigenvalues
    k = min(k, length(S))   # Ensure k is not greater than the amount of dimensions

    G_red = U[:, 1:k] * Diagonal(S[1:k]) # Get top k

    return Zonotope(c, G_red)
end