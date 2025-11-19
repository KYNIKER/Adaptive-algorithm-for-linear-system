using MAT, SparseArrays

function mna5_model()
    file = matopen("models/MNA5/mna5.mat")

    # system matrix
    A = sparse(read(file, "A"))

    # affine term
    b = sparsevec(19:27, [fill(-0.1, 5); fill(-0.2, 4)], 10913)


    return A, b
end