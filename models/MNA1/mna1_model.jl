using MAT, SparseArrays

function mna1_model()
    file = matopen("models/MNA1/mna1.mat")

    # system matrix
    A = sparse(read(file, "A"))

    # affine term
    b = sparsevec(570:578, [fill(-0.1, 5); fill(-0.2, 4)], 578)

    return A, b
end