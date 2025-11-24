using  LazySets, MAT

function fom_model()
    file = matopen("models/FOM/fom.mat")

    # system matrix
    A = float(read(file, "A"))  # the matrix has Int entries

    # input matrix
    B = read(file, "B")

    # input domain
    U = BallInf([0.0], 1.0)

    return A, B, U
end