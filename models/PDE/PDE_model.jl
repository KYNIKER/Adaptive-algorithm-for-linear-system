using LazySets, MAT

function pde_model()
    file = matopen(@current_path "pde.mat")

    # system matrix
    A = float(read(file, "A"))  # the matrix has Int entries

    # input matrix
    B = read(file, "B")

    # input domain
    U = BallInf([0.75], 0.25)

    return A, B, U
end