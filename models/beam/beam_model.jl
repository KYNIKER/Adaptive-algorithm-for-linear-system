using LazySets, MAT

function beam_model()
    file = matopen("models/beam/beam.mat")

    # system matrix
    A = read(file, "A")

    # input matrix
    B = read(file, "B")

    # input domain
    U = BallInf([0.5], 0.3)

    return A, B, U
end