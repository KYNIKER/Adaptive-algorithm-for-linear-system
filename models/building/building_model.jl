using LazySets, MAT

function building_model()
    file = matopen("models/building/building.mat")

    A = read(file, "A")

    B = read(file, "B")

    U = BallInf([0.9], 0.1)

    return A, B, U
end