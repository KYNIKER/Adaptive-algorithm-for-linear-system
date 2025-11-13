# ==================================
# Heat
#
# system type: continuous LTI system
# state dimension: 200
# input dimension: 1
# ==================================
using LazySets, MAT

function heat_model()
    file = matopen("heat/heat.mat")

    # system matrix
    A = Matrix(read(file, "A"))

    # input matrix
    B = sparse([67], [1], [1.0], size(A, 1), 1)

    # input domain
    U = BallInf([0.0], 0.5)

    return A, B, U
end