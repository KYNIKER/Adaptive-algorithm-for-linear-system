# ==================================
# Heat
#
# system type: continuous LTI system
# state dimension: 200
# input dimension: 1
# ==================================
using ReachabilityBenchmarks, MathematicalSystems, LazySets, MAT, SparseArraysfunction

function heat_model()
    file = matopen(@current_path "heat.mat")

    # system matrix
    A = read(file, "A")

    # input matrix
    B = sparse([67], [1], [1.0], size(A, 1), 1)

    # input domain
    U = BallInf([0.0], 0.5)

    return A, B, U
end