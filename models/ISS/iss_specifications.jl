using LazySets, MAT

function iss_specification()
    # initial set: xᵢ ∈ [-1e-4, 1e-4] for all i
    X0 = BallInf(zeros(270), 1e-4)

    # safety property: y ≤ 7e-4 for linear combination y (defined in out.mat)
    y = read(matopen("models/ISS/out.mat"), "M")[1, :]
    property = LazySets.HalfSpace(y, 7e-4)

    # time horizon: 20 time units
    time_horizon = 20.0

    return X0, time_horizon, property
end