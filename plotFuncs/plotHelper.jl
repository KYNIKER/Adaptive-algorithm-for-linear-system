using Plots, LazySets, LinearAlgebra, BenchmarkTools, Profile

include("../models/heat/heat_load.jl")
include("../models/motor/motor_load.jl")
include("../models/building/building_load.jl")
include("../models/PDE/pde_load.jl")
include("../models/ISS/iss_load.jl")
include("../models/beam/beam_load.jl")
include("../models/MNA1/mna1_load.jl")


function plotSupportFlowpipe(boxes2, timesteps, pos, neg)
    corners2 = Vector(undef, size(boxes2, 1))
    maxval = -Inf
    minval = Inf
    time = 0.0

    for i in 1:(size(boxes2, 1))
        t = timesteps[i]
        r = boxes2[i]

        maxcor2 = r[pos]
        mincor2 = -r[neg]
        corners2[i] = Shape([time, time + t, time + t, time], [mincor2, mincor2, maxcor2, maxcor2])
        time += t
        maxval = max(maxval, maxcor2)
        minval = min(minval, mincor2)
    end

    return corners2, maxval, minval
end