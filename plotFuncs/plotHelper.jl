using Plots, LazySets, LinearAlgebra, BenchmarkTools, Profile

include("../helperfunctions.jl")
include("../models/heat/heat_load.jl")
include("../models/motor/motor_load.jl")
include("../models/building/building_load.jl")
include("../models/PDE/pde_load.jl")
include("../models/ISS/iss_load.jl")
include("../models/beam/beam_load.jl")
include("../models/MNA1/mna1_load.jl")

function plotProjectedFlowpipe(boxes2, timesteps, dim1, dim2; approx=true)
    corners2 = Vector(undef, size(boxes2, 1))
    maxval = -Inf
    minval = Inf
    k = 0
    if approx
        if dim1 != 0
            dimSize = size(genmat(boxes2[1]), 1)
            projectionMatrix = zeros(Float64, dimSize, dimSize)
            projectionMatrix[dim1, dim1] = 1.0
            projectionMatrix[dim2, dim2] = 1.0


            for i in 1:(size(boxes2, 1))
                r = boxes2[i]
                G = genmat(r)
                c = r.center

                projectedG = projectionMatrix * G
                projectGDim1s = mapreduce(x -> sign(x[dim1]) * x, +, eachcol(projectedG))
                projectGDim2s = mapreduce(x -> sign(x[dim2]) * x, +, eachcol(projectedG))

                maxcor1s = c + projectGDim1s
                mincor1s = c - projectGDim1s
                maxcor2s = c + projectGDim2s
                mincor2s = c - projectGDim2s

                projectGDim1 = reduce(+, G[dim1, :])
                projectGDim2 = reduce(+, G[dim2, :])
                maxcor1 = c[dim1] + projectGDim1
                mincor1 = c[dim1] - projectGDim1
                maxcor2 = c[dim2] + projectGDim2
                mincor2 = c[dim2] - projectGDim2

                corners2[i] = Shape([(mincor1s[dim1], mincor1s[dim2]), (mincor2s[dim1], mincor2s[dim2]), (maxcor1s[dim1], maxcor1s[dim2]), (maxcor2s[dim1], maxcor2s[dim2])])
                maxval = max(maxval, maxcor1, maxcor2)
                minval = min(minval, mincor1, mincor2)
            end



        else

            time = 0.0
            for i in 1:(size(boxes2, 1))
                t = timesteps[i]
                r = boxes2[i]
                G = genmat(r)
                c = r.center
                G = abs.(genmat(r))
                c = r.center
                projectGDim2 = reduce(+, G[dim2, :])
                maxcor2 = c[dim2] + projectGDim2
                mincor2 = c[dim2] - projectGDim2
                corners2[i] = Shape([time, time + t, time + t, time], [mincor2, mincor2, maxcor2, maxcor2])
                time += t
                maxval = max(maxval, maxcor2)
                minval = min(minval, mincor2)
            end



        end
    else
        dimSize = size(genmat(boxes2[1]), 1)
        projectionMatrix = zeros(Float64, dimSize, dimSize)
        projectionMatrix[dim1, dim1] = 1.0
        projectionMatrix[dim2, dim2] = 1.0
        for i in 1:(size(boxes2, 1))
            r = boxes2[i]
            Go = genmat(r)
            c = r.center
            G = abs.(genmat(r))
            c = r.center
            projectGDim2 = reduce(+, G[dim2, :])
            maxcor2 = c[dim2] + projectGDim2

            mincor2 = c[dim2] - projectGDim2
            cProj = projectionMatrix * c
            GProj = projectionMatrix * Go
            GProjS = unique(Go[[dim1, dim2], :]; dims=2)
            corners2[i] = reduce_order(Zonotope([cProj[dim1], cProj[dim2]], GProjS), 8)

            maxval = max(maxval, maxcor2)
            minval = min(minval, mincor2)
        end
    end

    return corners2, maxval, minval
end