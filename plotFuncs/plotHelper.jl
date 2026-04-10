# Based on the paper JuliaReach: a Toolbox for Set-Based Reachability
using Plots, LazySets, LinearAlgebra, BenchmarkTools, Profile, PProf


include("../helperfunctions.jl")
include("../models/heat/heat_load.jl")
include("../models/motor/motor_load.jl")
include("../models/motor/motor_load.jl")
include("../models/building/building_load.jl")
include("../models/PDE/pde_load.jl")
include("../models/ISS/iss_load.jl")
include("../models/beam/beam_load.jl")
include("../models/FOM/fom_load.jl")
include("../models/MNA1/mna1_load.jl")
include("../models/MNA5/mna5_load.jl")
#include("../ReACT.jl")


function getShapes(boxes2, timesteps; dimToPlot=1, bounds=missing)

    corners2 = Vector(undef, size(boxes2, 1))

    begin
        for i in 1:(size(boxes2, 1))
            H = box_approximation(boxes2[i])
            H_proj = LazySets.project(H, [dimToPlot]) # Only get the dimension we care about
            if ismissing(bounds)

                corners2[i] = vertices_list(H_proj) #[[low(H_proj, 1), high(H_proj, 1)]] #
            else
                l = map(x -> x < bounds[1] ? bounds[1] : x, low(H_proj))
                h = map(x -> x > bounds[2] ? bounds[2] : x, high(H_proj))
                H_proj = Hyperrectangle(low=l, high=h)
                corners2[i] = vertices_list(H_proj)
            end
        end
    end

    # p = nothing
    # if plotConstraint
    #     cornersToSearch = [value[1] for box in corners2 for value in box]
    #     maxVal = maximum(cornersToSearch)
    #     minVal = minimum(cornersToSearch)
    #     # Take min and max with constraint
    #     constraintValAdjusted = constraint[1].b * 1.1
    #     maxVal = max(maxVal, constraintValAdjusted) 
    #     minVal = min(minVal, constraintValAdjusted)

    #     p = plot(dpi=300, thickness_scaling=1, ylims=(minVal, maxVal))
    # else
    #     p = plot(dpi=300, thickness_scaling=1)
    # end


    shapes2 = Vector{Shape}(undef, size(boxes2, 1))

    rectangleFromHBoxWithTimestepArray(shapes2, corners2, timesteps, minimum(T), 1)

    cornersToSearch = [value[1] for box in corners2 for value in box]
    maxVal = maximum(cornersToSearch)
    minVal = minimum(cornersToSearch)
    return shapes2, maxVal, minVal
end

function plotProjectedFlowpipe(boxes2, timesteps, dim1, dim2; approx=true)
    corners2 = Vector(undef, size(boxes2, 1))
    maxval = -∞
    minval = ∞
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
                #projectGDim1 = reduce(+, reduce(+, G, dims=dim1))
                projectGDim2 = reduce(+, G[dim2, :])
                #maxcor1 = c[dim1] + projectGDim1
                #mincor1 = c[dim1] - projectGDim1
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
            GProjS = unique(Go[[dim1, dim2], :]; dims=2)  #unique([permutedims(Go[dim1, :]) permutedims(Go[dim2, :])]) #[permutedims(GProj[dim1, :]); permutedims(GProj[dim2, :])]
            #println(eachcol(GProjS))
            corners2[i] = reduce_order(Zonotope([cProj[dim1], cProj[dim2]], GProjS), 8)

            maxval = max(maxval, maxcor2)
            minval = min(minval, mincor2)
            #=
            if k >= 10
                println(size(G))
                println(size(GProjS))
                println(c, " ", [cProj[dim1], cProj[dim2]])
                println(unique([permutedims(Go[dim1, :]); permutedims(Go[dim2, :])]))
                #println(GProj)
                #println(GProjS)

                return corners2, maxVal, minVal
            else
                k += 1
            end
            =#
        end
    end

    return corners2, maxval, minval
end