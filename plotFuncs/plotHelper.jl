# Based on the paper JuliaReach: a Toolbox for Set-Based Reachability
using Plots, LazySets, LinearAlgebra, BenchmarkTools, FastExpm, Profile, PProf


include("../helperfunctions.jl")
include("../models.jl")
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
#include("CegarFunctions.jl")
include("../CegarInhomogenous.jl")


function getShapes(boxes2, timesteps)
    
    corners2 = Vector(undef, size(boxes2, 1))

    begin
        for i in 1:(size(boxes2, 1))
            H = box_approximation(boxes2[i])
            H_proj = LazySets.project(H, [dimToPlot]) # Only get the dimension we care about
            corners2[i] = vertices_list(H_proj)
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