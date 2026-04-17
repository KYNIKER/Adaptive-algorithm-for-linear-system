using LazySets, LinearAlgebra, Plots

function rectangleFromHBoxWithTimestepArray(res::AbstractVector{Shape}, cornerss, timesteps::Vector{Float64}, startTime, dim::Int64)

    currentTime = startTime
    for i in 1:size(cornerss, 1)
        timestep::Float64 = timesteps[i]::Float64

        tope::Vector{Vector{Float64}} = getindex(cornerss, i)::Vector{Vector{Float64}}
        dimCoords::Vector{Float64} = getindex.(tope, dim)
        maxcor::Float64 = maximum(dimCoords)::Float64
        mincor::Float64 = minimum(dimCoords)::Float64
        res[i] = Shape([currentTime, currentTime + timestep, currentTime + timestep, currentTime], [mincor, mincor, maxcor, maxcor])
        currentTime::Float64 = currentTime + timestep::Float64
    end
    return res
end

Base.:+(z1::Zonotope, z2::Zonotope) = Zonotope(z1.center + z2.center, z1.generators + z2.generators)