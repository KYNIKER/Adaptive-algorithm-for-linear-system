# Based on the paper JuliaReach: a Toolbox for Set-Based Reachability
using Plots, LazySets, LinearAlgebra
include("helperfunctions.jl")
const tΔ = 0.01
const r = 1.2

#A = [-1. 0.; 0. -1.]#reshape([-1.0], 1, 1) #UniformScaling(-1.0) #implement with double and add

#A = [cos(r) sin(r); -sin(r) cos(r)]
const A = [0. 1.; -2.5 0.]
const μ = 0.001

P₁ = Zonotope([0., 1.5], [[0.0; 0.05]]) #[0.0 0.0; 0.0 0.5])

#const T = 8
T = 8

p = plot(dpi=300, thickness_scaling=1)


timesteps = []

values = [0.05, 0.20, 0.05]
amounts = [20, 5, 20]

for (value, amount) in zip(values, amounts)
    global timesteps = vcat(timesteps, fill(value, amount))
end


###
@time boxes2 = reachSetsForTimesteps(A, timesteps, [0, T], P₁, μ)


corners2 = Vector(undef, size(boxes2, 1))

@time begin
    for i in 1:(size(boxes2, 1))
        corners2[i] = vertices_list(box_approximation(boxes2[i]))
    end
end

shapes2 = Vector{Shape}(undef, size(boxes2, 1))
#@time rectangleFromHBox!(shapes2, corners2, 2*tΔ, 2)
@time rectangleFromHBoxWithTimestepArray(shapes2, corners2, timesteps, 0, 2)

for i in eachindex(shapes2)
    plot!(p, shapes2[i], vars=(1,0), c=:red, alpha=:0.5, lab="")
end

##
plot(p)