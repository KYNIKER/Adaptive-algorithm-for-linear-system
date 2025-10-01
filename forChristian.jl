using LazySets, LinearAlgebra, Plots


function reachSetsForTimesteps(A, timesteps, interval, X₀, μ)
    startTime = minimum(interval)
    
    ANorm = norm(A, Inf)

    timestep = timesteps[1] # Use first element

    α = (exp(ANorm*timestep)-1-timestep*ANorm)/norm(X₀, Inf)
    β = (exp(ANorm*timestep)-1)*μ/ANorm

    ϕ = exp(A*timestep)

    ϕp = (I+ϕ)/2
    ϕm = (I-ϕ)/2
    gens = hcat(ϕp*X₀.generators,ϕm*X₀.center, ϕm*X₀.generators)


    R₁ = minkowski_sum(Zonotope(ϕp*X₀.center, gens), Zonotope(zeros(2), (α+β)*I(2)))

    R = [R₁]

    totalLength = length(timesteps)

    for i in 2:totalLength
        timestep = timesteps[i]
        α = (exp(ANorm*timestep)-1-timestep*ANorm)/norm(X₀, Inf)
        β = (exp(ANorm*timestep)-1)*μ/ANorm

        ϕ = exp(A*timestep)

        
        ballβ = Zonotope(zeros(2), β*I(2))
        push!(R, minkowski_sum(linear_map(ϕ, R[i-1]), ballβ))

    #res = [R₁]

    end
    return R
end

function rectangleFromHBoxWithTimestepArray(res::AbstractVector{Shape}, cornerss, timesteps, startTime, dim)

    currentTime = startTime
    for i in 1:size(cornerss, 1)#ændrer til getindex -> max min 
        timestep = timesteps[i]

        tope = getindex(cornerss, i)
        dimCoords = getindex.(tope, dim)
        maxcor = maximum(dimCoords)
        mincor = minimum(dimCoords)
        res[i] = Shape([currentTime, currentTime + timestep, currentTime + timestep, currentTime], [mincor, mincor, maxcor, maxcor])
        #res[i] = Shape([deltt*(i-1), deltt*i, deltt*i, deltt*(i-1)], [mincor, mincor, maxcor, maxcor])
        currentTime = currentTime + timestep
    end
    return res
end



# Part 2 originally a different file


# Based on the paper JuliaReach: a Toolbox for Set-Based Reachability
#using Plots, LazySets, LinearAlgebra
# include("helperfunctions.jl")
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

values = [0.05, 0.2, 0.05]
amounts = [25, 4, 25]

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