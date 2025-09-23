# Based on the paper JuliaReach: a Toolbox for Set-Based Reachability
using Plots, LazySets, LinearAlgebra
include("helperfunctions.jl")
const tΔ = 0.05
const r = 1.2

#A = [-1. 0.; 0. -1.]#reshape([-1.0], 1, 1)#UniformScaling(-1.0)#implement with double and add

#A = [cos(r) sin(r); -sin(r) cos(r)]
const A = [0. 1.; -2. 0.]
const μ = 0.001

P₁ = Zonotope([0., 1.5], [[0.0; 0.05]])#[0.0 0.0; 0.0 0.5])


#=

T = 4
N = floor(Int, T/tΔ)
#P₁ = Zonotope([0., 2], [0.0 0.0; 0.0 1])
ANorm = norm(A, Inf)
α = (exp(ANorm*tΔ)-1-tΔ*ANorm)/norm(P₁, Inf)
β = (exp(ANorm*tΔ)-1)*μ/ANorm

ϕ = exp(A*tΔ)
println(dim(P₁))

#=As = Matrix(undef, ceil(Int, log2(N)), 1)
As[1] = ϕ

for j in 2:size(As, 1)
    As[j] = As[j-1]*As[j-1]
end=#

ϕp = (I+ϕ)/2
ϕm = (I-ϕ)/2
gens = hcat(ϕp*P₁.generators,ϕm*P₁.center, ϕm*P₁.generators)


R₁ = minkowski_sum(Zonotope(ϕp*P₁.center, gens), Zonotope(zeros(2), (α+β)*I(2)))
R = R₁
boxes = []
push!(boxes, R₁)


ballβ = Zonotope(zeros(2), β*I(2))

for i in 2:N
    push!(boxes, minkowski_sum(linear_map(ϕ, boxes[i-1]), ballβ))
    global R = R∪box_approximation(boxes[i])
end
xs = range(0, T, length=N)
box = []=#

const T = 4

p = plot(dpi=300, thickness_scaling=1)



@time boxes = reachsets(A, tΔ, [0, T], P₁, μ)

proj = [0. 0.; 0. 1.]
corners = Vector(undef, size(boxes, 1))

@time begin
    for i in 1:(size(boxes, 1))
        corners[i] = vertices_list(box_approximation(boxes[i]))
        #plot!(p, Shape(getindex.(corners, 1), getindex.(corners, 2)), vars=(0, 1), c=:blue,lab="")
    end
end

shapes = Vector{Shape}(undef, size(boxes, 1))
@time rectangleFromHBox!(shapes, corners, tΔ, 2)

for i in eachindex(shapes)
    plot!(p, shapes[i], vars=(1,0), c=:blue, lab="")
end

#plot!(p, rectangleFromHBox(corners, tΔ, 2), vars=(1,0), c=:blue,lab="")

plot(p)

#savefig("boxplot_-1_0_0_-1.png")