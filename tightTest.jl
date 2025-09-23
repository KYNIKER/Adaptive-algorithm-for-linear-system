# Based on the paper JuliaReach: a Toolbox for Set-Based Reachability
using Plots, LazySets, LinearAlgebra
include("helperfunctions.jl")

tΔ = 0.4



r = 1.5

T = 4
N = floor(Int, T/tΔ)
μ = 0.001
A = [cos(r) sin(r); -sin(r) cos(r)]
#A = [-1. 0.; 0. -1.]#reshape([-1.0], 1, 1)#UniformScaling(-1.0)#implement with double and add
P₁ = Zonotope([0., 1.5], [[0.0; 0.05]])#[0.0 0.0; 0.0 0.5])
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
gens = hcat(ϕp*P₁.generators,ϕm*P₁.center, ϕm*P₁.generators )


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
box = []

p = plot(dpi=300, thickness_scaling=1)

for i in 1:(N)
    approxBox = box_approximation(boxes[i])
    corners = vertices_list(boxes[i])
    plot!(p, rectangleFromHBox(corners), vars=(0,1), c=:blue,lab="")
    #plot!(p, Shape(getindex.(corners, 1), getindex.(corners, 2)), vars=(0, 1), c=:blue,lab="")
end

plot(p)

#savefig("boxplot_-1_0_0_-1.png")