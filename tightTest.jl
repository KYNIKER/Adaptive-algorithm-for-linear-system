# Based on the paper JuliaReach: a Toolbox for Set-Based Reachability
using Plots, LazySets, LinearAlgebra
tΔ = 0.02

function rectangleFromHBox(corners, offset)
    Shape((getindex.(corners, 1)).+(offset*tΔ), getindex.(corners, 2))
end

T = 4
N = floor(Int, T/tΔ )
μ = 0.001
A = [-1. 0.; 0. -1.]#UniformScaling(-1.0)
P₁ = Zonotope([0., 1.5], [0.0 0.0; 0.0 0.5])
ANorm = norm(A, Inf)
α = (exp(ANorm*tΔ)-1-tΔ*ANorm)/norm(P₁, Inf)
β = (exp(ANorm*tΔ)-1)*μ/ANorm

ϕ = exp(A*tΔ)

ϕp = (I+ϕ)/2
ϕm = (I-ϕ)/2
gens = hcat(ϕp*P₁.generators, ϕm*P₁.center, ϕm*P₁.generators)


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
maxes = []
mines = []
box = []

plot(size=(500,400), dpi=300, thickness_scaling=1)

for i in 1:(N)
    approxBox = box_approximation(boxes[i])
    corners = vertices_list(boxes[i])
    plot!(rectangleFromHBox(corners, i - 1), c=:blue,lab="")
    push!(maxes, norm(approxBox, Inf))
    push!(mines, low(approxBox, 2))
end

plot!(xs, 1.0 * exp.(-xs), vars=(0, 1), c=:magenta, lab="")
plot!(xs, 2.0 * exp.(-xs), vars=(0, 1), c=:magenta, lab="")
