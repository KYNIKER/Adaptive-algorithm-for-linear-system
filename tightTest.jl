# Based on the paper JuliaReach: a Toolbox for Set-Based Reachability
using Plots, LazySets, LinearAlgebra

T = 4
tΔ = 0.1
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
for i in 1:(N)
    approxBox = box_approximation(boxes[i])
    push!(maxes, norm(approxBox, Inf))
    push!(mines, low(approxBox, 2))
end
plot(xs, maxes, lab="max")
plot!(xs, mines, lab="min")
plot!(xs, 1.0 * exp.(-xs), vars=(0, 1), c=:magenta, lab="")
plot!(xs, 2.0 * exp.(-xs), vars=(0, 1), c=:magenta, lab="")
