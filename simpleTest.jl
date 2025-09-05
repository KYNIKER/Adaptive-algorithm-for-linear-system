using Plots, LazySets

K = 2
T = 2
tΔ = 0.05
N = T/tΔ    
μ = 0.05
ϵ = BallInf([0.0, 0.0], (μ/K)*(ℯ ^ (K*tΔ) - 1))
Ω = [1.0 0.0; 0.0 1.0-tΔ]
b = [tΔ, 0.0]

map(θ) = Ω*θ + b 

P₁ = Zonotope([tΔ/2, 1.5], [tΔ/2 0.0; 0.0 0.5])
Q₁ = P₁ ⊕ ϵ
R₁ = Q₁
Q = Q₁
R = R₁
boxes = [ ]
i = 1

while i < N
    P = map(Q)
    global Q = P ⊕ ϵ
    global R = R ∪ Q
    global i += 1
end

#println(R)
plot(R)

