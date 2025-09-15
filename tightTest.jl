using Plots, LazySets

K = 2
T = 2
tΔ = 0.01
N = T/tΔ    
μ = 0.001
U = BallInf([0.0, 0.0], (μ/K)*(ℯ ^ (K*tΔ) - 1))
#Ω = [1.0 0.0; 0.0 1.0-tΔ]
b = [tΔ, 0.0]
M(Θ) = [cos(Θ) sin(Θ); -sin(Θ) cos(Θ)]
circ = M(π/N)
map(θ) = Ω*θ + b 
A1 = [-1.0 -4.0; 4.0 -1.0]
A2 = [1.0 4.0; -4.0 -1.0]
ϕ = ℯ^(A1*tΔ)

P₁ = Zonotope([1.0, 0.0], [0.1 0.0; 0.0 0.1])
Ω₀= convex_hull(P₁,ϕ*P₁)
X = Ω₀
V = U
S = U
Q₁ = P₁ #⊕ ϵ
R₁ = Q₁
Q = Q₁
R = X
boxes = [ ]
i = 1

while i <= N
    if i < N/2
        global ϕ = ℯ^(A1*tΔ)
    else 
        global ϕ = ℯ^(A2*tΔ)
    end
    global X = ϕ*X
    if i > 1
        global S = S ⊕ V
    end
    global V = ϕ*V
    global Ω = X⊕S
    global R = R∪Ω
    global i += 1
end

#println(R)
plot(R)
savefig("myplotA1A2.png")