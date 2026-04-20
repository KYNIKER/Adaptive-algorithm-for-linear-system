using Plots, LazySets, LinearAlgebra, BenchmarkTools, Profile, Plots.PlotMeasures, LaTeXStrings


include("../helperfunctions.jl")
include("../models/heat/heat_load.jl")
include("../models/motor/motor_load.jl")
include("../models/building/building_load.jl")
include("../models/PDE/pde_load.jl")
include("../models/ISS/iss_load.jl")
include("../models/beam/beam_load.jl")
include("../models/MNA1/mna1_load.jl")
include("../ReACT.jl")
include("plotHelper.jl")
include("../ReACTDiscretize.jl")

function ReACTDiscretizeLazy(A, B, X0::Zonotope{N,Vector{N},Matrix{N}}, U::Zonotope, δ⁻, δ⁺, alg::ReachabilityAnalysis.Exponentiation.AbstractExpAlg=ReachabilityAnalysis.Exponentiation.BaseExp) where {N}
    phiDict = Dict{Float64,Matrix{Float64}}()
    discritezationDict = Dict()
    inputDiscritezationDict = Dict()

    U = B * U

    let ϕ::Matrix{Float64} = ReachabilityAnalysis.Exponentiation._exp(A, δ⁻, alg)
        tempM = similar(ϕ)
        d = δ⁻
        isInvA = isinvertible(A)

        A_abs = ReachabilityAnalysis.Exponentiation.elementwise_abs(A)
        Φcache = sum(A) == abs(sum(A)) ? ϕ : nothing
        P2A_abs = ReachabilityAnalysis.Exponentiation.Φ₂(A_abs, δ⁻, alg, isInvA, nothing)


        dU = d * U
        E_ψ = symmetric_interval_hull(P2A_abs * symmetric_interval_hull(A * U))
        E⁺ = symmetric_interval_hull(P2A_abs * symmetric_interval_hull((A * A) * X0))
        lt = ϕ * X0 ⊕ dU
        rt = E_ψ ⊕ E⁺
        f = lt ⊕ rt
        disc = CH(X0, f)
        Φ₁ = ReachabilityAnalysis.Exponentiation.Φ₁(A, d, alg, isInvA, nothing)
        P = Φ₁ * U
        while d < δ⁺
            discritezationDict[d] = copy(disc)
            inputDiscritezationDict[d] = copy(P)
            phiDict[d] = copy(ϕ)
            disc = copy(CH(disc, P ⊕ (ϕ * disc)))
            P = P ⊕ ϕ * P
            mul!(tempM, ϕ, ϕ)
            copy!(ϕ, tempM)
            d = d * 2
        end
        discritezationDict[d] = copy(disc)
        phiDict[d] = copy(ϕ)
        inputDiscritezationDict[d] = P

    end
    return discritezationDict, inputDiscritezationDict, phiDict
end

# We load a simple coswave
STRATEGY = 1

initialTimeStep = 0.8
m = initialTimeStep / 2
steps = 10

palette = Plots.palette(:fes10)
LazySets.Comparison.set_tolerance(Float64)
LazySets.Comparison.set_ztol(Float64, 1e-10)
alp = 1.0

A = [0. 1.;
    -2.5 0.]
P₁ = Zonotope([0.5, 1.5], [0.0 0.2 0.05; 0.45 -0.17 0.0]) #Zonotope([0., 1.5], [[0.0; 0.05]])

dimToPlot = 2
# No input
B = diagm([0.0, 0.0])
U = Zonotope(zeros(dim(P₁)), [zeros(dim(P₁))])

discritezationDict, _, phiDict = ReACTDiscretize(A, B, P₁, U, m, initialTimeStep, ReachabilityAnalysis.Exponentiation.BaseExp, 0, 5)
naivediscritezationDict, _, _ = ReACTDiscretize(A, B, P₁, U, initialTimeStep, initialTimeStep, ReachabilityAnalysis.Exponentiation.BaseExp, 0, 5)

p = plot(dpi=1200, thickness_scaling=1, guidefontsize=25, minorgrid=false, grid=false,
    legendfont=font(12, "Times"),
    tickfont=font(8, "Times"),
    xguidefont=font(12, "Times"),
    yguidefont=font(12, "Times"),
    xtick=([-1, 0, 1, 2], ["-1", "0", "1", "2"]),
    ytick=([-2, 0, 2, 4], ["-2", "0", "2", "4"]),
    bottom_margin=2mm,
    left_margin=2mm,
    right_margin=2mm,
    top_margin=2mm)
xlabel!(p, L"x")
ylabel!(p, L"y")
plot!(p, naivediscritezationDict[initialTimeStep], c=palette[9], alpha=alp, annotations=(2.3, -2.2, L"\Omega_{\left[0, \delta_1 \right]}"))
plot!(p, discritezationDict[initialTimeStep], c=palette[6], alpha=alp, annotations=(-0.1, 0.3, L"\Omega_{\left[0, \delta_1 \right]}"))
plot!(p, phiDict[m] * discritezationDict[m], c=palette[5], alpha=alp, annotations=(0.9, -1.0, L"\phi_{\delta_0} \Omega_{\left[0, \delta_0 \right]}"))
plot!(p, discritezationDict[m], c=palette[5], alpha=alp, annotations=(1.1, 1.0, L"\Omega_{\left[0, \delta_0 \right]}"))
plot!(p, P₁, c=:white, alpha=1., annotations=(0.5, 1.5, L"X_0"))

cX = copy(P₁)
tP = exp(A * (1 / steps))
for i in 0:steps-1
    global cX = tP * cX
    plot!(p, cX, c=:white, fa=0.0, la=0.3, ls=:dash, lc=:white)
end


savefig(p, "plots/" * "TD.pdf")

discritezationDict, _, phiDict = ReACTDiscretizeLazy(A, B, P₁, U, m, initialTimeStep, ReachabilityAnalysis.Exponentiation.BaseExp)
naivediscritezationDict, _, _ = ReACTDiscretizeLazy(A, B, P₁, U, initialTimeStep, initialTimeStep, ReachabilityAnalysis.Exponentiation.BaseExp)

p = plot(dpi=1200, thickness_scaling=1, guidefontsize=25, minorgrid=false, grid=false,
    legendfont=font(12, "Times"),
    tickfont=font(8, "Times"),
    xguidefont=font(12, "Times"),
    yguidefont=font(12, "Times"),
    xtick=([-1, 0, 1, 2], ["-1", "0", "1", "2"]),
    ytick=([-2, 0, 2, 4], ["-2", "0", "2", "4"]),
    bottom_margin=2mm,
    left_margin=2mm,
    right_margin=2mm,
    top_margin=2mm)
xlabel!(p, L"x")
ylabel!(p, L"y")

plot!(p, naivediscritezationDict[initialTimeStep], c=palette[9], alpha=alp, annotations=(2., -2.1, L"\Omega_{[0, \delta_1]}"))
plot!(p, discritezationDict[initialTimeStep], c=palette[6], alpha=alp, annotations=(0.1, 0.6, L"\Omega_{[0, \delta_1]}"))
plot!(p, phiDict[m] * discritezationDict[m], c=palette[5], alpha=alp, annotations=(0.9, -1.0, L"\phi_{\delta_0} \Omega_{[0, \delta_0]}"))
plot!(p, discritezationDict[m], c=palette[5], alpha=alp, annotations=(1.1, 1.0, L"\Omega_{[0, \delta_0]}"))
plot!(p, P₁, c=:white, alpha=1., annotations=(0.5, 1.5, L"X_0"))

cX = copy(P₁)
for i in 0:steps-1
    global cX = tP * cX
    plot!(p, cX, c=:white, fa=0.0, la=0.3, ls=:dash, lc=:white)
end

savefig(p, "plots/" * "TDL.pdf")

