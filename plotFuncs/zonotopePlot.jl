# Originally made by Christian Schilling
using LazySets, Plots, Plots.PlotMeasures

c = [2, 1.5]
g1 = [-0.5, 0.3]
g2 = [0.4, 0.8]
g3 = [0.3, 0]
g4 = [0, 0.3]
G = hcat(g1, g2, g3, g4)
Z = Zonotope(c, G)
println("\nzonotope: c = $c;   G = $G")

fig = plot(
    xlims = (0, 4),
    ylims = (-0.3, 3),
    ratio = 1,
    xlab = "X",
    ylab = "Y",
    legendfontsize = 25,
    tickfont = font(25, "Times"),
    guidefontsize = 25,
    xguidefont = font(25, "Times"),
    yguidefont = font(25, "Times"),
    xtick = ([0, 1, 2, 3, 4], ["0", "1", "2", "3", "4"]),
    ytick = ([0, 1, 2, 3], ["0", "1", "2", "3"]),
    bottom_margin = 0mm,
    left_margin = 0mm,
    right_margin = 0mm,
    top_margin = 0mm,
    size = (500, 450)
)

plot!(fig, Z)

plot!(
    fig, [c[1], c[1] + g1[1]], [c[2], c[2] + g1[2]],
    linecolor = :red, arrow = :arrow, linestyle = :dot, width = 2,
    annotations = (1.4, 2.0, text("g₁", 20)),
    lab = ""
)

plot!(
    fig, [c[1], c[1] + g2[1]], [c[2], c[2] + g2[2]],
    linecolor = :red, arrow = :arrow, linestyle = :dot, width = 2,
    annotations = (2.5, 2.5, text("g₂", 20)),
    lab = ""
)

plot!(
    fig, [c[1], c[1] + g3[1]], [c[2], c[2] + g3[2]],
    linecolor = :red, arrow = :arrow, linestyle = :dot, width = 2,
    annotations = (2.5, 1.5, text("g₃", 20)),
    lab = ""
)

plot!(
    fig, [c[1], c[1] + g4[1]], [c[2], c[2] + g4[2]],
    linecolor = :red, arrow = :arrow, linestyle = :dot, width = 2,
    annotations = (2.0, 2.0, text("g₄", 20)),
    lab = ""
)

plot!(
    fig, Singleton(c),
    annotations = (2, 1.3, text("c", 20)),
    color = :red, alpha = 1.
)

savefig(fig, "zonotope.pdf")