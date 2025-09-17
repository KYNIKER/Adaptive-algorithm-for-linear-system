using ReachabilityAnalysis, Plots

# initial-value problem specification
p = @ivp(x' = -x, x(0) ∈ Interval(1, 2))

# flowpipe computation
sol = solve(p, T=5)

plot(sol, vars=(0, 1), xlab="t", ylab="x(t)")