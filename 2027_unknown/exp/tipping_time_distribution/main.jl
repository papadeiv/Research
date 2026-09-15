"""
    Exit time distributions

Numerical check of the tipping time distributions of the saddle-node normal form:
the exponential law of the frozen problem, and the small- and large-noise laws
for the slowly ramped problem.  Produces `exit_times.pdf`.

Run with:  julia -t auto exit.jl
"""

using Statistics, Printf, Random
using StochasticDiffEq, SpecialFunctions, CairoMakie

Random.seed!(1)

# ---------------------------------------------------------------- Experiments
stationary = (μ0 = -1.0, ε = 0.0,  D = 0.25,   dt = 2e-2, T = 4.0e3)
smallnoise = (μ0 = -1.0, ε = 1e-2, D = 1e-3,   dt = 5e-3, T = 1.25e2)
largenoise = (μ0 = -1.0, ε = 1e-3, D = 0.0625, dt = 1e-2, T = 1.1e3)

Ne     = 10_000       # Ensemble size
npaths = 100          # Trajectories drawn per panel
xtip = -20.0          # Escape threshold, far past the unstable branch

# ------------------------------------------------------------------ Simulation
function simulate(P)
        σ = sqrt(2*P.D)
        f(x, p, t) = -(P.μ0 + P.ε*t) - x^2      # Drift
        η(x, p, t) = σ                          # Diffusion

        # Terminate a trajectory once it has escaped
        tipped = DiscreteCallback((x, t, integrator) -> x < xtip, terminate!)

        # Solve the ensemble problem starting from the stable equilibrium
        problem  = SDEProblem(f, η, sqrt(-P.μ0), (0.0, P.T))
        ensemble = solve(EnsembleProblem(problem), EM(), EnsembleThreads();
                         trajectories = Ne, dt = P.dt, adaptive = false,
                         saveat = P.T/500, callback = tipped)

        # Tipping times of every escaped trajectory
        ttip = [s.t[end] for s in ensemble.u if s.t[end] < P.T - P.dt]

        # Keep only npaths trajectories for plotting, whatever Ne is
        sols = [(t = collect(s.t), x = collect(s.u))
                for s in ensemble.u[1:max(1, Ne ÷ npaths):end]]

        return sols, ttip
end

# ------------------------------------------------------------------- Theory
ΔV(μ)   = 4/3*(-μ)^(3/2)                        # Potential barrier
k(μ, D) = sqrt(-μ)/π*exp(-ΔV(μ)/D)              # Kramers escape rate

# Exponential law of the frozen problem
p_stat(t, P) = k(P.μ0, P.D)*exp(-k(P.μ0, P.D)*t)

# Tail integral of the fourth power of Ai, by the trapezoidal rule
function G(w; wmax = 10.0, n = 3000)
        u = range(w, wmax, length = n)
        y = airyai.(u).^4
        return (sum(y) - (y[1] + y[end])/2)*step(u)
end

# Small-noise law: valid between the first two zeros of Bi
function p_small(t, P)
        τ = -P.μ0/P.ε
        w = P.ε^(1/3)*(τ - t)                   # = -μ(t)ε^{-2/3}
        -3.2711 < w < -1.1737 || return 0.0
        A, B, g = airyai(w), airybi(w), G(w)
        return P.ε^(5/6)/(sqrt(4π^5*P.D)*B^2*sqrt(g))*exp(-P.ε*A^2/(4π^2*P.D*B^2*g))
end

# Gaussian approximation of the small-noise law about its maximum
function p_gauss(t, P)
        τpb = -P.μ0/P.ε + 2.338107*P.ε^(-1/3)
        v   = P.D*P.ε^(-5/3)
        return exp(-(t - τpb)^2/(2v))/sqrt(2π*v)
end

# Large-noise law: the frozen rate modulated by the survival probability
function p_large(t, P)
        μ = P.μ0 + P.ε*t
        μ < 0 || return 0.0
        return k(μ, P.D)*exp(-P.D/(2π*P.ε)*exp(-ΔV(μ)/P.D))
end

# Deterministic pullback attractor, integrated from x0 at t = 0 by RK4
function pullback(P; x0 = 2.0)
        g(t, x) = -(P.μ0 + P.ε*t) - x^2
        h  = P.dt
        ts, xs = [0.0], [x0]
        for i in 1:round(Int, P.T/h)
                t, x = ts[end], xs[end]
                k1 = g(t, x)
                k2 = g(t + h/2, x + h*k1/2)
                k3 = g(t + h/2, x + h*k2/2)
                k4 = g(t + h,   x + h*k3)
                xn = x + h*(k1 + 2k2 + 2k3 + k4)/6
                xn < -2.5 && break                      # stop once it has run away
                push!(ts, t + h)
                push!(xs, xn)
        end
        return ts, xs
end

# --------------------------------------------------------------------- Run
runs    = (stationary, smallnoise, largenoise)
results = [simulate(P) for P in runs]

# ------------------------------------------------------------------- Figure
# Evenly thin a vector down to at most n points
sub(v, n) = v[round.(Int, range(1, length(v), length = min(n, length(v))))]

# Number of histogram bins by Scott's rule (1979): h = 3.49 σ n^(-1/3)
scott(x) = max(1, ceil(Int, (maximum(x) - minimum(x))/(3.49*std(x)*length(x)^(-1/3))))

fig = Figure(size = (1300, 700))

# Ticks pointing inwards, no grid lines
style = (xgridvisible = false, ygridvisible = false, xtickalign = 1, ytickalign = 1)

# Top row: sample trajectories with the equilibria (or bifurcation diagram)
for (col, P) in enumerate(runs)
        sols, _ = results[col]
        ax = Axis(fig[1, col]; xlabel = L"t",
                  ylabel = col == 1 ? L"x(t)" : "",
                  yticklabelsvisible = col == 1, style...)

        for s in sols
                lines!(ax, sub(s.t, 400), sub(s.x, 400), color = (:black, 0.5))
        end

        # Stable (blue) and unstable (red) equilibria, up to the breaking time
        τ  = P.ε > 0 ? -P.μ0/P.ε : P.T
        tb = range(0, min(τ, P.T), length = 400)
        xb = sqrt.(max.(-(P.μ0 .+ P.ε.*tb), 0.0))
        lines!(ax, tb,  xb, color = :blue, linewidth = 4.0)
        lines!(ax, tb, -xb, color = :red,  linewidth = 4.0)

        # Window of opportunity, between the breaking and the tipping time
        P.ε > 0 && vspan!(ax, τ, τ + 2.338107*P.ε^(-1/3), color = (:black, 0.15))

        # Pullback attractor, and the breaking (dashed) and tipping (solid) times
        if P.ε > 0
                vlines!(ax, τ, color = :black, linestyle = :dash, linewidth = 2)
                vlines!(ax, τ + 2.338107*P.ε^(-1/3), color = :black, linewidth = 2)
                lines!(ax, pullback(P)..., color = :gold3, linewidth = 2.0)
        end

        xlims!(ax, 0, P.T)
        ylims!(ax, -2.0, 2.0)
end

# Bottom row: histogram of the tipping times with the predicted density
for (col, P) in enumerate(runs)
        _, ttip = results[col]
        ax = Axis(fig[2, col]; xlabel = L"t",
                  ylabel = col == 1 ? "density" : "",
                  ytickformat = v -> [@sprintf("%.1e", y) for y in v],
                  style...)

        # Window of opportunity
        if P.ε > 0
                τ, τpb = -P.μ0/P.ε, -P.μ0/P.ε + 2.338107*P.ε^(-1/3)
                vspan!(ax, τ, τpb, color = (:black, 0.15))
        end

        # Empirical distribution
        hist!(ax, ttip, bins = scott(ttip), normalization = :pdf, color = (:hotpink1, 1.00))

        # Plotting range: the escapes, widened to include τ and τ_pb when ramped
        lo, hi = extrema(ttip)
        pad    = 0.05*(hi - lo)
        lo, hi = lo - pad, hi + pad
        if P.ε > 0
                τ, τpb = -P.μ0/P.ε, -P.μ0/P.ε + 2.338107*P.ε^(-1/3)
                lo, hi = min(lo, τ - pad), max(hi, τpb + pad)
        elseif col == 1
                lo = 0.0
        end
        tt = range(lo, hi, length = 1000)

        # Predicted densities
        color = :teal
        if col == 1
                lines!(ax, tt, p_stat.(tt, Ref(P)), color = color, linewidth = 4)
        elseif col == 2
                lines!(ax, tt, p_small.(tt, Ref(P)), color = color, linewidth = 4)
        else
                lines!(ax, tt, p_large.(tt, Ref(P)), color = color, linewidth = 4)
        end

        # Breaking time (dashed) and tipping time of the pullback attractor (solid)
        if P.ε > 0
                vlines!(ax, -P.μ0/P.ε, color = :black, linestyle = :dash, linewidth = 2)
                vlines!(ax, -P.μ0/P.ε + 2.338107*P.ε^(-1/3), color = :black, linewidth = 2)
        end

        xlims!(ax, lo, hi)
        ylims!(ax, 0, nothing)
end

save("exit_times.pdf", fig)
