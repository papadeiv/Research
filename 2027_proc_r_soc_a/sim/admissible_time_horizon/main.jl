"""
    Admissible time horizon

Top: slowly ramped saddle-node with its frozen branches, the pullback attractor and
a noisy timeseries.  Bottom: frozen and pullback-attractor return rates.
Dashed lines bound the admissible region [tlo, Tc].  Produces `horizon.pdf`.
"""

using CairoMakie, Random
Random.seed!(1)

# ------------------------------------------------------------------ Parameters
μ0, ε = -0.5, 1e-2              # Initial parameter and ramp speed
σ     = 0.6*sqrt(ε)             # Noise, small-noise regime σ ≪ √ε
x0    = 2.0                     # Initial condition
dt    = 1e-3                    # Time step

τ   = -μ0/ε                     # Breaking time
τpb = τ + 2.338107*ε^(-1/3)     # Tipping time of the pullback attractor
tlo = 10/(2*sqrt(-μ0))          # End of the transient: five relaxation times at t = 0
Tc  = τ - 2*ε^(-1/3)            # Admissible time horizon
T   = τpb# + 5                   # End of the plot

# -------------------------------------------------------------- Integration
# Euler–Maruyama for dX = (-μ(t) - X²)dt + s dW, stopped once X has run away
function integrate(s)
        t, x = [0.0], [x0]
        while t[end] < T && x[end] > -2.5
                μ = μ0 + ε*t[end]
                push!(x, x[end] + (-μ - x[end]^2)*dt + s*sqrt(dt)*randn())
                push!(t, t[end] + dt)
        end
        return t, x
end

tpb, xpb = integrate(0.0)       # Pullback attractor (deterministic)
tX,  X   = integrate(σ)         # Noisy timeseries

# Frozen branches and the two return rates
tb  = range(0, τ, length = 500)
xp  = sqrt.(-(μ0 .+ ε.*tb))     # Stable branch x₊(μ(t))
α   = 4 .* xp.^2                # Frozen return rate, 4ε(τ - t)
αpb = 4 .* xpb.^2               # Return rate of the pullback attractor

# Draw y(t) in red inside [tlo, Tc] and in grey outside
function segments!(ax, t, y; kw...)
        for (m, c) in ((t .<= tlo, :grey), (tlo .<= t .<= Tc, :hotpink1), (t .>= Tc, :grey))
                lines!(ax, t[m], y[m]; color = c, kw...)
        end
end

# ------------------------------------------------------------------ Figure
fig   = Figure(size = (800, 650))
style = (xgridvisible = false, ygridvisible = false, xtickalign = 1, ytickalign = 1)

# Top: branches, pullback attractor, and the noisy timeseries on top
ax1 = Axis(fig[1, 1]; ylabel = L"x(t)", xticklabelsvisible = false, style...)
lines!(ax1, tb,  xp, color = :blue, linewidth = 3)
lines!(ax1, tb, -xp, color = :red, linewidth = 3)
lines!(ax1, tpb, xpb, color = :orange, linewidth = 2)
segments!(ax1, tX, X, linewidth = 1)
ylims!(ax1, -sqrt(-μ0), 1)

# Bottom: frozen (solid) and pullback-attractor (dotted) return rates
ax2 = Axis(fig[2, 1]; xlabel = L"t", ylabel = "return rate", style...)
segments!(ax2, tb,  α,   linewidth = 2.5)
segments!(ax2, tpb, αpb, linewidth = 2.5, linestyle = :dash)
ylims!(ax2, 0, 2)

# Admissible region, shared time axis
for ax in (ax1, ax2)
        vlines!(ax, [tlo, Tc], color = :black, linestyle = :dash)
        xlims!(ax, 0, T)
end
linkxaxes!(ax1, ax2)

save("horizon.pdf", fig)
