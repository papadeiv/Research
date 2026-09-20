"""
    Main script

Run this file to execute the simulation, analyse and plot the results.
"""

using LinearAlgebra, Statistics, StochasticDiffEq
using CairoMakie, ProgressMeter

# Import the simulation's scripts
include("./scripts/sim.jl")
include("./scripts/proc.jl")
#include("./scripts/figs.jl")

# Generate the ensemble of timeseries
sols = generate_samples()

# Estimate the return rates across rolling windows
t, y1, v1 = analyse(sols[1])

# Y[window, member, method]
Y = fill(NaN, length(t), Ne, 3)
V = similar(Y)
Y[:,1,:] = y1
V[:,1,:] = v1

for n in 2:Ne
        _, yn, vn = analyse(sols[n])
        Y[:,n,:] = yn
        V[:,n,:] = vn
end

# Extrapolated tipping time, one per member per method
Ttip = fill(NaN, Ne, 3)
for n in 1:Ne, m in 1:3
        Ttip[n,m] = estimate_tipping_time(t, Y[:,n,m], V[:,n,m])
end

#@printf("fold τ = %.1f, tipping T = %.1f, %d windows of %d samples\n", τ, T, length(t), ΔT)
for m in 1:3
        c = filter(isfinite, Ttip[:,m])
        #@printf("  %s: median %7.1f (%+5.1f%%), IQR [%7.1f, %7.1f]\n", algs[m], median(c), 100*(median(c)/τ - 1), quantile(c,0.25), quantile(c,0.75))
end

# --- Figure 1: timeseries (top) and return rates (bottom) ---------------
fig1 = Figure(size = (900, 700))
npaths, npts = 15, 1500                         # ≈ 22k line vertices ⇒ well under 3 MB

# Evenly thin a vector down to at most n points (keeps the PDF small).
sub(v, n) = v[round.(Int, range(1, length(v), length = min(n, length(v))))]

ax1 = Axis(fig1[1,1], xlabel = L"t", ylabel = L"x(t)")
for n in 1:(Ne ÷ npaths):Ne
        lines!(ax1, sub(sols[n].t, npts), sub(sols[n].x, npts), color = (:gray, 0.4))
end
lines!(ax1, sub(sols[1].t, npts), sub(sols[1].x, npts), color = :black)
vlines!(ax1, T, color = :black, linestyle = :dash)
ylims!(ax1, 1.15xesc, -1.15xesc)

# Midpoint of the window
med(v) = (w = filter(isfinite, v); isempty(w) ? NaN : median(w))

ax2 = Axis(fig1[2,1], xlabel = L"t", ylabel = L"\lambda^2(t)")
for m in 1:3
        ym = [med(Y[j,:,m]) for j in axes(Y,1)]
        vm = [med(V[j,:,m])/Ne for j in axes(V,1)]  # variance of the ensemble median
        lines!(ax2, t, ym, color = cols[m], linewidth = 2.5, label = algs[m])
        Tm, c = estimate_tipping_time(t, ym, vm)
        te = LinRange(t[end], Tm, 50)
        lines!(ax2, te, c[1] .+ c[2].*te, color = cols[m], linestyle = :dash, linewidth = 2.5)
        scatter!(ax2, Tm, 0, color = cols[m], markersize = 14)
end
lines!(ax2, t, λ2_exact.(t), color = :black, linestyle = :dot, linewidth = 2)
hlines!(ax2, 0, color = (:black, 0.5)); vlines!(ax2, T, color = :black, linestyle = :dash)
xlims!(ax1, 0, Nt); xlims!(ax2, 0, Nt)
ylims!(ax2, -0.05λ2_exact(0), 1.05λ2_exact(0))
axislegend(ax2)
save("figure1.pdf", fig1)

# --- Figure 2: distribution of extrapolated tipping times ---------------
fig2 = Figure(size = (1250, 400))
#Label(fig2[0,1:3], @sprintf("Extrapolated tipping time, %s sampler (truth T = %.1f)", sampler, T), fontsize = 20)
for m in 1:3
        ax = Axis(fig2[1,m], title = algs[m], xlabel = L"T_{tip}", ylabel = "count")
        #Nb1 = convert(Int64, ceil(range/(3.49*std(u)*(Nt)^(-1.0/3.0))))   # Scott's rule (1985)
        hist!(ax, filter(isfinite, Ttip[:,m]), bins = 25, color = (cols[m], 0.65))
        vlines!(ax, T, color = :black, linestyle = :dash, linewidth = 2)
end
save("figure2.pdf", fig2)
