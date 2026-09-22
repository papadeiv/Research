"""
    Plotting script

Functions used to create the plots in each figure.
"""

# Plot the timeseries stationary distribution and nonlinear solutions 
function plot_solutions(timestep, timeseries, trend, residual, solution, idx)
        #--------------------------------------#
        #        Timeseries and residual       #
        #--------------------------------------#
        fig = Figure()

        ax1 = Axis(fig[1,1:2], title="Original signal w/ trend")
        lines!(ax1, timestep, timeseries, color = :black, linewidth = 1.0)
        lines!(ax1, timestep, trend, color = :red, linewidth = 3.0)

        ax2 = Axis(fig[2,1:2], title="Detrended signal")
        lines!(ax2, timestep, residual, color = idx, colormap = :managua, colorrange = (1,49), linewidth = 1.0)


        #--------------------------------------#
        #      Distribution and potential      #
        #--------------------------------------#
        ax3 = Axis(fig[1,3], title="Empirical distribution")
        Nb = convert(Int64, ceil(abs(maximum(residual)-minimum(residual))/(3.49*std(residual)*(length(residual))^(-1.0/3.0))))
        bins, pdf = fit_distribution(residual, n_bins=Nb+1)
        barplot!(ax3, bins, pdf, color = pdf, colormap = [:white,CtpGray], strokecolor = :black, strokewidth = 1)

        D = (std(residual)^2)/2
        p(x, c) = exp(-V(x,c)/D)
        N = normalise(p, solution) 
        ρ(x) = N*p(x, solution)
        domain = LinRange(minimum(bins), maximum(bins), 1000)
        lines!(ax3, domain, [ρ(x) for x in domain], color = idx, colormap = :managua, colorrange = (1,49), linewidth = 3.0)

        ax4 = Axis(fig[2,3], title="Reconstructed potential")
        domain = LinRange(-4, 4, 1000)
        lines!(ax4, domain, [V(x, solution) for x in domain], color = idx, colormap = :managua, colorrange = (1,49), linewidth = 3.0)

        savefig("amoc/residuals/$idx.png", fig)
end
