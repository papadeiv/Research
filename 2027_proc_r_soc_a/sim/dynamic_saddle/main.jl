"""
    Main script

Run this file to execute the simulation, analyse and plot the results.
"""

# Import the necessary packages and local modules
include("inc.jl")

# Import the simulation's scripts
include("./scripts/sim.jl")
include("./scripts/proc.jl")
include("./scripts/figs.jl")

# Define the main algorithm
function main()
        # Solve the ensemble problem 
        x0 = [2.0, μ0]
        ensemble = evolve(f, η, Λ, x0, steps=Nt, stepsize=Δt, saveat=dt, particles=Ne)

        # Loop over the ensemble's sample paths
        for (idx_sol, solution) in enumerate(ensemble.state)
                # Truncate the solution and timesteps to the tipping
                N_tip = length(solution)
                x = solution[1:N_tip]
                t = ensemble.time[1:N_tip]

                # Plot the timeseries
                lines!(ax1, t, x, color = :black, linewidth = 1.0)
        end

        # Loop over the parameter values
        xs = Vector{Real}(undef, length(ensemble.parameter))
        xu = Vector{Real}(undef, length(ensemble.parameter))
        α = Vector{Real}(undef, length(ensemble.parameter))
        for (idx_μ, μ) in enumerate(ensemble.parameter)
                # Compute the stable and unstable equilibria
                equilibria = get_equilibria(f, μ, domain=[-10,10])

                # Extract the stable equilibrium and the return rate
                if length(equilibria.stable) > 0
                        xs[idx_μ] = equilibria.stable[1]
                        α[idx_μ] = (J(xs[idx_μ]))^2
                else
                        xs[idx_μ] = NaN 
                        α[idx_μ] = NaN 
                end

                # Extract the stable equilibrium
                if length(equilibria.unstable) > 0
                        xu[idx_μ] = equilibria.unstable[1]
                else
                        xu[idx_μ] = NaN 
                end
        end

        # Plot the window-of-oppurtunity
        poly!(ax1, Point2f[(τ,-2), (T,-2), (T,2), (τ,2)], color = (:black,0.10), strokecolor = :black, strokewidth = 2.0)
        poly!(ax2, Point2f[(τ,-1), (T,-1), (T,5), (τ,5)], color = (:black,0.10), strokecolor = :black, strokewidth = 2.0)

        # Plot the bifurcation diagram and the return rate of the frozen system
        lines!(ax1, ensemble.time, xs, color = :blue, linewidth = 3.0)
        lines!(ax1, ensemble.time, xu, color = :red, linewidth = 3.0)
        scatter!(ax1, ensemble.time[findfirst(>(0), ensemble.parameter)], 0.0, markersize = 20, color = :yellow, strokewidth = 2.0)
        lines!(ax2, ensemble.time, α, color = :black, linewidth = 3.0, linestyle = :dash)
        hlines!(ax2, 0.0, color = (:black,0.15), linewidth = 2.0)

        # Solve the deterministic problem
        ηs(t) = 0.0
        x0 = [2.0, μ0]
        solution = evolve(f, ηs, Λ, x0, steps=Nt, stepsize=Δt, saveat=dt)

        # Plot the local pullback attractor and its return rate
        pb_attractor = solution.state[1] 
        lines!(ax1, solution.time[1:length(pb_attractor)], pb_attractor, color = :orange, linewidth = 3.0)
        lines!(ax2, solution.time, [(pullback_rate(t))^2 for t in solution.time], color = :black, linewidth = 3.0)

        # Plot the relevant timesteps 
        scatter!(ax2, τ, 0, markersize = 20, color = :yellow, strokewidth = 2.0)
        scatter!(ax2, T, 0, markersize = 20, color = :red, strokewidth = 2.0)

        # Format and export the figure
        ax1.limits = (ensemble.time[1], ensemble.time[end], -2, +2)
        ax2.limits = (ensemble.time[1], ensemble.time[end], -1, +5)
        savefig("dynamic_saddle.pdf", fig)
end

# Execute the main
main()
