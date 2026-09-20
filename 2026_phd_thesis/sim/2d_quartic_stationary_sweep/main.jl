"""
    Main script

Run this file to execute the simulation, analyse and plot the results.
"""

# Import the necessary packages and local modules
include("inc.jl")

# Import the simulation's scripts
include("./scripts/sim.jl")
include("./scripts/proc.jl")
include("./scripts/plot.jl")

# Define the main algorithm
function main()
        # Parameter sweep 
        for m in 1:length(μ)
                # Create empty layouts for the figures
                include("./scripts/figs.jl")

                # Define the initial condition 
                equilibria = get_equilibria(f1, f2, μ[m])
                z0 = [equilibria.stable[2]; μ[m]]

                # Solve the ensemble problem 
                ensemble = evolve(f, η, Λ, z0, stepsize=dt, steps=Nt, particles=Ne)

                # Loop over the ensemble's sample paths
                printstyled("μ=$(μ[m]): solving the weighted least-squares problems\n"; bold=true, underline=true, color=:light_blue)
                @showprogress for n in 1:length(ensemble.state)
                        # Extract the x- and y-components of the solution
                        x = (ensemble.state[n])[:, 1]
                        y = (ensemble.state[n])[:, 2]

                        # Plot the trajectory in the potential landscape
                        if n == 1 
                                plot_trajectory(x, y, m)
                        end

                        # Check for tipping
                        tipping = find_tipping(x, check = 0.010, criterion = 0.050, verbose=false)
                        if tipping.check
                                # Trajectory has tipped => discrad it
                        else
                                # Extract the trajectory until the tipping and center it 
                                u = x[1:tipping.index] .- z0[1] 
                                t = ensemble.time[1:tipping.index]

                                # Solve the nonlinear least-squares problem to fit a cubic potential
                                solution = fit_potential(u, noise=σ, transformation=[0.0,1.0,8.0], optimiser=β, attempts=Na)
                                push!(solutions, solution.fit)

                                # Perform postprocessing analysis on the solutions
                                Vs = analyse(solution.fit, μ[m])

                                # Plot the timeseries and the reconstructed potential 
                                if n < 20
                                        plot_solutions(t, u, Vs, m)
                                end
                        end
                end

                # Export the postprocessed data and vacate the analysis arrays for the next iteration
                writeout(results, "ews/$m.csv")
                empty!(results)
                empty!(solutions)

                # Export the potential reconstruction figure
                savefig("weighted/μ=$(μ[m]).png", fig1)
        end

        # Plot and export the early-warning signal figure
        plot_ews()
end

# Execute the main
main()
