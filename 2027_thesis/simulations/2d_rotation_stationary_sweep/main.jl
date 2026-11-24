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

                # Compute the equilibria in the original coordinates system 
                equilibria = get_equilibria(f1, f2, μ[m])
                z0 = [equilibria.stable[1]; μ[m]]

                # Solve the ensemble problem 
                ensemble = evolve(f, η, Λ, z0, stepsize=dt, steps=Nt, particles=Ne)

                # Loop over the ensemble's sample paths
                ews_obs = zeros(Float64, length(ensemble.state), 3) 
                printstyled("μ=$(μ[m]): solving the weighted least-squares problems\n"; bold=true, underline=true, color=:light_blue)
                @showprogress for n in 1:length(ensemble.state)
                        # Extract the x- and y-components of the solution and compute its observables
                        x = (ensemble.state[n])[:, 1]
                        y = (ensemble.state[n])[:, 2]
                        observables = compute_observables(hcat(x, y), μ[m])

                        # Plot the trajectory in the potential landscape
                        if n == 1
                                plot_trajectory(x, y, equilibria, m)
                        end

                        # Loop over the observables
                        for observable in values(observables)[1:3]
                                # Solve the nonlinear least-squares problem to fit a cubic potential
                                solution = fit_potential(observable, noise=σ, transformation=[0.0,1.0,8.0], optimiser=β, attempts=Na)
                                push!(solutions, solution.fit)

                                # Perform postprocessing analysis on the solutions
                                ews = analyse(solution.fit, μ[m])
                                push!(results, ews)
                        end

                        # Update the ews matrix
                        ews_obs[n,:] = results

                        # Plot the reconstructed potential 
                        if n < 20
                                plot_solutions(observables, z0[1:2], m)
                        end

                        # Vacate the utility arrays
                        empty!(results)
                        empty!(solutions)
                end

                # Export the postprocessed data and vacate the analysis arrays for the next iteration
                writeout(ews_obs, "ews/$m.csv")

                # Export the figures
                savefig("weighted/μ=$(μ[m]).png", fig1)
                savefig("2d_potential/μ=$(μ[m]).png", fig)
                savefig("projections/μ=$(μ[m]).png", Fig)
        end

        # Plot and export the early-warning signal figure
        plot_ews()
end

# Execute the main
main()
