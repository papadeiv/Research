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
        # Loop over the initial parameter values
        for h0 ∈ H0[1]
                # Plot the nullclines
                Sn_domain = LinRange(-1.0,1.0,1000)
                St_domain = LinRange(-1.0,1.0,1000)
                SN = [f1(Sn,St,h0) for Sn in Sn_domain, St in St_domain]
                ST = [f2(Sn,St,h0) for Sn in Sn_domain, St in St_domain]
                fig = Figure()
                ax = Axis(fig[1,1])
                contour!(ax, Sn_domain, St_domain, SN, levels=[0], color=:red, linewidth=3)
                contour!(ax, Sn_domain, St_domain, ST, levels=[0], color=:blue, linewidth=3)
                save("../../res/fig/wood/nullclines/$glb_idx.png", fig)

                # Compute the initial condition
                equilibria = get_equilibria(f1, f2, h0, guesses=[[Sn0,St0], [-0.5,0.5], [-0.1,0.6]])
                display(equilibria)
                z0 = [equilibria.stable[1]; h0]

                # Solve the ensemble slow-fast SDEs
                ensemble = evolve(f, η, Λ, z0, steps=Nt, stepsize=δt, particles=Ne)

                # Compute the drift of the quasi-steady equilibrium
                #qse = [(get_equilibria(f1, f2, h)).stable[1] for h in ensemble.parameter]

                n_bins = 20
                # Loop over the ensemble's sample paths
                printstyled("μ=$(h0): solving the weighted least-squares problems\n"; bold=true, underline=true, color=:light_blue)
                @showprogress for n in 1:length(ensemble.state)
                        # Extract the x-component of the solution (North-Atlantic salinity) 
                        x = (ensemble.state[n])[:, 1]# .- getindex.(qse, 1)
                        t = ensemble.time

                        # Check for tipping
                        tipping = find_tipping(x, check = 0.010, criterion = 0.050, verbose=true)
                        if tipping.check
                                # Trajectory has tipped => discard it
                        else
                                # Solve the nonlinear least-squares problem
                                solution = fit_potential(x, n_bins=n_bins)
                                push!(solutions, [ensemble.parameter[end]; solution.fit])

                                if n == 1
                                        fig = Figure()
                                        ax = Axis(fig[1,1])
                                        bins, pdf = fit_distribution(x, n_bins=n_bins)
                                        barplot!(ax, bins, pdf, strokecolor = :black, strokewidth = 1.0)

                                        V(x, μ) = μ[1]*x + μ[2]*(x^2) + μ[3]*(x^3)
                                        D = (std(x)^2)/2
                                        ρ(x, μ) = exp(-V(x, μ)/D)

                                        μ = solution.fit
                                        println("μ=$μ, var(x)= $(var(x)), D=$D")
                                        display(length(t))
                                        display(t[end])

                                        xs = +(1/(3*μ[3]))*(sqrt((μ[2])^2 - 3*μ[1]*μ[3]) - μ[2])
                                        xu = -(1/(3*μ[3]))*(sqrt((μ[2])^2 - 3*μ[1]*μ[3]) + μ[2])
                                        I = (-Inf,Inf)
                                        if xs > xu
                                                I = (xu, +Inf) 
                                        else
                                                I = (-Inf, xu) 
                                        end
                                        integral = IntegralProblem(ρ, (minimum(bins), maximum(bins)), μ)
                                        quadrature = solve(integral, QuadGKJL(; order=20000); maxiters=10000)
                                        p(x, μ) = ρ(x, μ)/(quadrature.u)

                                        domain = LinRange(minimum(bins), maximum(bins), 1000)
                                        lines!(ax, domain, [p(s, μ) for s in domain], color = :red, linewidth = 3.0)
                                        savefig("wood/distribution/$glb_idx.png", fig)
                                end

                                # Perform postprocessing analysis on the solutions and plot the results
                                analyse(solution.fit)
                        end
                end

                # Plot the reconstruction of the shifted potential
                include("./scripts/figs.jl")
                plot_solution(ensemble.parameter[end])

                # Export the postprocessed data and vacate the analysis arrays for the next iteration
                #export_data()
                empty!(results)
                empty!(solutions)

                # Export the figure
                savefig("wood/potential/$glb_idx.png", fig1)

                # Update the global index
                println()
                global glb_idx = glb_idx + 1
        end

        # Plot and export the early-warning signal figure
        #plot_ews()
end

# Execute the main
main()
