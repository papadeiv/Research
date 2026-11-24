"""
    Plotting script

Functions used to create the plots in each figure.
"""

# Create empty layouts for the figures
include("./figs.jl")

# Plot the trajectory of the particle inside the double-well potential
function plot_trajectory(x_component, y_component, equilibria, parameter_idx)
        # Define the domain of the potential
        X = collect(LinRange(-2.5, 2.5, 1000))
        Y = collect(LinRange(-2.5, 2.5, 1000))

        # Plot the contour of the potential
        U(x, y) = x^4 - 2*x^2 + μ[parameter_idx]*x + y^2 + 1 
        W(x, y) = U(a11*x + a12*y, a21*x + a22*y) 
        contourf!(ax, X, Y, [W(x, y) for x in X, y in Y], levels = 375)

        # Plot the projecting subspaces
        α = compute_subspace(μ[parameter_idx])
        lines!(ax, X, [α*x for x in X], color = CtpPink, linewidth = 3.0)
        α = compute_optimal(hcat(x_component, y_component))
        lines!(ax, X, [α*x for x in X], color = :red, linestyle = :dash, linewidth = 3.0)
        α = -8.0 
        lines!(ax, X, [α*x for x in X], color = CtpTeal, linewidth = 3.0)

        # Plot the trajectory in phase space
        lines!(ax, x_component, y_component, color = (:white, 0.75), linewidth = 0.125)

        # Plot the equilibria
        scatter!(ax, [eq[1] for eq in equilibria.stable], [eq[2] for eq in equilibria.stable], markersize = 15.0, color = :green1, strokewidth = 1.0, strokecolor = :black)
        scatter!(ax, [eq[1] for eq in equilibria.unstable], [eq[2] for eq in equilibria.unstable], markersize = 25.0, marker = :star5, color = :yellow, strokewidth = 1.0, strokecolor = :black)

        # Compute and plot the observables' timeseries
        observables = compute_observables(hcat(x_component, y_component), μ[parameter_idx])
        lines!(Ax1, 1:1000, observables.optimal[1:1000], color = :red)
        lines!(Ax2, 1:1000, observables.perturbed[1:1000], color = CtpPink)
        lines!(Ax3, 1:1000, observables.suboptimal[1:1000], color = CtpTeal)
end

# Plot the timeseries and the solution of the nonlinear reconstruction of the potential
function plot_solutions(observables, equilibrium, parameter_idx)
        # Loop over the observables
        for n in 1:length(solutions)
                # Extract NLLS solution and the projecting subspace of the current observable
                solution = solutions[n]
                subspace = observables.subspaces[n]

                # Compute the (1-dimensional) restriction of the 2-dimensional potential onto the subspace
                U_obs(x) = U0((1/sqrt(1 + subspace^2))*x, (1/sqrt(1 + subspace^2))*subspace*x, μ[parameter_idx])

                # Reconstruct a shifted potential to match the ground truth
                #xs, Vs = shift_potential(U, μ[parameter_idx], subspace, solution, equilibrium)

                # Plot the ground truth and the reconstructed shifted potential
                if n == 1           # Optimal
                        domain = LinRange(-5,5,5000)
                        lines!(ax1, domain, [U_obs(x) for x in domain], color = (:black,1.0), linewidth = 3.0)
                        lines!(ax1, domain, [V(x, solution) for x in domain], color = (:red,0.25), linewidth = 2.0)
                elseif n == 2       # Perturbed
                        domain = LinRange(-5,5,5000)
                        lines!(ax2, domain, [U_obs(x) for x in domain], color = (:black,1.0), linewidth = 3.0)
                        lines!(ax2, domain, [V(x, solution) for x in domain], color = (CtpPink,0.25), linewidth = 2.0)
                else                # Suboptimal
                        domain = LinRange(-5,5,5000)
                        lines!(ax3, domain, [U_obs(x) for x in domain], color = (:black,1.0), linewidth = 3.0)
                        lines!(ax3, domain, [V(x, solution) for x in domain], color = (CtpTeal,0.25), linewidth = 2.0)
                end
        end
end

# Plot the escape rate early-warning signal 
function plot_ews()
        # Plot the ews of the ground truth
        μ_range = LinRange(-0.5, 1.35, 1000)
        ews_gt = zeros(Float64, length(μ_range))
        for m in 1:length(μ_range)
                equilibria = get_equilibria(f1, f2, μ_range[m])
                Ps = equilibria.stable[1]
                Pu = equilibria.unstable[1]
                ΔU = U0(Pu[1], Pu[2], μ_range[m]) - U0(Ps[1], Ps[2], μ_range[m])
                ews_gt[m] = exp(-ΔU)
        end
        lines!(ax4, μ_range, ews_gt, color = CtpRed, linewidth = 5.0)
        lines!(ax5, μ_range, ews_gt, color = CtpRed, linewidth = 5.0)
        lines!(ax6, μ_range, ews_gt, color = CtpRed, linewidth = 5.0)

        # Loop over the parameter values
        m = 1::Integer
        for v in μ
                # Import the ews 
                ews_obs = readin("ews/$m.csv")
                m = m + 1

                # Optimal observable
                ews = ews_obs[:,1]
                lower_cutoff = quantile(filter(isfinite, ews), 0.05)
                upper_cutoff = quantile(filter(isfinite, ews), 0.95)
                filtered_ews = filter(x -> lower_cutoff ≤ x ≤ upper_cutoff, ews)
                errorbars!(ax4, [v], [median(filtered_ews)], [median(filtered_ews) - quantile(filtered_ews,0.25)], [quantile(filtered_ews,0.75) - median(filtered_ews)], color = :black, whiskerwidth = 15, linewidth = 3.0)
                scatter!(ax4, v, median(filtered_ews), color = :red, markersize = 25, strokewidth = 2.0)

                # Perturbed observable
                ews = ews_obs[:,2]
                lower_cutoff = quantile(filter(isfinite, ews), 0.05)
                upper_cutoff = quantile(filter(isfinite, ews), 0.95)
                filtered_ews = filter(x -> lower_cutoff ≤ x ≤ upper_cutoff, ews)
                errorbars!(ax5, [v], [median(filtered_ews)], [median(filtered_ews) - quantile(filtered_ews,0.25)], [quantile(filtered_ews,0.75) - median(filtered_ews)], color = :black, whiskerwidth = 15, linewidth = 3.0)
                scatter!(ax5, v, median(filtered_ews), color = CtpPink, markersize = 25, strokewidth = 2.0)

                # Subptimal observable
                ews = ews_obs[:,3]
                lower_cutoff = quantile(filter(isfinite, ews), 0.05)
                upper_cutoff = quantile(filter(isfinite, ews), 0.95)
                filtered_ews = filter(x -> lower_cutoff ≤ x ≤ upper_cutoff, ews)
                errorbars!(ax6, [v], [median(filtered_ews)], [median(filtered_ews) - quantile(filtered_ews,0.25)], [quantile(filtered_ews,0.75) - median(filtered_ews)], color = :black, whiskerwidth = 15, linewidth = 3.0)
                scatter!(ax6, v, median(filtered_ews), color = CtpTeal, markersize = 25, strokewidth = 2.0)
        end

        # Export the figures
        savefig("ews_optimal.png", fig4)
        savefig("ews_perturbed.png", fig5)
        savefig("ews_suboptimal.png", fig6)
end
