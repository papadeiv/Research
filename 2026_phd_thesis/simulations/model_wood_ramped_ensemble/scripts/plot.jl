"""
    Plotting script

Functions used to create the plots in each figure.
"""

# Create empty layouts for the figures
include("./figs.jl")

# Plot the reconstruction of the shifted potential
function plot_solution(μf)
        # Plot the shifted ground truth
        domain = collect(LinRange(-1, 1, 1000))

        # Loop over the solutions of the nonlinear least-squares reconstruction
        for solution in solutions
                # Extract the relevant data 
                μ = solution[1]
                c = solution[2:4]

                # Plot the reconstruction
                lines!(ax1, domain, [V(x, c) for x in domain], color = (:black, 0.5), linewidth = 3.0)
        end
end

# Plot the escape rate early-warning signal 
function plot_ews()
        # Loop over the parameter values
        for n in 1:length(H0)
                # Import the random variables 
                data = import_data(n)
                coefficients = data.coefficients
                ews = data.ews

                # Filter out the outliers and plot the ews 
                lower_cutoff = quantile(filter(isfinite, ews), 0.05)
                upper_cutoff = quantile(filter(isfinite, ews), 0.95)
                filtered_ews = filter(x -> lower_cutoff ≤ x ≤ upper_cutoff, ews)
                errorbars!(ax3, [H0[n]], [median(filtered_ews)], [median(filtered_ews) - quantile(filtered_ews,0.25)], [quantile(filtered_ews,0.75) - median(filtered_ews)], color = :black, whiskerwidth = 15, linewidth = 3.0)
                scatter!(ax3, H0[n], median(filtered_ews), color = CtpBlue, markersize = 25, strokewidth = 2.0)

                # Print the percentile distribution
                println("μ = ", H0[n])
                println("Lower cutoff = ", lower_cutoff, ". Upper cutoff = ", upper_cutoff)
                println("Q1 = ", quantile(filtered_ews,0.25), ", Q2 = ", median(filtered_ews), ", Q3 = ", quantile(filtered_ews,0.75))
                println("------------------------------------------------------------------")
        end

        # Export the figure
        savefig("wood/ews.png", fig3)
end
