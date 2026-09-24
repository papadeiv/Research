"""
    Postprocessing script

Collection of quantities and functions used to postprocess and analyse the results of a simulation.
"""

# Flag for tipped states
tipped = false

# Initialize data structures for the analysis
lower_envelope = Matrix{Real}(undef, Nμ, 3)
first_quartile = Matrix{Real}(undef, Nμ, 3)
median_estimate = Matrix{Real}(undef, Nμ, 3)
third_quartile = Matrix{Real}(undef, Nμ, 3)
upper_envelope = Matrix{Real}(undef, Nμ, 3)

# Perform the statistical analysis of the generated samples
function analysis()
        # Loop over the parameter values
        printstyled("Analyzing the samples\n"; bold=true, underline=true, color=:light_blue)
        @showprogress for (idx_μ, μ) in enumerate(μ_set)
                # Import the data 
                data = readin("solutions/$idx_μ.csv")

                # Compute the lower envelope of the three estimates 
                lower_envelope[idx_μ,1] = (extrema(data[:,1]))[1]
                lower_envelope[idx_μ,2] = (extrema(data[:,2]))[1]
                lower_envelope[idx_μ,3] = (extrema(data[:,3]))[1]

                # Compute the 1st quartile of the three estimates
                first_quartile[idx_μ,1] = quantile(data[:,1], 0.25)
                first_quartile[idx_μ,2] = quantile(data[:,2], 0.25)
                first_quartile[idx_μ,3] = quantile(data[:,3], 0.25)

                # Compute the median (2nd quartile) of the three estimates
                median_estimate[idx_μ,1] = median(data[:,1])
                median_estimate[idx_μ,2] = median(data[:,2])
                median_estimate[idx_μ,3] = median(data[:,3])

                # Compute the 3rd quartile of the three estimates
                third_quartile[idx_μ,1] = quantile(data[:,1], 0.75)
                third_quartile[idx_μ,2] = quantile(data[:,2], 0.75)
                third_quartile[idx_μ,3] = quantile(data[:,3], 0.75)

                # Compute the upper envelope of the three estimates
                upper_envelope[idx_μ,1] = (extrema(data[:,1]))[2]
                upper_envelope[idx_μ,2] = (extrema(data[:,2]))[2]
                upper_envelope[idx_μ,3] = (extrema(data[:,3]))[2]
        end

        # Plot the range 
        band!(ax1, μ_set, lower_envelope[:,1], upper_envelope[:,1], color = (:steelblue, 0.25))
        band!(ax2, μ_set, lower_envelope[:,2], upper_envelope[:,2], color = (:steelblue, 0.25))
        band!(ax3, μ_set, lower_envelope[:,3], upper_envelope[:,3], color = (:steelblue, 0.25))

        # Plot the envelopes
        lines!(ax1, μ_set, lower_envelope[:,1], color = :steelblue, linewidth = 2.0)
        lines!(ax1, μ_set, upper_envelope[:,1], color = :steelblue, linewidth = 2.0)
        lines!(ax2, μ_set, lower_envelope[:,2], color = :steelblue, linewidth = 2.0)
        lines!(ax2, μ_set, upper_envelope[:,2], color = :steelblue, linewidth = 2.0)
        lines!(ax3, μ_set, lower_envelope[:,3], color = :steelblue, linewidth = 2.0)
        lines!(ax3, μ_set, upper_envelope[:,3], color = :steelblue, linewidth = 2.0)

        # Plot the interquartile range
        band!(ax1, μ_set, first_quartile[:,1], third_quartile[:,1], color = :palegreen)
        band!(ax2, μ_set, first_quartile[:,2], third_quartile[:,2], color = :palegreen)
        band!(ax3, μ_set, first_quartile[:,3], third_quartile[:,3], color = :palegreen)

        # Plot the first and third quantiles 
        lines!(ax1, μ_set, first_quartile[:,1], color = :darkgreen, linewidth = 1.0)
        lines!(ax1, μ_set, third_quartile[:,1], color = :darkgreen, linewidth = 1.0)
        lines!(ax2, μ_set, first_quartile[:,2], color = :darkgreen, linewidth = 1.0)
        lines!(ax2, μ_set, third_quartile[:,2], color = :darkgreen, linewidth = 1.0)
        lines!(ax3, μ_set, first_quartile[:,3], color = :darkgreen, linewidth = 1.0)
        lines!(ax3, μ_set, third_quartile[:,3], color = :darkgreen, linewidth = 1.0)

        # Plot the median estimates
        lines!(ax1, μ_set, median_estimate[:,1], color = :brown1, linewidth = 2.0)
        lines!(ax2, μ_set, median_estimate[:,2], color = :brown1, linewidth = 2.0)
        lines!(ax3, μ_set, median_estimate[:,3], color = :brown1, linewidth = 2.0)
        
        # Compute the ground truth on an expanded parameter domain (up to SN bifurcation)
        exp_μ_set = LinRange(μ_set[1], 1.5396, 10000)
        ground_truth = Vector{Real}(undef, length(exp_μ_set))
        for (idx_μ, μ) in enumerate(exp_μ_set)
                # Compute the stable equilibrium
                xs = (get_equilibria(f, μ, domain=[-10,10])).stable[2]
                # Compute the leading eigenvalue
                ground_truth[idx_μ] = J(xs)
        end

        # Plot the ground truth
        lines!(ax1, exp_μ_set, ground_truth, color = :black, linewidth = 2.0)
        lines!(ax2, exp_μ_set, ground_truth, color = :black, linewidth = 2.0)
        lines!(ax3, exp_μ_set, ground_truth, color = :black, linewidth = 2.0)

        # Export the figure
        savefig("comparison.pdf", fig)
end
