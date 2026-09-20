"""
    Postprocessing script

Collection of quantities and functions used to postprocess and analyse the results of a simulation.
"""

# Perform the statistical analysis of the generated samples
function analysis()
        # Compute the ground truth of the leading eigenvalue
        xs = (get_equilibria(f, μ, domain=[-10,10])).stable[2]
        λ = J(xs)

        # Simulation parameters loop
        relative_errors = Matrix{Real}(undef, Ns, 3)
        for (idx_sim, (stepsize, steps)) in enumerate(Iterators.product(dt,Nt))
                # Import the data
                data = readin("accuracy/$stepsize/$steps.csv")

                # Compute the median of the three estimates
                c = median(data[:,1])
                α = median(data[:,2])
                θ = median(data[:,3])

                # Compute the relative errors between the median estimates and the ground truth
                relative_errors[idx_sim,1] = abs(λ-c)/abs(λ)
                relative_errors[idx_sim,2] = abs(λ-α)/abs(λ)
                relative_errors[idx_sim,3] = abs(λ-θ)/abs(λ)
        end

        # Extract the minimum and maximum relative errors
        err_min, err_max = minimum(relative_errors), maximum(relative_errors)
        println("------- Err. min                  Err.max")
        println("LLS:    $(minimum(relative_errors[:,1]))     $(maximum(relative_errors[:,1]))")
        println("LLS:    $(minimum(relative_errors[:,2]))     $(maximum(relative_errors[:,2]))")
        println("LLS:    $(minimum(relative_errors[:,3]))     $(maximum(relative_errors[:,3]))")

        # Transformation of the coordinate axes in logscale so that it does not interpolate negative values
        logedges(v) = (l = log10.(v); exp10.([l[1] - (l[2]-l[1])/2;
                      (l[1:end-1] .+ l[2:end]) ./ 2;
                      l[end] + (l[end]-l[end-1])/2]))
        x, y = logedges(dt), logedges(Nt)
        cmap = :thermal
        
        # LLS -----------------------------------------------------------------------------------------------------
        error_map = reshape(relative_errors[:,1], length(dt), length(Nt))
        heatmap!(ax1, x, y, error_map, colorrange = (err_min, err_max), colorscale = log10, colormap = cmap)

        # AR(1) ---------------------------------------------------------------------------------------------------
        error_map = reshape(relative_errors[:,2], length(dt), length(Nt))
        hm = heatmap!(ax2, logedges(dt), logedges(Nt), error_map, colorrange = (err_min, err_max), colorscale = log10, colormap = cmap)

        # OUP -----------------------------------------------------------------------------------------------------
        error_map = reshape(relative_errors[:,3], length(dt), length(Nt))
        heatmap!(ax3, logedges(dt), logedges(Nt), error_map, colorrange = (err_min, err_max), colorscale = log10, colormap = cmap)

        # Plot the location of the setup of Figure 1 in the error map
        scatter!(ax1, 1e-1, 1e+3, color = :brown1, markersize = 20, strokewidth = 1.5)
        scatter!(ax2, 1e-1, 1e+3, color = :brown1, markersize = 20, strokewidth = 1.5)
        scatter!(ax3, 1e-1, 1e+3, color = :brown1, markersize = 20, strokewidth = 1.5)

        # Add the colorbar 
        Colorbar(fig[1,4], 
                 size = 25,
                 label = "relative error",
                 spinewidth = border,
                 scale = log10,
                 hm,
                 #limits = (trunc(err_min,digits=1), trunc(err_max,digits=1)), 
                 #ticks = [trunc(err_min,digits=1), trunc(err_max,digits=1)], 
                )

        # Export the figure
        savefig("accuracy.pdf", fig)
end
