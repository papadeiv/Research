include("../../../../inc/IO.jl")
include("../../../../inc/PlottingTools.jl")

# Import the data
N_observations = readin("../data/N_observations.csv")
E_c_hat = readin("../data/expected_exact_error.csv")
estimate = readin("../data/estimate.csv")

d = 2

# Create figure 
fig, axs = rowfig(2, 1,
                  bg_out = :transparent,
                  size = [3000,1000],
                  toggle_y_lab = false,
                  toggle_ticks_lab = [true, false],
                  x_ticks_lab_trunc = zeros(Int64, d+1),
                  y_ticks_lab_trunc = 3
                  )

# Loop over the axis
for n in 1:(d)
        #set_ticks(axs[n], N_observations, estimate[:,n+1])
        axs[n].xticks = [10,1000]
        axs[n].yticks = [0.009,0.077]
        # Plot the exact error from the numerics
        lines!(axs[n], N_observations, E_c_hat[:,n+1], linewidth = 5, color = (:red, 1.0))
        # Plot the analytical estimate (20)
        lines!(axs[n], N_observations, estimate[:,n+1], linewidth = 5, color = (:black, 1.0))
end

# Export the figure
save("../fig/estimate_20.png", fig)
