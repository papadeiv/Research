include("../../../../inc/IO.jl")
include("../../../../inc/PlottingTools.jl")

# Import the data
N_observations = readin("../data/N_observations.csv")
estimate_8_1 = readin("../data/estimate_8_1.csv")
estimate_8_2 = readin("../data/estimate_8_2.csv")
estimate_9 = readin("../data/estimate_9.csv")
estimate_11 = readin("../data/estimate_11.csv")
c_hat = readin("../data/exact_error.csv")

d = 2

# Create figure for estimates (8)
fig, axsL = rowfig(2, 1,
                   bg_out = :transparent,
                   size = [3000,1000],
                   toggle_y_lab = false,
                   x_ticks_lab_trunc = zeros(Int64, d+1),
                   y_ticks_lab_trunc = 3
                  )
# Create the mirrored axis for the exact errors
axsR = [mirror_axis(fig, box_position = [1,n], toggle_lab = false, ticks_lab_trunc = 3) for n in 1:(d)]

# Loop over the axis
for n in 1:(d)
        axsL[n].xticks = [10,1000]
        axsL[n].yticks = [0,maximum(estimate_8_1[:,n+1]),maximum(estimate_8_2[:,n+1])]
        # Customise the ticks
        set_ticks(axsR[n], N_observations, c_hat[:,n+1])
        # Plot the exact error from the numerics
        lines!(axsR[n], N_observations, c_hat[:,n+1], linewidth = 5, color = (:red, 1.0))
        # Plot the analytical estimate (8.1)
        lines!(axsL[n], N_observations, estimate_8_1[:,n+1], linewidth = 5, color = (:black, 1.0))
        # Plot the analytical estimate (8.2)
        lines!(axsL[n], N_observations, estimate_8_2[:,n+1], linewidth = 5, color = (:black, 1.0), linestyle = :dash)
end

# Export the figure
save("../fig/estimates_8.png", fig)

# Create figure for estimates (8)
fig, axsL = rowfig(2, 1,
                   bg_out = :transparent,
                   size = [3000,1000],
                   toggle_y_lab = false,
                   toggle_ticks_lab = [true, false],
                   x_ticks_lab_trunc = zeros(Int64, d+1),
                   y_ticks_lab_trunc = 2
                  )
# Create the mirrored axis for the exact errors
#axsR = [mirror_axis(fig, box_position = [1,n], toggle_lab = false, ticks_lab_trunc = 3) for n in 1:(d)]

# Loop over the axis
for n in 1:(d)
        axsL[n].xticks = [10,1000]
        axsL[n].yticks = [0,0.5,1.0,1.5,1.85]
        # Customise the ticks
        #set_ticks(axsL[n], N_observations, c_hat[:,n+1])
        #set_ticks(axsL[n], N_observations, estimate_8_1[:,n+1])
        # Plot the exact error from the numerics
        lines!(axsL[n], N_observations, c_hat[:,n+1], linewidth = 5, color = (:red, 1.0))
        # Plot the analytical estimate (8.1)
        lines!(axsL[n], N_observations, estimate_9[:,n+1], linewidth = 5, color = (:black, 1.0))
        # Plot the analytical estimate (8.2)
        lines!(axsL[n], N_observations, estimate_11[:,n+1], linewidth = 5, color = (:black, 1.0), linestyle = :dash)
end

# Export the figure
save("../fig/estimates_9&11.png", fig)
