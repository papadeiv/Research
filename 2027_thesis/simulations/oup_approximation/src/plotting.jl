"""
    Plotting

Routines to format and render the plot results of a simulation.
"""

# Axes formats
const border = 2.0
const labels = 20

# Generate and format a panel of the figure
function panel(position)
  return Axis(position,
              spinewidth=border,
              xgridvisible=false, ygridvisible=false,
              xtickalign=1, ytickalign=1,
              xtickwidth=border, ytickwidth=border,
              xticklabelsize=labels, yticklabelsize=labels,
              xlabelvisible=false, ylabelvisible=false,
              xticklabelsvisible=false, yticklabelsvisible=false,
              xlabelsize=labels, ylabelsize=labels
             )
end

# Create and format the figure
fig = Figure(size=(1800, 1350))
ax_paths = [panel(fig[1, j]) for j in (1, 3)]
ax_hists = [panel(fig[1, j]) for j in (2, 4)]
ax_means = [panel(fig[2, j]) for j in (1, 3)]
ax_variances = [panel(fig[3, j]) for j in (1, 3)]

# Set the 2:1 timeseries boxes and the square histogram boxes
for j in (1, 3)
  colsize!(fig.layout, j, Aspect(1, 2.0))
end
for j in (2, 4)
  colsize!(fig.layout, j, Aspect(1, 1.0))
end

# Set the axes ranges
for ax in ax_paths
  ax.limits = (0.0, n_path*dt, 0.8, 1.2)
end
for ax in ax_hists
  ax.limits = (0.0, nothing, 0.8, 1.2)
end
for ax in ax_means
  ax.limits = (0.0, T, 0.99, 1.01)
end
for ax in ax_variances
  ax.limits = (0.0, T, 0.002, 0.003)
end

# Link the vertical axes of each row
linkyaxes!(ax_paths[1], ax_hists[1], ax_paths[2], ax_hists[2])
linkyaxes!(ax_means...)
linkyaxes!(ax_variances...)

# Generate the plot
function plot(paths, t, sample_means, sample_variances, true_means, true_variances)
  # Collect the stationary distributions of the saddle-node and of the OUP
  distributions = [saddle_node_distribution(μ0), x -> oup_distribution(x, μ0)]

  for j in 1:2
    # Plot the first timesteps of the sample path
    lines!(ax_paths[j], paths[j].t[1:n_path+1], paths[j].u[1:n_path+1], color=:black)

    # Plot the asymptotic mean and one standard deviation either side of it
    hlines!(ax_paths[j], true_means[j], color=:red, linewidth=3)
    hlines!(ax_paths[j], true_means[j] .+ [-1, 1]*sqrt(true_variances[j]), color=:black, linestyle=:dash, linewidth=2)

    # Plot the histogram of the whole path with Scott's number of bins
    x = paths[j].u
    bins = ceil(Int, (maximum(x) - minimum(x))/(3.49*std(x)*length(x)^(-1/3)))
    hist!(ax_hists[j], x, bins=bins, direction=:x, color=(:grey, 0.35), strokewidth = 1.0, strokecolor = :black)

    # Overlay the stationary distribution, scaled from density to counts per bin
    y = range(minimum(x), maximum(x), length=1000)
    lines!(ax_hists[j], length(x)*(maximum(x) - minimum(x))/bins*distributions[j].(y), y, color=:red, linewidth=3)

    # Plot the sample mean and its asymptotic value
    lines!(ax_means[j], t, sample_means[j], color=:black, linewidth=2)
    hlines!(ax_means[j], true_means[j], color=:red, linewidth=3)

    # Plot the sample variance and its asymptotic value
    lines!(ax_variances[j], t, sample_variances[j], color=:black, linewidth=2)
    hlines!(ax_variances[j], true_variances[j], color=:red, linewidth=3)
  end

  # Fit the figure to the panels
  resize_to_layout!(fig)
end
