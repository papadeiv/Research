"""
    Plotting

Routines to format and render the plot results of a simulation.
"""

# Axes formats
const border = 2.0
const labels = 20

# Create and format the figure
fig = Figure(size=(900, 450))
ax = Axis(fig[1, 1],
          spinewidth=border,
          xgridvisible=false, ygridvisible=false,
          xtickalign=1, ytickalign=1,
          xtickwidth=border, ytickwidth=border,
          xticklabelsize=labels, yticklabelsize=labels,
          xlabelvisible=false, ylabelvisible=false,
          xlabelsize=labels, ylabelsize=labels
         )

ax.limits = (μ_in, μ_f, 0.0, nothing)

# Generate the plot
function plot(ground_truth, sample_mean, lower_envelope, upper_envelope)
  # Collect the parameter space
  μ = collect(μ_range)

  band!(ax, μ, lower_envelope, upper_envelope, color=(:grey, 0.35))
  lines!(ax, μ, ground_truth, color=:red, linewidth=3)
  lines!(ax, μ, sample_mean, color=:black, linewidth=2)
end
