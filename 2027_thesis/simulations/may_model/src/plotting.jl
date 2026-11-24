"""
    Plotting

Routines to format and render the plot results of a simulation.
"""

# Axes formats
const border = 2.0
const labels = 20

# Generate and format the panels in the figure
function panel(fig, row)
  return Axis(fig[row, 1],
              spinewidth=border,
              xgridvisible=false, ygridvisible=false,
              xtickalign=1, ytickalign=1,
              xtickwidth=border, ytickwidth=border,
              xticklabelsize=labels, yticklabelsize=labels)
end

# Generate the plot 
function plot(sample_path)
  # Extract the state and timestep
  x = sample_path.u
  t = sample_path.t

  # Compute the parameter shift
  p = [μ(s) for s in t]

  # Find the tipping times
  T1_idx = findfirst(x -> x > 2.604, p) 
  T2_idx = length(p) - findfirst(x -> x > 1.787, reverse(p)) 

  # Create and format the figure
  fig = Figure(size=(900, 450))
  ax_top = panel(fig, 1)
  ax_top.limits = (0.0, T, nothing, nothing)
  ax_top.xticks = []
  ax_bottom = panel(fig, 2)
  ax_bottom.limits = (0.0, T, nothing, nothing)
  linkxaxes!(ax_top, ax_bottom)

  # Plot the tipping times and saddle-node values 
  hlines!(ax_bottom, p[T1_idx], color=:grey, linewidth=1)
  hlines!(ax_bottom, p[T2_idx], color=:grey, linewidth=1)
  vlines!(ax_top, [t[T1_idx], t[T2_idx]], color=:black, linestyle=:dash, linewidth=2)
  vlines!(ax_bottom, [t[T1_idx], t[T2_idx]], color=:black, linestyle=:dash, linewidth=2)

  # Plot the timeseries 
  lines!(ax_top, t, x, color=:black, linewidth=1)
  lines!(ax_bottom, t, p, color=:black, linewidth=2)

  save("may_model.pdf", fig)
end
