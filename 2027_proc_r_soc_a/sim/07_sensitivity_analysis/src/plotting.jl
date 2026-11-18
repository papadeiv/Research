"""
    Plotting

Routines to format and render the plot results.
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

# Plot the forecasted tipping time against the window length
function plot_sensitivity(window, forecast)
  # Create and format the figure
  fig = Figure(size=(1000, 500))
  ax = panel(fig, 1)
  ax.limits = (window[1], window[end], 0, 5000)

  # Plot the forecasted tipping time
  band!(ax, window, mean(forecast) - std(forecast), mean(forecast) + std(forecast), color=(:grey, 0.25))
  lines!(ax, window, forecast, color=:hotpink, linewidth=3)
  hlines!(ax, mean(forecast), color=:red, linewidth=3)
  hlines!(ax, τ_tip, color=:black, linewidth=3)
  scatter!(ax, 25, 1755.13, color=:yellow, markersize=16, strokecolor=:black, strokewidth=1, marker = :star5)

  # Export the figure
  save("sensitivity_analysis.pdf", fig)
end
