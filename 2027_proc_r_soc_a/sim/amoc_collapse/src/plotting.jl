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

# Plot the timeseries and its forecast
function plot_forecast(time, flow, rate, forecast)
  # Create and format the figure
  fig = Figure(size=(1000, 750))
  ax_top = panel(fig, 1)
  ax_top.limits = (0.0, time[end], nothing, nothing)
  ax_top.xticks = []
  ax_bottom = panel(fig, 2)
  ax_bottom.limits = (0.0, time[end], nothing, nothing)
  linkxaxes!(ax_top, ax_bottom)

  # Plot the observation window and the tipping time
  vlines!(ax_top, [rate.t[1], rate.t[end]], color=:gray, linestyle=:dash, linewidth=2)
  vlines!(ax_bottom, [rate.t[1], rate.t[end]], color=:gray, linestyle=:dash, linewidth=2)
  vlines!(ax_top, τ_tip, color=:black, linewidth=2)
  vlines!(ax_bottom, τ_tip, color=:black, linewidth=2)

  # Plot the observable timeseries
  T1_idx = searchsortedfirst(time, rate.t[1])
  T2_idx = searchsortedlast(time, rate.t[end])
  lines!(ax_top, time[1:T1_idx], flow[1:T1_idx], color=:grey)
  lines!(ax_top, time[T1_idx:T2_idx], flow[T1_idx:T2_idx], color=:hotpink1)
  lines!(ax_top, time[T2_idx:end], flow[T2_idx:end], color=:grey)

  # Plot the estimate return rate sequence and the extrapolated breaking time
  hlines!(ax_bottom, 0.0, color=(:black, 0.5))
  lines!(ax_bottom, [rate.t[1], forecast.τ], [forecast.q + forecast.m*rate.t[1], 0.0], color=:red, linewidth=2.5)
  scatter!(ax_bottom, rate.t, rate.α, color=:hotpink1, markersize=16, strokecolor=:black, strokewidth=1)
  scatter!(ax_bottom, [forecast.τ], [0.0], color=:red, markersize=16, strokecolor=:black, strokewidth=1)

  # Export the figure
  save("amoc_collapse.pdf", fig)
end
