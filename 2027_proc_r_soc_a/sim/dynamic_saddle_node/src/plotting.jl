"""
    Plotting

Routines to format and render the plot results.
"""

# Axes formats
const border = 2.0
const labels = 20

# Generate and format the panels in the figure
function panel(fig, row, col)
  return Axis(fig[row, col],
              spinewidth=border,
              xgridvisible=false, ygridvisible=false,
              xtickalign=1, ytickalign=1,
              xtickwidth=border, ytickwidth=border,
              xticklabelsize=labels, yticklabelsize=labels)
end

# Generate the plot 
function plot(sample, rate, forecast, tipping_times, breaking_times)
  # Create and format the figure
  fig = Figure(size=(1500, 750))
  ax_top = panel(fig, 1, 1:2)
  ax_top.limits = (0.0, T, -0.5, 1.2)
  ax_top.xticks = []
  ax_bottom = panel(fig, 2, 1:2)
  ax_bottom.limits = (0.0, T, -0.2, 4.2)
  linkxaxes!(ax_top, ax_bottom)
  ax_tip = panel(fig, 1, 3)
  ax_breaking = panel(fig, 2, 3)

  # Plot the observation window, the breaking time and the tipping time
  vlines!(ax_top, [rate.t[1], rate.t[end]], color=:gray, linestyle=:dash, linewidth=2)
  vlines!(ax_bottom, [rate.t[1], rate.t[end]], color=:gray, linestyle=:dash, linewidth=2)
  vlines!(ax_top, τ, color=:black, linestyle=:dash, linewidth=2)
  vlines!(ax_bottom, τ, color=:black, linestyle=:dash, linewidth=2)
  vlines!(ax_breaking, τ, color=:black, linestyle=:dash, linewidth=2)
  vlines!(ax_top, τpb, color=:black, linewidth=2)
  vlines!(ax_bottom, τpb, color=:black, linewidth=2)
  vlines!(ax_tip, τpb, color=:black, linewidth=2)
  vlines!(ax_breaking, τpb, color=:black, linewidth=2)

  # Plot the sample path
  T1_idx = searchsortedfirst(sample.t, rate.t[1])
  T2_idx = searchsortedlast(sample.t, rate.t[end])
  lines!(ax_top, sample.t[1:T1_idx], sample.x[1:T1_idx], color=:grey)
  lines!(ax_top, sample.t[T1_idx:T2_idx], sample.x[T1_idx:T2_idx], color=:hotpink1)
  lines!(ax_top, sample.t[T2_idx:end], sample.x[T2_idx:end], color=:grey)

  # Plot the analytical and the median return rate and the extrapolated breaking time
  hlines!(ax_bottom, 0.0, color=(:black, 0.5))
  lines!(ax_bottom, [0.0, τ], [α(0.0), 0.0], color=:black, linewidth=2)
  lines!(ax_bottom, [rate.t[1], forecast.τ], [forecast.q + forecast.m*rate.t[1], 0.0], color=:red, linewidth=2.5)
  scatter!(ax_bottom, rate.t, rate.α, color=:hotpink1, markersize=16, strokecolor=:black, strokewidth=1)
  scatter!(ax_bottom, [forecast.τ], [0.0], color=:red, markersize=16, strokecolor=:black, strokewidth=1)

  # Plot the empirical and the analytical distribution of the tipping times
  low = minimum(tipping_times)
  high = maximum(tipping_times)
  ax_tip.limits = (low, high, 0.0, nothing)
  bins = ceil(Int, (high - low)/(3.49*std(tipping_times)*length(tipping_times)^(-1/3)))
  hist!(ax_tip, tipping_times, bins=bins, normalization=:pdf, color=(:hotpink1, 0.75))
  t = range(low, high, length=1000)
  lines!(ax_tip, t, tipping_distribution.(t), color=:teal, linewidth=4)

  # Plot the distribution of the breaking times, fading the bins past the tipping time
  low = quantile(breaking_times, 0.01)
  high = quantile(breaking_times, 0.99)
  ax_breaking.limits = (low, high, 0.0, nothing)
  w = (high - low)/25
  bins = τpb .+ w.*(floor((low - τpb)/w):ceil((high - τpb)/w))
  hist!(ax_breaking, breaking_times[low .<= breaking_times .< τpb], bins=bins, color=:red)
  hist!(ax_breaking, breaking_times[τpb .<= breaking_times .<= high], bins=bins, color=(:red, 0.5))

  # Make the histogram panels square
  colsize!(fig.layout, 3, Aspect(1, 1.0))
  resize_to_layout!(fig)

  # Export the figure
  save("dynamic_saddle_node.pdf", fig)
end
