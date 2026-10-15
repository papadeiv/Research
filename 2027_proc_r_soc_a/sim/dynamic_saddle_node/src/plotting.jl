"""
    Plotting

Routines to format and render the plot results of a simulation.
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

# Mark the admissible horizon, the breaking time and the tipping time
function horizon!(ax)
  vlines!(ax, [tw[1], tw[end]], color=:gray, linestyle=:dash, linewidth=2)
  vlines!(ax, τ, color=:black, linestyle=:dash, linewidth=2)
  vlines!(ax, τpb, color=:black, linewidth=2)
end

# Draw a sample in colour inside the admissible horizon and in grey outside
function admissible!(ax, t, y)
  masks = split_admissible(t)
  lines!(ax, t[masks.before], y[masks.before], color=:grey)
  lines!(ax, t[masks.admissible], y[masks.admissible], color=:hotpink1)
  lines!(ax, t[masks.after], y[masks.after], color=:grey)
end

# Create and format the figure
function plot(sample, estimates, uncertainties, breaking_times, tipping_times)
  # Create the figure
  fig = Figure(size=(1500, 750))

  # Plot the sample in the top panel
  ax_top = panel(fig, 1, 1:2)
  ax_top.limits = (0.0, T, -0.5, 1.2)
  ax_top.xticks = []
  admissible!(ax_top, sample.t, sample.x)
  horizon!(ax_top)

  # Extrapolate the median return rate from the first window to its zero crossing
  ym = compute_median(estimates[:, 1:Nm])
  fit = estimate_breaking_time(tw, ym, compute_median(uncertainties[:, 1:Nm])./Nm)
  te = collect(range(tw[1], fit.τ, length=50))

  # Plot the analytical and the estimated return rates in the bottom panel
  ax_bottom = panel(fig, 2, 1:2)
  ax_bottom.limits = (0.0, T, -0.2, 4.2)
  lines!(ax_bottom, [0.0, τ], [α(0.0), 0.0], color=:black, linewidth=2)
  hlines!(ax_bottom, 0.0, color=(:black, 0.5))
  lines!(ax_bottom, te, fit.q .+ fit.m .* te, color=:red, linewidth=2.5)
  scatter!(ax_bottom, tw, ym, color=:hotpink1, markersize=16, strokecolor=:black, strokewidth=1)
  scatter!(ax_bottom, [fit.τ], [0.0], color=:red, markersize=16, strokecolor=:black, strokewidth=1)
  horizon!(ax_bottom)

  # Plot the empirical and the analytical distribution of the tipping times in the top right panel
  low, high = extrema(tipping_times)
  ax_tip = panel(fig, 1, 3)
  ax_tip.limits = (low, high, 0.0, nothing)
  hist!(ax_tip, tipping_times, bins=scott(tipping_times), normalization=:pdf, color=(:hotpink1, 0.75))
  t = collect(range(low, high, length=1000))
  lines!(ax_tip, t, p_small_noise.(t), color=:teal, linewidth=4)
  vlines!(ax_tip, τpb, color=:black, linewidth=2)

  # Plot the distribution of the breaking times between the 1st and 99th percentiles in the bottom right panel
  τ̂ = filter(isfinite, breaking_times)
  low, high = quantile(τ̂, 0.01), quantile(τ̂, 0.99)
  ax_breaking = panel(fig, 2, 3)
  ax_breaking.limits = (low, high, 0.0, nothing)

  # Align the bin edges with the tipping time and fade the bins beyond it
  w = (high - low)/25
  edges = τpb .+ w .* (floor((low - τpb)/w):ceil((high - τpb)/w))
  hist!(ax_breaking, filter(x -> low <= x < τpb, τ̂), bins=edges, color=:red)
  hist!(ax_breaking, filter(x -> τpb <= x <= high, τ̂), bins=edges, color=(:red, 0.5))
  vlines!(ax_breaking, τ, color=:black, linestyle=:dash, linewidth=2)
  vlines!(ax_breaking, τpb, color=:black, linewidth=2)

  # Make the histogram panels square and shrink the figure around the layout
  colsize!(fig.layout, 3, Aspect(1, 1.0))
  resize_to_layout!(fig)

  save("dynamic_saddle_node.pdf", fig)
  return fig
end
