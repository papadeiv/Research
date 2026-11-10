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

# Plot the samples and the distributions of the tipping times
function plot(Ω_stationary, Ω_small, Ω_large, pullback_small, pullback_large)
  # Create the figure
  fig = Figure(size=(1300, 700))

  # Stationary problem: format the panels
  ax_top = panel(fig, 1, 1)
  ax_top.limits = (0.0, stationary.T, -2.0, 2.0)
  ax_top.xticks = []
  ax_bottom = panel(fig, 2, 1)
  ax_bottom.limits = (0.0, maximum(Ω_stationary.tip), 0.0, nothing)
  ax_bottom.ytickformat = v -> [@sprintf("%.1e", y) for y in v]

  # Stationary problem: plot the samples and the branches
  for path in Ω_stationary.paths
    lines!(ax_top, path.t, path.u, color=(:black, 0.5))
  end
  hlines!(ax_top, xs(0.0, stationary), color=:blue, linewidth=4)
  hlines!(ax_top, -xs(0.0, stationary), color=:red, linewidth=4)

  # Stationary problem: plot the empirical and the analytical distribution
  tip = Ω_stationary.tip
  t = range(0.0, maximum(tip), length=1000)
  hist!(ax_bottom, tip, bins=ceil(Int, (maximum(tip) - minimum(tip))/(3.49*std(tip)*length(tip)^(-1/3))),
        normalization=:pdf, color=:hotpink1)
  lines!(ax_bottom, t, p_stationary.(t, Ref(stationary)), color=:teal, linewidth=4)

  # Small noise: format the panels
  ax_top = panel(fig, 1, 2)
  ax_top.limits = (0.0, small_noise.T, -2.0, 2.0)
  ax_top.xticks = []
  low = min(minimum(Ω_small.tip), τ(small_noise))
  high = max(maximum(Ω_small.tip), τpb(small_noise))
  ax_bottom = panel(fig, 2, 2)
  ax_bottom.limits = (low, high, 0.0, nothing)
  ax_bottom.ytickformat = v -> [@sprintf("%.1e", y) for y in v]

  # Small noise: mark the window of opportunity
  vspan!(ax_top, τ(small_noise), τpb(small_noise), color=(:black, 0.15))
  vspan!(ax_bottom, τ(small_noise), τpb(small_noise), color=(:black, 0.15))
  vlines!(ax_top, τ(small_noise), color=:black, linestyle=:dash, linewidth=2)
  vlines!(ax_bottom, τ(small_noise), color=:black, linestyle=:dash, linewidth=2)
  vlines!(ax_top, τpb(small_noise), color=:black, linewidth=2)
  vlines!(ax_bottom, τpb(small_noise), color=:black, linewidth=2)

  # Small noise: plot the samples, the branches and the pullback attractor
  for path in Ω_small.paths
    lines!(ax_top, path.t, path.u, color=(:black, 0.5))
  end
  t = range(0.0, τ(small_noise), length=400)
  lines!(ax_top, t, xs.(t, Ref(small_noise)), color=:blue, linewidth=4)
  lines!(ax_top, t, -xs.(t, Ref(small_noise)), color=:red, linewidth=4)
  lines!(ax_top, pullback_small.t, pullback_small.x, color=:gold3, linewidth=2)

  # Small noise: plot the empirical and the analytical distribution
  tip = Ω_small.tip
  t = range(low, high, length=1000)
  hist!(ax_bottom, tip, bins=ceil(Int, (maximum(tip) - minimum(tip))/(3.49*std(tip)*length(tip)^(-1/3))),
        normalization=:pdf, color=:hotpink1)
  lines!(ax_bottom, t, p_small_noise.(t, Ref(small_noise)), color=:teal, linewidth=4)

  # Large noise: format the panels
  ax_top = panel(fig, 1, 3)
  ax_top.limits = (0.0, large_noise.T, -2.0, 2.0)
  ax_top.xticks = []
  low = min(minimum(Ω_large.tip), τ(large_noise))
  high = max(maximum(Ω_large.tip), τpb(large_noise))
  ax_bottom = panel(fig, 2, 3)
  ax_bottom.limits = (low, high, 0.0, nothing)
  ax_bottom.ytickformat = v -> [@sprintf("%.1e", y) for y in v]

  # Large noise: mark the window of opportunity
  vspan!(ax_top, τ(large_noise), τpb(large_noise), color=(:black, 0.15))
  vspan!(ax_bottom, τ(large_noise), τpb(large_noise), color=(:black, 0.15))
  vlines!(ax_top, τ(large_noise), color=:black, linestyle=:dash, linewidth=2)
  vlines!(ax_bottom, τ(large_noise), color=:black, linestyle=:dash, linewidth=2)
  vlines!(ax_top, τpb(large_noise), color=:black, linewidth=2)
  vlines!(ax_bottom, τpb(large_noise), color=:black, linewidth=2)

  # Large noise: plot the samples, the branches and the pullback attractor
  for path in Ω_large.paths
    lines!(ax_top, path.t, path.u, color=(:black, 0.5))
  end
  t = range(0.0, τ(large_noise), length=400)
  lines!(ax_top, t, xs.(t, Ref(large_noise)), color=:blue, linewidth=4)
  lines!(ax_top, t, -xs.(t, Ref(large_noise)), color=:red, linewidth=4)
  lines!(ax_top, pullback_large.t, pullback_large.x, color=:gold3, linewidth=2)

  # Large noise: plot the empirical and the analytical distribution
  tip = Ω_large.tip
  t = range(low, high, length=1000)
  hist!(ax_bottom, tip, bins=ceil(Int, (maximum(tip) - minimum(tip))/(3.49*std(tip)*length(tip)^(-1/3))),
        normalization=:pdf, color=:hotpink1)
  lines!(ax_bottom, t, p_large_noise.(t, Ref(large_noise)), color=:teal, linewidth=4)

  # Export the figure
  save("tipping_time_distribution.pdf", fig)
end
