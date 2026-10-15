"""
    Plotting

Routines to format and render the plot results of a simulation.
"""

# Axes formats
const border = 2.0
const labels = 20

# Generate and format the panels in the figure
function panel(fig, row, col, xlabel, ylabel, title)
  return Axis(fig[row, col],
              title=title,
              spinewidth=border,
              xgridvisible=false, ygridvisible=false,
              xtickalign=1, ytickalign=1,
              xtickwidth=border, ytickwidth=border,
              xlabel=xlabel, ylabel=ylabel, ylabelrotation=0.0,
              xlabelvisible=(row != 1), xticklabelsvisible=(row != 1),
              ylabelvisible=(first(col) == 1),
              xlabelsize=labels, ylabelsize=labels,
              xticklabelsize=labels, yticklabelsize=labels)
end

# Mark the window of opportunity between the breaking and the tipping time
function opportunity!(ax)
  vspan!(ax, τ, τpb, color=(:black, 0.15))
  vlines!(ax, τ, color=:black, linestyle=:dash, linewidth=2)
  vlines!(ax, τpb, color=:black, linewidth=2)
end

# Create and format the figure
function plot(solutions, estimates, uncertainties, breaking_times)
  # Create the figure and initialize the estimators
  fig = Figure(size=(1250, 1100))
  panels = ((:OUP, "OUP", :blue),
            (:AC1, "AC1", :green),
            (:OLS, "LLS", :red))

  # Compute the ground truth for the return rate
  t_fine = collect(range(0.0, τ, length=500))
  ground_truth = [α(t) for t in t_fine]

  # Format the top panel (solutions)
  ax_top = panel(fig, 1, 1:3, L"t", L"x", "")
  ax_top.limits = (0.0, T, -1.2, 1.2)
  ax_top.xticks = []

  # Plot the solutions of the ensemble and highlight the first one
  for path in solutions
    lines!(ax_top, path.t, path.x, color=(:gray, 0.4))
  end
  lines!(ax_top, solutions[1].t, solutions[1].x, color=:black)
  opportunity!(ax_top)

  # Format the middle panel (return rates)
  ax_middle = panel(fig, 2, 1:3, L"t", L"\hat\alpha", "")
  ax_middle.limits = (0.0, T, -0.2, 4.2)

  # Loop over the estimators
  for (key, label, color) in panels
    # Compute the ensemble median of the estimates and the variance of the median
    ym = compute_median(getproperty(estimates, key))
    um = compute_median(getproperty(uncertainties, key))./Ne

    # Extrapolate the median return rate to its zero crossing
    fit = estimate_breaking_time(tw, ym, um)
    te = collect(range(tw[end], fit.τ, length=50))

    # Plot the median return rate, its extrapolation and the breaking time
    lines!(ax_middle, tw, ym, color=color, linewidth=2.5)
    lines!(ax_middle, te, fit.q .+ fit.m .* te, color=color, linestyle=:dash, linewidth=2.5)
    scatter!(ax_middle, [fit.τ], [0.0], color=color, markersize=14)
  end

  # Plot the ground truth and mark the breaking time
  lines!(ax_middle, t_fine, ground_truth, color=:black, linestyle=:dot, linewidth=2)
  hlines!(ax_middle, 0.0, color=(:black, 0.5))
  opportunity!(ax_middle)

  # Loop over the bottom panels (breaking times)
  for (col, (key, label, color)) in enumerate(panels)
    τ̂ = filter(isfinite, getproperty(breaking_times, key))

    # Format the panel
    ax_bottom = panel(fig, 3, col, L"\hat\tau", "count", label)
    low, high = compute_range(τ̂)
    ax_bottom.limits = (low, high, 0.0, nothing)

    # Plot the distribution of the breaking times within the range and mark the true one
    hist!(ax_bottom, filter(x -> low <= x <= high, τ̂), bins=25, color=(color, 0.65))
    vlines!(ax_bottom, τ, color=:black, linestyle=:dash, linewidth=2)
  end

  save("dynamic_saddle_node.pdf", fig)
  return fig
end
