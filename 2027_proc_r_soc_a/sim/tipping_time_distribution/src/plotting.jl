"""
    Plotting

Routines to format and render the plot results of a simulation.
"""

# Axes formats
const border = 2.0
const labels = 20

# Generate and format the panels in the figure
function panel(fig, row, col, label, title)
  return Axis(fig[row, col],
              title=title,
              spinewidth=border,
              xgridvisible=false, ygridvisible=false,
              xtickalign=1, ytickalign=1,
              xtickwidth=border, ytickwidth=border,
              xlabel=L"t", ylabel=label, ylabelrotation=0.0,
              xlabelvisible=(row == 2), xticklabelsvisible=(row == 2),
              ylabelvisible=(col == 1),
              xlabelsize=labels, ylabelsize=labels,
              xticklabelsize=labels, yticklabelsize=labels)
end

# Mark the window of opportunity between the breaking and the tipping time
function opportunity!(ax, P)
  vspan!(ax, τ(P), τpb(P), color=(:black, 0.15))
  vlines!(ax, τ(P), color=:black, linestyle=:dash, linewidth=2)
  vlines!(ax, τpb(P), color=:black, linewidth=2)
end

# Create and format the figure
function plot(results)
  # Create the figure
  fig = Figure(size=(1300, 700))

  # Loop over the experiments
  for (col, P) in enumerate(runs)
    result = results[col]

    # Format the top panel (solutions)
    ax_top = panel(fig, 1, col, L"x", "")
    ax_top.limits = (0.0, P.T, -2.0, 2.0)
    ax_top.xticks = []

    # Plot the solutions of the ensemble
    for path in result.paths
      lines!(ax_top, thin(path.t, 400), thin(path.x, 400), color=(:black, 0.5))
    end

    # Plot the stable and the unstable branch up to the breaking time
    tb = collect(range(0.0, min(τ(P), P.T), length=400))
    xb = sqrt.(max.(-(P.μ0 .+ P.ε .* tb), 0.0))
    lines!(ax_top, tb, xb, color=:blue, linewidth=4.0)
    lines!(ax_top, tb, -xb, color=:red, linewidth=4.0)

    # Plot the pullback attractor and mark the window of opportunity
    if P.ε > 0
      opportunity!(ax_top, P)
      lines!(ax_top, result.pullback..., color=:gold3, linewidth=2.0)
    end

    # Format the bottom panel (tipping times)
    low, high = compute_range(result.ttip, P)
    ax_bottom = panel(fig, 2, col, "density", "")
    ax_bottom.limits = (low, high, 0.0, nothing)
    ax_bottom.ytickformat = v -> [@sprintf("%.1e", y) for y in v]

    # Mark the window of opportunity and plot the empirical distribution
    P.ε > 0 && opportunity!(ax_bottom, P)
    hist!(ax_bottom, result.ttip, bins=scott(result.ttip),
          normalization=:pdf, color=(:hotpink1, 1.00))

    # Plot the predicted density
    t = collect(range(low, high, length=1000))
    lines!(ax_bottom, t, compute_density(t, P, col), color=:teal, linewidth=4)
  end

  save("tipping_time_distribution.pdf", fig)
  return fig
end
