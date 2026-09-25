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
              ylabelvisible=(col == 1), yticklabelsvisible=(col != 2),
              xlabelsize=labels, ylabelsize=labels,
              xticklabelsize=labels, yticklabelsize=labels)
end

# Draw a sample in colour inside the admissible horizon and in grey outside
function admissible!(ax, t, y; kw...)
  masks = split_admissible(t)
  lines!(ax, t[masks.transient], y[masks.transient]; color=:grey, kw...)
  lines!(ax, t[masks.admissible], y[masks.admissible]; color=:hotpink1, kw...)
  lines!(ax, t[masks.inadmissible], y[masks.inadmissible]; color=:grey, kw...)
end

# Create and format the figure
function plot(branches, solutions)
  # Create the figure
  fig = Figure(size=(800, 650))

  # Format the top panel (solutions)
  ax_top = panel(fig, 1, 1, L"x", "")
  ax_top.limits = (0.0, T, -sqrt(-μ0), 1.0)
  ax_top.xticks = []

  # Plot the frozen branches, the pullback attractor and the noisy solution
  lines!(ax_top, branches.t, branches.stable, color=:blue, linewidth=3)
  lines!(ax_top, branches.t, branches.unstable, color=:red, linewidth=3)
  lines!(ax_top, solutions.pullback.t, solutions.pullback.x, color=:orange, linewidth=2)
  admissible!(ax_top, solutions.noisy.t, solutions.noisy.x, linewidth=1)

  # Format the bottom panel (return rates)
  ax_bottom = panel(fig, 2, 1, L"\alpha", "")
  ax_bottom.limits = (0.0, T, 0.0, 2.0)

  # Plot the frozen and the pullback attractor return rates
  admissible!(ax_bottom, branches.t, branches.rate, linewidth=2.5)
  admissible!(ax_bottom, solutions.pullback.t, solutions.pullback.rate,
              linewidth=2.5, linestyle=:dash)

  # Mark the end of the transient and the admissible time horizon
  for ax in (ax_top, ax_bottom)
    vlines!(ax, [tlo, Tc], color=:black, linestyle=:dash)
  end

  save("admissible_time_horizon.pdf", fig)
  return fig
end
