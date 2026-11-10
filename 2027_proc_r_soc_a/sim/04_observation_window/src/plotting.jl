"""
    Plotting

Routines to format and render the plot results of a simulation.
"""

# Axes formats
const border = 2.0
const labels = 20

# Generate and format the panels in the figure
function panel(fig, row, col, xlabel, ylabel, yscale)
  return Axis(fig[row, col],
              yscale=yscale,
              spinewidth=border,
              xgridvisible=false, ygridvisible=false,
              xtickalign=1, ytickalign=1,
              xtickwidth=border, ytickwidth=border,
              xlabel=xlabel, ylabel=ylabel, ylabelrotation=0.0,
              ylabelvisible=(first(col) == 1), yticklabelsvisible=(first(col) == 1),
              xlabelsize=labels, ylabelsize=labels,
              xticklabelsize=labels, yticklabelsize=labels)
end

# Generate the plot 
function plot(samples)
  # Create the figure
  fig = Figure(size=(1200, 1200))

  # Compute the timesteps before and after the observation time
  t_before = collect(range(0.0, tn, length=500))
  t_after = τ .- (τ - tn).*exp10.(range(0.0, -6.0, length=500))

  # Format the top panel (window lengths)
  ax_top = panel(fig, 1, 1:3, "", "", log10)
  ax_top.limits = (0.0, 1.02*τ, 0.5*L(0.0), 2*U(0.0))

  # Plot the breaking time and the admissible time horizon
  vlines!(ax_top, τ, color=:gray, linewidth=2.0)
  vlines!(ax_top, Tc, color=:gray, linestyle=:dash, linewidth=2.0)

  # Plot the lower bound, the upper bound and the fixed window length
  lines!(ax_top, t_before, L.(t_before), color=:blue, linewidth=3.0)
  lines!(ax_top, t_after, L.(t_after), color=(:blue, 0.5), linewidth=3.0)
  lines!(ax_top, t_before, U.(t_before), color=:red, linewidth=3.0)
  lines!(ax_top, t_after, U.(t_after), color=(:red, 0.5), linewidth=3.0)
  lines!(ax_top, [0.0, tn], [Tw, Tw], color=:purple, linewidth=3.0)
  lines!(ax_top, [tn, Tx], [Tw, Tw], color=(:purple, 0.5), linewidth=3.0)

  # Plot the window lengths at the observation time
  scatter!(ax_top, [tn], [L(tn)], color=:blue, markersize=20, strokecolor=:black, strokewidth=1)
  scatter!(ax_top, [tn], [Tw], color=:purple, markersize=20, strokecolor=:black, strokewidth=1)
  scatter!(ax_top, [tn], [U(tn)], color=:red, markersize=20, strokecolor=:black, strokewidth=1)

  # Format the middle panel (autocorrelation)
  s = Δ.*(0:(length(samples.L) - 1))
  ax_middle = panel(fig, 2, 1, "", "", identity)
  ax_middle.limits = (0.0, s[end], 0.0, 1.0)

  # Plot the correlation function across the window
  lines!(ax_middle, s, exp.(-s./L(tn)), color=:blue, linewidth=3.0)

  # Format the bottom panel (eigenvalue change)
  ax_bottom = panel(fig, 3, 1, "", "", log10)
  ax_bottom.limits = (0.0, tn, 1e-3, 3.0)

  # Compute and plot the relative change of the leading eigenvalue
  hlines!(ax_bottom, 1.0, color=:gray, linewidth=2.0)
  lines!(ax_bottom, t_before, [(λ(t) - λ(t - L(t)))/abs(λ(t)) for t in t_before], color=:blue, linewidth=3.0)

  # Format the middle panel (autocorrelation)
  s = Δ.*(0:(length(samples.Tw) - 1))
  ax_middle = panel(fig, 2, 2, "", "", identity)
  ax_middle.limits = (0.0, s[end], 0.0, 1.0)

  # Plot the correlation function across the window
  lines!(ax_middle, s, exp.(-s./L(tn)), color=:purple, linewidth=3.0)

  # Format the bottom panel (eigenvalue change)
  ax_bottom = panel(fig, 3, 2, "", "", log10)
  ax_bottom.limits = (0.0, tn, 1e-3, 3.0)

  # Compute and plot the relative change of the leading eigenvalue
  hlines!(ax_bottom, 1.0, color=:gray, linewidth=2.0)
  lines!(ax_bottom, t_before, [(λ(t) - λ(t - Tw))/abs(λ(t)) for t in t_before], color=:purple, linewidth=3.0)

  # Format the middle panel (autocorrelation)
  s = Δ.*(0:(length(samples.U) - 1))
  ax_middle = panel(fig, 2, 3, "", "", identity)
  ax_middle.limits = (0.0, s[end], 0.0, 1.0)

  # Plot the correlation function across the window
  lines!(ax_middle, s, exp.(-s./L(tn)), color=:red, linewidth=3.0)

  # Format the bottom panel (eigenvalue change)
  ax_bottom = panel(fig, 3, 3, "", "", log10)
  ax_bottom.limits = (0.0, tn, 1e-3, 3.0)

  # Compute and plot the relative change of the leading eigenvalue
  hlines!(ax_bottom, 1.0, color=:gray, linewidth=2.0)
  lines!(ax_bottom, t_before, [(λ(t) - λ(t - U(t)))/abs(λ(t)) for t in t_before], color=:red, linewidth=3.0)

  save("observation_window.pdf", fig)
  #return fig
end
