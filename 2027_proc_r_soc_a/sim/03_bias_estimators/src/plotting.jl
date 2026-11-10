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
              xlabel=L"\mu", ylabel=label, ylabelrotation=0.0,
              xlabelvisible=(row == 2), xticklabelsvisible=(row == 2),
              ylabelvisible=(col == 1), yticklabelsvisible=(col != 2),
              xlabelsize=labels, ylabelsize=labels,
              xticklabelsize=labels, yticklabelsize=labels)
end

# Generate the plot 
function plot(estimates)
  # Create and format the figure
  fig = Figure(size=(1200, 700))
  ax_top = panel(fig, 1, 1, L"b_{\hat\alpha}", "OUP")
  ax_top.limits = (μ_fine[1], μ_fine[end], 0, 0.6)
  ax_top.xticks = []
  ax_top.yticks = [0, 0.6]
  ax_bottom = panel(fig, 2, 1, L"\hat\alpha", "")
  ax_bottom.limits = (-1.0, -0.05, 0.0, 9.0) 
  ax_bottom.xticks = [-1.0, -0.05] 
  ax_bottom.yticks = [0.0, 9.0] 

  # Compute the ground truth for the return rate
  μ_fine = collect(range(μ_set[1], μ_set[end], length=500))
  ground_truth = [α(μ) for μ in μ_fine]

  # Compute and plot the biases in the return rate
  lines!(ax_top, μ_fine, [compute_bias(μ).OUP for μ in μ_fine], color=:black, linewidth=2.0)

  # Compute and plot the summary statistics
  stats = compute_summary_statistics(estimates.OUP)
  band!(ax_bottom, μ_set, stats.low_bound, stats.high_bound, color=(:steelblue, 0.25))
  lines!(ax_bottom, μ_set, stats.low_bound, color=:steelblue, linewidth=2.0)
  lines!(ax_bottom, μ_set, stats.high_bound, color=:steelblue, linewidth=2.0)
  band!(ax_bottom, μ_set, stats.first_quant, stats.third_quant, color=:palegreen)
  lines!(ax_bottom, μ_set, stats.first_quant, color=:darkgreen, linewidth=1.0)
  lines!(ax_bottom, μ_set, stats.third_quant, color=:darkgreen, linewidth=1.0)
  lines!(ax_bottom, μ_set, stats.second_quant, color=:brown1, linewidth=2.0)
  lines!(ax_bottom, μ_fine, ground_truth, color=:black, linewidth=2.0)

  # Format the top panel (bias)
  ax_top = panel(fig, 1, 2, L"b_{\hat\alpha}", "AR1")
  ax_top.limits = (μ_fine[1], μ_fine[end], 0, 0.6)
  ax_top.xticks = []
  ax_top.yticks = [0, 0.6]

  # Compute and plot the biases in the return rate
  lines!(ax_top, μ_fine, [compute_bias(μ).AC1 for μ in μ_fine], color=:black, linewidth=2.0)

  # Format the bottom panel (estimators)
  ax_bottom = panel(fig, 2, 2, L"\hat\alpha", "")
  ax_bottom.limits = (-1.0, -0.05, 0.0, 9.0) 
  ax_bottom.xticks = [-1.0, -0.05] 
  ax_bottom.yticks = [0.0, 9.0] 

  # Compute and plot the summary statistics
  stats = compute_summary_statistics(estimates.AC1)
  band!(ax_bottom, μ_set, stats.low_bound, stats.high_bound, color=(:steelblue, 0.25))
  lines!(ax_bottom, μ_set, stats.low_bound, color=:steelblue, linewidth=2.0)
  lines!(ax_bottom, μ_set, stats.high_bound, color=:steelblue, linewidth=2.0)
  band!(ax_bottom, μ_set, stats.first_quant, stats.third_quant, color=:palegreen)
  lines!(ax_bottom, μ_set, stats.first_quant, color=:darkgreen, linewidth=1.0)
  lines!(ax_bottom, μ_set, stats.third_quant, color=:darkgreen, linewidth=1.0)
  lines!(ax_bottom, μ_set, stats.second_quant, color=:brown1, linewidth=2.0)
  lines!(ax_bottom, μ_fine, ground_truth, color=:black, linewidth=2.0)

  # Format the top panel (bias)
  ax_top = panel(fig, 1, 3, L"b_{\hat\alpha}", "OLS")
  ax_top.limits = (μ_fine[1], μ_fine[end], -0.4, 0.2)
  ax_top.xticks = []
  ax_top.yticks = [-0.4, 0.2]

  # Compute and plot the biases in the return rate
  hlines!(ax_top, 0.0, color=:gray, linewidth=1.0)
  lines!(ax_top, μ_fine, [compute_bias(μ).OLS for μ in μ_fine], color=:black, linewidth=2.0)

  # Format the bottom panel (estimators)
  ax_bottom = panel(fig, 2, 3, L"\hat\alpha", "")
  ax_bottom.limits = (-1.0, -0.05, 0.0, 9.0) 
  ax_bottom.xticks = [-1.0, -0.05] 
  ax_bottom.yticks = [0.0, 9.0] 

  # Compute and plot the summary statistics
  stats = compute_summary_statistics(estimates.OLS)
  band!(ax_bottom, μ_set, stats.low_bound, stats.high_bound, color=(:steelblue, 0.25))
  lines!(ax_bottom, μ_set, stats.low_bound, color=:steelblue, linewidth=2.0)
  lines!(ax_bottom, μ_set, stats.high_bound, color=:steelblue, linewidth=2.0)
  band!(ax_bottom, μ_set, stats.first_quant, stats.third_quant, color=:palegreen)
  lines!(ax_bottom, μ_set, stats.first_quant, color=:darkgreen, linewidth=1.0)
  lines!(ax_bottom, μ_set, stats.third_quant, color=:darkgreen, linewidth=1.0)
  lines!(ax_bottom, μ_set, stats.second_quant, color=:brown1, linewidth=2.0)
  lines!(ax_bottom, μ_fine, ground_truth, color=:black, linewidth=2.0)

  save("bias_estimators.pdf", fig)
  return fig
end
