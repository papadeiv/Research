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

# Create and format the figure
function plot(estimates)
  # Create the figure and initialize the panels
  fig = Figure(size=(1200, 700))
  panels = ((:OUP, "From OUP variance"),
            (:AC1, "From AR(1) fit"),
            (:OLS, "From LLS solutions"))

  # Compute the ground truth for the return rate
  μ_fine = collect(range(μ_set[1], μ_set[end], length=500))
  ground_truth = [α(μ) for μ in μ_fine]

  # Loop over the panels
  for (col, (key, title)) in enumerate(panels)

    # Format the top panel (bias)
    ax_top = panel(fig, 1, col, L"b_{\hat\alpha}", title)
    col == 3 ? ax_top.limits = (μ_fine[1], μ_fine[end], -0.4, 0.2) : ax_top.limits = (μ_fine[1], μ_fine[end], 0, 0.6)
    ax_top.xticks = []
    col == 3 ? ax_top.yticks = [-0.4, 0.2] : ax_top.yticks = [0, 0.6]

    # Compute and plot the biases in the return rate
    col == 3 && hlines!(ax_top, 0.0, color=:gray, linewidth=1.0)
    lines!(ax_top, μ_fine, [getproperty(compute_bias(μ), key) for μ in μ_fine], color=:black, linewidth=2.0)

    # Format the bottom panel (estimators)
    ax_bottom = panel(fig, 2, col, L"\hat\alpha", "")
    ax_bottom.limits = (-1.0, -0.05, 0.0, 9.0) 
    ax_bottom.xticks = [-1.0, -0.05] 
    ax_bottom.yticks = [0.0, 9.0] 

    # Compute and plot the summary statistics
    stats = compute_summary_statistics(getproperty(estimates, key))
    band!(ax_bottom, μ_set, stats.low_bound, stats.high_bound, color=(:steelblue, 0.25))
    lines!(ax_bottom, μ_set, stats.low_bound, color=:steelblue, linewidth=2.0)
    lines!(ax_bottom, μ_set, stats.high_bound, color=:steelblue, linewidth=2.0)
    band!(ax_bottom, μ_set, stats.first_quant, stats.third_quant, color=:palegreen)
    lines!(ax_bottom, μ_set, stats.first_quant, color=:darkgreen, linewidth=1.0)
    lines!(ax_bottom, μ_set, stats.third_quant, color=:darkgreen, linewidth=1.0)
    lines!(ax_bottom, μ_set, stats.second_quant, color=:brown1, linewidth=2.0)
    lines!(ax_bottom, μ_fine, ground_truth, color=:black, linewidth=2.0)
  end

  save("bias_estimators.pdf", fig)
  return fig
end
