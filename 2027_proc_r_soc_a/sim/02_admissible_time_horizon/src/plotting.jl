"""
    Plotting

Routines to format and render the plot results of a simulation.
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

# Generate the plot 
function plot(t_eq, x_eq, α_eq, t_pb, x_pb, t, x)
  # Create and format the figure
  fig = Figure(size=(800, 650))
  ax_top = panel(fig, 1)
  ax_top.limits = (0.0, T, -sqrt(-μ0), 1.0)
  ax_top.xticks = []
  ax_bottom = panel(fig, 2)
  ax_bottom.limits = (0.0, T, 0.0, 2.0)
  linkxaxes!(ax_top, ax_bottom)

  # Plot the admissible time horizon
  vlines!(ax_top, [Tb, Tc], color=:black, linestyle=:dash)
  vlines!(ax_bottom, [Tb, Tc], color=:black, linestyle=:dash)

  # Plot the quasi-dteady equilibria, the pullback attractor and the sample path 
  Tb_idx = searchsortedfirst(t, Tb)
  Tc_idx = searchsortedlast(t, Tc)
  lines!(ax_top, t_eq, x_eq, color=:blue, linewidth=3)
  lines!(ax_top, t_eq, -x_eq, color=:red, linewidth=3)
  lines!(ax_top, t_pb, x_pb, color=:orange, linewidth=2)
  lines!(ax_top, t[1:Tb_idx], x[1:Tb_idx], color=:grey, linewidth=1)
  lines!(ax_top, t[Tb_idx:Tc_idx], x[Tb_idx:Tc_idx], color=:hotpink1, linewidth=1)
  lines!(ax_top, t[Tc_idx:end], x[Tc_idx:end], color=:grey, linewidth=1)

  # Plot the return rates
  Tb_eq_idx = searchsortedfirst(t_eq, Tb)
  Tc_eq_idx = searchsortedlast(t_eq, Tc)
  lines!(ax_bottom, t_eq[1:Tb_eq_idx], α_eq[1:Tb_eq_idx], color=:grey, linewidth=2.5)
  lines!(ax_bottom, t_eq[Tb_eq_idx:Tc_eq_idx], α_eq[Tb_eq_idx:Tc_eq_idx], color=:hotpink1, linewidth=2.5)
  lines!(ax_bottom, t_eq[Tc_eq_idx:end], α_eq[Tc_eq_idx:end], color=:grey, linewidth=2.5)
  lines!(ax_bottom, t_pb[1:Tb_idx], α_pb(x_pb[1:Tb_idx]), color=:grey, linewidth=2.5, linestyle=:dash)
  lines!(ax_bottom, t_pb[Tb_idx:Tc_idx], α_pb(x_pb[Tb_idx:Tc_idx]), color=:hotpink1, linewidth=2.5, linestyle=:dash)
  lines!(ax_bottom, t_pb[Tc_idx:end], α_pb(x_pb[Tc_idx:end]), color=:grey, linewidth=2.5, linestyle=:dash)

  save("admissible_time_horizon.pdf", fig)
end
