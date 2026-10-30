"""
    Plotting

Routines to format and render the plot results of a simulation.
"""

# Axes formats
const border = 2.0
const labels = 20

# Create and format the figure
fig = Figure(size=(900, 450))
ax_left = Axis(fig[1,1],
               spinewidth=border,
               xgridvisible=false, ygridvisible=false,
               xtickalign=1, ytickalign=1,
               xtickwidth=border, ytickwidth=border,
               xticklabelsize=labels, yticklabelsize=labels,
              )

ax_left.limits = (-1, 0.3, -1, 1.5)

ax_right = Axis(fig[1,2],
               spinewidth=border,
               xgridvisible=false, ygridvisible=false,
               xtickalign=1, ytickalign=1,
               xtickwidth=border, ytickwidth=border,
               xticklabelsize=labels, yticklabelsize=labels,
               ylabelvisible = false
              )

ax_right.limits = (-1, 0.3, -1, 1.5)
ax_right.yticks = []

linkyaxes!(ax_left, ax_right)
