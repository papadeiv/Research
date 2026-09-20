"""
    Figures layout

Generation of the layouts and formats of the figures.
"""

# Specify the figure dimensions and style
fig = Figure(; size = (900, 900), figure_padding = (30,30,30,30))

# Set thickness and size of axes elements 
border = 2.0
labels = 20
ticks = 22

# Timeseries and bifurcation diagram 
ax1 = Axis(fig[1,1], 
           spinewidth = border,
           xgridvisible = false,
           ygridvisible = false,
           xlabelvisible = false,
           ylabelvisible = true,
           xticklabelsvisible = false,
           yticklabelsvisible = true,
           xtickalign = 1,
           ytickalign = 1,
           xtickwidth = border,
           ytickwidth = border,
           xlabel = L"t",
           ylabel = L"X(t)",
           xlabelsize = labels,
           ylabelsize = labels,
           xticklabelsize = labels,
           yticklabelsize = labels,
          )

# Early-warning signal 
ax2 = Axis(fig[2,1], 
           spinewidth = border,
           xgridvisible = false,
           ygridvisible = false,
           xlabelvisible = true,
           ylabelvisible = true,
           xticklabelsvisible = true,
           yticklabelsvisible = true,
           xtickalign = 1,
           ytickalign = 1,
           xtickwidth = border,
           ytickwidth = border,
           xlabel = L"t",
           ylabel = L"\alpha(t)",
           xlabelsize = labels,
           ylabelsize = labels,
           xticklabelsize = labels,
           yticklabelsize = labels,
          )
