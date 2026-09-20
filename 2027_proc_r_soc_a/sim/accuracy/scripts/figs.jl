"""
    Figures layout

Generation of the layouts and formats of the figures.
"""

# Specify the figure dimensions and style
fig = Figure(; size = (1200, 400))

# Set thickness and size of axes elements 
border = 2.0
labels = 20
ticks = 22

ax1 = Axis(fig[1,1], 
           title = "From LLS solutions",
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
           xlabel = L"\Delta t", 
           ylabel = L"N_t", 
           xlabelsize = labels,
           ylabelsize = labels,
           xticklabelsize = labels,
           yticklabelsize = labels,
           xscale = log10, 
           yscale = log10
          )

ax2 = Axis(fig[1,2], 
           title = "From AR(1) fit",
           spinewidth = border,
           xgridvisible = false,
           ygridvisible = false,
           xlabelvisible = true,
           ylabelvisible = false,
           xticklabelsvisible = true,
           yticklabelsvisible = false,
           xtickalign = 1,
           ytickalign = 1,
           xtickwidth = border,
           ytickwidth = border,
           xlabel = L"\Delta t", 
           ylabel = L"N_t",
           xlabelsize = labels,
           ylabelsize = labels,
           xticklabelsize = labels,
           yticklabelsize = labels,
           xscale = log10, 
           yscale = log10
          )

ax3 = Axis(fig[1,3], 
           title = "From OUP variance",
           spinewidth = border,
           xgridvisible = false,
           ygridvisible = false,
           xlabelvisible = true,
           ylabelvisible = false,
           xticklabelsvisible = true,
           yticklabelsvisible = false,
           xtickalign = 1,
           ytickalign = 1,
           xtickwidth = border,
           ytickwidth = border,
           xlabel = L"\Delta t", 
           ylabel = L"N_t",
           xlabelsize = labels,
           ylabelsize = labels,
           xticklabelsize = labels,
           yticklabelsize = labels,
           xscale = log10, 
           yscale = log10
          )
