"""
    Figures layout

Generation of the layouts and formats of the figures.
"""

#----------------#
#    Figure 1    #
#----------------#

fig1, ax1 = makefig(size = [1600,600],
                    pad = (20,40,20,20), # Order is: left, right, bottom, top 
                    bg_out = :white,
                    box_position = [1,1:2],
                    lab = [L"\mathbf{t}", L"\mathbf{x_t}"],
                    ticks_lab_trunc = [0,1]
                   )

fig2, ax2 = makefig(fig = fig1,
                    box_position = [1,3],
                    limits = ((-1.5, 1.5), (-2, 2)),
                    lab = [L"\mathbf{x}", L"\mathbf{V(x)}"],
                    x_ticks = [-1.5,0,1.5],
                    y_ticks = [-2,0,2],
                    ticks_lab_trunc = [1,0]
                   )

#----------------#
#    Figure 2    #
#----------------#

fig3, ax3 = makefig(size = [1500,1000],
                    pad = (30,35,20,20), # Order is: left, right, bottom, top 
                    bg_out = :white,
                    box_position = [1,1],
                    limits = ((-0.05,1.35),(-0.05,1.05)),
                    lab = [L"\mathbf{\mu}", L"\textbf{exp}\mathbf{(-ΔV)}"],
                    x_ticks = μ,
                    y_ticks = [0,0.25,0.5,0.75,1],
                    ticks_lab_trunc = [1,2]
                   )
