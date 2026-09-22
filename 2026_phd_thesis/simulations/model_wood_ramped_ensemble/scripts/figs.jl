"""
    Figures layout

Generation of the layouts and formats of the figures.
"""

#----------------#
#    Figure 1    #
#----------------#

fig1, ax1 = makefig(size = [900,900],
                    pad = (20,60,20,20), # Order is: left, right, bottom, top 
                    bg_out = :white,
                    box_position = [1,1],
                    lab = [L"\mathbf{x}", L"\mathbf{V(x\,;\,\,\mu(t))}"],
                    limits = ((-0.5, 0.5), (-0.5, 0.5)),
                    x_ticks = [-0.5, 0, 0.5],
                    y_ticks = [-0.5, 0, 0.5],
                    ticks_lab_trunc = [0,0]
                   )

#----------------#
#    Figure 2    #
#----------------#

fig3, ax3 = makefig(size = [1500,1000],
                    pad = (30,35,20,20), # Order is: left, right, bottom, top 
                    bg_out = :white,
                    box_position = [1,1],
                    limits = ((-0.375,0.385),(-0.05,1.05)),
                    lab = [L"\mathbf{\mu}", L"\textbf{exp}\mathbf{(-ΔV)}"],
                    x_ticks = H0,
                    y_ticks = [0,0.25,0.5,0.75,1],
                    ticks_lab_trunc = [2,2]
                   )
