"""
    Figures layout

Generation of the layouts and formats of the figures.
"""

#----------------#
#    Figure 1    #
#----------------#

fig1, ax1 = makefig(size = [2000,1000],
                    pad = (20,35,20,20), # Order is: left, right, bottom, top 
                    bg_out = :white,
                    box_position = [1,1:2],
                    limits = ((0,2200),(-7,25)),
                    lab = [L"\textbf{time (millenia)}", L"\textbf{observable (Sv)}"],
                    toggle_lab = [false,true],
                    x_ticks = [0,1000,2000],
                    y_ticks = [-7,25],
                    toggle_ticks_lab = [false,true],
                    ticks_lab_trunc = [0,0]
                   )

fig2, ax2 = makefig(fig = fig1,
                    box_position = [1:2,3:4],
                    lab = [L"\mathbf{x}", L"\mathbf{V(x)}"],
                    limits = ((-10, 10), (-1, 10)),
                    x_ticks = [-10, 0, 10],
                    y_ticks = [0, 10],
                    ticks_lab_trunc = [0,0]
                   )

fig3, ax3 = makefig(fig = fig1,
                    box_position = [2,1:2],
                    limits = ((0,2200), (-0.05,1.05)),
                    lab = [L"\textbf{time (millenia)}", L"\textbf{exp}\mathbf{(-ΔV)}"],
                    x_ticks = [0,1000,2000],
                    y_ticks = [0,0.25,0.5,0.75,1],
                    ticks_lab_trunc = [0,2]
                   )

#----------------#
#    Figure 2    #
#----------------#

fig4, ax4 = makefig(size = [2000,1000],
                    pad = (20,35,20,20), # Order is: left, right, bottom, top 
                    bg_out = :white,
                    box_position = [1,1:2],
                    lab = [L"\mathbf{t}", L"\mathbf{x_t}"],
                    toggle_lab = [false,true],
                    ticks_lab_trunc = [0,0],
                    toggle_ticks_lab = [false,true],
                   )

fig5, ax5 = makefig(fig = fig4,
                    box_position = [1,3],
                    limits = ((nothing,nothing),(-0.01,0.22)),
                    lab = [L"\mathbf{x_t}", L"\textbf{distribution}"],
                    toggle_lab = [false,true],
                    y_ticks = [0,0.2],
                    ticks_lab_trunc = [0,1]
                   )

fig6, ax6 = makefig(fig = fig4,
                    box_position = [2,1:2],
                    limits = ((nothing,nothing),(-8,8)),
                    lab = [L"\textbf{time (years)}", L"\mathbf{x_t}"],
                    y_ticks = [-7,0,7],
                    ticks_lab_trunc = [0,0]
                   )

fig7, ax7 = makefig(fig = fig4,
                    box_position = [2,3],
                    limits = ((-8,8),(-0.01,0.22)),
                    lab = [L"\mathbf{x_t}", L"\textbf{distribution}"],
                    x_ticks = [-8,0,8],
                    y_ticks = [0,0.2],
                    ticks_lab_trunc = [0,1]
                   )
