function show(L::Lattice, states::Vector{Vector{Float64}}, eigenvalues::Vector{Float64}, eigenvectors::Array{Float64})
        # Define sliding window parameters
        width = length(states) - length(eigenvalues)
        time = 0.0
        dt = 1e-2
        # Define vector storing the proposed spatial ews
        ews = Vector{Float64}() 
        # Define vector storing the parameter values
        params = Vector{Float64}()
        # Define counter for eigenvectors 
        ctr = 1
        CairoMakie.activate!()
        for j=width:(length(states)-1)
                # Get the current system's state values
                snapshot = states[j][1:end-3]
                snapshot = reshape(snapshot, (L.rows, L.cols))
                # Get the current eigenvector's state values
                eigenvec = eigenvectors[:, ctr]
                eigenvec = reshape(eigenvec, (L.rows, L.cols))
                # Get the current leading eigenvalue
                push!(ews, eigenvalues[ctr])
                # Get the model's parameters values
                parameter = round(states[j][end-2], digits=3)
                push!(params, parameter)
                noise = states[j][end-1]
                drift = states[j][end]

                # Reshape the state and eignevector for the heatmap
                state = Array{Float64}(undef, L.rows, L.cols)
                mode = Array{Float64}(undef, L.rows, L.cols)
                for k=0:(L.rows-1)
                        state[k+1,:] = snapshot[L.rows-k,:]
                        mode[k+1,:] = eigenvec[L.rows-k,:]
                end

                # Plot the lattice dynamical system evolution
                fig1 = Figure(; size = (1000, 1000), backgroundcolor = "#eeeeeeff")
                ax1 = Axis(fig1[1,1],
                          #title = L"\mu = %$parameter",
                          titlesize = 25,
                          backgroundcolor = :white,
                          spinewidth = 5.0,
                          xticksvisible = false,
                          yticksvisible = false,
                          xticklabelsvisible = false,
                          yticklabelsvisible = false,
                          aspect = DataAspect()
                )
                heatmap!(ax1, state', colorrange = (0,1))

                #=
                # Plot the leading eigenvector evolution 
                fig2 = Figure(; size = (1000, 1000), backgroundcolor = :transparent)
                ax2 = Axis(fig2[1,1],
                          #title = L"window\; size = %$width",
                          backgroundcolor = :transparent,
                          xticksvisible = false,
                          yticksvisible = false,
                          xticklabelsvisible = false,
                          yticklabelsvisible = false,
                          aspect = DataAspect()
                )
                heatmap!(ax2, mode', colorrange = ((findmin(mode))[1],(findmax(mode))[1]))

                # Plot the leading eigenvalue evolution 
                fig3 = Figure(; size = (1000, 400))
                ax3 = Axis(fig3[1,1],
                           limits = ((0, dt*(length(eigenvalues)-width)),((findmin(eigenvalues))[1],(findmax(eigenvalues))[1])) #(0.98, 1.0))
                )
                lines!(ax3, LinRange(0.0, time, length(ews)), ews)
                =#

                # Export the figures
                save("../results/lattice/lung_ventilation/snapshots/$ctr.png", fig1)

                # Update the time variable
                time += dt
                # Update the image counter
                ctr += 1
        end
        #=
        # Plot the leading eigenvalue evolution 
        CairoMakie.activate!(; px_per_unit = 3)
        fig3 = Figure(size = (1200, 600), backgroundcolor = :transparent)
        ax3 = Axis(fig3[1,1],
                # Background
                backgroundcolor = :transparent,
                xgridvisible = false,
                ygridvisible = false,
                limits = ((params[1], params[end]), (0.9985, 1.0015)),
                # Title
                #title = (L"\mu=%$(round(Y; digits=3))"),
                titlevisible = false,
                titlesize = 22,
                titlealign = :center,
                titlegap = 4.0,
                # x-axis
                xlabel = L"\mu",
                xlabelvisible = true,
                xlabelsize = 24,
                xlabelcolor = :black,
                xlabelpadding = -20.0,
                xticks = [params[1], 0.7, 0.9, 1.08, 1.3, params[end]],
                xticksvisible = true,
                xticksize = 6,
                xticklabelsvisible = true,
                xticklabelsize = 18,
                xtickformat = "{:.1f}",
                xscale = identity, #log10,
                xaxisposition = :bottom,
                # y-axis
                ylabel = L"\text{leading eigenvalue}",
                ylabelvisible = true,
                ylabelsize = 22,
                ylabelcolor = :black,
                ylabelpadding = -20.0,
                yticks = [0.9985, 1.0015],
                yticksvisible = true,
                yticksize = 6,
                yticklabelsvisible = true,
                yticklabelsize = 18,
                ytickformat = "{:.3f}",
                yscale = identity,
                yaxisposition = :left,
        )
        rowsize!(fig3.layout, 1, Aspect(1, 0.25))
        resize_to_layout!(fig3)

        # This is for the lung ventilation model
        lines!(ax3, LinRange(params[1],params[end],length(eigenvalues)), eigenvalues, linewidth = 2)
        lines!(ax3, 1.08.*ones(2), [0.9985,1.0015], color = :black, linestyle = :dash)
        
        # This is for the vegetation turbidity model
        #lines!(ax3, LinRange(3.0,3.1,length(eigenvalues)), eigenvalues.*1e4, linewidth = 2)
        #lines!(ax3, 3.04.*ones(2), [eigenvalues[argmin(eigenvalues)]*1e4, eigenvalues[argmax(eigenvalues)]*1e4], color = :black, linestyle = :dash)
        #text!(3.01, eigenvalues[argmax(eigenvalues)]*1e4, text = L"\times 10^{-4}", align = [:right, :top], color = :black, fontsize = 28)
        save("../results/lattice/lung_ventilation/LungVentilationEWS.png", fig3)
        =#
end
