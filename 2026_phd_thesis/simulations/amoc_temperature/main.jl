"""
    Main script

Run this file to execute the simulation, analyse and plot the results.
"""

# Import the necessary packages and local modules
include("inc.jl")

# Import the simulation's scripts
include("./scripts/proc.jl")
include("./scripts/plot.jl")

# Define the main algorithm
function main()
        # Data structures
        solutions = Vector{Vector{Float64}}()           # Solutions of the inference method 
        time_idx = Vector{Float64}()                    # Index of timesteps 
        ews = Vector{Float64}()                         # Statistical distribution over the ensemble

        # Import the AMOC timeseries
        data = NCDataset("../../res/data/amoc/Atlantic_250-500m_year_0001-2200.nc")
        time = data["time"][:]
        latitude = data["lat"][:]
        observable = data["TEMP"]
        N = size(observable, 1)
        
        # Create timeseries figure
        fig = Figure()
        ax = Axis(fig[1,1], xlabel = L"\textbf{years}", ylabel = L"\textbf{temperature (250-500m)}")

        # Loop over the ensemble's sample paths
        ews = Float64[]
        @showprogress for n in 1:N
                # Extract windowed subseries
                t = copy(time[1:idx])
                u = copy(observable[n,1:idx]) 

                # Plot the timeseries
                lines!(ax, t, u, color = n, colormap = :managua, colorrange = (1,N), linewidth = 1.0, label = "$(trunc(Int, latitude[n]))")
                lines!(ax, time[(idx+1):end], observable[n,(idx+1):end], color = (:black,0.35), linewidth = 1.0)

                # Detrend the timeseries 
                detrended_solution = detrend(u, alg = "emd", n_modes=1)
                trend = detrended_solution.trend
                residual = detrended_solution.residuals

                # Solve the nonlinear least-squares problem to fit a cubic potential
                solution = fit_potential(residual, transformation=[0.0,1.0,8.0], optimiser=β, attempts=Na, verbose=true).fit

                # Perform postprocessing analysis on the solutions
                push!(ews, analyse(solution))
                
                # Plot the timeseries's stationary distribution and the nonlinear solution
                plot_solutions(t, u, trend, residual, solution, n)
        end

        # Plot the location of the tipping point
        lines!(ax, [time[idx], time[idx]], [3, 17], color = :red, linestyle = :dash, linewidth = 2.0)

        # Plot the legend
        grid = GridLayout(fig[2,1], tellwidth = false)
        Legend(grid[1,1], ax, "Latitude (deg N)", orientation = :horizontal, nbanks = 5, halign = :right, valign = :bottom, tellwidth = true, marker = :circle, markersize = 0.1, labelsize = 10.0)

        # Plot the ews
        fig2 = Figure()
        ax2 = Axis(fig2[1,1], xlabel = L"\textbf{latitude (deg N)}", ylabel = L"\textbf{ews}", xticks = [-30,-20,-20,-10,0,10,20,30,40,50,60])
        for n in 1:length(ews)
                scatter!(ax2, latitude[n], ews[n], color = n, colormap = :managua, colorrange = (1,N), markersize = 10.0, strokecolor = :black, strokewidth = 1.0)
        end

        # Export the figure
        savefig("amoc/temperature.png", fig)
        savefig("amoc/ews.png", fig2)
end

# Execute the main
main()
