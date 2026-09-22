"""
    Main script

Run this file to execute the simulation, analyse and plot the results.
"""

# Import the necessary packages and local modules
include("inc.jl")

# Import the simulation's scripts
include("./scripts/proc.jl")

# Define the main algorithm
function main()
        # Import the AMOC timeseries
        data = NCDataset("../../res/data/amoc/Atlantic_250-500m_year_0001-2200.nc")
        time = data["time"][:]
        latitude = data["lat"][:]
        temperature = data["TEMP"]
        salinity = data["SALT"]

        # Loop over the different latitudes
        for N in 1:length(latitude)
                # Extract subseries until tipping
                t = copy(time[1:idx])
                T = copy(temperature[N,1:idx]) 
                S = copy(salinity[N,1:idx]) 

                # Define the observable (flow)
                ϕ = T - S
                detrended_solution = detrend(ϕ, alg = "emd", n_modes=1)
                trend = detrended_solution.trend
                residual = detrended_solution.residuals

                ensemble = preprocess_solution(t, residual, window_size)
                Ne = length(ensemble.trajectories)
                Nt = length(ensemble.trajectories[1])

                # Generate figure and axis
                fig = Figure(size = (1200,800))
                ax1 = Axis(fig[1,1], xlabel = L"\textbf{years}", ylabel = L"\textbf{flow (250-500m)}", limits = (time[1], time[end], -32,-20), title = "Latitude = $(latitude[N])")
                lines!(ax1, t, ϕ, color = :red, linewidth = 1.0)
                lines!(ax1, time[(idx+1):end], temperature[N,(idx+1):end] - salinity[N,(idx+1):end], color = (:black,0.35), linewidth = 1.0)
                lines!(ax1, [time[idx], time[idx]], [-32, -20], color = :red, linestyle = :dash, linewidth = 2.0)
                ax2 = Axis(fig[2,1], xlabel = L"\textbf{years}", ylabel = L"\textbf{early-warning signal}", limits = (time[1], time[end], -0.1, 1.1))

                lines!(ax2, [time[idx], time[idx]], [-0.1, 1.1], color = :red, linestyle = :dash, linewidth = 2.0)

                # Loop over the ensemble
                printstyled("Latitude = $(latitude[N]): computing the least-squares solutions across the ensemble\n"; bold=true, underline=true, color=:light_blue)
                @showprogress for n in 1:Ne
                        # Extract windowed subseries
                        tw = copy(ensemble.timesteps[n])
                        ϕw = copy(ensemble.trajectories[n]) 

                        # Assemble histogram of the windowed subseries
                        σ = std(ϕw)
                        n_bins = convert(Int64, ceil(abs(maximum(ϕw)-minimum(ϕw))/(3.49*σ*(length(ϕw))^(-1.0/3.0))))
                        display(n_bins)
                        bins, pdf = fit_distribution(ϕw, n_bins=20)

                        # Solve the nonlinear least-squares problem and compute the early-warning signal
                        solution = fit_potential(ϕw, transformation=[0.0,0.0,16.0], n_bins=n_bins).fit
                        ews = analyse(solution)

                        # Plot location of the endpoint of the window 
                        lines!(ax1, [tw[end], tw[end]], [-32, -20], color = :black, linewidth = 1.0)
                        lines!(ax2, [tw[end], tw[end]], [-0.1, 1.1], color = :black, linewidth = 1.0)

                        # Plot the ews
                        scatter!(ax2, tw[end], ews, color = :red, markersize = 20.0, strokewidth = 1.0, strokecolor = :black)
                end

                # Export the figure
                savefig("amoc/$(N).png", fig)
        end
end

# Execute the main
main()
