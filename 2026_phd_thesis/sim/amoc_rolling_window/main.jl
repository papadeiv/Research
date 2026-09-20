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

        # Loop over the latitude subsets
        for M in 1:length(subsets)
                # Extract the latitudes subset
                latitudes = subsets[M] 

                # Generate figure and axis
                fig = Figure(size = (1200,800))
                ax1 = Axis(fig[1,1], ylabel = L"\textbf{flow (250-500m)}", limits = (time[1], time[end], -32,-20), title = "Latitude range = [$(latitude[(subsets[M])[1]]), $(latitude[(subsets[M])[end]])]")
                ax2 = Axis(fig[2,1], xlabel = L"\textbf{years}", ylabel = L"\textbf{early-warning signal}", limits = (time[1], time[end], -0.1, 1.1))

                # Loop over the different latitudes
                printstyled("Batch $M, latitudes = [$(latitude[latitudes[1]]), $(latitude[latitudes[end]])]: computing the least-squares solutions across the ensemble\n"; bold=true, underline=true, color=:light_blue)
                @showprogress for N in latitudes
                        # Extract subseries until tipping
                        t = copy(time[1:idx])
                        T = copy(temperature[N,1:idx]) 
                        S = copy(salinity[N,1:idx]) 

                        # Define the observable (flow)
                        ϕ = T - S
                        detrended_solution = detrend(ϕ, alg = "emd", n_modes=1)
                        trend = detrended_solution.trend
                        residual = detrended_solution.residuals

                        # Plot the observable timeseries
                        if N ≥ N_emph 
                                lines!(ax1, t, ϕ, color = N, colormap = (:managua, 0.75), colorrange = (1,(subsets[end])[end]), linewidth = 2.0)
                        else
                                lines!(ax1, t, ϕ, color = N, colormap = (:managua, 0.5), colorrange = (1,(subsets[end])[end]), linewidth = 1.0)
                        end
                        lines!(ax1, time[(idx+1):end], temperature[N,(idx+1):end] - salinity[N,(idx+1):end], color = (:black,0.35), linewidth = 1.0)

                        # Convert the sliding window problem into an ensemble problem
                        ensemble = preprocess_solution(t, residual, window_size)
                        Ne = length(ensemble.trajectories)
                        Nt = length(ensemble.trajectories[1])

                        # Loop over the ensemble
                        ews = Vector{Float64}(undef, Ne)
                        for n in 1:Ne
                                # Extract windowed subseries
                                tw = copy(ensemble.timesteps[n])
                                ϕw = copy(ensemble.trajectories[n]) 

                                # Assemble histogram of the windowed subseries
                                σ = std(ϕw)
                                n_bins = convert(Int64, ceil(abs(maximum(ϕw)-minimum(ϕw))/(3.49*σ*(length(ϕw))^(-1.0/3.0))))
                                bins, pdf = fit_distribution(ϕw, n_bins=n_bins)

                                # Solve the nonlinear least-squares problem and compute the early-warning signal
                                solution = fit_potential(ϕw, transformation=[0.0,0.0,16.0], n_bins=n_bins).fit
                                ews[n] = analyse(solution)
                        end

                        # Plot the early-warning signal
                        if N ≥ N_emph 
                                lines!(ax2, t[(Nt-1):(Nt+Ne-2)], ews, color = N, colormap = :managua, colorrange = (1,(subsets[end])[end]), linewidth = 3.0)
                        else
                                lines!(ax2, t[(Nt-1):(Nt+Ne-2)], ews, color = N, colormap = (:managua, 0.5), colorrange = (1,(subsets[end])[end]), linewidth = 1.0)
                        end
                end

                # Plot the location of the tipping point
                lines!(ax1, [time[idx], time[idx]], [-32, -20], color = :red, linestyle = :dash, linewidth = 3.0)
                lines!(ax2, [time[idx], time[idx]], [-0.1, 1.1], color = :red, linestyle = :dash, linewidth = 3.0)

                # Export the figure
                savefig("amoc/$(window_size*100)%/$M.png", fig)
        end
end

# Execute the main
main()
