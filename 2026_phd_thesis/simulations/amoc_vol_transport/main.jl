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
        # Data structures
        solutions = Vector{Vector{Float64}}()           # Solutions of the inference method 
        time_idx = Vector{Float64}()                    # Index of timesteps 
        ews = Vector{Float64}()                         # Statistical distribution over the ensemble

        # Import the AMOC timeseries
        data = NCDataset("../../res/data/amoc/AMOC_transport_depth_0-1000m_monthly.nc")
        time = data["time"][:]
        observable = data["Transport"][:]

        # Convert the analysis of the non-autonomous drift into an ensemble problem
        ensemble = preprocess_solution(time, observable, window_size)
        tipping = ensemble.tipping_point
        Ne = length(ensemble.trajectories)
        Nt = length(ensemble.trajectories[1])

        # Loop over the ensemble's sample paths
        residuals = Matrix{Float64}(undef, Nt, Ne)
        trends = Matrix{Float64}(undef, Nt, Ne)
        print_lock = ReentrantLock()
        counter = Base.Threads.Atomic{Int}(0)
        #=Threads.@threads=#for n in 1:Ne
                # Display the counter status
                c = atomic_add!(counter, 1)
                if c % 1 == 0
                        lock(print_lock) do
                                printstyled("\rComputing the least-squares solutions across the ensemble using $(Threads.nthreads()) threads: $c/$Ne"; bold=true, underline=true, color=:light_blue)
                                flush(stdout)
                        end
                end

                # Extract windowed subseries
                tw = copy(ensemble.timesteps[n])
                uw = copy(ensemble.trajectories[n]) 

                # Detrend it 
                detrended_solution = detrend(uw, alg = "emd", n_modes=2)
                trend = detrended_solution.trend
                residual = detrended_solution.residuals
                trends[:,n] = trend 
                residuals[:,n] = residual

                # Solve the nonlinear least-squares problem to fit a cubic potential
                solution = fit_potential(residual, transformation=[0.0,1.0,8.0], optimiser=β, attempts=Na, verbose=true)
                push!(solutions, solution.fit)

                # Perform postprocessing analysis on the solutions
                push!(ews, analyse(solution.fit))
                push!(time_idx, tw[end])
        end

        # Loop over the solutions
        @showprogress for n in 1:Ne
                if n == 1
                        # Sort the ensemble solutions by temporal order
                        sort_idx = sortperm(time_idx)
                        time_idx = time_idx[sort_idx]
                        solutions = solutions[sort_idx]
                        ews = ews[sort_idx]
                        writeout(ews, "amoc/ews.csv")
                end

                # Create empty layouts for the figures
                include("./scripts/figs.jl")

                # Plot the full timeseries
                lines!(ax1, time, observable, color = (:black,0.5), linewidth = 0.5)
                labels = ["0", "1", "2"]
                ax1.xtickformat = values -> ["$(label)" for label in labels]
                ax3.xtickformat = values -> ["$(label)" for label in labels]

                # Plot the sliding window
                poly!(ax1, Point2f[(time[n], -7), 
                                   (time_idx[n], -7),
                                   (time_idx[n], 25),
                                   (time[n], 25),
                                  ], color = (CtpGray, 0.1), strokecolor = :black, strokewidth = 1.0)

                # Plot the subseries
                lines!(ax1, time[n:(n+Nt)], observable[n:(Nt+n)], color = (:black,1.0), linewidth = 0.5)

                # Plot the location of the tipping
                lines!(ax1, [time[tipping], time[tipping]], [-7, 25], color = (:red, 1.0), linewidth = 3.0, linestyle = :dash)
                lines!(ax3, [time[tipping], time[tipping]], [-0.5, 1.5], color = (:red, 1.0), linewidth = 3.0, linestyle = :dash)

                # Plot the ews
                lines!(ax3, time_idx[1:n], ews[1:n], color = CtpBlue, linewidth = 3.0)
                scatter!(ax3, time_idx[n], ews[n], color = CtpBlue, strokewidth = 3.0, markersize = 25)

                # Plot the potential reconstruction
                domain = collect(LinRange(-15, 15, 3000))
                lines!(ax2, domain, [V(x, solutions[n]) for x in domain], color = (:black, 1.00), linewidth = 5.0)

                # Plot the windowed subseries, its trend and its detrended residuals
                lines!(ax4, time[n:(n+Nt)], observable[n:(Nt+n)], color = :black, linewidth = 1.0) 
                lines!(ax4, time[n:(n+Nt-1)], trends[:,n], color = :red, linewidth = 3.0) 
                lines!(ax6, time[n:(n+Nt-1)], residuals[:,n], color = :black, linewidth = 1.0)
                
                # Compute and plot the empirical distribution of the above
                Nb = convert(Int64, ceil(abs(maximum(observable[n:(Nt+n)])-minimum(observable[n:(Nt+n)]))/(3.49*std(observable[n:(Nt+n)])*(length(observable[n:(Nt+n)]))^(-1.0/3.0))))
                bins, pdf = fit_distribution(observable[n:(Nt+n)], n_bins=Nb+1)
                barplot!(ax5, bins, pdf, color = pdf, colormap = [CtpWhite,CtpGray], strokecolor = :black, strokewidth = 2)
                Nb = convert(Int64, ceil(abs(maximum(residuals[:,n])-minimum(residuals[:,n]))/(3.49*std(residuals[:,n])*(length(residuals[:,n]))^(-1.0/3.0))))
                bins, pdf = fit_distribution(residuals[:,n], n_bins=Nb+1)
                barplot!(ax7, bins, pdf, color = pdf, colormap = [CtpWhite,CtpGray], strokecolor = :black, strokewidth = 2)

                # Export the figures
                savefig("amoc/ews/$n.png", fig1)
                savefig("amoc/residuals/$n.png", fig4)
        end
end

# Execute the main
main()
