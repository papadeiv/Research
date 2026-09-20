"""
Analyse the solutions of stochastic processes whose deterministic term is nonautonomus (i.e. time-dependent).

Author: Davide Papapicco
Affil: U. of Auckland
Date: 03-09-2025
"""

function detrend(timeseries; alg = "exact", timestamps = Float64[], qse = Float64[], n_modes=0::Integer)
        # Initialise arrays for the trend and the residuals
        trend = Float64[]
        residuals = Float64[]

        # Detrend the timeseries
        if alg == "mean"
                # Compute the mean of the timeseries to define the trend
                trend = mean(timeseries).*ones(length(timeseries))
                # Remove the trend to find the residuals
                residuals = timeseries - trend 

        elseif alg == "linear"
                # Assemble the model matrix
                X = hcat(ones(length(timestamps)), timestamps)
                # Solve the least-squares problem
                c = X\timeseries
                # Compute the linear trend and the residuals
                trend = X*c
                residuals = timeseries - trend

        elseif alg == "exact"
                trend = qse 
                residuals = timeseries - trend

        else alg == "emd"
                emd = PyEMD.EMD()
                imfs = Array(emd(timeseries))
                trend = sum(imfs[(end-n_modes):end,:], dims=1)[:]
                residuals = timeseries .- trend
        end

        return (
                trend = trend, 
                residuals = residuals
               )
end

# Converts the sliding window problem into an ensemble one 
function build_window(timestamps, timeseries, width)
        # Compute the number of strides
        strides = (length(timestamps) - width) + 1

        timesteps = [@view timestamps[n:(n+width-1)] for n in 1:strides]
        ensemble = [@view timeseries[n:(n+width-1)] for n in 1:strides]

        # Export the ensemble problem 
        return (
                timesteps = timesteps,
                trajectories = ensemble 
               ) 
end

# Truncate a timeseries once it crosses a threshold
function find_tipping(timeseries, threshold)
        tipping_index = findfirst(u -> u < threshold, timeseries)
        if isnothing(tipping_index)
                return length(timeseries)
        else
                return tipping_index
        end 
end
