"""
Statistical methods to estimate parameters of a sample.

Author: Davide Papapicco
Affil: U. of Auckland
Date: 03-09-2026
"""

# Estimate the return rate of a sample path 
function estimate_return_rate(timeseries, dt; alg = "LLS", σ = nothing)
        # Initialize the return rate
        α = 0.0

        # Fit a linear model 
        if alg == "LLS"
                α = (estimate_parameters(timeseries, dt))[2]

        # Fit a lag-1 autocorrelation model
        elseif alg == "AC1"
                ρ = sum(timeseries[1:end-1].*timeseries[2:end])/sum(abs2, timeseries)
                α = log(ρ)/dt

        # Fit an OUP
        elseif alg == "OUP"
                α = -σ^2/(2*var(timeseries))
        end

        return α 
end

# Solve a linear least-square regression on the increments
function estimate_parameters(timeseries, dt; α = 0.0)
        # Define the observation and data vectors 
        Xn = timeseries[1:end-1]
        Y  = (timeseries[2:end] .- timeseries[1:end-1])./dt

        # Assemble the model matrix
        A = hcat(ones(length(Xn)), Xn)

        # Solve the regularized (linear) least-squares problem
        θ = (A'*A + α.*I(size(A,2)))\(A'*Y)
        return θ 
end

# Estimate the tipping time of a system by linear extrapolation of the return rate
function estimate_tipping_time(timestamps, timeseries)
        # Assemble the model matrix
        X = hcat(ones(length(timestamps)), timestamps)

        # Solve the least-squares problem
        c = X\timeseries
        #trend = X*c

        # Return the solution
        #return trend 
        return c 
end
