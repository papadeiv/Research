"""
    Postprocessing script

Collection of quantities and functions used to postprocess and analyse the results of a simulation.
"""

algs = ("LLS", "AC1", "OUP")
cols = (:red, :green, :blue)

# Linear detrending of a timeseries 
function detrend(t, x)
    A = hcat(ones(length(t)), t .- mean(t))
    return x .- A*(A\x)
end

# Estimation of the return rate of a linearised system
function estimate_return_rate(x, dt, alg)
        # Compute timeseries length and initialise the estimate
        Nt = length(x)
        α = 0.0

        # Ordinary least-squares regression on the first-order increments
        if alg == "LLS"
                A = hcat(ones(Nt-1), x[1:end-1])
                Y = (x[2:end] .- x[1:end-1])./dt
                α = ((A'*A)\(A'*Y))[2]

        # Lag-1 autocorrelation regression
        elseif alg == "AC1"
                ρ = sum(x[1:end-1].*x[2:end])/sum(abs2, x)
                α = log(ρ)/dt

        # Inversion of the OUP stationary variance
        else
                α = -σ^2/(2*var(x))
        end

        # Compute the bias in the estimate and remove it
        bias = 2*abs(α)/(Nt*dt)
        uncertainty = 4*α^2*bias
        ews = α^2 - bias 

        return ews, uncertainty 
end

# Estimate the return rates across rolling strides of the inference window
function analyse(solutions)
        # Extract indices of the inference window
        idx_T1, idx_T2 = searchsortedfirst(solutions.t, T1), searchsortedfirst(solutions.t, T2)

        # Extract the timesteps of the left edges of the strides rolling window
        strides = idx_T1:Δs:(idx_T2 - ΔT + 1)
        ts = [solutions.t[idx + ΔT÷2] for idx in strides]

        # Loop over the strides
        ews = fill(NaN, length(strides), 3)               # λ̂²
        uncertainty = fill(NaN, length(strides), 3)       # Var(λ̂²)
        for (idx_str, stride) in enumerate(strides)
                # Identify the bounding indices of the current window
                indices = stride:(stride + ΔT - 1)

                # Compute quasi-stationary residuals 
                z = detrend(solutions.t[indices], solutions.x[indices])

                # Loop over the estimating algorithms
                for (idx_est, estimator) in enumerate(estimators)
                        α, δα = estimate_return_rate(z, dt, estimator)
                        ews[idx_str, idx_est] = α 
                        uncertainty[idx_str, idx_est] = δα 
                end
        end

        return ts, ews, uncertainty 
end

# Extrapolate the estimated return rates via weighted least-squares to estimate the tipping time
function estimate_tipping_time(t, α, u)
        # Assemble the model matrix
        A = hcat(ones(length(t)), t)

        # Assemble the weight matrix
        W = Diagonal(1./u)

        # Solve the weighted normal equations
        c = (A'*W*A) \ (A'*W*α)

        # Compute the zero-crossing time of the return rate
        T_tip = -c[1]/c[2]
        return T_tip 
end
