include("../../../../inc/IO.jl")
include("../../../../inc/PotentialLearning.jl")

# Define the domain of the observations
a = -1.000 
b = -a 
domain = [a, b]

# Define the deterministic function
f(x) = x^2

# define the degree of the polynomial least-squares regression
d = 2

# Vector of different number of observations per sample 
N_observations = collect(StepRange(10 ,Int64(10), 1000))
N = length(N_observations)

# Standard deviation of the observations 
σ = 0.100

# Number of samples
N_samples = convert(Int64,5e3)

# Matrix storing the expected exact and estimated errors
estimate = Matrix{Float64}(undef, N, d+1)
E_c_hat = Matrix{Float64}(undef, N, d+1)

# Loop over the different number of observations 
using ProgressMeter
@showprogress for n in 1:N
        # Generate the deterministic data 
        x = LinRange(domain[1], domain[end], N_observations[n])
        y = [f(xn) for xn in x]

        # Get the least-squares solution of the deterministic data
        c_bar = get_coefficients(x, y, degree=d)[3]

        # Matrix storing the exact errors
        c_hat = Matrix{Float64}(undef, N_samples, d+1)

        # Loop over the different samples
        for s in 1:N_samples
                # Generate normally distributed random samples
                ξ = σ.*randn(N_observations[n])
        
                # Generate the noisy data
                z = y + ξ

                # Get the least-squares solution of the noisy data
                c = get_coefficients(x, z, degree=d)[3]
                # Compute the exact error
                c_hat[s,:] = abs.(c - c_bar)
        end

        # Compute the expected value of the exact error
        E_c_hat[n,:] = [mean(c_hat[:,m]) for m in 1:(d+1)]

        # Compute the expected error estimate (20)
        estimate[n,:] = estimate_error(x, y, σ, d, estimate=5)
end

# Export the data
writeout(N_observations, "../data/N_observations.csv")
writeout(E_c_hat, "../data/expected_exact_error.csv")
writeout(estimate, "../data/estimate.csv")

# Execute the additional scripts for this simulation
include("../plotting/03_plotting.jl")
