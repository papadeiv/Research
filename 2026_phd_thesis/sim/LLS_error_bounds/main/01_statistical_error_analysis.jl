include("../../../../inc/IO.jl")
include("../../../../inc/PotentialLearning.jl")

# Define the domain of the observations
a = -1.000 
b = -a 
domain = [a, b]

# Define the deterministic function
f(x) = x^2

# Define the degree of the polynomial least-squares regression
d = 2

# Vector of different number of observations per sample 
N_observations = collect(StepRange(10 ,Int64(10), 1000))
N = length(N_observations)

# Vectors storing the errors
c_hat = Matrix{Float64}(undef, N, d+1)
estimate_8_1 = Matrix{Float64}(undef, N, d+1)
estimate_8_2 = Matrix{Float64}(undef, N, d+1)
estimate_9 = Matrix{Float64}(undef, N, d+1)
estimate_11 = Matrix{Float64}(undef, N, d+1)

# Standard deviation of the observations 
σ = 0.100

# Loop over the different number of observations
for n in 1:N
        # Generate normally distributed random samples
        ξ = σ.*randn(N_observations[n])
        
        # Generate the observation functions
        x = LinRange(domain[1], domain[end], N_observations[n])
        y = [f(xn) for xn in x] # Deterministic data
        z = y + ξ               # Noisy data

        # Get the least-squares solution of the noisy data
        c = get_coefficients(x, z, degree=d)[3]
        # Get the least-squares solution of the deterministic data
        c_bar = get_coefficients(x, y, degree=d)[3]
        # Compute the exact error between the mean fit and the noisy fit
        c_hat[n,:] = abs.(c - c_bar)
        # Compute the estimated errors between the mean fit and the noisy fit
        estimate_8_1[n,:] = estimate_error(x, z, σ, d, estimate=1)
        estimate_8_2[n,:] = estimate_error(x, z, σ, d, estimate=2)
        estimate_9[n,:] = estimate_error(x, z, σ, d, estimate=3)
        estimate_11[n,:] = estimate_error(x, z, σ, d, estimate=4)
end

# Export the data
writeout(N_observations, "../data/N_observations.csv")
writeout(estimate_8_1, "../data/estimate_8_1.csv")
writeout(estimate_8_2, "../data/estimate_8_2.csv")
writeout(estimate_9, "../data/estimate_9.csv")
writeout(estimate_11, "../data/estimate_11.csv")
writeout(c_hat , "../data/exact_error.csv")

# Execute the additional scripts for this simulation
#include("../postprocessing/00_postprocessing.jl")
include("../plotting/01_plotting.jl")
