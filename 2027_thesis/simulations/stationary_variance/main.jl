"""
    Main script

Run this file to execute the simulation, analyse and plot the results.
"""

using CairoMakie
using StochasticDiffEq
using ProgressMeter 
using QuadGK
using Statistics

# Import local functions
include("src/setup.jl")
include("src/processing.jl")
include("src/plotting.jl")

# Main algorithm
function main()
  # Initialize data arrays 
  ground_truth = zeros(Nμ)
  sample_mean = zeros(Nμ)
  lower_envelope = zeros(Nμ)
  upper_envelope = zeros(Nμ)

  # Sweep the parameter space
  @showprogress for (i, μ) in enumerate(μ_range)
    # Compute the variance of the stationary distribution
    ground_truth[i] = stationary_variance(μ)

    # Compute the sample variance of each trajectory
    sample_variance = [var(path.u[path.t .≥ T0]) for path in generate_samples(μ) if path.t[end] ≥ T]

    # Collect the ensemble statistics
    sample_mean[i] = mean(sample_variance)
    lower_envelope[i] = quantile(sample_variance, 0.05)
    upper_envelope[i] = quantile(sample_variance, 0.95)
  end

  # Plot the results
  plot(ground_truth, sample_mean, lower_envelope, upper_envelope)
  save("stationary_variance.pdf", fig)
end

# Execute main
main()
