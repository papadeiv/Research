"""
    Main script

Run this file to execute the simulation, analyse and plot the results.
"""

using StochasticDiffEq
using SpecialFunctions
using Statistics
using Random
using Printf
using CairoMakie

Random.seed!(1)

# Import local functions
include("src/setup.jl")
include("src/processing.jl")
include("src/plotting.jl")

# Main algorithm
function main()
  # Solve the ensemble problems
  Ω_stationary = generate_samples(stationary)
  Ω_small = generate_samples(small_noise)
  Ω_large = generate_samples(large_noise)

  # Solve the pullback attractors
  pullback_small = generate_pullback(small_noise)
  pullback_large = generate_pullback(large_noise)

  # Plot the samples and the distributions of the tipping times
  plot(Ω_stationary, Ω_small, Ω_large, pullback_small, pullback_large)
end

# Execute main
main()
