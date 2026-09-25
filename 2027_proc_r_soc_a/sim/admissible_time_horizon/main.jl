"""
    Main script

Run this file to execute the simulation, analyse and plot the results.
"""

using Random
using CairoMakie

Random.seed!(1)

# Import local functions
include("src/setup.jl")
include("src/processing.jl")
include("src/plotting.jl")

# Main algorithm
function main()

  # Integrate the pullback attractor and a single noisy solution
  t_pb, x_pb = generate_sample(0.0)
  t_X, X = generate_sample(σ)

  # Collect the solutions and their return rates
  solutions = (pullback = (t = t_pb, x = x_pb, rate = compute_return_rate(x_pb)),
               noisy = (t = t_X, x = X))

  # Compute the frozen branches and their return rate
  branches = compute_frozen_branches()

  # Analyse and plot the solutions
  plot(branches, solutions)
end

# Execute main
main()
