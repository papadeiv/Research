"""
    Main script

Run this file to execute the simulation, analyse and plot the results.
"""

using CairoMakie
using StochasticDiffEq
using QuadGK
using Statistics

# Import local functions
include("src/setup.jl")
include("src/processing.jl")
include("src/plotting.jl")

# Main algorithm
function main()
  # Generate the sample paths of the saddle-node and of the approximating OUP
  paths = [generate_sample(f), generate_sample(g)]

  # Expand the window from the first n0 timesteps to the whole timeseries
  window = n0+1:length(paths[1].t)

  # Compute the sample mean and variance over each window
  sample_means = [[mean(@view path.u[1:n]) for n in window] for path in paths]
  sample_variances = [[var(@view path.u[1:n]) for n in window] for path in paths]

  # Compute the asymptotic mean and variance of each process
  saddle_node_mean, saddle_node_variance = stationary_moments(μ0)
  true_means = [saddle_node_mean, sqrt(-μ0)]
  true_variances = [saddle_node_variance, (σ^2/2)/(2*sqrt(-μ0))]

  # Plot the results
  plot(paths, paths[1].t[window], sample_means, sample_variances, true_means, true_variances)
  save("oup_approximation.pdf", fig)
end

# Execute main
main()
