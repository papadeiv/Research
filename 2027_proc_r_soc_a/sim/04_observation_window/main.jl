"""
    Main script

Run this file to execute the simulation, analyse and plot the results.
"""

using StochasticDiffEq 
using Random
using CairoMakie

Random.seed!(1)

# Import local functions
include("src/setup.jl")
include("src/plotting.jl")

# Main algorithm
function main()

  # Solve the sample paths over the three window lengths
  samples = (L = generate_sample(L(tn)),
             Tw = generate_sample(Tw),
             U = generate_sample(U(tn)))

  # Analyse and plot the samples
  plot(samples)
end

# Execute main
main()
