"""
    Main script

Run this file to execute the simulation, analyse and plot the results.
"""

using CairoMakie
using StochasticDiffEq

# Import local functions
include("src/setup.jl")
include("src/plotting.jl")

# Main algorithm
function main()
  # Solve May's tipping model
  sample_path = generate_samples()

  # Plot the solution
  plot(sample_path)
end

# Execute main
main()
