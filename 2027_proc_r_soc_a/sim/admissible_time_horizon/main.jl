"""
    Main script

Run this file to execute the simulation, analyse and plot the results.
"""

using Random
using CairoMakie

Random.seed!(1)

# Import local functions
include("src/setup.jl")
include("src/plotting.jl")

# Main algorithm
function main()

  # Solve the dynamic saddle-node normal form
  t_pb, x_pb = generate_sample(0.0)     # Pullback attractor
  t, x = generate_sample(σ)             # Sample path

  # Compute the return rate of the quasi-steady equilibrium 
  t_eq = collect(range(0.0, τ, length=500))
  x_eq = xs.(t_eq)
  α_eq = α.(t_eq)

  # Analyse and plot the solutions
  plot(t_eq, x_eq, α_eq, t_pb, x_pb, t, x)
end

# Execute main
main()
