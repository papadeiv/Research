"""
    Main script

Run this file to execute the simulation, analyse and plot the results.
"""

using StochasticDiffEq, Statistics, Random
using CairoMakie, SpecialFunctions, Printf, ProgressMeter

Random.seed!(1)

# Import local functions
include("src/setup.jl")
include("src/processing.jl")
include("src/plotting.jl")

# Main algorithm
function main()

  # Initialize empty data structure
  results = []

  # Loop over the three experiments
  printstyled("Simulating the experiments\n"; bold = true, color = :light_blue)
  @showprogress for P in runs
    # Solve the ensemble problem
    Ω = generate_samples(P)

    # Collect the tipping times and the solutions to be drawn
    push!(results, (paths = collect_solutions(Ω),
                    pullback = generate_pullback(P),
                    ttip = collect_tipping_times(Ω, P)))
  end

  # Analyse and plot the results
  plot(results)
end

# Execute main
main()
