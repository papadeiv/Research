"""
    Main script

Run this file to execute the simulation, analyse and plot the results.
"""

using StochasticDiffEq 
using Statistics
using Random
using CairoMakie
using ProgressMeter

Random.seed!(1)

# Import local functions
include("src/setup.jl")
include("src/processing.jl")
include("src/plotting.jl")

# Main algorithm
function main()

  # Initialize empty data structure
  estimates = (OUP = fill(NaN, Nμ, Ne),
               AC1 = fill(NaN, Nμ, Ne),
               OLS = fill(NaN, Nμ, Ne))

  # Loop over the parameter set
  printstyled("Performing the parameter sweep\n"; bold = true, color = :light_blue)
  @showprogress for (idx_μ, μ) in enumerate(μ_set)
    # Solve the ensemble problem
    Ω = generate_samples(μ)

    # Loop over the samples
    for (idx_x, x) in enumerate(Ω)
      # Estimate the return rate of the sample
      α = estimate_return_rate(x)
      estimates.OUP[idx_μ,idx_x] = α.OUP
      estimates.AC1[idx_μ,idx_x] = α.AC1
      estimates.OLS[idx_μ,idx_x] = α.OLS
    end
  end

  # Analyse and plot the estimates
  plot(estimates)
end

# Execute main
main()
