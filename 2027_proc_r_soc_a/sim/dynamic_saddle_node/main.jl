"""
    Main script

Run this file to execute the simulation, analyse and plot the results.
"""

using StochasticDiffEq, LinearAlgebra, Statistics, Random
using CairoMakie, SpecialFunctions, Printf

Random.seed!(1)

# Import local functions
include("src/setup.jl")
include("src/processing.jl")
include("src/plotting.jl")

# Main algorithm
function main()

  # Solve the ensemble problem and estimate the return rates of every sample
  printstyled("Simulating the ensemble and estimating the return rates\n"; bold = true, color = :light_blue)
  Ω = generate_samples()

  # Collect the estimates in matrices of size (windows, samples)
  estimates = reduce(hcat, [sample.α for sample in Ω])
  uncertainties = reduce(hcat, [sample.u for sample in Ω])
  tipping_times = filter(isfinite, [sample.tip for sample in Ω])

  # Extrapolate the breaking time of every sample
  breaking_times = [estimate_breaking_time(tw, estimates[:,j], uncertainties[:,j]).τ for j in 1:Ne]
  τ̂ = filter(isfinite, breaking_times)
  @printf("τ = %.1f, τpb = %.1f, median breaking time %.1f\n", τ, τpb, median(τ̂))
  @printf("forecasts before the tipping time of the pullback attractor: %.1f%%\n", 100*count(τ̂ .< τpb)/length(τ̂))

  # Plot the results
  plot(Ω[1], estimates, uncertainties, breaking_times, tipping_times)
end

# Execute main
main()
