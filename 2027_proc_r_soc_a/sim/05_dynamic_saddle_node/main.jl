"""
    Main script

Run this file to execute the simulation, analyse and plot the results.
"""

using StochasticDiffEq
using SpecialFunctions
using LinearAlgebra
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
  # Solve the ensemble problem
  Ω = generate_samples()

  # Compute the ensemble tipping time 
  tipping_times = filter(isfinite, [ω.tip for ω in Ω])

  # Estimate the ensemble breaking time
  breaking_times = filter(isfinite, [extrapolate_breaking_time(ω.rate).τ for ω in Ω])

  @printf("Median breaking time: %.1f\n", median(breaking_times))
  @printf("Forecasts before the tipping time: %.1f%%\n", 100*count(breaking_times .< τpb)/length(breaking_times))

  # Extract the ensemble median
  α = median(reduce(hcat, [ω.rate.α for ω in Ω[1:Nm]]), dims=2)[:]
  u = median(reduce(hcat, [ω.rate.u for ω in Ω[1:Nm]]), dims=2)[:]./Nm
  rate = (t = Ω[1].rate.t, α = α, u = u)

  # Extrapolate the breaking time of the median return rate
  forecast = extrapolate_breaking_time(rate)

  # Plot the results
  plot(Ω[1], rate, forecast, tipping_times, breaking_times)
end

# Execute main
main()
