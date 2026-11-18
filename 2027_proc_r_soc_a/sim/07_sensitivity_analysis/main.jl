"""
    Main script

Run this file to run the analysis and plot the results.
"""

using LinearAlgebra
using NCDatasets
using CairoMakie
using Statistics
using Printf

# Import local functions
include("src/setup.jl")
include("src/processing.jl")
include("src/plotting.jl")

# Main algorithm
function main()
  # Import the AMOC timeseries
  data = NCDataset("data/Atlantic_250-500m_year_0001-2200.nc")
  t = Float64.(data["time"][:])
  L = data["lat"][:]
  T = data["TEMP"]
  S = data["SALT"]

  # Compute the observable at the 34S latitude
  idx = findfirst(==(-34.0), L)
  ϕ = Float64.(T[idx,:] - S[idx,:])

  # Forecast the tipping time for each window length
  τ = Float64[]
  for W in Tw
    # Estimate the return rate
    α = estimate_return_rate(t, ϕ, round(Int, W/Δ))

    # Extrapolate the breaking time
    push!(τ, extrapolate_breaking_time(α).τ)
  end

  @printf("Mean forecasted breaking time: %.1f\n", mean(τ))
  @printf("Estimates within one standard deviation from the mean: %.1f%%\n", 100*count(t -> mean(τ) - std(τ) < t < mean(τ) + std(τ), τ)/length(τ))

  # Plot the forecasted tipping time against the window length
  plot_sensitivity(Tw, τ)
end

# Execute main
main()
