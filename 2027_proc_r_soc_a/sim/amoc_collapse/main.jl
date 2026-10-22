"""
    Main script

Run this file to run the analysis and plot the results.
"""

using NCDatasets
using CairoMakie
using Statistics
using LinearAlgebra

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

  # Estimate the return rate
  α = estimate_return_rate(t, ϕ)

  # Extrapolate the breaking time
  τ = extrapolate_breaking_time(α)
  display(τ)

  # Plot the timeseries and the forecast
  plot_forecast(t, ϕ, α, τ)
end

# Execute main
main()
