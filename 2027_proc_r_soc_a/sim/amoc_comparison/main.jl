"""
    Main script

Run this file to execute the simulation, analyse and plot the results.
"""

using StochasticDiffEq, LinearAlgebra, Statistics, Random
using CairoMakie, Printf

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
  estimates = (OUP = reduce(hcat, [sample.estimates.OUP for sample in Ω]),
               AC1 = reduce(hcat, [sample.estimates.AC1 for sample in Ω]),
               OLS = reduce(hcat, [sample.estimates.OLS for sample in Ω]))
  uncertainties = (OUP = reduce(hcat, [sample.uncertainties.OUP for sample in Ω]),
                   AC1 = reduce(hcat, [sample.uncertainties.AC1 for sample in Ω]),
                   OLS = reduce(hcat, [sample.uncertainties.OLS for sample in Ω]))

  # Extrapolate the breaking time of every sample
  breaking_times = (OUP = [estimate_breaking_time(tw, estimates.OUP[:,j], uncertainties.OUP[:,j]).τ for j in 1:Ne],
                    AC1 = [estimate_breaking_time(tw, estimates.AC1[:,j], uncertainties.AC1[:,j]).τ for j in 1:Ne],
                    OLS = [estimate_breaking_time(tw, estimates.OLS[:,j], uncertainties.OLS[:,j]).τ for j in 1:Ne])

  # Analyse and plot the results
  print_breaking_times(breaking_times)
  plot(collect_solutions(Ω), estimates, uncertainties, breaking_times)
end

# Execute main
main()
