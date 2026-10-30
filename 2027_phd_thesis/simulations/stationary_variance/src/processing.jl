"""
    Processing

Routines to postprocess and analyse the results of a simulation.
"""

# Stationary distribution
ρ(x, μ) = exp(2*quadgk(y -> f(y, μ, 0.0), sqrt(-μ), x)[1]/σ^2)

# Compute the variance of the stationary distribution by quadrature
function stationary_variance(μ)
  # Normalise the distribution over I
  N = 1/quadgk(x -> ρ(x, μ), x_u(μ), sqrt(-μ), Inf)[1]
  p(x) = N*ρ(x, μ)

  # Compute the expectation of the distribution
  E = quadgk(x -> x*p(x), x_u(μ), sqrt(-μ), Inf)[1]

  # Return the central variance
  return quadgk(x -> ((x - E)^2)*p(x), x_u(μ), sqrt(-μ), Inf)[1]
end
