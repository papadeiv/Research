"""
    Processing

Routines to postprocess and analyse the results of a simulation.
"""

# Stationary distribution
ρ(x, μ) = exp(2*quadgk(y -> f(y, μ, 0.0), sqrt(-μ), x)[1]/σ^2)

# Compute the mean and variance of the stationary distribution by quadrature
function stationary_moments(μ)
  # Normalise the distribution over I
  N = 1/quadgk(x -> ρ(x, μ), x_u(μ), sqrt(-μ), Inf)[1]
  p(x) = N*ρ(x, μ)

  # Compute the expectation of the distribution
  E = quadgk(x -> x*p(x), x_u(μ), sqrt(-μ), Inf)[1]

  # Return the expectation and the central variance
  return E, quadgk(x -> ((x - E)^2)*p(x), x_u(μ), sqrt(-μ), Inf)[1]
end

# Stationary distribution of the saddle-node normalised over I = [x_u, ∞)
function saddle_node_distribution(μ)
  N = 1/quadgk(x -> ρ(x, μ), x_u(μ), sqrt(-μ), Inf)[1]
  return x -> N*ρ(x, μ)
end

# Stationary distribution of the OUP, with ω = 2√(-μ) and D = σ²/2
oup_distribution(x, μ) = sqrt(2*sqrt(-μ)/(π*σ^2))*exp(-2*sqrt(-μ)*(x - sqrt(-μ))^2/σ^2)
