"""
    Processing

Routines to postprocess and analyse the results of a simulation.
"""

# Collect the tipping time of every solution that escaped
function collect_tipping_times(Ω, P)
  return [path.t[end] for path in Ω.u if path.t[end] < P.T - P.dt]
end

# Thin the ensemble down to the solutions that are drawn
function collect_solutions(Ω)
  return [(t = path.t, x = path.u) for path in Ω.u[1:max(1, Ne ÷ Np):end]]
end

# Potential barrier and Kramers escape rate
ΔV(μ) = 4/3*(-μ)^(3/2)
k(μ, D) = sqrt(-μ)/π*exp(-ΔV(μ)/D)

# Compute the exponential law of the frozen problem
p_stationary(t, P) = k(P.μ0, P.D)*exp(-k(P.μ0, P.D)*t)

# Compute the tail integral of the fourth power of Ai by the trapezoidal rule
function airy_tail(w; wmax=10.0, n=3000)
  u = range(w, wmax, length=n)
  y = airyai.(u).^4

  return (sum(y) - (y[1] + y[end])/2)*step(u)
end

# Compute the small noise law, supported between the first two zeros of Bi
function p_small_noise(t, P)
  # Rescale the time and reject the points outside the support
  w = P.ε^(1/3)*(τ(P) - t)
  -3.2711 < w < -1.1737 || return 0.0

  # Evaluate the Airy functions and their tail integral
  A, B, g = airyai(w), airybi(w), airy_tail(w)

  return P.ε^(5/6)/(sqrt(4π^5*P.D)*B^2*sqrt(g))*exp(-P.ε*A^2/(4π^2*P.D*B^2*g))
end

# Compute the large noise law, the frozen rate modulated by the survival probability
function p_large_noise(t, P)
  # Reject the points beyond the breaking time
  μ = P.μ0 + P.ε*t
  μ < 0 || return 0.0

  return k(μ, P.D)*exp(-P.D/(2π*P.ε)*exp(-ΔV(μ)/P.D))
end

# Compute the predicted density of the tipping times of an experiment
function compute_density(t, P, col)
  col == 1 && return p_stationary.(t, Ref(P))
  col == 2 && return p_small_noise.(t, Ref(P))

  return p_large_noise.(t, Ref(P))
end

# Compute the plotting range of the tipping times, widened to include τ and τpb
function compute_range(ttip, P)
  low, high = extrema(ttip)
  pad = 0.05*(high - low)
  low, high = low - pad, high + pad

  P.ε > 0 && return (min(low, τ(P) - pad), max(high, τpb(P) + pad))

  return (0.0, high)
end

# Compute the number of histogram bins by Scott's rule
scott(x) = max(1, ceil(Int, (maximum(x) - minimum(x))/(3.49*std(x)*length(x)^(-1/3))))

# Thin a vector down to at most n points
thin(v, n) = v[round.(Int, range(1, length(v), length=min(n, length(v))))]
