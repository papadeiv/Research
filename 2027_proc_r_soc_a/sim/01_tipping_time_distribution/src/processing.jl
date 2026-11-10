"""
    Processing

Routines to compute the analytical densities of the tipping time.
"""

# Exponential density of the stationary problem
function p_stationary(t, P)
  # Compute the Kramers escape rate
  k = sqrt(-μ0)/π*exp(-4/3*(-μ0)^(3/2)/P.D)

  return k*exp(-k*t)
end

# Small noise density of the tipping time
function p_small_noise(t, P)
  # Rescale the time and reject the points outside the support
  w = P.ε^(1/3)*(τ(P) - t)
  -3.2711 < w < -1.1737 || return 0.0

  # Integrate the fourth power of Ai over the tail by the trapezoidal rule
  s = range(w, 10.0, length=3000)
  y = airyai.(s).^4
  g = (sum(y) - (y[1] + y[end])/2)*step(s)

  # Evaluate the density
  A = airyai(w)
  B = airybi(w)
  return P.ε^(5/6)/(sqrt(4π^5*P.D)*B^2*sqrt(g))*exp(-P.ε*A^2/(4π^2*P.D*B^2*g))
end

# Large noise density of the tipping time
function p_large_noise(t, P)
  # Reject the points past the breaking time
  μ(t, P) < 0 || return 0.0

  # Compute the Kramers escape rate
  V = 4/3*(-μ(t, P))^(3/2)
  k = sqrt(-μ(t, P))/π*exp(-V/P.D)

  return k*exp(-P.D/(2π*P.ε)*exp(-V/P.D))
end
