"""
    Processing

Routines to postprocess and analyse the results of a simulation.
"""

# Compute the frozen branches and their return rate over the ramp
function compute_frozen_branches()
  # Sample the ramp up to the breaking time
  t = collect(range(0.0, τ, length=500))

  return (t = t,
          stable = xs.(t),
          unstable = -xs.(t),
          rate = α.(t))
end

# Compute the return rate along a solution of the normal form
compute_return_rate(x) = 4 .* x.^2

# Split a sample into the parts lying outside and inside the admissible horizon
function split_admissible(t)
  return (transient = t .<= tlo,
          admissible = tlo .<= t .<= Tc,
          inadmissible = t .>= Tc)
end
