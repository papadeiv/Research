"""
    Setup

System's parameter and simulation settings.
"""

# Settings
const μ0 = -0.5                                      # Initial parameter value
const ε = 1e-2                                       # Ramp speed
const σ = 0.6*sqrt(ε)                                # Noise intensity
const x0 = 2.0                                       # Initial condition
const dt = 1e-3                                      # Timestep size
const xesc = -2.5                                    # Escape threshold

# Parameters
const τ = -μ0/ε                                      # Breaking time
const τpb = τ + 2.338107*ε^(-1/3)                    # Tipping time
const Tb = 10/(2*sqrt(-μ0))                          # Initial time (end of the transient)
const Tc = τ - 2*ε^(-1/3)                            # Admissible time horizon
const T = τpb                                        # Termination condition

# System
μ(t) = μ0 + ε*t                                      # Linear ramp
f(x, t) = -μ(t) - x^2                                # Normal form
xs(t) = sqrt(-μ(t))                                  # Stable branch
λ(t) = -2*sqrt(-μ(t))                                # Leading eigenvalue
α(t) = λ(t)^2                                        # Return rate

# Compute the return rate along a solution of the normal form
α_pb(x) = 4 .* x.^2

# Forward solution of the initial value problem 
function generate_sample(η)
  # Initialize the sample path
  t, x = [0.0], [x0]

  # Euler-Maruyama time stepping 
  while t[end] < T && x[end] > xesc
    push!(x, x[end] + f(x[end], t[end])*dt + η*sqrt(dt)*randn())
    push!(t, t[end] + dt)
  end

  return t, x
end
