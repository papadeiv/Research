"""
    Setup

System's parameter and simulation settings.
"""

# Settings
const μ0 = -0.5                                      # Initial parameter
const ε = 1e-2                                       # Ramp speed
const σ = 0.6*sqrt(ε)                                # Noise intensity
const x0 = 2.0                                       # Initial condition
const dt = 1e-3                                      # Timestep size
const xesc = -2.5                                    # Escape threshold

# Parameters
const τ = -μ0/ε                                      # Breaking time
const τpb = τ + 2.338107*ε^(-1/3)                    # Tipping time
const tlo = 10/(2*sqrt(-μ0))                         # End of the transient
const Tc = τ - 2*ε^(-1/3)                            # Admissible time horizon
const T = τpb                                        # Termination condition

# System
μ(t) = μ0 + ε*t                                      # Linear ramp
f(x, t) = -μ(t) - x^2                                # Normal form
xs(t) = sqrt(-μ(t))                                  # Stable branch
λ(t) = -2*sqrt(-μ(t))                                # Leading eigenvalue
α(t) = λ(t)^2                                        # Return rate

# Solve the initial value problem by Euler-Maruyama
function generate_sample(s)
  # Initialize the sample at the initial condition
  t, x = [0.0], [x0]

  # Propagate until the termination condition or the escape threshold
  while t[end] < T && x[end] > xesc
    push!(x, x[end] + f(x[end], t[end])*dt + s*sqrt(dt)*randn())
    push!(t, t[end] + dt)
  end

  return t, x
end
