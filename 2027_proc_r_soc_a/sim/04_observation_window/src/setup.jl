"""
    Setup

System's parameter and simulation settings.
"""

# Parameters
const μ0 = -0.12                                     # Initial parameter value
const ε = 1e-3                                       # Ramp speed
const σ = 0.014                                      # Noise intensity 
const Δ = 0.1                                        # Sampling frequency 
const n = 10                                         # Timestep size 

# Observation protocol
const τ = -μ0/ε                                      # Breaking time
const Tc = τ - ε^(-1/3)                              # Admissible time horizon
const Tx = τ - (16*ε)^(-1/3)                         # Intersection of the bounds
const Tw = 2*(τ - Tx)                                # Fixed window length
const tn = 2*τ/3                                     # Observation time

# System 
μ(t) = μ0 + ε*t                                      # Linear ramp
f(x, t) = -μ(t) - x^2                                # Normal form
xs(t) = sqrt(-μ(t))                                  # Stable branch
λ(t) = -2*sqrt(-μ(t))                                # Leading eigenvalue
L(t) = 1/abs(λ(t))                                   # Lower bound of the window length
U(t) = λ(t)^2/(2*ε)                                  # Upper bound of the window length

# Solve the initial value problem over a window of length T ending at the observation time
function generate_sample(T)
  # Define the sampling times
  Nw = round(Int, T/Δ)
  timesteps = tn .- (Nw:-1:0).*Δ

  # Define the slow-fast SDE
  drift(x, p, t) = f(x, t)
  diffusion(x, p, t) = σ

  # Propagate the initial condition 
  problem = SDEProblem(drift, diffusion, xs(timesteps[1]), (timesteps[1], timesteps[end]))
  path = solve(problem, SOSRA(); dt=Δ/n, adaptive=false, saveat=timesteps)

  return path.u
end
