"""
    Setup

System's parameter and simulation settings.
"""

# Settings
const σ = 1e-1                                       # Noise intensity
const dt = 1e-1                                      # Timestep size
const T = 5e3                                        # Time horizon
const n0 = 100                                       # Initial window size
const n_path = 1000                                  # Displayed timesteps
const μ0 = -1.0                                      # Parameter value
const x0 = 1.25                                      # Initial condition
x_u(μ) = -sqrt(-μ)                                   # Unstable branch

# Dynamical systems
f(x, p, t) = -p - x^2                                # Saddle-node drift
g(x, p, t) = 2*sqrt(-p)*(sqrt(-p) - x)               # OUP drift
η(x) = σ                                             # Diffusion

# Solve the initial-valued problem
function generate_sample(drift)
  # Define the differential equation
  diffusion(x, p, t) = η(x)

  # Define and solve the problem
  problem = SDEProblem(drift, diffusion, x0, (0.0, T), μ0)
  path = solve(problem, SOSRA(); dt=dt, adaptive=false, saveat=dt)

  # Return the sample path
  return path
end
