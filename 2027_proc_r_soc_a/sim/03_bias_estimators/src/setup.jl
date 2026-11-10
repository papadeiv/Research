"""
    Setup

System's parameter and simulation settings.
"""

# Settings
const Nμ = 200                                       # Parameter sweeps
const Ne = 200                                       # Ensemble size 
const Nt = 400                                       # Sample size 

# Parameters
const μ_set = collect(range(-1.0, -1e-3, length=Nμ)) # Parameter set 
const σ = 0.02                                       # Noise intensity 
const Δ = 0.1                                        # Sampling frequency 
const n = 10                                         # Timestep size 
const T = (Nt - 1)*Δ                                 # Termination condition 

# System 
f(x, μ) = -μ - x^2                                   # Normal form
xs(μ) = sqrt(-μ)                                     # Stable branch
λ(μ) = -2*sqrt(-μ)                                   # Leading eigenvalue
α(μ) = λ(μ)^2                                        # Return rate

# Solve the ensemble problem
function generate_samples(μ)
  # Define the escape condition
  escaped = DiscreteCallback((x, t, integrator) -> x < -xs(μ), terminate!)

  # Define the slow-fast SDE
  drift(x, p, t) = f(x, μ)
  diffusion(x, p, t) = σ

  # Propagate the initial condition 
  problem = SDEProblem(drift, diffusion, xs(μ), (0.0, T))
  paths = solve(EnsembleProblem(problem), SOSRA(), EnsembleThreads();
                trajectories=Ne, dt=Δ/n, adaptive=false, saveat=Δ, callback=escaped)

  return [path.u[1:Nt] for path in paths.u if length(path.u) >= Nt]
end
