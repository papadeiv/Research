"""
    Setup

System's parameter and simulation settings.
"""

# Settings
const σ = 1e-1                                       # Noise intensity
const dt = 1e-1                                      # Timestep size
const T = 2e2                                        # Time horizon
const T0 = 5e1                                       # Transient treshold
const Ne = 100                                       # Ensemble size 
const μ_in = -1.0                                    # Initial parameter value
const μ_f = -1e-16                                   # Final parameter value
const Nμ = 1000                                      # Number of sweep points 
const μ_range = range(μ_in, μ_f, length=Nμ)          # Parameter space
x_u(μ) = -sqrt(-μ)                                   # Unstable branch 

# Dynamical system
f(x, p, t) = -p - x^2                                # Drift
η(x) = σ                                             # Diffusion

# Solve the initial-valued problem
function generate_samples(μ)
  # Define the reflecting boundary at the unstable branch
  reflect!(integrator) = integrator.u = max(integrator.u, 2*x_u(μ) - integrator.u)
  reflected = DiscreteCallback((x, t, integrator) -> true, reflect!; save_positions=(false, true))

  # Define the differential equation
  drift(x, p, t) = f(x, p, t)
  diffusion(x, p, t) = η(x)

  # Define and solve the ensemble problem
  problem = SDEProblem(drift, diffusion, sqrt(-μ), (0.0, T), μ)
  ensemble = solve(EnsembleProblem(problem), SOSRA(), EnsembleThreads();
                   trajectories=Ne, dt=dt, adaptive=false, save_everystep=false, callback=reflected)

  # Return the ensemble of sample paths
  return ensemble.u
end
