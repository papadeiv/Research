"""
    Setup

System's parameter and simulation settings.
"""

# Dynamical system  
f(μ0, t, ε, x) = -(μ0 + ε*t) - x^2                   # Drift
η(x) = σ                                             # Diffusion

# Solve the initial-valued problem
function generate_samples(IC, rate, noise, termination)
  # Define the differential equation
  drift(x, p, t) = f(IC[2], t, rate, x)
  diffusion(x, p, t) = noise

  # Define and solve the ensemble problem
  problem = SDEProblem(drift, diffusion, IC[1], (0.0, 10000.0))
  sample = solve(EnsembleProblem(problem), SOSRA(), EnsembleThreads();
                 trajectories=1, dt=1e-1, adaptive=false, saveat=1e-1, callback=termination)

  # Return the solution tuple (state and timestep)
  return sample.u[1]
end
