"""
    Setup

System's parameter and simulation settings.
"""

# Settings
const μ0 = 1.20                                      # Initial parameter value
const μf = 2.80                                      # Final parameter value
const ε = 1e-3                                       # Ramp speed
const σ = 1e-1                                       # Noise intensity
const dt = 1e-1                                      # Timestep size
const r = 1.0                                        # Growth rate
const k = 10.0                                       # Carrying capacity
const h = 1.0                                        # Half-grazing biomass
const T = 0.5*(2*π/ε)                                # Time horizon 
const x_in = 8.0                                     # Initial condition
const x_esc = -1.0                                   # Termination condition

# Dynamical system  
μ(t) = μ0 + (μf - μ0)*sin(ε*t)                       # Shift
f(x, t) = r*x*(1-x/k)-μ(t)*((x^2)/((x^2)+(h^2)))     # Drift
η(x) = σ                                             # Diffusion

# Solve the initial-valued problem
function generate_samples()
  # Define the escape condition
  escaped = DiscreteCallback((x, t, integrator) -> x < x_esc, terminate!)

  # Define the SDE
  drift(x, p, t) = f(x, t)
  diffusion(x, p, t) = η(x) 

  # Compute a sample path 
  problem = SDEProblem(drift, diffusion, x_in, (0.0, T))
  sample = solve(EnsembleProblem(problem), SOSRA(), EnsembleThreads();
                 trajectories=1, dt=dt, adaptive=false, saveat=dt, callback=escaped)

  # Return the solution tuple (state and timestep)
  return sample.u[1]
end
