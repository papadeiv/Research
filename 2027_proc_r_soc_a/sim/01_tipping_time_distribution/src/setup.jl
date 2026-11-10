"""
    Setup

System's parameter and simulation settings.
"""

# Settings
const Ne = 10_000                                    # Ensemble size
const Np = 100                                       # Samples drawn per panel
const x0 = 2.0                                       # Initial condition of the pullback attractor
const xesc = -2.5                                    # Escape threshold of the pullback attractor
const xtip = -20.0                                   # Tipping threshold
const μ0 = -1.0                                      # Initial parameter

# Experiments
const stationary = (ε = 0.0, D = 0.25, dt = 2e-2, T = 4000.0)
const small_noise = (ε = 1e-2, D = 1e-3, dt = 5e-3, T = 125.0)
const large_noise = (ε = 1e-3, D = 0.0625, dt = 1e-2, T = 1100.0)

# System
μ(t, P) = μ0 + P.ε*t                                 # Linear ramp
f(x, t, P) = -μ(t, P) - x^2                          # Normal form
xs(t, P) = sqrt(-μ(t, P))                            # Stable branch
τ(P) = -μ0/P.ε                                       # Breaking time
τpb(P) = τ(P) + 2.338107*P.ε^(-1/3)                  # Tipping time

# Solve the ensemble problem
function generate_samples(P)
  # Define the tipping condition
  tipped = DiscreteCallback((x, t, integrator) -> x < xtip, terminate!)

  # Define the ramped SDE
  drift(x, p, t) = f(x, t, P)
  diffusion(x, p, t) = sqrt(2*P.D)

  # Propagate the initial condition
  problem = SDEProblem(drift, diffusion, xs(0.0, P), (0.0, P.T))
  paths = solve(EnsembleProblem(problem), EM(), EnsembleThreads();
                trajectories=Ne, dt=P.dt, adaptive=false, saveat=P.T/500, callback=tipped)

  # Keep the samples to be drawn and the tipping times
  return (paths = paths.u[1:Np],
          tip = [path.t[end] for path in paths.u if path.t[end] < P.T])
end

# Solve the deterministic problem
function generate_pullback(P)
  # Initialize the solution
  t, x = [0.0], [x0]

  # Euler time stepping
  while t[end] < P.T && x[end] > xesc
    push!(x, x[end] + f(x[end], t[end], P)*P.dt)
    push!(t, t[end] + P.dt)
  end

  return (t = t, x = x)
end
