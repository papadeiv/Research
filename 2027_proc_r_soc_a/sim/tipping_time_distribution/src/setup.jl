"""
    Setup

System's parameter and simulation settings.
"""

# Settings
const Ne = 10_000                                    # Ensemble size
const Np = 100                                       # Solutions drawn per panel
const x0 = 2.0                                       # Initial condition
const xesc = -2.5                                    # Escape threshold
const xtip = -20.0                                   # Tipping threshold

# Parameters: stationary, small noise and large noise experiments
const runs = ((μ0 = -1.0, ε = 0.0, D = 0.25, dt = 2e-2, T = 4.0e3),
              (μ0 = -1.0, ε = 1e-2, D = 1e-3, dt = 5e-3, T = 1.25e2),
              (μ0 = -1.0, ε = 1e-3, D = 0.0625, dt = 1e-2, T = 1.1e3))

# System
f(x, t, P) = -(P.μ0 + P.ε*t) - x^2                   # Normal form
xs(P) = sqrt(-P.μ0)                                  # Stable branch
τ(P) = P.ε > 0 ? -P.μ0/P.ε : P.T                     # Breaking time
τpb(P) = τ(P) + 2.338107*P.ε^(-1/3)                  # Tipping time

# Solve the ensemble problem
function generate_samples(P)
  # Define the escape condition
  tipped = DiscreteCallback((x, t, integrator) -> x < xtip, terminate!)

  # Define the ramped SDE
  drift(x, p, t) = f(x, t, P)
  diffusion(x, p, t) = sqrt(2*P.D)

  # Propagate the initial condition
  problem = SDEProblem(drift, diffusion, xs(P), (0.0, P.T))
  paths = solve(EnsembleProblem(problem), EM(), EnsembleThreads();
                trajectories=Ne, dt=P.dt, adaptive=false,
                saveat=P.T/500, callback=tipped)

  return paths
end

# Solve the deterministic problem by Runge-Kutta
function generate_pullback(P)
  # Initialize the solution at the initial condition
  t, x = [0.0], [x0]

  # Propagate until the termination condition or the escape threshold
  for i in 1:round(Int, P.T/P.dt)
    k1 = f(x[end], t[end], P)
    k2 = f(x[end] + P.dt*k1/2, t[end] + P.dt/2, P)
    k3 = f(x[end] + P.dt*k2/2, t[end] + P.dt/2, P)
    k4 = f(x[end] + P.dt*k3, t[end] + P.dt, P)
    xn = x[end] + P.dt*(k1 + 2*k2 + 2*k3 + k4)/6

    xn < xesc && break
    push!(x, xn)
    push!(t, t[end] + P.dt)
  end

  return t, x
end
