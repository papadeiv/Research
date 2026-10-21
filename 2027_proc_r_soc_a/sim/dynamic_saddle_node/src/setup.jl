"""
    Setup

System's parameter and simulation settings.
"""

# Settings
const Ne = 1000                                      # Ensemble size
const Nm = 20                                        # Samples in the drawn median
const Ns = 100                                       # Stride of the drawn sample
const n = 10                                         # Timestep size

# Parameters
const μ0 = -1.0                                      # Initial parameter
const ε = 5e-4                                       # Ramp speed
const σ = 0.01                                       # Noise intensity
const D = σ^2/2                                      # Diffusion intensity
const τ = -μ0/ε                                      # Breaking time
const τpb = τ + 2.338107*ε^(-1/3)                    # Tipping time
const T = 1.1*τpb                                    # Termination condition
const xtip = -20.0                                   # Tipping threshold

# Admissible time horizon
const T1 = 0.2*τ
const T2 = 0.8*τ

# Rolling window
const Δ = 0.01
const Tw = 2*(τ - T2)/10
const Nw = round(Int, Tw/Δ)

# System
μ(t) = μ0 + ε*t                                      # Linear ramp
f(x, t) = -μ(t) - x^2                                # Normal form
xs(t) = sqrt(-μ(t))                                  # Stable branch
λ(t) = -2*sqrt(-μ(t))                                # Leading eigenvalue
α(t) = λ(t)^2                                        # Return rate

# Solve the ensemble problem
function generate_samples()
  # Define the tipping condition
  tipped = DiscreteCallback((x, t, integrator) -> x < xtip, terminate!)

  # Define the ramped SDE
  drift(x, p, t) = f(x, t)
  diffusion(x, p, t) = σ

  # Reduce every sample path to its return rate, its tipping time and a thinned copy
  function reduction(path, i)
    rate = estimate_return_rate(path.t, path.u)
    tip = path.t[end] < T ? path.t[end] : NaN
    return ((t = path.t[1:Ns:end], x = path.u[1:Ns:end], rate = rate, tip = tip), false)
  end

  # Propagate the initial condition
  problem = SDEProblem(drift, diffusion, xs(0.0), (0.0, T))
  paths = solve(EnsembleProblem(problem, output_func=reduction), SOSRA(), EnsembleThreads();
                trajectories=Ne, dt=Δ/n, adaptive=false, saveat=Δ, callback=tipped)

  return paths.u
end
