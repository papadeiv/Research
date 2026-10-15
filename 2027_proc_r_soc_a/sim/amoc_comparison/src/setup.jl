"""
    Setup

System's parameter and simulation settings.
"""

# Settings
const Ne = 1000                                      # Ensemble size
const Np = 15                                        # Solutions drawn
const Nt = 1500                                      # Points drawn per solution
const n = 10                                         # Timestep size
const mrg = 10.0                                     # Stationarity margin

# Parameters
const μ0 = -1.0                                      # Initial parameter
const ε = 5e-4                                       # Ramp speed
const σ = 0.02                                       # Noise intensity
const Δ = 0.01                                       # Sampling frequency
const τ = -μ0/ε                                      # Breaking time
const τpb = τ + 2.338107*ε^(-1/3)                    # Tipping time
const T = 1.1*τpb                                    # Termination condition
const xesc = -sqrt(-μ0)                              # Escape threshold

# Rolling windows
const T1 = 0.2*τ                                     # Start of the inference window
const T2 = 0.8*τ                                     # End of the inference window
const Tw = 4*(-(μ0 + ε*T2))/(2*ε*mrg)                # Window length
const Nw = round(Int, Tw/Δ)                          # Window size
const starts = collect((round(Int, T1/Δ) + 1):Nw:(round(Int, T2/Δ) + 2 - Nw))
const tw = [(s - 1 + Nw ÷ 2)*Δ for s in starts]      # Window midpoints

# System
μ(t) = μ0 + ε*t                                      # Linear ramp
f(x, t) = -μ(t) - x^2                                # Normal form
xs(t) = sqrt(-μ(t))                                  # Stable branch
λ(t) = -2*sqrt(-μ(t))                                # Leading eigenvalue
α(t) = λ(t)^2                                        # Return rate

# Solve the ensemble problem
function generate_samples()
  # Define the escape condition
  escaped = DiscreteCallback((x, t, integrator) -> x < xesc, terminate!)

  # Define the ramped SDE
  drift(x, p, t) = f(x, t)
  diffusion(x, p, t) = σ

  # Reduce every solution to its summary as soon as it is computed
  reduction(path, i) = (summarise_sample(path), false)

  # Propagate the initial condition
  problem = SDEProblem(drift, diffusion, xs(0.0), (0.0, T))
  paths = solve(EnsembleProblem(problem, output_func=reduction), SOSRA(), EnsembleThreads();
                trajectories=Ne, dt=Δ/n, adaptive=false, saveat=Δ, callback=escaped)

  return paths.u
end
