"""
    Main script

Run this file to execute the simulation, analyse and plot the results.
"""

using CairoMakie
using StochasticDiffEq

# Import local functions
include("src/setup.jl")
include("src/plotting.jl")

# Main algorithm
function main()
  # Compute the equilibria branches 
  μ = LinRange(-1,-1e-16,500)
  a = [sqrt(-p) for p in μ]
  b = [-sqrt(-p) for p in μ]

  lines!(ax_left, μ, a, color = :blue, linewidth = 2)
  lines!(ax_left, μ, b, color = :red, linewidth = 2)
  lines!(ax_right, μ, a, color = :blue, linewidth = 2)
  lines!(ax_right, μ, b, color = :red, linewidth = 2)

  # Compute the sample paths
  # ------------------------ Case 1: singular limit (fast)
  x0, μ0 = 1.5, -0.8
  ε, σ = 0.0, 0.0
  xf = DiscreteCallback((x, t, integrator) -> x < sqrt(-μ0) + 1e-3, terminate!)
  path = generate_samples([x0,μ0], ε, σ, xf)
  t = path.t
  x = path.u
  μ = [μ0 for s in t]

  lines!(ax_left, μ, x, color = :grey, linewidth = 2.0)

  # ------------------------ Case 2: singular limit (slow)
  μ = LinRange(μ0,-1e-16,500) 
  x = [sqrt(-p) for p in μ] 
  lines!(ax_left, μ, x, color = :grey, linewidth = 2.0)
  display(x[end])
  display(μ[end])

  # ------------------------ Case 3: singular limit (fast)
  x0, μ0 = -1e-3, 0.0 
  ε, σ = 0.0, 0.0
  xf = DiscreteCallback((x, t, integrator) -> x < -1.0, terminate!)
  path = generate_samples([x0,μ0], ε, σ, xf)
  t = path.t
  x = path.u
  μ = [μ0 for s in t]
  display(x)

  lines!(ax_left, μ, x, color = :grey, linewidth = 2.0)

  # ------------------------ Case 4: pullback attractor 
  x0, μ0 = 1.5, -0.8
  ε, σ = 3e-2, 0.0
  xf = DiscreteCallback((x, t, integrator) -> x < -1.0, terminate!)
  path = generate_samples([x0,μ0], ε, σ, xf)
  t = path.t
  x = path.u
  μ = [μ0 + ε*s for s in t]

  lines!(ax_right, μ, x, color = :orange, linewidth = 3.0)

  # ------------------------ Case 5: slow-fast sample path 
  x0, μ0 = 1.5, -0.8
  ε, σ = 3e-2, 5e-2
  xf = DiscreteCallback((x, t, integrator) -> x < -1.0, terminate!)
  path = generate_samples([x0,μ0], ε, σ, xf)
  t = path.t
  x = path.u
  μ = [μ0 + ε*s for s in t]

  lines!(ax_right, μ, x, color = :grey, linewidth = 2.0)

  #display(fig)
  save("slow_fast_model.pdf", fig)
end

# Execute main
main()
