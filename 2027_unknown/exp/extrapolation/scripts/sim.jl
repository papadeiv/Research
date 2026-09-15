"""
    Simulation script

Storage of the definitions of the system alongside all the settings of the problem.
"""

# Simulation parameters
μ0 = -4.0                          # Initial bifurcation parameter
ε  = 1e-3                          # Parameter ramp rate
σ  = 5e-2                          # Noise level
dt = 1e-1                          # Sampling frequency
Ne = 1000                          # Number of samples
τ  = -μ0/ε                         # Breaking time
t  = τ + 2.338107*ε^(-1/3)         # Tipping time
T = 1.1*t                          # Total time 
T1, T2 = 0.20*τ, 0.80*τ            # Inference window 
mrg = 5.0                          # Stationarity margin
Tw  = 4*(-(μ0 + ε*T2))/(2ε*mrg)    # Rolling window length T_w = λ²(T2)/(2ε·mrg)
ΔT  = round(Int, Tw/dt)            
Δs  = ΔT-1#ΔT ÷ 2                  # Rolling window stride
xtip = -sqrt(-μ0)                  # Termination condition

# Dynamical system
f(x, p, t) = -(μ0 + ε*t) - x^2     # Drift
η(x, p, t) = σ                     # Diffusion
λ2_exact(t) = 4*ε*(τ - t)           # Return rate

# Generate solutions of the ensemble problem
tipped = DiscreteCallback((x, t, integrator) -> x < xtip, terminate!)
function generate_samples()
        # Define the ensemble problem
        problem = SDEProblem(f, η, sqrt(-μ0), (0.0, T))

        # Solve the ensemble problem 
        ensemble = solve(EnsembleProblem(problem), EM(), EnsembleThreads();
                         trajectories = Ne, dt = dt, adaptive = false, saveat = dt, callback = tipped)

        # Format and export the solutions
        return [(t = collect(solution.t), x = collect(solution.u)) for solution in ensemble.u]
end
