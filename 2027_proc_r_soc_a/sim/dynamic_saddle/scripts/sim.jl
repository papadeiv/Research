"""
    Simulation script

Storage of the definitions of the system alongside all the settings of the problem.
"""

# System parameters
μ0 = -1.0                                         # Initial value of the bifurcation parameter
μf = 0.5                                          # Final value of the bifurcation parameter
ε = 3e-2                                          # Timescale separation
σ = 0.100                                         # Noise level (additive)
D = (σ^2)/2.0                                     # Diffusion level (additive) 
τ = -μ0/ε                                         # Breaking time
T = -μ0/ε + 1.018793*ε^(-1/3)                     # Tipping time (point of no-return)

# Dynamical system  
f(x, μ) = -μ - x^2                                # Drift
J(x) = -2*x                                       # Jacobian 
Λ(t) = ε                                          # Shift/Ramp
η(x) = σ                                          # Diffusion

# Simulation parameters
Δt = 1e-2                                         # Timestep of the solver
dt = 1e-1                                         # Timestep of the solution
Nt = 1e+4                                         # Number of steps 
Ne = 1e+1                                         # Number of particles in the ensemble 
