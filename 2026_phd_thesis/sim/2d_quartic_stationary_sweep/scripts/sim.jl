"""
    Simulation script

Storage of the definitions of the system alongside all the settings of the problem.
"""

# System parameters
μ = collect(range(0.0, stop=1.3, step=0.1))   # Bifurcation parameter value
ε = 0.0                                       # Timescale separation
σ = 0.100                                     # Noise level (additive)
D = (σ^2)/2.0                                 # Diffusion level (additive) 

# Dynamical system  
f1(x, y, μ) = -μ + 4*x - 4*x^3                # Drift for x
f2(x, y, μ) = -2*y                            # Drift for y
f = [f1, f2]                                  # Drift vector field
Λ(t) = ε                                      # Shift/Ramp
g(x, y) = σ                                   # Diffusion
η = [g, g]                                    # Diagonal diffusion

# Simulation parameters
dt = 5e-2                                     # Timestep
Nt = 2e5                                      # Total number of steps
Ne = 2e2                                      # Number of particles in the ensemble 
