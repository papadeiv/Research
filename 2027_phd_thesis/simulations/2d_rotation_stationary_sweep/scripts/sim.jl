"""
    Simulation script

Storage of the definitions of the system alongside all the settings of the problem.
"""

# System parameters
μ = collect(range(0.0, stop=1.3, step=0.1))   # Bifurcation parameter value
ε = 0.0                                       # Timescale separation
σ = 0.100                                     # Noise level (additive)
D = (σ^2)/2.0                                 # Diffusion level (additive) 

# Rotation coefficients
α = pi/4 
β = pi/2 
a11 = cos(α) 
a12 = sin(α) 
a21 = -sin(β) 
a22 = cos(β) 
A = [a11 a12
     a21 a22]

# Dynamical system: drift 
f1(x, y, μ) = a11*(-μ + 4*(a11*x+a12*y) - 4*(a11*x+a12*y)^3) - 2*a21*(a21*x+a22*y)
f2(x, y, μ) = a12*(-μ + 4*(a11*x+a12*y) - 4*(a11*x+a12*y)^3) - 2*a22*(a21*x+a22*y)
f = [f1, f2]

# Dynamical system: shift 
Λ(t) = ε

# Dynamical system: diffusion 
g(x, y) = σ
η = [g, g]

# Simulation parameters
dt = 5e-2                                     # Timestep
Nt = 1e6                                      # Total number of steps
Ne = 1e2                                      # Number of particles in the ensemble 
