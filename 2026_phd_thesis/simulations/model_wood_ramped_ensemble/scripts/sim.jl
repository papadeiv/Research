"""
    Simulation script

Storage of the definitions of the system alongside all the settings of the problem.
"""

# Model parameters
S0 = 0.035

α = 0.12
β = 790.0
Υ = 100*3.15*1e7

Sn_eq = 0.034912
St_eq = 0.035435
Ss_eq = 0.034427
Sip_eq = 0.034668
Sb_eq = 0.034538

Vn = 3.683*1e16
Vt = 5.418*1e16
Vs = 6.097*1e16
Vip = 1.486*1e17
Vb = 9.925*1e17

Ts = 7.919
T0 = 3.87

C = Vn*Sn_eq + Vt*St_eq + Vip*Sip_eq + Vs*Ss_eq + Vb*Sb_eq

λ = 1.62*1e7
γ = 0.36                                           # Control parameter (location of HB_sub)
μ = 22*1e-8

Kn = 1.762*1e6
Ks = 1.872*1e6

Fn(H) = (0.486 + 0.131*H)*1e6
Ft(H) = (-0.997 + 0.696*H)*1e6

Sn0 = (Sn_eq - S0)*100.0
St0 = (St_eq - S0)*100.0
Ss = (Ss_eq - S0)*100.0 
Sb = (Sb_eq - S0)*100.0 
Sip(Sn, St) = 100.0*(C - (Vn*Sn + Vs*Ss + Vb*Sb + Vt*St)/100.0 - S0*(Vb + Vn + Vt + Vip + Vs))/Vip

q(Sn) = λ*(α*(Ts - T0) + (β/100.0)*(Sn - Ss))/(1.0 + λ*α*μ)

# System parameters
H0 = collect(range(-0.38, stop=0.38, step=0.1))    # Collection of initial bifurcation parameters
ε = 1e-6                                           # Slow timescale
σ = 1e-3                                           # Noise level (additive)

# Dynamical system: drift 
f1_pos(Sn, St, H) = (Υ/Vn)*(q(Sn)*(St/100.0 - Sn/100.0) + Kn*(St/100.0 - Sn/100.0) - Fn(H)*S0)
f1_neg(Sn, St, H) = (Υ/Vn)*(abs(q(Sn))*(Sb/100.0 - Sn/100.0) + Kn*(St/100.0 - Sn/100.0) - Fn(H)*S0)
f1(Sn, St, H) = q(Sn) >= 0 ? f1_pos(Sn, St, H) : f1_neg(Sn, St, H)
f2_pos(Sn, St, H) = (Υ/Vt)*(q(Sn)*(γ*Ss/100.0 + (1 - γ)*Sip(Sn, St)/100.0 - St/100.0) + Ks*(Ss/100.0 - St/100.0) + Kn*(Sn/100.0 - St/100.0) - Ft(H)*S0)
f2_neg(Sn, St, H) = (Υ/Vt)*(abs(q(Sn))*(Sn/100.0 - St/100.0) + Ks*(Ss/100.0 - St/100.0) + Kn*(Sn/100.0 - St/100.0) - Ft(H)*S0)
f2(Sn, St, H) = q(Sn) >= 0 ? f2_pos(Sn, St, H) : f2_neg(Sn, St, H)
f = [f1, f2]

# Dynamical system: shift
Λ(t) = ε

# Dynamical system: diffusion
g(x, y) = σ
η = [g, g]

# Simulation parameters
δt = 5e-2                                     # Timestep
Nt = 1e5                                      # Total number of steps
Ne = 1e0                                      # Number of particles in the ensemble 
