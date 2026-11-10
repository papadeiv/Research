"""
    Setup

Settings for the observation protocol.
"""

# Admissible time horizon
const T1 = 200.0
const T2 = 1600.0

# Tipping time (ground truth)
const τ_tip = 1758.0

# Rolling window
const Δ = 1.0
const Tw = 25.0
const Nw = round(Int, Tw/Δ)
