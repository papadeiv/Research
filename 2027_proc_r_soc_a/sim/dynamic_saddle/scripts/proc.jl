"""
    Postprocessing script

Collection of quantities and functions used to postprocess and analyse the results of a simulation.
"""

# Compute the return rate of the pullback attractor
function pullback_rate(t)
        z = -(μ0 + ε*t)/ε^(2/3)
        airy_ratio(z) = z > 0 ? airyaiprimex(z)/airyaix(z) : airyaiprime(z)/airyai(z)
        return 2*ε^(1/3)*airy_ratio(z)
end
