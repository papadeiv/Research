"""
    Postprocessing script

Collection of quantities and functions used to postprocess and analyse the results of a simulation.
"""

# Parameters of the scalar potential method
window_size = 0.25                              # Relative width of the slinding window
idx = 1800                                      # Time index of the tipping point 
Na = convert(Int64, 1e4)                        # Number of attempts per guess 
β = 1e-2                                        # Std of the guess perturbation 

# Arbitrary cubic potential
V(x, c) = c[1]*x + c[2]*(x^2) + c[3]*(x^3)

function analyse(solution)
        # Compute estimated stable and unstable equilibria of the cubic
        xs = +(1/(3*solution[3]))*(sqrt((solution[2])^2 - 3*solution[1]*solution[3]) - solution[2])
        xu = -(1/(3*solution[3]))*(sqrt((solution[2])^2 - 3*solution[1]*solution[3]) + solution[2])

        # Compute the ews
        ΔV = abs(V(xu, solution) - V(xs, solution))
        return exp(-ΔV)
end
