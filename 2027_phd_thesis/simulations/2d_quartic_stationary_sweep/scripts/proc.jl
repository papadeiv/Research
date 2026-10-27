"""
    Postprocessing script

Collection of quantities and functions used to postprocess and analyse the results of a simulation.
"""

# Parameters of the scalar potential method
Nc = convert(Int64, 3e0)                                                # Solution space dim. of the method 
Na = convert(Int64, 1e4)                                                # Number of attempts per guess 
Nx = 0::Integer                                                         # Number of escapes prior to the end of the simulation
β = 1e-3                                                                # Std of the guess perturbation 

# Scalar potential of the conservative system 
U0(x, μ) = 1 + μ*x - 2*(x^2) + x^4              # Ground truth...
V(x, c) = c[1]*x + c[2]*(x^2) + c[3]*(x^3)      # Arbitrary cubic potential...
 
# Data structures
solutions = Vector{Vector{Float64}}()           # Solutions of the inference method 
results = Vector{Float64}()                     # Statistical distribution over the ensemble

# Shift the reconstructed potential to match the local minimum of the ground truth
function shift_potential(U0::Function, parameter, c)
        # Compute the stable equilibrium (center of the shift)
        xs = +(1/(3*c[3]))*(sqrt((c[2])^2 - 3*c[1]*c[3]) - c[2])

        # Compute the shifts
        f0(x, μ) = f1(x, 0, μ)
        x0 = (get_equilibria(f0, parameter, domain=[-10,10]).stable)[2]
        δx = x0 - xs 
        δy = U0(x0, parameter) - (Polynomial([0.0; c]))(xs)

        # Define the shifted potential
        Vs(x) = δy + c[1]*(x - δx) + c[2]*(x - δx)^2 + c[3]*(x - δx)^3

        return xs, Vs
end

function analyse(solutions, parameter)
        # Reconstruct a shifted potential to match the ground truth (for error and plotting purposes)
        xs, Vs = shift_potential(U0, parameter, solutions)

        # Compute estimated stable and unstable equilibria of the cubic
        xs = +(1/(3*solutions[3]))*(sqrt((solutions[2])^2 - 3*solutions[1]*solutions[3]) - solutions[2])
        xu = -(1/(3*solutions[3]))*(sqrt((solutions[2])^2 - 3*solutions[1]*solutions[3]) + solutions[2])

        # Compute the escape rate's components
        ΔV = abs(V(xu, solutions) - V(xs, solutions))
        decay = exp(-ΔV)

        # Update the results vector
        push!(results, decay)

        # Return the shifted potential (for plotting)
        return Vs 
end
