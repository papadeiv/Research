"""
    Postprocessing script

Collection of quantities and functions used to postprocess and analyse the results of a simulation.
"""

# Parameters of the scalar potential method
window_size = 0.195                             # Relative width of the slinding window
idx = 1800                                      # Time index of the tipping point 
Na = convert(Int64, 1e4)                        # Number of attempts per guess 
β = 1e-2                                        # Std of the guess perturbation 

# Arbitrary cubic potential
V(x, c) = c[1]*x + c[2]*(x^2) + c[3]*(x^3)
xs(C) = +(1/(3*C[3]))*(sqrt((C[2])^2 - 3*C[1]*C[3]) - C[2])
xu(C) = -(1/(3*C[3]))*(sqrt((C[2])^2 - 3*C[1]*C[3]) + C[2])

# Converts the non-stationary timeseries into an ensemble of subseries associated to the strides of a sliding window 
function preprocess_solution(timestamps, timeseries, width)
        # Extract the subseries up to the tipping
        t = timestamps[1:idx] 
        u = timeseries[1:idx]
        Nt = length(u)

        # Assemble the sliding window
        window = build_window(Nt, width)
        Nw = window.size 
        Ns, Nr = divrem(Nt, Nw)

        # Convert the sliding window subseries into an ensemble of timeseries
        timesteps = [@view t[(n+(n-1)*Nw):((n+(n-1)*Nw)+Nw-1)] for n in 1:Ns] 
        ensemble = [@view u[(n+(n-1)*Nw):((n+(n-1)*Nw)+Nw-1)] for n in 1:Ns] 

        # Export the parameters of the ensemble problem 
        return (
                timesteps = timesteps,
                trajectories = ensemble 
               ) 
end

function analyse(solution)
        # Compute estimated stable and unstable equilibria of the cubic
        xs = +(1/(3*solution[3]))*(sqrt((solution[2])^2 - 3*solution[1]*solution[3]) - solution[2])
        xu = -(1/(3*solution[3]))*(sqrt((solution[2])^2 - 3*solution[1]*solution[3]) + solution[2])

        # Compute the ews
        ΔV = abs(V(xu, solution) - V(xs, solution))
        return exp(-ΔV)
end
