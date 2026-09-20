"""
    Simulation script

Storage of the definitions of the system alongside all the settings of the problem.
"""

# System parameters
μ = 1.20                                                   # Set of bifurcation parameter values 
ε = 0.0                                                    # Timescale separation
σ = 0.050                                                  # Noise level (additive)
D = (σ^2)/2.0                                              # Diffusion level (additive) 

# Dynamical system  
f(x, μ) = -μ + 4*x - 4*x^3                                 # Drift
J(x) = + 4 - 12*x^2                                        # Jacobian 
Λ(t) = ε                                                   # Shift/Ramp
η(x) = σ                                                   # Diffusion

# Simulation parameters
dt = exp10.(range(log10(0.2), log10(0.001), length = 150)) # Timestep size
Nt = round.(Integer, exp10.(range(2, 5, length = 150)))    # Number of timesteps
Ns = length(dt)*length(Nt)                                 # Number of simulations
Ne = 1e+2                                                  # Number of particles

# Compute solutions and error approximations of the parameter estimation problems
function generate_samples(stepsize, steps)
        # Compute the relevant equilibria
        equilibria = get_equilibria(f, μ, domain=[-10,10])
        xu = equilibria.unstable[1] 
        xs = equilibria.stable[2]

        # Solve the ensemble problem
        x0 = [xs, μ]
        ensemble = evolve(f, η, Λ, x0, stepsize=stepsize, steps=steps, particles=Ne)

        # Ensemble loop
        data = Matrix{Float64}(undef, convert(Integer, Ne), 3)
        for (idx_sol, solution) in enumerate(ensemble.state)
                # Truncate and detrend the solution
                idx_tip = find_tipping(solution, xu)
                u = solution[1:idx_tip]
                u = u .- x0[1]

                # Estimate the return rate of the samples
                α1 = estimate_return_rate(u, stepsize)
                α2 = estimate_return_rate(u, stepsize, alg = "AC1")
                α3 = estimate_return_rate(u, stepsize, alg = "OUP", σ = σ)

                # Update the data matrix
                data[idx_sol, :] = [α1, α2, α3]
        end

        # Export the data
        writeout(data, "accuracy/$stepsize/$steps.csv")
end
