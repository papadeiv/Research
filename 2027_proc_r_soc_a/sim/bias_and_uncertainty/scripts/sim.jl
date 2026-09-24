"""
    Simulation script

Storage of the definitions of the system alongside all the settings of the problem.
"""

# System parameters
Nμ = 100                                                  # Number of parameter values in the sweep
μ_set = collect(LinRange(0.0, 1.45, Nμ))                  # Set of bifurcation parameter values 
ε = 0.0                                                   # Timescale separation
σ = 0.050                                                 # Noise level (additive)
D = (σ^2)/2.0                                             # Diffusion level (additive) 

# Dynamical system  
f(x, μ) = -μ + 4*x - 4*x^3                                # Drift
J(x) = + 4 - 12*x^2                                       # Jacobian 
Λ(t) = ε                                                  # Shift/Ramp
η(x) = σ                                                  # Diffusion

# Simulation parameters
dt = 1e-1                                                 # Timestep size
Nt = 1e+3                                                 # Number of timesteps
Ne = 1e+2                                                 # Number of particles

# Compute solutions and error approximations of the parameter estimation problems
function generate_samples()
        # Initialize stop flag
        stop_flag = false 

        # Parameter sweep loop
        printstyled("Generating the samples\n"; bold=true, underline=true, color=:light_blue)
        @showprogress for (idx_μ, μ) in enumerate(μ_set)
                # Compute the relevant equilibria
                equilibria = get_equilibria(f, μ, domain=[-10,10])
                xu = equilibria.unstable[1] 
                xs = equilibria.stable[2]

                # Solve the ensemble problem
                x0 = [xs, μ]
                ensemble = evolve(f, η, Λ, x0, stepsize=dt, steps=Nt, particles=Ne)

                # Ensemble loop
                data = Matrix{Float64}(undef, convert(Integer, Ne), 3)
                for (idx_sol, solution) in enumerate(ensemble.state)
                        # Truncate and detrend the solution
                        idx_tip = find_tipping(solution, xu)
                        u = solution[1:idx_tip]
                        u = u .- x0[1]

                        # Stop the simulation if a single sample path has tipped
                        if idx_tip < Nt
                                stop_flag = true 
                                @goto tipped 
                        end

                        # Estimate the return rate of the samples
                        α1 = estimate_return_rate(u, dt)
                        α2 = estimate_return_rate(u, dt, alg = "AC1")
                        α3 = estimate_return_rate(u, dt, alg = "OUP", σ = σ)

                        # Update the data matrix
                        data[idx_sol, :] = [α1, α2, α3]
                end

                # Export the data
                writeout(data, "solutions/$idx_μ.csv")
        end

        # Return the flag
        @label tipped
        return stop_flag
end
