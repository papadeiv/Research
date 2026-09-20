# Solves the high-dimensional, slow-fast, lattice dynamical system

function evolve_lattice(lattice::Lattice, f::Function, g::Function, h::Function, u0::Vector{Float64}; δt=5e-2, Nt=1000::Int64, saveat=δt, μf=nothing)
        # Compute the temporal quantities according to the user's choice
        if μf==nothing
                # Specify the simulation endtime according to the number of steps
                global T = δt*Nt
        else
                # Specify the simulation endtime according to the parameter range
                global T = (μf-u0[end])/h(0)
        end

        # Get the number of nodes (lattice sites)
        Ns = lattice.rows*lattice.cols

        # Deterministic dynamics 
        function drift!(dudt, u, μ, t)
                # Vector field for the lattice nodes
                for j in 1:Ns
                        dudt[j] = f(u, j, lattice)
                end
                # Parameter ramp
                dudt[Ns+1] = h(u[end])

        end
        # Stochastic dynamics
        function diffusion!(dudt, u, μ, t)
                # Vector field for the lattice nodes
                for j in 1:Ns
                        dudt[j] = g(u, j, lattice)
                end
                # Parameter ramp
                dudt[Ns+1] = 0.0::Float64
        end
        
        # Define the differential equation
        dynamics = SDEProblem(drift!, diffusion!, u0, (0.0, T))

        # Solve the SDE forward in time
        sol = solve(dynamics, EM(), dt=δt, verbose=false, saveat=saveat)

        # Extract the timestamps of the (saved) states
        time = sol.t
        Nt = length(time)

        # Extract the states of the lattice solution
        states = [(sol.u)[t][1:(end-1)] for t in 1:Nt]
        
        # Extract the states of the parameter ram[]
        parameter = [(sol.u)[t][end] for t in 1:Nt]

        # Return the solutions and their timestamps
        return time, parameter, states
end
