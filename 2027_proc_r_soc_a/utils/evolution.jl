"""
Propagating solutions of a system of nonautonomous, stochastic differential equations forward in time.

Author: Davide Papapicco
Affil: U. of Auckland
Date: 21-08-2025
"""

#___________________________#
#                           #
#                           #
#   1-dimensional systems   # 
#                           #
#___________________________#

function evolve(f::Function, η::Function, Λ::Function, u0::Vector; endparameter=nothing, timerange=nothing, stepsize=1e-2, saveat=stepsize, steps=nothing, particles=nothing)
        # Initialise internal quantities
        T_start = 0.0::Real
        T_end = 1.0::Real
        Nt = 0::Integer
        δt = stepsize
        Ne = 1::Int64
        boolean_placeholder = particles == nothing ? Ne = 1 : Ne = convert(Int64, particles)

        # Check input's sanity (solution's range)
        if timerange == nothing
                if endparameter == nothing
                        if steps == nothing
                                throw(ArgumentError("You have to specify either the time interval or the parameter range of the solution!"))
                        elseif steps isa AbstractFloat
                                Nt = convert(Int64, steps)
                                T_start = 0.0::Real
                                T_end = δt*Nt
                        elseif steps isa Integer
                                Nt = steps
                                T_start = 0.0::Real
                                T_end = δt*Nt
                        else
                                throw(ArgumentError("The total number of steps can either be an AbstractFloat or an Integer, got $(typeof(steps))"))
                        end
                else
                        T_start = 0.0::Real 
                        T_end = (endparameter-u0[2])/Λ(0)
                end
        else
                # Compute the time range from the parameter range assuming a linear ramp
                if timerange isa AbstractFloat
                        T_start = 0.0::Real
                        T_end = timerange 
                elseif timerange isa Vector{<:Real}
                        T_start = minimum(timerange)
                        T_end = maximum(timerange)
                else
                        throw(ArgumentError("The time interval is either a Real or a Vector{<:Real}, got $(typeof(timerange))"))
                end
        end

        # Compute the total number of timesteps
        if steps == nothing
                Nt = (T_end - T_start)/δt
        elseif steps isa AbstractFloat
                Nt = convert(Int64, steps)
        else
                Nt = steps
        end

        # Deterministic dynamics 
        function drift!(dudt, u, μ, t)
                dudt[1] = f(u[1], u[2])
                dudt[2] = Λ(u[2])
        end
        # Stochastic dynamics
        function diffusion!(dudt, u, μ, t)
                dudt[1] = η(u[1])
                dudt[2] = 0 
        end
        
        # Define the differential equation
        dynamics = SDEProblem(drift!, diffusion!, u0, (T_start, T_end))

        # Convert it into an ensemble problem
        ensemble = EnsembleProblem(dynamics)

        # Solve the problem forward in time
        solutions = solve(ensemble, EM(), dt=δt, EnsembleDistributed(), trajectories=Ne, verbose = None())

        # Extract the trajectories of the state 
        state = Vector{Vector{Float64}}() 
        parameter_length = Vector{Int64}()
        time_length = Vector{Int64}()
        for n in 1:Ne
                time = solutions.u[n].t 
                push!(time_length, length(time))

                parameter = [(solutions.u[n].u)[t][2] for t in 1:length(time)]
                push!(parameter_length, length(parameter))

                push!(state, [(solutions.u[n].u)[t][1] for t in 1:length(time)])
        end

        # Extract the timestamps and parameter realisations
        longest_run = findmax(time_length)[2]
        time = solutions.u[longest_run].t
        parameter = [(solutions.u[longest_run].u)[t][2] for t in 1:length(solutions.u[longest_run].u)]

        return (
                time = time,
                parameter = parameter,
                state = state
               )
end

#___________________________#
#                           #
#                           #
#   2-dimensional systems   # 
#                           #
#___________________________#

function evolve(f::Vector{<:Function}, η::Vector{<:Function}, Λ::Function, u0::Vector{Float64}; endparameter=nothing, timerange=nothing, stepsize=1e-2, saveat=stepsize, steps=nothing, particles=nothing)
        # Initialise internal quantities
        T_start = 0.0::Real
        T_end = 1.0::Real
        Nt = 0::Integer
        δt = stepsize
        Ne = 1::Int64
        boolean_placeholder = particles == nothing ? Ne = 1 : Ne = convert(Int64, particles)

        # Check input's sanity (solution's range)
        if timerange == nothing
                if endparameter == nothing
                        if steps == nothing
                                throw(ArgumentError("You have to specify either the time interval or the parameter range of the solution!"))
                        elseif steps isa AbstractFloat
                                Nt = convert(Int64, steps)
                                T_start = 0.0::Real
                                T_end = δt*Nt
                        elseif steps isa Integer
                                Nt = steps
                                T_start = 0.0::Real
                                T_end = δt*Nt
                        else
                                throw(ArgumentError("The total number of steps can either be an AbstractFloat or an Integer, got $(typeof(steps))"))
                        end
                else
                        T_start = 0.0::Real 
                        T_end = (endparameter-u0[3])/Λ(0)
                end
        else
                # Compute the time range from the parameter range assuming a linear ramp
                if timerange isa AbstractFloat
                        T_start = 0.0::Real
                        T_end = timerange 
                elseif timerange isa Vector{<:Real}
                        T_start = minimum(timerange)
                        T_end = maximum(timerange)
                else
                        throw(ArgumentError("The time interval is either a Real or a Vector{<:Real}, got $(typeof(timerange))"))
                end
        end

        # Compute the total number of timesteps
        if steps == nothing
                Nt = (T_end - T_start)/δt
        elseif steps isa AbstractFloat
                Nt = convert(Int64, steps)
        else
                Nt = steps
        end

        # Deterministic dynamics 
        function drift!(dudt, u, μ, t)
                dudt[1] = f[1](u[1], u[2], u[3])
                dudt[2] = f[2](u[1], u[2], u[3])
                dudt[3] = Λ(u[3])
        end
        # Stochastic dynamics
        function diffusion!(dudt, u, μ, t)
                dudt[1] = η[1](u[1], u[2])
                dudt[2] = η[2](u[1], u[2])
                dudt[3] = 0 
        end
      
        # Define the differential equation
        dynamics = SDEProblem(drift!, diffusion!, u0, (T_start, T_end), Λ)

        # Convert it into an ensemble problem
        ensemble = EnsembleProblem(dynamics)

        # Solve the problem forward in time
        solutions = solve(ensemble, EM(), dt=δt, EnsembleDistributed(), trajectories=Ne)

        # Extract the trajectories of the state 
        state = Vector{Matrix{Float64}}() 
        parameter_length = Vector{Int64}()
        time_length = Vector{Int64}()
        for n in 1:Ne
                time = solutions[n].t 
                push!(time_length, length(time))

                parameter = [(solutions[n].u)[t][3] for t in 1:length(solutions[n].u)]
                push!(parameter_length, length(parameter))

                x = [(solutions[n].u)[t][1] for t in 1:length(solutions[n].u)]
                y = [(solutions[n].u)[t][2] for t in 1:length(solutions[n].u)]
                push!(state, hcat(x, y))
        end

        # Extract the timestamps and parameter realisations
        longest_run = findmax(time_length)[2]
        time = solutions[longest_run].t
        parameter = [(solutions[longest_run].u)[t][3] for t in 1:length(solutions[longest_run].u)]

        return (
                time = time,
                parameter = parameter,
                state = state
               )
end
